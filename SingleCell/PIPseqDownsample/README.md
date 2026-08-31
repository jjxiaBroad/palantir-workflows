# PIPseq Downsampling Pipelines

Nextflow DSL2 pipelines for downsampling already-processed PIPseq single-cell RNA-seq samples for cross-sample or cross-depth comparison.

Both entrypoints operate on already-published DRAGEN scRNA output + gRNA AnnData (from `../PIPseqPipeline`, or any other pipeline producing the same shapes) — neither runs DRAGEN itself, so both are cheap to rerun repeatedly at different depths.

There are two entrypoints, sharing a common per-sample downsampling engine (`workflows/downsample_core.nf`):
- **`downsample.nf`** — normalizes a *batch* of samples to one common sequencing depth per modality, then combines them into one AnnData for cross-sample comparison. Validated against `nextflow_schema.json`.
- **`downsample_single_sample.nf`** — downsamples *one* sample to several caller-specified target depths (e.g. to build a within-sample depth-response comparison), then combines the results across depths into one AnnData. Also supports an optional GEX sequencing-saturation sweep. Validated against `nextflow_schema_single_sample.json`.

The pipelines are designed to run on Illumina Connected Analytics (ICA) but are plain Nextflow DSL2 and can run anywhere a compatible executor and the required container image are available.

## Overview

### Batch downsample-and-combine (`downsample.nf`)

Given a `--samplesheet` (columns: `sample_id, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv, crispr_h5ad`, plus any per-sample metadata columns) describing a batch of samples:

1. Resolves one common GEX target depth and one common gRNA target depth for the whole batch — each independently defaults to the minimum "mean reads per cell" observed across the batch (you can only downsample down, never up), or an explicit `--gex_target_depth`/`--crispr_target_depth` override. Published as `gex_downsample_summary.csv`/`grna_downsample_summary.csv` (useful on their own, as a batch depth-landscape report).
2. Downsamples every sample's GEX molecule-info (raw and/or filtered matrix layout per `--gex_matrix_format`, via `bin/downsample_molecule_info.py`'s IMI-based downsampler) and/or gRNA AnnData (via `bin/downsample_crispr_anndata.py`'s binomial-thinning downsampler) to those targets. If `--run_saturation` is set, also runs a GEX sequencing-saturation sweep per sample (see below).
3. Combines every sample's downsampled GEX + gRNA data into one AnnData (`bin/combine_downsampled_batch.py`) — only if both branches ran (`--run_combine`, default `true`). Barcodes between a sample's filtered GEX and gRNA are expected to already match; a mismatch is reported as a warning (only the shared barcodes are kept) rather than failing the run.

### Single-sample multi-depth (`downsample_single_sample.nf`)

Given one sample's `--input_files` (molecule_info_h5/filtered_barcodes_tsv/scrna_metrics_csv and/or crispr_h5ad, role auto-detected per file — see below), plus one or more `--gex_target_depths`/`--crispr_target_depths`:

1. Downsamples the sample's GEX molecule-info and/or gRNA AnnData to *every* requested depth (one downsampled matrix/AnnData per depth), via the same two scripts `downsample.nf` uses. If `--run_saturation` is set, also runs a GEX sequencing-saturation sweep (see below).
2. Combines the downsampled GEX + gRNA data *across depths* into one AnnData (`bin/combine_downsampled_depths.py`), annotated with integer `depth` and `modality` columns — only if both branches ran (`--run_combine`, default `true`). Every depth either modality produced is included: `--gex_target_depths` and `--crispr_target_depths` are independent and need not line up. A depth carrying both modalities is merged on the feature axis (`modality == "gex+grna"`); a depth carrying only one has the other's features zero-filled (`modality == "gex"` / `"grna"`), so filter on `.obs['modality']` before treating those zeros as measured counts.

### GEX sequencing-saturation sweep (`--run_saturation`, both entrypoints)

`bin/downsample_molecule_info.py` can optionally also run a sweep over a built-in reads-per-cell ladder (plus any `--saturation_extra_depths`), writing `saturation.csv` (median transcripts/genes per cell and % sequencing saturation at each depth) alongside the sample's downsampled matrices, without materializing full matrices at every ladder depth. This is GEX-only — gRNA/CRISPR Guide Capture reads have no meaningful molecule deduplication (see `bin/downsample_crispr_anndata.py`'s docstring), so "sequencing saturation" isn't a meaningful concept for that modality. `--run_saturation` requires `--run_gex_downsample`.

## Pipeline Structure

```
SingleCell/PIPseqDownsample/
├── downsample.nf                     # Batch downsample-and-combine entrypoint (--samplesheet)
├── downsample_single_sample.nf       # Single-sample multi-depth entrypoint
├── workflows/
│   └── downsample_core.nf             # Shared per-sample downsampling engine called by both entrypoints
├── nextflow.config                   # Pipeline configuration (default params, resources, reports)
├── nextflow_schema.json              # Parameter schema for downsample.nf
├── nextflow_schema_single_sample.json # Parameter schema for downsample_single_sample.nf
├── modules/
│   ├── summarize_gex_downsample_targets.nf   # Resolve the batch's common GEX target depth (downsample.nf only)
│   ├── summarize_grna_downsample_targets.nf  # Resolve the batch's common gRNA target depth (downsample.nf only)
│   ├── downsample_molecule_info.nf   # Downsample one sample's GEX molecule-info to one or more depths, +optional saturation sweep
│   ├── downsample_crispr_anndata.nf  # Downsample one sample's gRNA AnnData to one or more depths
│   ├── combine_downsampled_batch.nf  # Combine a batch's downsampled GEX+gRNA into one AnnData (downsample.nf only)
│   └── combine_downsampled_depths.nf # Combine one sample's downsampled GEX+gRNA across depths into one AnnData (downsample_single_sample.nf only)
├── bin/
│   ├── downsample_molecule_info.py   # GEX molecule-info downsampler (IMI-based) + saturation sweep
│   ├── downsample_crispr_anndata.py  # gRNA AnnData downsampler (binomial thinning)
│   ├── summarize_gex_downsample_targets.py   # GEX batch-target summary/resolution script
│   ├── summarize_grna_downsample_targets.py  # gRNA batch-target summary/resolution script
│   ├── downsample_combine_common.py  # Shared GEX+gRNA merge helpers used by both combine scripts below
│   ├── combine_downsampled_batch.py  # GEX+gRNA batch-combine script (across samples, one depth)
│   └── combine_downsampled_depths.py # GEX+gRNA depth-combine script (across depths, one sample)
├── docker/                          # Dockerfile/build scripts for the qc_container image
├── stub_test/                       # Example inputs + `-stub-run` test setup
├── test/                            # Flat pipeline-inputs JSON for each entrypoint, used by ica_tools/start_analysis.py
├── ica_tools/                       # Scripts for publishing/running these pipelines on ICA
│   ├── export_pipeline_to_ica.py     # Imports the current commit into ICA as a git-backed pipeline
│   ├── start_analysis.py             # Interactively starts an ICA analysis run
│   ├── ica_common.py                 # Shared helpers (API key, project list, prompt_choice())
│   └── inputforms/<downsample|single_sample>/inputForm.json  # Hand-maintained ICA launch-form definitions
└── README.md                        # This file
```

## Quick Start

### Prerequisites

- Nextflow (>= 25.10.0 — required by the pinned `nf-schema` validation plugin; requires network access to the Nextflow plugin registry on first run)
- A QC container image (`--qc_container`) — built from `docker/qc/Dockerfile`, see `docker/SETUP.md`
- An executor/environment that can run the `container` directive (e.g. Nextflow's k8s executor on ICA, or Docker/Singularity enabled locally via your own config)

There is no bundled Nextflow profile for local/container-less execution — every process declares a `container`, so running for real requires a container-capable executor.

### Stub run (no container required)

Every process has a `stub:` block that just touches placeholder output files, so you can validate a pipeline's wiring without running any real downsampling:

```bash
cd stub_test
nextflow run ../downsample.nf -stub-run -params-file stub_inputs/pipeline_input.json

# Single-sample multi-depth entrypoint
nextflow run ../downsample_single_sample.nf -stub-run -params-file stub_inputs/pipeline_input_single_sample.json
```

### Running downsample.nf (batch)

```bash
nextflow run downsample.nf \
  --samplesheet samplesheet.csv \
  --staged_files sample1/sample1.scRNA.moleculeInfo.h5,sample1/sample1.scRNA.filtered.barcodes.tsv.gz,sample1/sample1.scRNA_metrics.csv,sample1/adata/sample1.crispr.h5ad,sample2/sample2.scRNA.moleculeInfo.h5,sample2/sample2.scRNA.filtered.barcodes.tsv.gz,sample2/sample2.scRNA_metrics.csv,sample2/adata/sample2.crispr.h5ad \
  --batch_id "Batch_A" \
  --batch_basename "batch_a" \
  --qc_container <qc image> \
  --outdir results
```

**samplesheet format** — CSV with columns `sample_id, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv, crispr_h5ad`, plus any other columns you want (they're carried through as per-sample metadata):

```csv
sample_id,molecule_info_h5,filtered_barcodes_tsv,scrna_metrics_csv,features_tsv,crispr_h5ad,condition
sample1,/path/to/sample1/sample1.scRNA.moleculeInfo.h5,/path/to/sample1/sample1.scRNA.filtered.barcodes.tsv.gz,/path/to/sample1/sample1.scRNA_metrics.csv,,/path/to/sample1/adata/sample1.crispr.h5ad,treated
sample2,/path/to/sample2/sample2.scRNA.moleculeInfo.h5,/path/to/sample2/sample2.scRNA.filtered.barcodes.tsv.gz,/path/to/sample2/sample2.scRNA_metrics.csv,,/path/to/sample2/adata/sample2.crispr.h5ad,control
```

- `molecule_info_h5`, `filtered_barcodes_tsv`, `scrna_metrics_csv` are each required per row when `--run_gex_downsample` is `true` (the default).
- `features_tsv` is optional — leave it blank unless the row's `molecule_info_h5` is a combined GEX+CRISPR DRAGEN h5 (see `bin/downsample_molecule_info.py`).
- `crispr_h5ad` is required per row when `--run_crispr_downsample` is `true` (the default).
- Every other column (e.g. `condition` above) ends up attached to `.obs` on the final combined AnnData.

**Required:**
- `--samplesheet`: CSV described above
- `--staged_files`: Every `molecule_info_h5`/`filtered_barcodes_tsv`/`scrna_metrics_csv`/`features_tsv`/`crispr_h5ad` file referenced by `--samplesheet`'s CSV, all in one list. Not read by the pipeline itself — it's what makes ICA localize those files onto the compute node (ICA has no way to know the CSV references them otherwise); the samplesheet's CSV is what actually assigns each file its role, matched by filename once staged, so a file can be listed here in any order.
- `--batch_id`: Batch identifier
- `--batch_basename`: Batch basename for output organization
- `--qc_container`: Container image for QC/downsample processing

**Optional:**
- `--run_gex_downsample` / `--run_crispr_downsample` / `--run_combine`: toggle each stage (defaults: all `true`; `--run_combine` requires both downsample stages to be enabled)
- `--gex_target_depth` / `--crispr_target_depth`: override the auto-computed common target
- `--run_saturation` / `--saturation_extra_depths`: see [GEX sequencing-saturation sweep](#gex-sequencing-saturation-sweep---run_saturation-both-entrypoints) above (default: `false`; requires `--run_gex_downsample`)
- `--gex_matrix_format`: which GEX matrix layout(s) to write per depth — `raw` (all barcodes), `filtered` (cells in the DRAGEN filtered-barcodes list), or `both` (default: `both`). `--run_combine` reads the filtered layout, so it can't be used with `raw`
- `--min_reads_per_cell`: floor below which a requested depth is dropped (default: `1`)
- `--random_seed`: seed for the downsampling draws (default: `42`)
- Resource params for each process — see `nextflow_schema.json` for defaults
- `--outdir`: Output directory (default: `out`)
- `--help`: Show help message

**Output**, under `${params.outdir}/${params.batch_basename}/`:
- **`gex_downsample_summary/`**: batch-wide GEX metrics summary + resolved common target depth (only if `--run_gex_downsample true`)
- **`grna_downsample_summary/`**: batch-wide gRNA metrics summary + resolved common target depth (only if `--run_crispr_downsample true`)
- **`<sample_id>/gex_downsample/`**: that sample's downsampled GEX matrix — `downsampled_matrix/<depth>rpc/` holding `raw_matrix/` and/or `filtered_matrix/` per `--gex_matrix_format` — plus `saturation.csv` if `--run_saturation true` (only if `--run_gex_downsample true`)
- **`<sample_id>/crispr_downsample/`**: that sample's downsampled gRNA AnnData (only if `--run_crispr_downsample true`)
- **`combined/<batch_basename>.combined.h5ad`**: the final combined, metadata-annotated AnnData (only if `--run_combine true`)
- **`pipeline_info/`**: Nextflow reports (`timeline.html`, `report.html`, `trace.txt`, `dag.svg`)

A `README.txt` describing this layout is written directly into `${params.outdir}/${params.batch_basename}/` when the run finishes successfully.

### Running downsample_single_sample.nf

```bash
nextflow run downsample_single_sample.nf \
  --sample_id "Sample_A" \
  --input_files sample_a/sample_a.scRNA.moleculeInfo.h5,sample_a/sample_a.scRNA.filtered.barcodes.tsv.gz,sample_a/sample_a.scRNA_metrics.csv,sample_a/adata/sample_a.crispr.h5ad \
  --gex_target_depths 5000 10000 20000 \
  --crispr_target_depths 500 1000 2000 \
  --qc_container <qc image> \
  --outdir results
```

**Required:**
- `--sample_id`: Sample identifier
- `--input_files`: All of the sample's files, in any order. Each file's role is auto-detected from its DRAGEN/CRISPR filename suffix — no per-role flags to fill in separately:
  - `*scRNA.moleculeInfo.h5` → molecule_info_h5, `*scRNA.filtered.barcodes.tsv.gz` → filtered_barcodes_tsv, `*scRNA_metrics.csv` → scrna_metrics_csv (all three required when `--run_gex_downsample` is `true`)
  - `*scRNA.features.tsv.gz` → features_tsv (optional — only needed when the molecule_info_h5 is a combined GEX+CRISPR DRAGEN h5)
  - `*.h5ad` → crispr_h5ad (required when `--run_crispr_downsample` is `true`)
- `--gex_target_depths`: one or more GEX target reads/cell depths (required when `--run_gex_downsample` is `true`)
- `--crispr_target_depths`: one or more gRNA target reads/cell depths (required when `--run_crispr_downsample` is `true`)
- `--qc_container`: Container image for QC/downsample processing

**Optional:** the same `--run_gex_downsample` / `--run_crispr_downsample` / `--run_combine` / `--run_saturation` / `--saturation_extra_depths` / `--gex_matrix_format` / `--min_reads_per_cell` / `--random_seed` / `--outdir` / `--help` params as `downsample.nf` — see `nextflow_schema_single_sample.json` for defaults.

**Output**, under `${params.outdir}/${params.sample_id}/`:
- **`gex_downsample/`**: one `<depth>rpc/` subdir per requested GEX depth, each holding `raw_matrix/` and/or `filtered_matrix/` per `--gex_matrix_format`, plus `saturation.csv` if `--run_saturation true` (only if `--run_gex_downsample true`)
- **`crispr_downsample/`**: one `<depth>rpc.h5ad` file per requested gRNA depth (only if `--run_crispr_downsample true`)
- **`combined/<sample_id>.depths_combined.h5ad`**: the combined AnnData across every depth either modality produced, annotated with `depth`, `modality` and `sample_id` (only if `--run_combine true`)
- **`pipeline_info/`**: Nextflow reports, nested under `batch_basename`'s config-level default rather than `sample_id` — a cosmetic-only quirk (see the comment in `nextflow.config`); every actual data output above already publishes under `sample_id` correctly

A `README.txt` describing this layout is written directly into `${params.outdir}/${params.sample_id}/` when the run finishes successfully.

## Parameter validation

Each entrypoint has its own schema (`downsample.nf` → `nextflow_schema.json`, `downsample_single_sample.nf` → `nextflow_schema_single_sample.json`). Required/typed params (presence, type, allowed range/pattern) are validated against the relevant schema via the [`nf-schema`](https://nextflow-io.github.io/nf-schema/) plugin as soon as the pipeline starts — a missing, mistyped, or out-of-range param fails immediately with a clear message rather than partway through the run.

Both schemas declare the same full `resource_options`/`batch_basename` set (marking whichever half a given entrypoint doesn't use as `hidden`), even though each entrypoint's Nextflow processes only read the subset relevant to it — this keeps `nf-schema` from warning about the config-level defaults in the shared `nextflow.config` that the *other* entrypoint's schema doesn't declare.

A few checks that can't be expressed in JSON Schema are enforced separately, right after schema validation:
- `downsample.nf`: every samplesheet row needs non-empty `molecule_info_h5`/`filtered_barcodes_tsv`/`scrna_metrics_csv` when `--run_gex_downsample` is true (and likewise `crispr_h5ad` for `--run_crispr_downsample`); `--run_combine` requires both downsample stages enabled; `--run_saturation` requires `--run_gex_downsample`.
- `downsample_single_sample.nf`: `--input_files` must include a molecule_info_h5/filtered_barcodes_tsv/scrna_metrics_csv (by filename suffix) and `--gex_target_depths` must be given when `--run_gex_downsample` is true (and likewise a crispr_h5ad in `--input_files`/`--crispr_target_depths` for `--run_crispr_downsample`); same `--run_combine`/`--run_saturation` cross-checks as above.

**Don't add hand-rolled `if (!params.x) { exit 1 }` checks for anything expressible in JSON Schema** — add/edit the corresponding property (and its `required` list) in the relevant schema file instead, so ICA's rendered form and the pipeline's own validation never drift apart.

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run downsample.nf ... -resume
```

## Troubleshooting

**Error**: `ModuleNotFoundError` / missing Python package inside a process
- **Solution**: The package needs to be added to `docker/qc/Dockerfile` (or `docker/qc/requirements.txt`) and the `qc_container` image rebuilt/pushed — pipeline processes run inside the container, not your local environment.

**Error**: Pipeline exits immediately with a parameter validation error
- **Solution**: Required/typed params are validated against the entrypoint's schema via the `nf-schema` plugin (see [Parameter validation](#parameter-validation)); a few additional business-logic checks run right after. Check the exact list of required flags above.

**Out of memory / timeout**: Adjust resource allocations in `nextflow.config` for specific processes, e.g.:
  ```groovy
  process {
      withName: DOWNSAMPLE_MOLECULE_INFO {
          memory = 64.GB
          time = 24.h
      }
  }
  ```

## ICA export

All ICA-facing tooling lives in `ica_tools/`: `export_pipeline_to_ica.py`, `start_analysis.py`, `ica_common.py` (shared helpers), and `ica_tools/inputforms/<downsample|single_sample>/inputForm.json` (hand-maintained ICA launch-form definitions, one per entrypoint).

`ica_tools/export_pipeline_to_ica.py` imports the current git commit into an ICA project as a git-backed Nextflow pipeline (reads an API key from `~/.icav2/api_key.txt`, prompts interactively for which ICA project **and which entrypoint** — `downsample.nf` or `downsample_single_sample.nf`, since they're registered as separate ICA pipelines sharing the same `nextflow.config`). Exported pipeline codes are `PIPseq_BCL_Downsample_<hash>` (batch) / `PIPseq_BCL_Downsample_SingleSample_<hash>` (single-sample). Keep the corresponding `ica_tools/inputforms/<entrypoint>/inputForm.json` in sync whenever a param is added/changed/removed in that entrypoint's schema — it's a separately hand-maintained file (ICA doesn't derive its launch form from the schema), so nothing enforces this at runtime the way `validateParameters()` does for the schema itself.

`ica_tools/start_analysis.py` (no CLI args other than `--dry-run`) interactively starts an actual analysis run for a pipeline already imported into ICA: pick the ICA project, pick which already-imported downsample pipeline to run (listed via `GET /projects/{projectId}/pipelines`, filtered to pipelines whose `code` starts with `PIPseq_BCL_Downsample`, newest first, flagging whichever one's `gitPipelineImportDto.commitId` matches the current git HEAD), then pick `test/test_inputs_downsample.json` or `test/test_inputs_single_sample.json` to submit. `--dry-run` resolves inputs and prints the request payload without submitting.

## Additional Resources

- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
- [Nextflow Container Documentation](https://www.nextflow.io/docs/latest/container.html)
