# PIPseq Batch Downsample-and-Combine Pipeline

A Nextflow DSL2 pipeline for normalizing a *batch* of already-processed PIPseq single-cell RNA-seq samples to a common sequencing depth per cell, then combining them into one AnnData for cross-sample comparison.

It operates on already-published DRAGEN scRNA output + gRNA AnnData (from `../PIPseqPipeline`, or any other pipeline producing the same shapes) — it never runs DRAGEN itself, so it's cheap to rerun repeatedly at different depths. There is a single entrypoint, `downsample.nf`, validated against `nextflow_schema.json`.

The pipeline is designed to run on Illumina Connected Analytics (ICA) but is plain Nextflow DSL2 and can run anywhere a compatible executor and the required container image are available.

## Overview

Given a `--samplesheet` (columns: `sample_id, dragen_results_dir, crispr_h5ad`, plus any per-sample metadata columns) describing a batch of samples, `downsample.nf`:

1. Resolves one common GEX target depth and one common gRNA target depth for the whole batch — each independently defaults to the minimum "mean reads per cell" observed across the batch (you can only downsample down, never up), or an explicit `--gex_target_depth`/`--crispr_target_depth` override. Published as `gex_downsample_summary.csv`/`grna_downsample_summary.csv` (useful on their own, as a batch depth-landscape report).
2. Downsamples every sample's GEX molecule-info (filtered-matrix layout, via `bin/downsample_molecule_info.py`'s IMI-based downsampler) and/or gRNA AnnData (via `bin/downsample_crispr_anndata.py`'s binomial-thinning downsampler) to those targets.
3. Combines every sample's downsampled GEX + gRNA data into one AnnData (`bin/combine_downsampled_batch.py`) — only if both branches ran (`--run_combine`, default `true`). Barcodes between a sample's filtered GEX and gRNA are expected to already match; a mismatch is reported as a warning (only the shared barcodes are kept) rather than failing the run.

## Pipeline Structure

```
SingleCell/PIPseqDownsample/
├── downsample.nf                    # The single entrypoint (--samplesheet)
├── nextflow.config                  # Pipeline configuration (default params, resources, reports)
├── nextflow_schema.json             # Parameter schema (drives the ICA-rendered input form)
├── modules/
│   ├── summarize_gex_downsample_targets.nf   # Resolve the batch's common GEX target depth
│   ├── summarize_grna_downsample_targets.nf  # Resolve the batch's common gRNA target depth
│   ├── downsample_molecule_info.nf   # Downsample one sample's GEX molecule-info
│   ├── downsample_crispr_anndata.nf  # Downsample one sample's gRNA AnnData
│   └── combine_downsampled_batch.nf  # Combine a batch's downsampled GEX+gRNA into one AnnData
├── bin/
│   ├── downsample_molecule_info.py   # GEX molecule-info downsampler (IMI-based)
│   ├── downsample_crispr_anndata.py  # gRNA AnnData downsampler (binomial thinning)
│   ├── summarize_gex_downsample_targets.py   # GEX batch-target summary/resolution script
│   ├── summarize_grna_downsample_targets.py  # gRNA batch-target summary/resolution script
│   └── combine_downsampled_batch.py  # GEX+gRNA batch-combine script
├── docker/                          # Dockerfile/build scripts for the qc_container image
├── stub_test/                       # Example inputs + `-stub-run` test setup
├── test/                            # Flat pipeline-inputs JSON used by ica_tools/start_analysis.py
├── ica_tools/                       # Scripts for publishing/running this pipeline on ICA
│   ├── export_pipeline_to_ica.py     # Imports the current commit into ICA as a git-backed pipeline
│   ├── start_analysis.py             # Interactively starts an ICA analysis run
│   ├── ica_common.py                 # Shared helpers (API key, project list, prompt_choice())
│   └── inputforms/inputForm.json     # Hand-maintained ICA launch-form definition
└── README.md                        # This file
```

## Quick Start

### Prerequisites

- Nextflow (>= 25.10.0 — required by the pinned `nf-schema` validation plugin; requires network access to the Nextflow plugin registry on first run)
- A QC container image (`--qc_container`) — built from `docker/qc/Dockerfile`, see `docker/SETUP.md`
- An executor/environment that can run the `container` directive (e.g. Nextflow's k8s executor on ICA, or Docker/Singularity enabled locally via your own config)

There is no bundled Nextflow profile for local/container-less execution — every process declares a `container`, so running for real requires a container-capable executor.

### Stub run (no container required)

Every process has a `stub:` block that just touches placeholder output files, so you can validate the pipeline's wiring without running any real downsampling:

```bash
cd stub_test
nextflow run ../downsample.nf -stub-run -params-file stub_inputs/pipeline_input.json
```

### Running the Pipeline

```bash
nextflow run downsample.nf \
  --samplesheet samplesheet.csv \
  --dragen_results_dirs sample1/dragen_output,sample2/dragen_output \
  --crispr_h5ads sample1/adata/sample1.crispr.h5ad,sample2/adata/sample2.crispr.h5ad \
  --batch_id "Batch_A" \
  --batch_basename "batch_a" \
  --qc_container <qc image> \
  --outdir results
```

### samplesheet format

CSV with columns `sample_id, dragen_results_dir, crispr_h5ad`, plus any other columns you want (they're carried through as per-sample metadata):

```csv
sample_id,dragen_results_dir,crispr_h5ad,condition
sample1,/path/to/sample1/dragen_output,/path/to/sample1/adata/sample1.crispr.h5ad,treated
sample2,/path/to/sample2/dragen_output,/path/to/sample2/adata/sample2.crispr.h5ad,control
```

- `dragen_results_dir` is required per row when `--run_gex_downsample` is `true` (the default).
- `crispr_h5ad` is required per row when `--run_crispr_downsample` is `true` (the default).
- Every other column (e.g. `condition` above) ends up attached to `.obs` on the final combined AnnData.

### Command-Line Options

**Required:**
- `--samplesheet`: CSV described above
- `--dragen_results_dirs`: All DRAGEN results directories referenced by `--samplesheet`'s CSV. Not read by the pipeline itself — it's what makes ICA localize those directories onto the compute node (ICA has no way to know the CSV references them otherwise).
- `--crispr_h5ads`: All `crispr_h5ad` files referenced by `--samplesheet`'s CSV — same reasoning as `--dragen_results_dirs`.
- `--batch_id`: Batch identifier
- `--batch_basename`: Batch basename for output organization
- `--qc_container`: Container image for QC/downsample processing

**Optional:**
- `--run_gex_downsample` / `--run_crispr_downsample` / `--run_combine`: toggle each stage (defaults: all `true`; `--run_combine` requires both downsample stages to be enabled)
- `--gex_target_depth` / `--crispr_target_depth`: override the auto-computed common target
- `--min_reads_per_cell`: floor below which a requested depth is dropped (default: `1`)
- `--random_seed`: seed for the downsampling draws (default: `42`)
- Resource params for each of the 5 processes (`summarize_gex_downsample_targets`, `summarize_grna_downsample_targets`, `downsample_molecule_info`, `downsample_crispr_anndata`, `combine_downsampled_batch`) — see `nextflow_schema.json` for defaults
- `--outdir`: Output directory (default: `out`)
- `--help`: Show help message

## Parameter validation

Required/typed params (presence, type, allowed range/pattern) are validated against `nextflow_schema.json` via the [`nf-schema`](https://nextflow-io.github.io/nf-schema/) plugin as soon as the pipeline starts — a missing, mistyped, or out-of-range param fails immediately with a clear message rather than partway through the run.

A few checks that can't be expressed in JSON Schema are enforced separately, right after schema validation: every samplesheet row needs a non-empty `dragen_results_dir` when `--run_gex_downsample` is true (and likewise `crispr_h5ad` for `--run_crispr_downsample`), and `--run_combine` requires both downsample stages to be enabled.

**Don't add hand-rolled `if (!params.x) { exit 1 }` checks for anything expressible in JSON Schema** — add/edit the corresponding property (and its `required` list) in `nextflow_schema.json` instead, so ICA's rendered form and the pipeline's own validation never drift apart.

## Output

Results are organized under `${params.outdir}/${params.batch_basename}/`:

- **`gex_downsample_summary/`**: batch-wide GEX metrics summary + resolved common target depth (only if `--run_gex_downsample true`)
- **`grna_downsample_summary/`**: batch-wide gRNA metrics summary + resolved common target depth (only if `--run_crispr_downsample true`)
- **`<sample_id>/gex_downsample/`**: that sample's downsampled GEX filtered matrix (only if `--run_gex_downsample true`)
- **`<sample_id>/crispr_downsample/`**: that sample's downsampled gRNA AnnData (only if `--run_crispr_downsample true`)
- **`combined/<batch_basename>.combined.h5ad`**: the final combined, metadata-annotated AnnData (only if `--run_combine true`)
- **`pipeline_info/`**: Nextflow reports (`timeline.html`, `report.html`, `trace.txt`, `dag.svg`)

A `README.txt` describing this layout is written directly into `${params.outdir}/${params.batch_basename}/` when the run finishes successfully.

## Resuming Failed Runs

Nextflow caches completed tasks. Resume a failed pipeline with `-resume`:
```bash
nextflow run downsample.nf ... -resume
```

## Troubleshooting

**Error**: `ModuleNotFoundError` / missing Python package inside a process
- **Solution**: The package needs to be added to `docker/qc/Dockerfile` (or `docker/qc/requirements.txt`) and the `qc_container` image rebuilt/pushed — pipeline processes run inside the container, not your local environment.

**Error**: Pipeline exits immediately with a parameter validation error
- **Solution**: Required/typed params are validated against `nextflow_schema.json` via the `nf-schema` plugin (see [Parameter validation](#parameter-validation)); a few additional business-logic checks (samplesheet content) run right after. Check the exact list of required flags above.

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

All ICA-facing tooling lives in `ica_tools/`: `export_pipeline_to_ica.py`, `start_analysis.py`, `ica_common.py` (shared helpers), and `ica_tools/inputforms/inputForm.json` (a hand-maintained ICA launch-form definition).

`ica_tools/export_pipeline_to_ica.py` imports the current git commit of this pipeline into an ICA project as a git-backed Nextflow pipeline (reads an API key from `~/.icav2/api_key.txt`, prompts interactively for which ICA project). Keep `ica_tools/inputforms/inputForm.json` in sync whenever a param is added/changed/removed in `nextflow_schema.json` — it's a separately hand-maintained file (ICA doesn't derive its launch form from the schema), so nothing enforces this at runtime the way `validateParameters()` does for the schema itself.

`ica_tools/start_analysis.py` (no CLI args other than `--dry-run`) interactively starts an actual analysis run for a pipeline already imported into ICA: pick the ICA project, pick which already-imported downsample pipeline to run (listed via `GET /projects/{projectId}/pipelines`, filtered to pipelines whose `code` starts with `PIPseq_BCL_Downsample`, newest first, flagging whichever one's `gitPipelineImportDto.commitId` matches the current git HEAD), then submits `test/test_inputs.json`. `--dry-run` resolves inputs and prints the request payload without submitting.

## Additional Resources

- [Nextflow Configuration Reference](https://www.nextflow.io/docs/latest/config.html)
- [Nextflow Container Documentation](https://www.nextflow.io/docs/latest/container.html)
