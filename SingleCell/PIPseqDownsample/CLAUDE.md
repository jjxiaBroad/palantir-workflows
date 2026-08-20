# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Nextflow DSL2 pipelines for downsampling already-processed PIPseq single-cell RNA-seq samples (GEX + gRNA/CRISPR Guide Capture) for cross-sample or cross-depth comparison. Sibling to, but fully independent of, `../PIPseqPipeline` (own entrypoints, schemas, `bin/`, `modules/`, Docker container, `ica_tools/`). See `README.md` for full CLI/params/output details — trust it over this file for day-to-day usage; this file is for things not obvious from the code/README alone.

Two entrypoints: `downsample.nf` (batch, `--samplesheet`) and `downsample_single_sample.nf` (one sample, multiple target depths). Deployment target is Illumina Connected Analytics (ICA); `ica_tools/` publishes/runs these on ICA.

## Commands

**Always set `NXF_VER=25.10.0`.** The machine's default Nextflow (26.04.6) fails to *compile* these scripts: its strict parser rejects the top-level `workflow.onComplete { ... }` blocks in both entrypoints ("Statements cannot be mixed with script declarations"). 25.10.0 is also the floor required by the pinned `nf-schema@2.7.3` plugin.

```bash
# Stub run — the only self-contained way to validate pipeline wiring (no container needed).
# Every process has a stub: block that just touches placeholder outputs.
cd stub_test
NXF_VER=25.10.0 nextflow run ../downsample.nf -stub-run -params-file stub_inputs/pipeline_input.json
NXF_VER=25.10.0 nextflow run ../downsample_single_sample.nf -stub-run -params-file stub_inputs/pipeline_input_single_sample.json
# Clean up afterwards: rm -rf work out .nextflow .nextflow.log*

# Any bin/ script standalone (they're plain argparse CLIs, no Nextflow needed) — needs the qc env
python bin/downsample_molecule_info.py --help

# Build/push the qc_container image
export ECR_REGISTRY=<account>.dkr.ecr.us-east-1.amazonaws.com
docker/build_and_push.sh [tag]        # repo name: pipseq-downsample-qc

# ICA (both interactive; read ~/.icav2/api_key.txt)
python ica_tools/export_pipeline_to_ica.py   # import current git commit as an ICA pipeline
python ica_tools/start_analysis.py [--dry-run]
```

There is **no unit-test suite and no CI** for this directory — the repo-root GitHub Actions workflows don't cover `SingleCell/`. `-stub-run` plus a real ICA run are the whole verification story, so run a stub run after any channel/tuple change: an input-tuple arity mismatch shows up there as `WARN: Input tuple does not match tuple declaration` followed by a hard failure.

There's no local/container-less profile — every process declares `container "${params.qc_container}"`, so a non-stub run needs a container-capable executor.

## Architecture

Both entrypoints are thin: they resolve *which* depths to use, then delegate all per-sample work to `workflows/downsample_core.nf` (`DOWNSAMPLE_SAMPLES`), which applies a given depth list uniformly to every sample in its input channels. The difference is entirely in how the depth list is produced:

- `downsample.nf` — one common depth per modality for the whole batch, computed by `SUMMARIZE_GEX_DOWNSAMPLE_TARGETS` / `SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS` (default: the minimum mean-reads-per-cell across the batch — you can only downsample down) or overridden by `--gex_target_depth`/`--crispr_target_depth`. Each summarize process writes a `*_target_depth.txt` that the workflow reads back with `.map { it.text.trim() as Integer }`, so the depth arrives as a channel value resolved *asynchronously* by an upstream process.
- `downsample_single_sample.nf` — the depth list comes straight from `--gex_target_depths`/`--crispr_target_depths` params.

Conventions inside `downsample_core.nf` worth knowing before editing it:

- **Depth lists are wrapped in a single-key map** (`[depths: [...]]`) rather than passed as bare `List`s. Nextflow's `combine()`/tuple machinery silently flattens a bare `List` value into separate tuple positions; a map value isn't flattened.
- **`run_basename` controls output nesting.** It's carried in each process's `meta` map and interpolated into `publishDir`: non-blank (batch, `batch_basename`) nests `${outdir}/${run_basename}/${sample_id}/...`, null/blank (single-sample) collapses to `${outdir}/${sample_id}/...`.
- **Aggregate combine steps use `stageAs`.** Every sample's `DOWNSAMPLE_MOLECULE_INFO` output dir has the same basename (`downsampled_matrix`), so `COMBINE_DOWNSAMPLED_BATCH` (and the summarize processes) stage their input lists as `gex_matrix_dir_*` etc. Don't "fix" this by embedding `sample_id` into the produced dir name.
- **GEX inputs are 4 explicit per-file paths, not a directory.** `DOWNSAMPLE_MOLECULE_INFO`'s input tuple is `(meta, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv)` — `features_tsv` is optional and passed as a `file('NO_FILE')` sentinel when absent (checked in the module via `features_tsv.name != 'NO_FILE'`), matching the same convention `../PIPseqPipeline/workflows/pipseq_core.nf` uses for its own optional guide-assignments file. This was a deliberate move away from a `--dragen_results_dir` directory param (see "Why GEX inputs are per-file" below).
- `params.x = ...` defaults are re-declared inside modules/`downsample_core.nf` purely to silence Nextflow's "access to undefined parameter" warning; the entrypoint is what actually sets them.
- In each entrypoint, **all `def` declarations must precede any top-level statement** (including `params.x = ...`) — hence `helpMessage()`/`writeOutputManifest()` sit at the top of the file.

Python side (`bin/`): every script is a standalone CLI usable outside Nextflow.

- `downsample_molecule_info.py` (~1700 lines) is the heavy one — an IMI-only molecule-info downsampler ported from the parent pipeline's `counting/minf.py`. DRAGEN 4.6+ HDF5 molecule-info only (needs `/binning_idx`); no UMI counting, no pre-4.6 CSV parsers, no CRISPR path. Features/barcodes are read from the h5 itself; `--features-tsv` is only needed for *combined* GEX+CRISPR molecule-info h5s. It also owns the optional saturation sweep (`--saturation`), which walks a built-in reads-per-cell ladder without materializing full matrices at every rung. It still accepts `--dragen-results-dir` and will glob-search it for the other files (`find_molecule_info_h5`/`find_filtered_barcodes`/`find_scrna_metrics_csv`/`find_features_tsv`/`detect_dragen_prefix`) — that path only matters for ad hoc standalone CLI use now; the pipeline itself always passes `--molecule-info-h5`/`--filtered-barcodes`/`--scrna-metrics-csv`/`--features-tsv` explicitly.
- `downsample_crispr_anndata.py` uses plain binomial thinning of raw counts — gRNA reads have no binning index and no molecule dedup, so there's no unique-molecule correction and **no meaningful "sequencing saturation"** for that modality (that's why `--run_saturation` is GEX-only).
- The two combine scripts share `downsample_combine_common.py` for barcode normalization (strips the trailing `-1`-style 10x suffix) and the feature-axis merge. A GEX/gRNA barcode mismatch is a *warning* (keep the intersection), not a failure.
- Both summarize scripts and `combine_downsampled_batch.py` each carry their own `SAMPLESHEET_PATH_COLUMNS` set — the samplesheet columns treated as file locators rather than per-sample metadata (`sample_id, molecule_info_h5, scrna_metrics_csv, filtered_barcodes_tsv, features_tsv, crispr_h5ad`). Everything *not* in that set is carried through into the summary CSVs and the combined AnnData's `.obs`. **These three sets must stay in sync** whenever the samplesheet's columns change.

### Adding or changing a param

Params are validated by `nf-schema`'s `validateParameters()` against a per-entrypoint schema (`downsample.nf` → `nextflow_schema.json`, `downsample_single_sample.nf` → `nextflow_schema_single_sample.json`). **Don't add hand-rolled `if (!params.x) { exit 1 }` checks for anything expressible in JSON Schema** — edit the schema instead, so ICA's rendered form and the pipeline's validation can't drift. Only genuinely non-schema-expressible checks (samplesheet *contents*, and cross-param rules like "`--run_combine` requires both downsample stages") live as inline `log.error` + `exit 1` right after `validateParameters()`.

Both schemas deliberately declare the *union* of both entrypoints' `resource_options`/`batch_basename` params (marking the unused half `hidden`) so `nf-schema` doesn't warn about the shared `nextflow.config` defaults the other entrypoint's schema wouldn't know about.

A param change touches, in practice: the schema(s), the entrypoint's `params.x` default + `helpMessage()`, `ica_tools/inputforms/<downsample|single_sample>/inputForm.json` (hand-maintained — ICA does *not* derive it from the schema, and nothing checks this at runtime), `test/test_inputs_*.json`, `stub_test/stub_inputs/*`, and `README.md`.

`nextflow.config`'s `timeline`/`report`/`trace`/`dag` paths are resolved at config-parse time, before `validateParameters()` applies schema defaults, so they can only key off `batch_basename` (the one param with a config-level default). That's why a single-sample run writes `pipeline_info/` under `batch_output/` instead of its `sample_id` — cosmetic only; every data output publishes under `sample_id` correctly via its module's own `publishDir`.

## Why GEX inputs are per-file, not a directory

`--dragen_results_dirs`/`--dragen_results_dir` used to take a *directory*, and `bin/downsample_molecule_info.py`'s `_files_with_suffix()` glob-searched inside it (non-recursively) for `*scRNA_metrics.csv`, `*scRNA.moleculeInfo.h5`, filtered barcodes, and `*scRNA.features.tsv.gz`. This was a problem on ICA specifically: every sample's DRAGEN output folder is conventionally named the same thing (e.g. `dragen_output`), and ICA's Nextflow samplesheet-file staging model stages selected data into one shared working directory **by basename only** (see the ICA reference memory / `reference_ica_samplesheet_path_convention` note) — so multiple same-named folders across samples risked colliding when staged together, unlike the individual DRAGEN files themselves, which have a unique per-sample prefix.

This has been refactored (2026-08-19): the directory param is now split into 4 explicit per-file params (`--molecule_info_h5(s)`, `--filtered_barcodes_tsv(s)`, `--scrna_metrics_csv(s)`, `--features_tsv(s)`, the last optional) on `downsample_single_sample.nf`, mirroring how `--crispr_h5ads`/`--crispr_h5ad` already worked. Touched: `workflows/downsample_core.nf`'s `DOWNSAMPLE_SAMPLES` input tuple, both entrypoints' params/samplesheet-parsing/`helpMessage()`, both `nextflow_schema*.json`, the samplesheet format itself (4 new columns replacing `dragen_results_dir`, carried into `stub_test/stub_inputs/` and `test/test_inputs_*.json`), both `ica_tools/inputforms/*/inputForm.json`, and `README.md`. Verified via both entrypoints' `-stub-run`. Not touched (deliberately): `bin/downsample_molecule_info.py`'s directory-globbing helpers (`find_*`, `detect_dragen_prefix`, `--dragen-results-dir`) — still there for standalone/ad hoc CLI use, just no longer exercised by the pipeline.

`downsample.nf` (batch) took this one step further (2026-08-20): its 5 per-role list params (`--molecule_info_h5s`, `--filtered_barcodes_tsvs`, `--scrna_metrics_csvs`, `--features_tsvs`, `--crispr_h5ads`) were collapsed into a single `--staged_files` list param. These were *never read by the pipeline* in the first place — only `--samplesheet`'s CSV assigns each file its role; the list params exist solely so ICA localizes the files (see "ICA notes" below). Since nothing checks which field a file came through, one field is sufficient and lets you multi-select every file for a sample (or the whole batch) in one browse action instead of one per role.

`downsample_single_sample.nf` has no samplesheet, so its 5 per-role params (`--molecule_info_h5`, `--filtered_barcodes_tsv`, `--scrna_metrics_csv`, `--features_tsv`, `--crispr_h5ad`) *are* real values the entrypoint reads directly — they couldn't be merged the same way without losing role information. Instead (also 2026-08-20) they were replaced with a single `--input_files` list param; `classifyInputFiles()` (defined at the top of `downsample_single_sample.nf`, alongside `helpMessage()`/`writeOutputManifest()` per the `def`-ordering rule above) assigns each file's role by matching its filename against the same DRAGEN suffixes `bin/downsample_molecule_info.py`'s standalone-CLI directory-globbing helpers use (`*scRNA.moleculeInfo.h5`, `*scRNA.filtered.barcodes.tsv.gz`, `*scRNA_metrics.csv`, `*scRNA.features.tsv.gz`) plus `*.h5ad` for `crispr_h5ad`, erroring out on an unmatched or ambiguous (multiple files matching one role) file. This reintroduces filename-convention inference deliberately avoided elsewhere in this pipeline (see the samplesheet docstrings' "no filename-convention guessing") — acceptable here specifically because there's exactly one sample, so there's no cross-sample basename-collision risk to guess wrong about.

## ICA notes (learned 2026-08-18, first real test run)

- **Samplesheet file references on ICA must be bare filenames, no path at all** (not even the ICA-project-relative path shown in the data browser) — see `reference_ica_samplesheet_path_convention`. Any file a samplesheet references must *also* be separately selected in a `"type": "data"` input-form field (`downsample.nf`: the single `staged_files` field; `downsample_single_sample.nf`: its own per-role fields) so ICA localizes it into the run's working directory; the samplesheet is just how the pipeline looks it up by name once staged.
- Real ICA run (2026-08-20) hit `downsample_crispr_anndata.py`: "No requested depths remain after filtering." Root cause: the batch gRNA/GEX target depth is the exact minimum sample's mean reads/cell, and `summarize_gex/grna_downsample_targets.py` used `int(round(...))` to write it — rounding *up* past that same sample's true depth, so its own downsample factor computed as `> 1` ("would be upsampling") and got rejected. Fixed by rounding the target down (`math.floor`) in both summarize scripts, and by relaxing `downsample_molecule_info.py`/`downsample_crispr_anndata.py`'s exclusion from `factor >= 1` to `factor > 1` (an exact factor of 1 is a legitimate no-op, not upsampling).
- `ica_tools/export_pipeline_to_ica.py` requires the commit to exist on a repo ICA can clone, and hardcodes the *upstream* URL (see the git-remotes reference memory). Testing so far has bypassed it by pointing ICA's GUI git import at the user's fork (`https://github.com/jjxiaBroad/palantir-workflows`), avoiding the need to merge/push upstream first.
- Exported pipeline codes: `PIPseq_BCL_Downsample_<hash>` (batch) / `PIPseq_BCL_Downsample_SingleSample_<hash>` (single-sample). The two entrypoints are registered as *separate* ICA pipelines sharing one `nextflow.config`.
- A first real ICA run runs slower than a warm local run: `nf-schema` plugin download, qc_container image pull, compute cold-start. Not pipeline-logic overhead.
- ICA's GUI "Rerun" restarts the whole analysis from scratch; there's no confirmed equivalent of Nextflow's own `-resume` task-level caching exposed through the ICA GUI/API.
