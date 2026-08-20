#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Batch downsample-and-combine workflow.
 *
 * Separate entrypoint from main.nf/main_simple.nf/workflows/pipseq_core.nf -- this operates on
 * already-published DRAGEN scRNA output + gRNA AnnData for a batch of independently processed
 * samples, rather than running DRAGEN itself. Given a --samplesheet describing the batch, it:
 *   1. Resolves a common GEX downsampling target (default: the minimum "Mean reads per cell"
 *      observed across the batch) and a common gRNA downsampling target (same idea, independently
 *      -- the two modalities are on unrelated numeric scales), and publishes both as summary CSVs.
 *   2. Downsamples every sample's GEX molecule-info and gRNA AnnData to those common targets.
 *   3. Combines every sample's downsampled GEX + gRNA data into one AnnData, annotated with the
 *      samplesheet's per-sample metadata columns, for cross-sample analysis.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { SUMMARIZE_GEX_DOWNSAMPLE_TARGETS } from './modules/summarize_gex_downsample_targets'
include { SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS } from './modules/summarize_grna_downsample_targets'
include { DOWNSAMPLE_SAMPLES } from './workflows/downsample_core'
include { COMBINE_DOWNSAMPLED_BATCH } from './modules/combine_downsampled_batch'

// All `def` declarations must come before any top-level statement (including `params.x = ...`
// assignments) in the same script -- the strict script parser rejects declarations that follow
// a statement. Hence helpMessage()/writeOutputManifest() are defined here, ahead of the param
// defaults and workflow blocks below.

def helpMessage() {
    log.info"""
    Usage:
      nextflow run downsample.nf --samplesheet <samplesheet.csv> --batch_id <id> --batch_basename <name> \\
          --qc_container <image> [options]

    Required arguments:
      --samplesheet              CSV describing the batch (columns: sample_id, molecule_info_h5,
                                  filtered_barcodes_tsv, scrna_metrics_csv, features_tsv, crispr_h5ad, +metadata)
      --molecule_info_h5s        All molecule_info_h5 files referenced by --samplesheet (needed so ICA localizes them onto the compute node; not read directly by the pipeline)
      --filtered_barcodes_tsvs   All filtered_barcodes_tsv files referenced by --samplesheet (same reasoning as --molecule_info_h5s)
      --scrna_metrics_csvs       All scrna_metrics_csv files referenced by --samplesheet (same reasoning as --molecule_info_h5s)
      --crispr_h5ads             All crispr_h5ad files referenced by --samplesheet (same reasoning as --molecule_info_h5s)
      --batch_id                 Batch identifier
      --batch_basename           Batch basename for output organization
      --qc_container             Container image for QC/downsample processing

    Optional data arguments:
      --features_tsvs            All features_tsv files referenced by --samplesheet, for any row that
                                  sets one -- only needed when a sample's molecule_info_h5 is a combined
                                  GEX+CRISPR DRAGEN h5 (see downsample_molecule_info.py); same
                                  localization reasoning as --molecule_info_h5s

    Samplesheet format:
      CSV file with columns: sample_id, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv, crispr_h5ad
      - molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv are required for every row when
        --run_gex_downsample is true
      - features_tsv is optional (blank unless the row's molecule_info_h5 is a combined GEX+CRISPR h5)
      - crispr_h5ad is required for every row when --run_crispr_downsample is true
      - any other column is treated as per-sample metadata and carried through into the
        downsample summary CSVs and the final combined AnnData's .obs

    Optional arguments:
      --run_gex_downsample       Whether to downsample GEX molecule-info (default: true)
      --run_crispr_downsample    Whether to downsample gRNA AnnData (default: true)
      --run_combine              Whether to combine downsampled GEX+gRNA into one AnnData -- requires both
                                  --run_gex_downsample and --run_crispr_downsample (default: true)
      --gex_target_depth         Override the common GEX target reads/cell (default: minimum observed across the batch)
      --crispr_target_depth      Override the common gRNA target reads/cell (default: minimum observed across the batch)
      --run_saturation           Also run a GEX sequencing-saturation sweep per sample (default: false; requires --run_gex_downsample)
      --saturation_extra_depths  Extra reads-per-cell depths to include in the saturation sweep, beyond the script's built-in ladder
      --min_reads_per_cell       Floor below which a requested depth is dropped (default: 1)
      --random_seed              Seed for the downsampling draws (default: 42)
      --outdir                   Output directory (default: out)
      --help                     Show this help message

    Behavior:
      - Resolves one common GEX target depth and one common gRNA target depth for the whole batch,
        publishing a summary CSV for each under gex_downsample_summary/ and grna_downsample_summary/
      - Downsamples every sample's GEX molecule-info (filtered-matrix layout only) and/or gRNA AnnData
        to those targets, under <sample_id>/gex_downsample/ and <sample_id>/crispr_downsample/
      - If --run_saturation is set, also writes a saturation.csv per sample alongside its downsampled
        GEX matrix (median transcripts/genes and % sequencing saturation across a depth ladder)
      - Combines every sample's downsampled GEX+gRNA into one AnnData under combined/, annotated with
        the samplesheet's per-sample metadata columns -- only if both branches ran
    """.stripIndent()
}

def writeOutputManifest() {
    def manifest = file("${params.outdir}/${params.batch_basename}/README.txt")
    manifest.text = """
        Output layout for batch '${params.batch_id}' (${params.batch_basename}):

          gex_downsample_summary/         Batch-wide GEX metrics summary + resolved common target
                                           depth (only if --run_gex_downsample true)
          grna_downsample_summary/        Batch-wide gRNA metrics summary + resolved common target
                                           depth (only if --run_crispr_downsample true)
          <sample_id>/gex_downsample/     That sample's downsampled GEX filtered matrix, plus
                                           saturation.csv if --run_saturation true (only if
                                           --run_gex_downsample true)
          <sample_id>/crispr_downsample/  That sample's downsampled gRNA AnnData (only if
                                           --run_crispr_downsample true)
          combined/                       Combined AnnData for the whole batch (<basename>.combined.h5ad),
                                           annotated with the samplesheet's per-sample metadata columns
                                           (only if --run_combine true)
          pipeline_info/                  Nextflow execution reports (timeline, report, trace, DAG)

        See README.md in the pipeline repository for parameter and output details.
        """.stripIndent()
}

// Define parameters
params.samplesheet = null              // CSV: sample_id, molecule_info_h5, filtered_barcodes_tsv,
                                        // scrna_metrics_csv, features_tsv, crispr_h5ad, +metadata columns
params.molecule_info_h5s = null        // Unused by the pipeline directly -- exists so ICA localizes the
                                        // files referenced by path inside --samplesheet's CSV.
params.filtered_barcodes_tsvs = null   // Same idea as --molecule_info_h5s, for the filtered_barcodes_tsv column.
params.scrna_metrics_csvs = null       // Same idea as --molecule_info_h5s, for the scrna_metrics_csv column.
params.features_tsvs = []              // Same idea as --molecule_info_h5s, for the (optional) features_tsv column.
params.crispr_h5ads = null             // Same idea as --molecule_info_h5s, for the crispr_h5ad column.
// batch_id/batch_basename are intentionally not declared here (same convention as
// main.nf's supersample_id/supersample_basename) -- batch_basename's default lives in
// nextflow.config (needed early, for the timeline/report/trace/dag file paths), and
// batch_id has no default since it's required by nextflow_schema.json.
params.outdir = "out"                  // Output directory
params.help = false
params.qc_container = null             // QC container image
params.run_gex_downsample = true       // Whether to downsample GEX molecule-info
params.run_crispr_downsample = true    // Whether to downsample gRNA AnnData
params.run_combine = true              // Whether to combine downsampled GEX+gRNA into one AnnData
params.gex_target_depth = null         // Optional override for the common GEX target (reads/cell)
params.crispr_target_depth = null      // Optional override for the common gRNA target (reads/cell)
params.run_saturation = false          // Also run a GEX sequencing-saturation sweep per sample
params.saturation_extra_depths = []    // Extra reads-per-cell depths for the saturation sweep
params.min_reads_per_cell = 1          // Floor below which a requested depth is dropped
params.random_seed = 42                // Seed for the downsampling draws

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required/typed params against nextflow_schema.json. Don't add
    // hand-rolled `if (!params.x) exit 1` checks here for anything the schema already declares.
    validateParameters(parameters_schema: 'nextflow_schema.json')
    log.info paramsSummaryLog(workflow)

    log.info "Reading samplesheet..."

    // --- Business-logic checks that can't be expressed in JSON Schema (samplesheet content) ---
    sample_info = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            [
                sample_id: row.sample_id,
                molecule_info_h5: row.molecule_info_h5,
                filtered_barcodes_tsv: row.filtered_barcodes_tsv,
                scrna_metrics_csv: row.scrna_metrics_csv,
                features_tsv: row.features_tsv,
                crispr_h5ad: row.crispr_h5ad
            ]
        }
        .toList()
        .flatMap { rows ->
            if (params.run_gex_downsample) {
                def missing = rows.findAll {
                    !it.molecule_info_h5?.trim() || !it.filtered_barcodes_tsv?.trim() || !it.scrna_metrics_csv?.trim()
                }
                if (missing) {
                    log.error "ERROR: --run_gex_downsample is true, but the following sample(s) in " +
                        "--samplesheet are missing molecule_info_h5/filtered_barcodes_tsv/scrna_metrics_csv: " +
                        missing.collect { it.sample_id }.join(', ')
                    exit 1
                }
            }
            if (params.run_crispr_downsample) {
                def missing = rows.findAll { !it.crispr_h5ad?.trim() }
                if (missing) {
                    log.error "ERROR: --run_crispr_downsample is true, but the following sample(s) in " +
                        "--samplesheet have no crispr_h5ad: " + missing.collect { it.sample_id }.join(', ')
                    exit 1
                }
            }
            if (params.run_combine && !(params.run_gex_downsample && params.run_crispr_downsample)) {
                log.error "ERROR: --run_combine requires both --run_gex_downsample and --run_crispr_downsample to be true."
                exit 1
            }
            if (params.run_saturation && !params.run_gex_downsample) {
                log.error "ERROR: --run_saturation requires --run_gex_downsample to be true."
                exit 1
            }
            rows
        }

    gex_target_depths_ch = Channel.value([depths: []])
    gex_input_ch = Channel.empty()
    if (params.run_gex_downsample) {
        log.info "Summarizing GEX downsample targets across the batch..."

        gex_summary_input = sample_info
            .map { s -> tuple(file(s.scrna_metrics_csv), s.sample_id) }
            .toList()
            .map { pairs -> tuple(pairs.collect { it[0] }, pairs.collect { it[1] }, file(params.samplesheet)) }

        SUMMARIZE_GEX_DOWNSAMPLE_TARGETS(gex_summary_input)

        gex_target_depths_ch = SUMMARIZE_GEX_DOWNSAMPLE_TARGETS.out.target_depth
            .map { [depths: [it.text.trim() as Integer]] }

        gex_input_ch = sample_info.map { s ->
            tuple(
                s.sample_id,
                file(s.molecule_info_h5),
                file(s.filtered_barcodes_tsv),
                file(s.scrna_metrics_csv),
                s.features_tsv?.trim() ? file(s.features_tsv) : file('NO_FILE')
            )
        }
    }

    crispr_target_depths_ch = Channel.value([depths: []])
    crispr_input_ch = Channel.empty()
    if (params.run_crispr_downsample) {
        log.info "Summarizing gRNA downsample targets across the batch..."

        grna_summary_input = sample_info
            .map { s -> tuple(file(s.crispr_h5ad), s.sample_id) }
            .toList()
            .map { pairs -> tuple(pairs.collect { it[0] }, pairs.collect { it[1] }, file(params.samplesheet)) }

        SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS(grna_summary_input)

        crispr_target_depths_ch = SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS.out.target_depth
            .map { [depths: [it.text.trim() as Integer]] }

        crispr_input_ch = sample_info.map { s -> tuple(s.sample_id, file(s.crispr_h5ad)) }
    }

    log.info "Downsampling GEX molecule-info and/or gRNA AnnData for each sample..."

    DOWNSAMPLE_SAMPLES(
        gex_input_ch,
        crispr_input_ch,
        params.batch_basename,
        gex_target_depths_ch,
        crispr_target_depths_ch
    )

    if (params.run_combine) {
        log.info "Combining downsampled GEX+gRNA for the batch..."

        combine_input = DOWNSAMPLE_SAMPLES.out.matrices
            .join(DOWNSAMPLE_SAMPLES.out.h5ads)
            .toList()
            .map { rows ->
                tuple(
                    rows.collect { it[1] },
                    rows.collect { it[2] },
                    rows.collect { it[0] },
                    file(params.samplesheet)
                )
            }

        COMBINE_DOWNSAMPLED_BATCH(combine_input)
    }
}

workflow.onComplete {
    if (workflow.success) writeOutputManifest()
}
