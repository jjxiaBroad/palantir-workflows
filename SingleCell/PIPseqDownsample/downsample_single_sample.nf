#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Single-sample, multiple-target-depths downsampling workflow.
 *
 * Separate entrypoint from downsample.nf -- rather than normalizing a batch of samples to
 * one common depth, this takes exactly one already-processed sample and downsamples it to
 * several caller-specified target depths (e.g. to build a within-sample depth-response /
 * saturation-style comparison), optionally also running a GEX sequencing-saturation sweep.
 * It:
 *   1. Downsamples the sample's GEX molecule-info and/or gRNA AnnData to every requested depth.
 *   2. Combines the downsampled GEX + gRNA data across depths into one AnnData annotated with
 *      integer `depth` and `modality` columns, for within-sample cross-depth analysis. Every
 *      depth either modality produced is included -- the two modalities' depth lists are
 *      independent and need not line up.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { DOWNSAMPLE_SAMPLES } from './workflows/downsample_core'
include { COMBINE_DOWNSAMPLED_DEPTHS } from './modules/combine_downsampled_depths'

// All `def` declarations must come before any top-level statement (including `params.x = ...`
// assignments) in the same script -- the strict script parser rejects declarations that follow
// a statement. Hence helpMessage()/writeOutputManifest() are defined here, ahead of the param
// defaults and workflow blocks below.

def helpMessage() {
    log.info"""
    Usage:
      nextflow run downsample_single_sample.nf --sample_id <id> \\
          --input_files <molecule_info_h5>,<filtered_barcodes_tsv>,<scrna_metrics_csv>,<crispr_h5ad> \\
          --gex_target_depths <ints> --crispr_target_depths <ints> --qc_container <image> [options]

    Required arguments:
      --sample_id                Sample identifier
      --input_files              All of the sample's DRAGEN scRNA output + gRNA AnnData files, in
                                  any order: molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv
                                  (required if --run_gex_downsample), the optional features_tsv (only
                                  for a combined GEX+CRISPR molecule_info_h5), and crispr_h5ad
                                  (required if --run_crispr_downsample). Each file's role is
                                  auto-detected from its DRAGEN filename suffix:
                                    *scRNA.moleculeInfo.h5        -> molecule_info_h5
                                    *scRNA.filtered.barcodes.tsv.gz -> filtered_barcodes_tsv
                                    *scRNA_metrics.csv             -> scrna_metrics_csv
                                    *scRNA.features.tsv.gz         -> features_tsv
                                    *.h5ad                         -> crispr_h5ad
      --gex_target_depths        One or more GEX target reads/cell depths (required if --run_gex_downsample)
      --crispr_target_depths     One or more gRNA target reads/cell depths (required if --run_crispr_downsample)
      --qc_container             Container image for QC/downsample processing

    Optional arguments:
      --run_gex_downsample       Whether to downsample GEX molecule-info (default: true)
      --run_crispr_downsample    Whether to downsample gRNA AnnData (default: true)
      --run_combine               Whether to combine downsampled GEX+gRNA across depths into one
                                  AnnData -- requires both --run_gex_downsample and
                                  --run_crispr_downsample (default: true)
      --run_saturation            Also run a GEX sequencing-saturation sweep (default: false;
                                  requires --run_gex_downsample)
      --saturation_extra_depths   Extra reads-per-cell depths for the saturation sweep, beyond
                                  the script's built-in ladder
      --gex_matrix_format         Which GEX matrix layout(s) to write per depth: raw (all
                                  barcodes), filtered (cells only), or both (default: both).
                                  --run_combine reads the filtered layout, so it can't be
                                  used with raw
      --min_reads_per_cell       Floor below which a requested depth is dropped (default: 1)
      --random_seed              Seed for the downsampling draws (default: 42)
      --outdir                   Output directory (default: out)
      --help                     Show this help message

    Behavior:
      - Downsamples the sample's GEX molecule-info (raw and/or filtered matrix layout, see
        --gex_matrix_format) and/or gRNA
        AnnData to every depth in --gex_target_depths / --crispr_target_depths, under
        gex_downsample/ and crispr_downsample/
      - If --run_saturation is set, also writes saturation.csv alongside the downsampled GEX
        matrix (median transcripts/genes and % sequencing saturation across a depth ladder)
      - Combines the downsampled GEX+gRNA across depths into one AnnData under combined/,
        annotated with integer depth and modality columns -- only if both branches ran.
        Every depth either modality produced is kept; --gex_target_depths and
        --crispr_target_depths need not match, and a depth present in only one modality
        gets the other's features zero-filled (.obs.modality says which)
    """.stripIndent()
}

def writeOutputManifest() {
    def manifest = file("${params.outdir}/${params.sample_id}/README.txt")
    manifest.text = """
        Output layout for sample '${params.sample_id}':

          gex_downsample/     Downsampled GEX matrix, one <depth>rpc/ subdir per requested
                               depth holding raw_matrix/ and/or filtered_matrix/ per
                               --gex_matrix_format, plus saturation.csv if --run_saturation
                               true (only if --run_gex_downsample true)
          crispr_downsample/  Downsampled gRNA AnnData, one <depth>rpc.h5ad file per requested
                               depth (only if --run_crispr_downsample true)
          combined/            Combined AnnData across depths (<sample_id>.depths_combined.h5ad),
                               annotated with integer depth and modality columns; a depth present
                               in only one modality has the other's features zero-filled, so
                               filter on .obs['modality'] before reading those zeros as measured
                               counts (only if --run_combine true)
          pipeline_info/       Nextflow execution reports (timeline, report, trace, DAG)

        See README.md in the pipeline repository for parameter and output details.
        """.stripIndent()
}

// Classify --input_files by DRAGEN filename suffix -- this is what lets ICA users select every
// file for the sample together in one field instead of one per role (there's no samplesheet here
// to assign roles explicitly, unlike downsample.nf's batch mode). Suffixes mirror
// bin/downsample_molecule_info.py's DRAGEN_*_SUFFIX constants. Order matters only in that
// features.tsv.gz must be checked before a hypothetical looser molecule-info match; none of these
// suffixes actually overlap.
def classifyInputFiles(input_files) {
    def result = [molecule_info_h5: null, filtered_barcodes_tsv: null, scrna_metrics_csv: null,
                   features_tsv: null, crispr_h5ad: null]
    def matchedBy = [:]
    def unmatched = []

    input_files.each { path ->
        def name = file(path).name
        def role =
            name.endsWith('scRNA.features.tsv.gz') ? 'features_tsv' :
            name.endsWith('scRNA.moleculeInfo.h5') ? 'molecule_info_h5' :
            name.endsWith('scRNA.filtered.barcodes.tsv.gz') ? 'filtered_barcodes_tsv' :
            name.endsWith('scRNA_metrics.csv') ? 'scrna_metrics_csv' :
            name.endsWith('.h5ad') ? 'crispr_h5ad' :
            null

        if (role == null) {
            unmatched << path
            return
        }
        if (result[role] != null) {
            log.error "ERROR: --input_files has more than one file matching '${role}': " +
                "'${matchedBy[role]}' and '${path}'."
            exit 1
        }
        result[role] = file(path)
        matchedBy[role] = path
    }

    if (unmatched) {
        log.error "ERROR: --input_files has file(s) that don't match any known DRAGEN/CRISPR " +
            "filename suffix (expected one of *scRNA.moleculeInfo.h5, " +
            "*scRNA.filtered.barcodes.tsv.gz, *scRNA_metrics.csv, *scRNA.features.tsv.gz, " +
            "*.h5ad): " + unmatched.join(', ')
        exit 1
    }
    return result
}

// Define parameters
params.sample_id = null                // Sample identifier
params.input_files = []                // All of the sample's DRAGEN scRNA + gRNA AnnData files;
                                        // role is auto-detected per file, see classifyInputFiles()
params.outdir = "out"                  // Output directory
params.help = false
params.qc_container = null             // QC container image
params.run_gex_downsample = true       // Whether to downsample GEX molecule-info
params.run_crispr_downsample = true    // Whether to downsample gRNA AnnData
params.run_combine = true              // Whether to combine downsampled GEX+gRNA across depths
params.gex_target_depths = []          // GEX target reads/cell depths
params.crispr_target_depths = []       // gRNA target reads/cell depths
params.run_saturation = false          // Also run a GEX sequencing-saturation sweep
params.saturation_extra_depths = []    // Extra reads-per-cell depths for the saturation sweep
params.gex_matrix_format = 'both'      // Which GEX matrix layout(s) to write: raw, filtered, or both
params.min_reads_per_cell = 1          // Floor below which a requested depth is dropped
params.random_seed = 42                // Seed for the downsampling draws

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required/typed params against nextflow_schema_single_sample.json. Don't add
    // hand-rolled `if (!params.x) exit 1` checks here for anything the schema already declares.
    validateParameters(parameters_schema: 'nextflow_schema_single_sample.json')
    log.info paramsSummaryLog(workflow)

    // --- Business-logic checks that can't be expressed in JSON Schema ---
    def input = classifyInputFiles(params.input_files)

    if (params.run_gex_downsample && !(input.molecule_info_h5 && input.filtered_barcodes_tsv && input.scrna_metrics_csv)) {
        log.error "ERROR: --run_gex_downsample is true, but --input_files did not include a " +
            "molecule_info_h5/filtered_barcodes_tsv/scrna_metrics_csv (all three required)."
        exit 1
    }
    if (params.run_gex_downsample && !params.gex_target_depths) {
        log.error "ERROR: --run_gex_downsample is true, but --gex_target_depths was not given."
        exit 1
    }
    if (params.run_crispr_downsample && !input.crispr_h5ad) {
        log.error "ERROR: --run_crispr_downsample is true, but --input_files did not include a crispr_h5ad."
        exit 1
    }
    if (params.run_crispr_downsample && !params.crispr_target_depths) {
        log.error "ERROR: --run_crispr_downsample is true, but --crispr_target_depths was not given."
        exit 1
    }
    if (params.run_combine && !(params.run_gex_downsample && params.run_crispr_downsample)) {
        log.error "ERROR: --run_combine requires both --run_gex_downsample and --run_crispr_downsample to be true."
        exit 1
    }
    if (params.run_saturation && !params.run_gex_downsample) {
        log.error "ERROR: --run_saturation requires --run_gex_downsample to be true."
        exit 1
    }
    if (params.run_combine && params.gex_matrix_format == 'raw') {
        log.error "ERROR: --run_combine reads each depth's filtered_matrix/ layout, so it is " +
            "incompatible with --gex_matrix_format raw. Use 'filtered' or 'both'."
        exit 1
    }

    gex_input_ch = Channel.empty()
    gex_target_depths_ch = Channel.value([depths: []])
    if (params.run_gex_downsample) {
        gex_input_ch = Channel.of(tuple(
            params.sample_id,
            input.molecule_info_h5,
            input.filtered_barcodes_tsv,
            input.scrna_metrics_csv,
            input.features_tsv ?: file('NO_FILE')
        ))
        gex_target_depths_ch = Channel.value([depths: params.gex_target_depths])
    }

    crispr_input_ch = Channel.empty()
    crispr_target_depths_ch = Channel.value([depths: []])
    if (params.run_crispr_downsample) {
        crispr_input_ch = Channel.of(tuple(params.sample_id, input.crispr_h5ad))
        crispr_target_depths_ch = Channel.value([depths: params.crispr_target_depths])
    }

    log.info "Downsampling GEX molecule-info and/or gRNA AnnData for ${params.sample_id}..."

    DOWNSAMPLE_SAMPLES(
        gex_input_ch,
        crispr_input_ch,
        null,  // run_basename: no separate outer grouping folder for a single-sample run
        gex_target_depths_ch,
        crispr_target_depths_ch
    )

    if (params.run_combine) {
        log.info "Combining downsampled GEX+gRNA across depths for ${params.sample_id}..."

        combine_input = DOWNSAMPLE_SAMPLES.out.matrices
            .join(DOWNSAMPLE_SAMPLES.out.h5ads)

        COMBINE_DOWNSAMPLED_DEPTHS(combine_input)
    }
}

workflow.onComplete {
    if (workflow.success) writeOutputManifest()
}
