/*
 * Process module wrapping bin/downsample_molecule_info.py: downsamples one sample's
 * DRAGEN scRNA molecule-info h5 to one or more requested GEX target depths, optionally
 * also running a sequencing-saturation sweep. Shared by both downsample.nf (batch) and
 * downsample_single_sample.nf via workflows/downsample_core.nf.
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint is what actually sets these; see pipseq_core.nf (PIPseqPipeline) for the same note.
params.min_reads_per_cell = 1
params.random_seed = 42
params.run_saturation = false
params.saturation_extra_depths = []

process DOWNSAMPLE_MOLECULE_INFO {
    tag "${meta.sample_id}"
    // run_basename is optional -- null/blank (single-sample callers, which have no separate
    // outer grouping folder) collapses to just "${outdir}/${sample_id}/...", non-blank (batch
    // callers) nests under "${outdir}/${run_basename}/${sample_id}/...".
    publishDir "${params.outdir}/${meta.run_basename ? meta.run_basename + '/' : ''}${meta.sample_id}/gex_downsample", mode: 'copy'
    container "${params.qc_container}"
    // cpus reserved here are also handed to the script's own --threads flag below, so
    // the multiprocessing pool it spins up actually matches what's been reserved.
    cpus params.cpu_downsample_molecule_info
    memory "${params.memory_gb_downsample_molecule_info}.GB"

    input:
    // meta: [sample_id, run_basename, target_depths] -- run_basename is whatever top-level
    // output folder the calling entrypoint organizes this run under (batch_basename for
    // downsample.nf, sample_id for downsample_single_sample.nf), target_depths a list of one
    // or more reads-per-cell depths. features_tsv is a 'NO_FILE' sentinel when not given --
    // only combined GEX+CRISPR molecule-info h5s need it (see downsample_molecule_info.py).
    tuple val(meta), path(molecule_info_h5), path(filtered_barcodes_tsv), path(scrna_metrics_csv), path(features_tsv)

    output:
    // Aggregate combine steps collect one of these per sample into a list -- since every
    // sample's output dir has the same basename ("downsampled_matrix"), those steps stage
    // this input with stageAs to avoid a name collision, rather than this process embedding
    // sample_id into the dir name itself. One <depth>rpc/filtered_matrix/ subdir is written
    // per requested depth in meta.target_depths.
    tuple val(meta.sample_id), path("downsampled_matrix"), emit: matrices
    tuple val(meta.sample_id), path("saturation.csv"), emit: saturation_csv, optional: true

    script:
    def saturation_args = params.run_saturation
        ? ' --saturation' + (params.saturation_extra_depths ? " --saturation-extra-depths ${params.saturation_extra_depths.join(' ')}" : '')
        : ''
    def features_args = features_tsv.name != 'NO_FILE' ? " --features-tsv ${features_tsv}" : ''
    """
    set -ex

    downsample_molecule_info.py \\
        --molecule-info-h5 ${molecule_info_h5} \\
        --filtered-barcodes ${filtered_barcodes_tsv} \\
        --scrna-metrics-csv ${scrna_metrics_csv} \\
        --output-dir . \\
        --matrix-depths ${meta.target_depths.join(' ')} \\
        --matrix-format filtered \\
        --min-reads-per-cell ${params.min_reads_per_cell} \\
        --random-seed ${params.random_seed} \\
        --threads ${task.cpus}${saturation_args}${features_args}
    """

    stub:
    def depth_dirs = meta.target_depths.collect { depth ->
        "mkdir -p downsampled_matrix/${depth}rpc/filtered_matrix\n" +
        "    touch downsampled_matrix/${depth}rpc/filtered_matrix/matrix.mtx.gz\n" +
        "    touch downsampled_matrix/${depth}rpc/filtered_matrix/barcodes.tsv.gz\n" +
        "    touch downsampled_matrix/${depth}rpc/filtered_matrix/features.tsv.gz"
    }.join('\n    ')
    def saturation_stub = params.run_saturation ? 'touch saturation.csv' : ''
    """
    echo "[STUB] Would downsample molecule info for ${meta.sample_id} to ${meta.target_depths.join(',')} reads/cell"

    ${depth_dirs}
    ${saturation_stub}
    """
}
