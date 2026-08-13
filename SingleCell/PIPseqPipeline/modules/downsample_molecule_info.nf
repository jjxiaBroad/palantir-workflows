/*
 * Process module wrapping bin/downsample_molecule_info.py: downsamples one sample's
 * DRAGEN scRNA molecule-info h5 to the batch's resolved GEX target depth.
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint (downsample.nf) is what actually sets these; see pipseq_core.nf for the same note.
params.min_reads_per_cell = 1
params.random_seed = 42

process DOWNSAMPLE_MOLECULE_INFO {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${params.batch_basename}/${meta.sample_id}/gex_downsample", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_downsample_molecule_info
    memory "${params.memory_gb_downsample_molecule_info}.GB"

    input:
    // meta: [sample_id, target_depth]
    tuple val(meta), path(dragen_results_dir)

    output:
    // Aggregate steps (COMBINE_DOWNSAMPLED_BATCH) collect one of these per sample into a
    // list -- since every sample's output dir has the same basename ("downsampled_matrix"),
    // that step stages this input with stageAs to avoid a name collision, rather than this
    // process embedding sample_id into the dir name itself.
    tuple val(meta.sample_id), path("downsampled_matrix"), emit: matrices

    script:
    """
    set -ex

    downsample_molecule_info.py \\
        --dragen-results-dir ${dragen_results_dir} \\
        --output-dir . \\
        --matrix-depths ${meta.target_depth} \\
        --matrix-format filtered \\
        --min-reads-per-cell ${params.min_reads_per_cell} \\
        --random-seed ${params.random_seed}
    """

    stub:
    """
    echo "[STUB] Would downsample molecule info for ${meta.sample_id} to ${meta.target_depth} reads/cell"

    mkdir -p downsampled_matrix/${meta.target_depth}rpc/filtered_matrix
    touch downsampled_matrix/${meta.target_depth}rpc/filtered_matrix/matrix.mtx.gz
    touch downsampled_matrix/${meta.target_depth}rpc/filtered_matrix/barcodes.tsv.gz
    touch downsampled_matrix/${meta.target_depth}rpc/filtered_matrix/features.tsv.gz
    """
}
