/*
 * Process module wrapping bin/downsample_crispr_anndata.py: downsamples one sample's
 * cell x gRNA AnnData to the batch's resolved gRNA target depth.
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint (downsample.nf) is what actually sets these; see pipseq_core.nf for the same note.
params.min_reads_per_cell = 1
params.random_seed = 42

process DOWNSAMPLE_CRISPR_ANNDATA {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/${params.batch_basename}/${meta.sample_id}/crispr_downsample", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_downsample_crispr_anndata
    memory "${params.memory_gb_downsample_crispr_anndata}.GB"

    input:
    // meta: [sample_id, target_depth]
    tuple val(meta), path(crispr_h5ad)

    output:
    // Aggregate steps (COMBINE_DOWNSAMPLED_BATCH) collect one of these per sample into a
    // list and stage it with stageAs to avoid a name collision (every sample's output dir
    // has the same basename "downsampled_crispr").
    tuple val(meta.sample_id), path("downsampled_crispr"), emit: h5ads

    script:
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}

    mkdir -p downsampled_crispr
    downsample_crispr_anndata.py \\
        --input-h5ad ${crispr_h5ad} \\
        --output-dir downsampled_crispr \\
        --matrix-depths ${meta.target_depth} \\
        --min-reads-per-cell ${params.min_reads_per_cell} \\
        --random-seed ${params.random_seed}
    """

    stub:
    """
    echo "[STUB] Would downsample gRNA AnnData for ${meta.sample_id} to ${meta.target_depth} reads/cell"

    mkdir -p downsampled_crispr
    touch downsampled_crispr/${meta.target_depth}rpc.h5ad
    """
}
