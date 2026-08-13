/*
 * Process module wrapping bin/combine_downsampled_batch.py: merges every sample's
 * downsampled GEX + gRNA data into one AnnData, annotated with samplesheet metadata.
 */

process COMBINE_DOWNSAMPLED_BATCH {
    tag "Combining ${sample_ids.size()} samples"
    publishDir "${params.outdir}/${params.batch_basename}/combined", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_combine_downsampled_batch
    memory "${params.memory_gb_combine_downsampled_batch}.GB"

    input:
    // Every sample's downsampled-matrix/downsampled-crispr dir shares the same basename,
    // so both lists are staged with stageAs to avoid a collision.
    tuple path(gex_matrix_dirs, stageAs: 'gex_matrix_dir_*'),
          path(crispr_h5ad_dirs, stageAs: 'crispr_h5ad_dir_*'),
          val(sample_ids),
          path(samplesheet)

    output:
    path "${params.batch_basename}.combined.h5ad", emit: combined_adata

    script:
    def gex_dirs = gex_matrix_dirs.join(' ')
    def crispr_dirs = crispr_h5ad_dirs.join(' ')
    def ids = sample_ids.join(',')
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}

    combine_downsampled_batch.py \\
        --gex-matrix-dirs ${gex_dirs} \\
        --crispr-h5ad-dirs ${crispr_dirs} \\
        --sample-ids ${ids} \\
        --samplesheet ${samplesheet} \\
        --output-h5ad ${params.batch_basename}.combined.h5ad
    """

    stub:
    """
    echo "[STUB] Would combine ${sample_ids.size()} samples: ${sample_ids.join(',')}"

    touch ${params.batch_basename}.combined.h5ad
    """
}
