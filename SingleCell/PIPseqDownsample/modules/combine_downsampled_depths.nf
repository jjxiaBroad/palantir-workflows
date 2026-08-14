/*
 * Process module wrapping bin/combine_downsampled_depths.py: merges one sample's
 * downsampled GEX + gRNA data across multiple depths into one AnnData annotated with
 * an integer depth column, for within-sample depth-response analysis.
 */

process COMBINE_DOWNSAMPLED_DEPTHS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/combined", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_combine_downsampled_depths
    memory "${params.memory_gb_combine_downsampled_depths}.GB"

    input:
    // gex_matrix_dir / crispr_h5ad_dir: this sample's DOWNSAMPLE_MOLECULE_INFO /
    // DOWNSAMPLE_CRISPR_ANNDATA output dirs, each holding one subdir/file per requested depth.
    tuple val(sample_id), path(gex_matrix_dir, stageAs: 'gex_matrix_dir'), path(crispr_h5ad_dir, stageAs: 'crispr_h5ad_dir')

    output:
    path "${sample_id}.depths_combined.h5ad", emit: combined_adata

    script:
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}

    combine_downsampled_depths.py \\
        --gex-matrix-dir ${gex_matrix_dir} \\
        --crispr-h5ad-dir ${crispr_h5ad_dir} \\
        --sample-id ${sample_id} \\
        --output-h5ad ${sample_id}.depths_combined.h5ad
    """

    stub:
    """
    echo "[STUB] Would combine downsampled depths for sample ${sample_id}"

    touch ${sample_id}.depths_combined.h5ad
    """
}
