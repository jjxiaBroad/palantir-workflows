/*
 * Process module for summarizing gRNA AnnData metrics across a batch and resolving
 * the common gRNA downsampling target (see bin/summarize_grna_downsample_targets.py).
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint (downsample.nf) is what actually sets these; see pipseq_core.nf for the same note.
params.crispr_target_depth = null

process SUMMARIZE_GRNA_DOWNSAMPLE_TARGETS {
    tag "Summarizing gRNA downsample targets for ${sample_ids.size()} samples"
    publishDir "${params.outdir}/${params.batch_basename}/grna_downsample_summary", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_summarize_grna_downsample_targets
    memory "${params.memory_gb_summarize_grna_downsample_targets}.GB"

    input:
    // stageAs avoids a staging collision when samples' crispr_h5ad files share a basename.
    tuple path(crispr_h5ads, stageAs: 'crispr_h5ad_*'), val(sample_ids), path(samplesheet)

    output:
    path "grna_downsample_summary.csv", emit: summary
    path "grna_target_depth.txt", emit: target_depth

    script:
    def h5ads = crispr_h5ads.join(' ')
    def ids = sample_ids.join(',')
    def target_override = params.crispr_target_depth ? "--target-reads-per-cell ${params.crispr_target_depth}" : ''
    """
    set -ex
    export NUMBA_CACHE_DIR=${workflow.launchDir}
    export MPLCONFIGDIR=${workflow.launchDir}

    summarize_grna_downsample_targets.py \\
        --crispr-h5ads ${h5ads} \\
        --sample-ids ${ids} \\
        --samplesheet ${samplesheet} \\
        --output grna_downsample_summary.csv \\
        --target-depth-output grna_target_depth.txt \\
        ${target_override}
    """

    stub:
    """
    echo "[STUB] Would summarize gRNA downsample targets for: ${sample_ids.join(',')}"

    echo "sample_id,grna_mean_reads_per_cell,grna_target_mean_reads_per_cell" > grna_downsample_summary.csv
    echo "500" > grna_target_depth.txt
    """
}
