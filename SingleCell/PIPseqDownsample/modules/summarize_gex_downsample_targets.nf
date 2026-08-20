/*
 * Process module for summarizing DRAGEN scRNA metrics across a batch and resolving
 * the common GEX downsampling target (see bin/summarize_gex_downsample_targets.py).
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint (downsample.nf) is what actually sets these; see pipseq_core.nf for the same note.
params.gex_target_depth = null

process SUMMARIZE_GEX_DOWNSAMPLE_TARGETS {
    tag "Summarizing GEX downsample targets for ${sample_ids.size()} samples"
    publishDir "${params.outdir}/${params.batch_basename}/gex_downsample_summary", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_summarize_gex_downsample_targets
    memory "${params.memory_gb_summarize_gex_downsample_targets}.GB"

    input:
    // stageAs avoids a staging collision if multiple samples' scRNA_metrics.csv happen to
    // share a basename.
    tuple path(scrna_metrics_csvs, stageAs: 'scrna_metrics_csv_*'), val(sample_ids), path(samplesheet)

    output:
    path "gex_downsample_summary.csv", emit: summary
    path "gex_target_depth.txt", emit: target_depth

    script:
    def csvs = scrna_metrics_csvs.join(' ')
    def ids = sample_ids.join(',')
    def target_override = params.gex_target_depth ? "--target-reads-per-cell ${params.gex_target_depth}" : ''
    """
    set -ex

    summarize_gex_downsample_targets.py \\
        --scrna-metrics-csvs ${csvs} \\
        --sample-ids ${ids} \\
        --samplesheet ${samplesheet} \\
        --output gex_downsample_summary.csv \\
        --target-depth-output gex_target_depth.txt \\
        ${target_override}
    """

    stub:
    """
    echo "[STUB] Would summarize GEX downsample targets for: ${sample_ids.join(',')}"

    echo "sample_id,mean_reads_per_cell,target_mean_reads_per_cell" > gex_downsample_summary.csv
    echo "10000" > gex_target_depth.txt
    """
}
