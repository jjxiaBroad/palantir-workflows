/*
 * Process module wrapping bin/downsample_crispr_anndata.py: downsamples one sample's
 * cell x gRNA AnnData to one or more requested gRNA target depths. Shared by both
 * downsample.nf (batch) and downsample_single_sample.nf via workflows/downsample_core.nf.
 */

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint is what actually sets these; see pipseq_core.nf (PIPseqPipeline) for the same note.
params.min_reads_per_cell = 1
params.random_seed = 42

process DOWNSAMPLE_CRISPR_ANNDATA {
    tag "${meta.sample_id}"
    // run_basename is optional -- null/blank (single-sample callers, which have no separate
    // outer grouping folder) collapses to just "${outdir}/${sample_id}/...", non-blank (batch
    // callers) nests under "${outdir}/${run_basename}/${sample_id}/...".
    publishDir "${params.outdir}/${meta.run_basename ? meta.run_basename + '/' : ''}${meta.sample_id}/crispr_downsample", mode: 'copy'
    container "${params.qc_container}"
    cpus params.cpu_downsample_crispr_anndata
    memory "${params.memory_gb_downsample_crispr_anndata}.GB"

    input:
    // meta: [sample_id, run_basename, target_depths] -- run_basename is whatever top-level
    // output folder the calling entrypoint organizes this run under (batch_basename for
    // downsample.nf, sample_id for downsample_single_sample.nf), target_depths a list of one
    // or more reads-per-cell depths.
    tuple val(meta), path(crispr_h5ad)

    output:
    // Aggregate combine steps collect one of these per sample into a list and stage it with
    // stageAs to avoid a name collision (every sample's output dir has the same basename
    // "downsampled_crispr"). One <depth>rpc.h5ad file is written per requested depth in
    // meta.target_depths.
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
        --matrix-depths ${meta.target_depths.join(' ')} \\
        --min-reads-per-cell ${params.min_reads_per_cell} \\
        --random-seed ${params.random_seed}
    """

    stub:
    def depth_touches = meta.target_depths.collect { "touch downsampled_crispr/${it}rpc.h5ad" }.join('\n    ')
    """
    echo "[STUB] Would downsample gRNA AnnData for ${meta.sample_id} to ${meta.target_depths.join(',')} reads/cell"

    mkdir -p downsampled_crispr
    ${depth_touches}
    """
}
