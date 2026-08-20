/*
 * Shared per-sample downsampling engine, called by both entrypoints (downsample.nf for a
 * batch normalized to one common depth per modality, downsample_single_sample.nf for one
 * sample swept across several depths). Callers resolve their own target-depth list(s) --
 * this workflow just applies whatever list it's given, uniformly, to every sample in its
 * input channels.
 */

include { DOWNSAMPLE_MOLECULE_INFO } from '../modules/downsample_molecule_info'
include { DOWNSAMPLE_CRISPR_ANNDATA } from '../modules/downsample_crispr_anndata'

// Declared here too so Nextflow doesn't warn about "access to undefined parameter" -- the
// entrypoint is what actually sets these; see pipseq_core.nf (PIPseqPipeline) for the same note.
params.run_gex_downsample = true
params.run_crispr_downsample = true

workflow DOWNSAMPLE_SAMPLES {
    take:
    gex_input            // channel of (sample_id, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv,
                          // features_tsv) -- only consumed if run_gex_downsample. features_tsv is a
                          // file('NO_FILE') sentinel when not given (see downsample_molecule_info.py).
    crispr_input         // channel of (sample_id, crispr_h5ad) -- only consumed if run_crispr_downsample
    run_basename         // String or null: optional top-level output folder nested above each
                          // sample's own subdirectory (batch_basename for downsample.nf; null for
                          // downsample_single_sample.nf, which has no separate outer grouping folder)
    gex_target_depths    // one-element channel carrying [depths: List<Integer>] -- reads-per-cell depths
                          // applied to every sample in gex_input; a channel (not a plain List) so callers
                          // can feed in a value resolved asynchronously by an upstream process (as
                          // downsample.nf does for its batch-wide common depth). Wrapped in a single-key
                          // map rather than passed as a bare List -- Nextflow's combine()/tuple machinery
                          // silently flattens a bare List value into separate tuple positions, which a
                          // map value is not subject to.
    crispr_target_depths // one-element channel carrying [depths: List<Integer>], same idea for crispr_input

    main:
    gex_matrices_ch = Channel.empty()
    saturation_csv_ch = Channel.empty()
    crispr_h5ads_ch = Channel.empty()

    if (params.run_gex_downsample) {
        downsample_gex_input = gex_input
            .combine(gex_target_depths)
            .map { sample_id, molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv, depths_map ->
                tuple(
                    [sample_id: sample_id, run_basename: run_basename, target_depths: depths_map.depths],
                    molecule_info_h5, filtered_barcodes_tsv, scrna_metrics_csv, features_tsv
                )
            }
        DOWNSAMPLE_MOLECULE_INFO(downsample_gex_input)
        gex_matrices_ch = DOWNSAMPLE_MOLECULE_INFO.out.matrices
        saturation_csv_ch = DOWNSAMPLE_MOLECULE_INFO.out.saturation_csv
    }

    if (params.run_crispr_downsample) {
        downsample_crispr_input = crispr_input
            .combine(crispr_target_depths)
            .map { sample_id, crispr_h5ad, depths_map ->
                tuple([sample_id: sample_id, run_basename: run_basename, target_depths: depths_map.depths], crispr_h5ad)
            }
        DOWNSAMPLE_CRISPR_ANNDATA(downsample_crispr_input)
        crispr_h5ads_ch = DOWNSAMPLE_CRISPR_ANNDATA.out.h5ads
    }

    emit:
    matrices = gex_matrices_ch          // (sample_id, downsampled_matrix dir)
    h5ads = crispr_h5ads_ch             // (sample_id, downsampled_crispr dir)
    saturation_csvs = saturation_csv_ch // (sample_id, saturation.csv) -- empty unless --run_saturation
}
