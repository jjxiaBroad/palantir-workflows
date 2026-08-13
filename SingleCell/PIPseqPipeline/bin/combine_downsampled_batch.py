#!/usr/bin/env python3
"""Combine a batch's per-sample downsampled GEX + gRNA (CRISPR Guide Capture)
data into a single AnnData for cross-sample comparison.

For each sample: read its downsampled GEX filtered matrix (10x-mtx format,
written by downsample_molecule_info.py under
<gex-matrix-dir>/<depth>rpc/filtered_matrix/) and its downsampled gRNA h5ad
(written by downsample_crispr_anndata.py under
<crispr-h5ad-dir>/<depth>rpc.h5ad), normalize the trailing 10x barcode suffix
(e.g. "-1") on both sides so barcodes line up, then merge the two modalities
along the feature axis (inner join on barcode -- GEX and gRNA are expected to
already share the same filtered-cell barcodes, but a mismatch is reported
rather than silently ignored).

All samples are then concatenated along the cell axis into one combined
AnnData, and every column from --samplesheet other than
sample_id/dragen_results_dir/crispr_h5ad is attached to the combined
AnnData's .obs (keyed by sample_id), so every cell carries its sample's full
metadata row.

Usage:
    python combine_downsampled_batch.py \\
        --gex-matrix-dirs sample1_gex/ sample2_gex/ \\
        --crispr-h5ad-dirs sample1_crispr/ sample2_crispr/ \\
        --sample-ids sample1,sample2 \\
        --samplesheet samplesheet.csv \\
        --output-h5ad batch.combined.h5ad
"""

import argparse
import glob
import os
import re
import sys

import anndata as ad
import pandas as pd
import scanpy as sc

BARCODE_SUFFIX_RE = re.compile(r"-\d+$")

# Samplesheet columns that locate files rather than describe sample metadata.
SAMPLESHEET_PATH_COLUMNS = {"sample_id", "dragen_results_dir", "crispr_h5ad"}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--gex-matrix-dirs",
        nargs="+",
        required=True,
        help="Each sample's DOWNSAMPLE_MOLECULE_INFO output dir, aligned with --sample-ids.",
    )
    parser.add_argument(
        "--crispr-h5ad-dirs",
        nargs="+",
        required=True,
        help="Each sample's DOWNSAMPLE_CRISPR_ANNDATA output dir, aligned with --sample-ids.",
    )
    parser.add_argument(
        "--sample-ids",
        required=True,
        help="Comma-separated sample IDs, aligned with --gex-matrix-dirs/--crispr-h5ad-dirs.",
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="Batch samplesheet CSV, used to attach per-sample metadata columns to .obs.",
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help="Path to write the combined AnnData to.",
    )
    return parser.parse_args(argv)


def find_one(pattern, what):
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise FileNotFoundError(f"No {what} found for pattern: {pattern}")
    if len(matches) > 1:
        print(f"WARNING: multiple {what} matches for {pattern}; using first: {matches}", file=sys.stderr)
    return matches[0]


def load_gex(gex_matrix_dir):
    matrix_dir = find_one(os.path.join(gex_matrix_dir, "*", "filtered_matrix"), "GEX filtered_matrix dir")
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    adata.var["feature_types"] = "Gene Expression"
    return adata


def load_downsampled_grna(crispr_h5ad_dir):
    h5ad_path = find_one(os.path.join(crispr_h5ad_dir, "*.h5ad"), "downsampled gRNA .h5ad")
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names_make_unique()
    adata.var["feature_types"] = "CRISPR Guide Capture"
    return adata


def normalize_barcodes(adata):
    adata.obs_names = adata.obs_names.str.replace(BARCODE_SUFFIX_RE, "", regex=True)
    return adata


def merge_sample(sample_id, gex_matrix_dir, crispr_h5ad_dir):
    gex = normalize_barcodes(load_gex(gex_matrix_dir))
    grna = normalize_barcodes(load_downsampled_grna(crispr_h5ad_dir))

    shared = gex.obs_names.intersection(grna.obs_names)
    if len(shared) != gex.n_obs or len(shared) != grna.n_obs:
        print(
            f"WARNING: sample '{sample_id}': GEX has {gex.n_obs} filtered cells, gRNA has "
            f"{grna.n_obs} filtered cells, only {len(shared)} barcodes are shared -- "
            "only shared barcodes will be kept for this sample.",
            file=sys.stderr,
        )

    combined = ad.concat([gex, grna], axis=1, join="inner", merge="unique")
    combined.obs["sample_id"] = sample_id
    print(f"{sample_id}: GEX {gex.shape} + downsampled gRNA {grna.shape} -> combined {combined.shape}")
    return combined


def attach_samplesheet_metadata(combined, samplesheet_path):
    df = pd.read_csv(samplesheet_path, dtype=str)
    if "sample_id" not in df.columns:
        raise ValueError(f"{samplesheet_path} has no 'sample_id' column.")
    metadata_cols = [c for c in df.columns if c not in SAMPLESHEET_PATH_COLUMNS and c != "sample_id"]
    metadata_df = df.set_index("sample_id")[metadata_cols]
    combined.obs = combined.obs.join(metadata_df, on="sample_id")
    return combined


def main(argv=None):
    args = parse_args(argv)
    sample_ids = [s.strip() for s in args.sample_ids.split(",") if s.strip()]

    if not (len(sample_ids) == len(args.gex_matrix_dirs) == len(args.crispr_h5ad_dirs)):
        raise ValueError(
            f"--sample-ids ({len(sample_ids)}), --gex-matrix-dirs ({len(args.gex_matrix_dirs)}), "
            f"and --crispr-h5ad-dirs ({len(args.crispr_h5ad_dirs)}) must all be aligned 1:1."
        )

    per_sample = [
        merge_sample(sample_id, gex_dir, crispr_dir)
        for sample_id, gex_dir, crispr_dir in zip(sample_ids, args.gex_matrix_dirs, args.crispr_h5ad_dirs)
    ]

    combined = ad.concat(
        per_sample,
        axis=0,
        join="outer",
        label="sample_id",
        keys=sample_ids,
        index_unique="_",
    )
    # ad.concat's label column takes precedence; each sample already set
    # obs["sample_id"] to its own ID in merge_sample, so this is a no-op
    # sanity-preserving assignment, not a second source of truth.
    combined = attach_samplesheet_metadata(combined, args.samplesheet)

    print(f"\nFinal combined AnnData: {combined}")
    combined.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"\nWrote combined AnnData to {args.output_h5ad}")


if __name__ == "__main__":
    main()
