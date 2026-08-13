#!/usr/bin/env python3
"""Summarize gRNA (CRISPR Guide Capture) AnnData metrics across a batch of
samples and compute the common gRNA downsampling target needed to normalize
mean reads-per-cell across the batch before cross-sample comparison.

Standalone -- does not touch GEX/DRAGEN scRNA_metrics.csv at all (see
summarize_gex_downsample_targets.py for that side).

Each sample's gRNA AnnData is located via its samplesheet's crispr_h5ad
column (an explicit, required path -- no filename-convention guessing).
``.X`` is assumed to hold raw, undeduplicated read counts (no IMI/IPM --
gRNA/CRISPR Guide Capture data has no meaningful binning index).

The common target -- grna_target_mean_reads_per_cell -- defaults to the
minimum mean reads/cell observed across the batch (or an explicit
--target-reads-per-cell override). This value (rounded to the nearest
integer) is exactly what downsample_crispr_anndata.py --matrix-depths should
be given for every sample in the batch.

Every samplesheet column other than sample_id/dragen_results_dir/crispr_h5ad
is carried through into the output summary verbatim.

Usage:
    python summarize_grna_downsample_targets.py \\
        --crispr-h5ads sample1.crispr.h5ad sample2.crispr.h5ad \\
        --sample-ids sample1,sample2 \\
        --samplesheet samplesheet.csv \\
        --output grna_downsample_summary.csv \\
        --target-depth-output grna_target_depth.txt
"""

import argparse
import sys

import numpy as np
import pandas as pd
import anndata as ad

# Samplesheet columns that locate files rather than describe sample metadata.
SAMPLESHEET_PATH_COLUMNS = {"sample_id", "dragen_results_dir", "crispr_h5ad"}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--crispr-h5ads",
        nargs="+",
        required=True,
        help="Cell x gRNA .h5ad files, one per sample, aligned with --sample-ids.",
    )
    parser.add_argument(
        "--sample-ids",
        required=True,
        help="Comma-separated sample IDs, aligned with --crispr-h5ads.",
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="Batch samplesheet CSV, used only to pull per-sample metadata columns.",
    )
    parser.add_argument(
        "--target-reads-per-cell",
        type=float,
        default=None,
        help="Target mean reads per cell for gRNA. Defaults to the minimum gRNA "
             "mean reads/cell observed across the batch.",
    )
    parser.add_argument(
        "--output",
        default="grna_downsample_summary.csv",
        help="Path to write the summary CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--target-depth-output",
        default="grna_target_depth.txt",
        help="Path to write the resolved integer target depth to (default: %(default)s)",
    )
    return parser.parse_args(argv)


def load_samplesheet_metadata(samplesheet_path):
    """Return {sample_id: {metadata column: value}} for every non-path column."""
    df = pd.read_csv(samplesheet_path, dtype=str)
    if "sample_id" not in df.columns:
        raise ValueError(f"{samplesheet_path} has no 'sample_id' column.")
    metadata_cols = [c for c in df.columns if c not in SAMPLESHEET_PATH_COLUMNS]
    return {
        row["sample_id"]: {c: row[c] for c in metadata_cols}
        for _, row in df.iterrows()
    }


def parse_grna_h5ad(h5ad_path):
    """Load a cell x gRNA AnnData and return its key summary stats."""
    adata = ad.read_h5ad(h5ad_path)
    mat = adata.X
    total_reads = float(mat.sum())
    n_cells = adata.n_obs
    n_guides = adata.n_vars
    guides_detected_per_cell = np.asarray((mat > 0).sum(axis=1)).flatten()
    return {
        "grna_total_reads": total_reads,
        "grna_n_cells": n_cells,
        "grna_n_guides": n_guides,
        "grna_mean_reads_per_cell": total_reads / n_cells if n_cells else float("nan"),
        "grna_median_guides_detected_per_cell": float(np.median(guides_detected_per_cell)) if n_cells else float("nan"),
    }


def build_summary(crispr_h5ads, sample_ids, samplesheet_path):
    if len(crispr_h5ads) != len(sample_ids):
        raise ValueError(
            f"Got {len(crispr_h5ads)} --crispr-h5ads but {len(sample_ids)} "
            f"--sample-ids; these must be aligned 1:1."
        )
    metadata_by_sample = load_samplesheet_metadata(samplesheet_path)

    rows = []
    for sample_id, h5ad_path in zip(sample_ids, crispr_h5ads):
        row = {"sample_id": sample_id}
        row.update(metadata_by_sample.get(sample_id, {}))
        row["grna_h5ad"] = h5ad_path
        row.update(parse_grna_h5ad(h5ad_path))
        rows.append(row)

    return pd.DataFrame(rows)


def add_downsampling_targets(df, target_reads_per_cell=None):
    """Add grna_target_mean_reads_per_cell and grna_downsample_fraction columns."""
    df = df.copy()

    bad = df["grna_mean_reads_per_cell"].isna() | (df["grna_mean_reads_per_cell"] <= 0)
    if bad.any():
        print(
            f"WARNING: {bad.sum()} sample(s) have missing/zero grna_mean_reads_per_cell: "
            f"{df.loc[bad, 'sample_id'].tolist()}",
            file=sys.stderr,
        )

    if target_reads_per_cell is None:
        target_reads_per_cell = df["grna_mean_reads_per_cell"].min()
    df["grna_target_mean_reads_per_cell"] = target_reads_per_cell
    df["grna_downsample_fraction"] = (
        target_reads_per_cell / df["grna_mean_reads_per_cell"]
    ).clip(upper=1.0)
    return df, target_reads_per_cell


def main(argv=None):
    args = parse_args(argv)
    sample_ids = [s.strip() for s in args.sample_ids.split(",") if s.strip()]

    summary = build_summary(args.crispr_h5ads, sample_ids, args.samplesheet)
    summary, target_reads_per_cell = add_downsampling_targets(
        summary, target_reads_per_cell=args.target_reads_per_cell
    )

    print(summary.to_string(index=False))
    summary.to_csv(args.output, index=False)
    print(f"\nWrote gRNA summary to {args.output}")

    target_depth = int(round(target_reads_per_cell))
    with open(args.target_depth_output, "w") as fh:
        fh.write(f"{target_depth}\n")
    print(f"Resolved gRNA target depth: {target_depth} reads/cell -> {args.target_depth_output}")


if __name__ == "__main__":
    main()
