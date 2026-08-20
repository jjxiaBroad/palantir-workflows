#!/usr/bin/env python3
"""Summarize DRAGEN scRNA metrics across a batch of samples and compute the
common GEX downsampling target needed to normalize sequencing depth per cell
across the batch before cross-sample comparison.

For each sample this reads its DRAGEN ``<prefix>.scRNA_metrics.csv`` (an
explicit, required path per sample -- no filename-convention guessing) and
pulls out a handful of QC metrics, keyed by the exact metric name DRAGEN
writes in column 3 of that CSV -- see METRIC_FIELDS below. "Mean reads per
cell (Total input reads / Passing cells)" is DRAGEN's own precomputed value,
so no separate read-count/cell-count arithmetic is needed here.

The common target -- ``target_mean_reads_per_cell`` -- defaults to the
minimum "Mean reads per cell" observed across the batch (you can only
downsample down, never up), or an explicit --target-reads-per-cell override.
This value (rounded to the nearest integer) is what
``downsample_molecule_info.py --matrix-depths`` should be given for every
sample in the batch.

Every samplesheet column other than sample_id/molecule_info_h5/
scrna_metrics_csv/filtered_barcodes_tsv/features_tsv/crispr_h5ad is carried
through into the output summary verbatim, whatever a given batch's
samplesheet happens to include.

Usage:
    python summarize_gex_downsample_targets.py \\
        --scrna-metrics-csvs sample1.scRNA_metrics.csv sample2.scRNA_metrics.csv \\
        --sample-ids sample1,sample2 \\
        --samplesheet samplesheet.csv \\
        --output gex_downsample_summary.csv \\
        --target-depth-output gex_target_depth.txt
"""

import argparse
import csv
import sys

import pandas as pd

# Metrics pulled out of each scRNA_metrics.csv, keyed by the exact metric
# name DRAGEN writes in column 3 of the CSV.
METRIC_FIELDS = {
    "total_input_reads": "Total input reads",
    "total_barcoded_reads": "Total barcoded reads",
    "mapped_reads": "Mapped reads",
    "total_molecules": "Total molecules",
    "sequencing_saturation": "Sequencing saturation",
    "passing_cells": "Passing cells",
    "fraction_reads_in_passing_cells": "Fraction of reads in passing cells",
    "mean_reads_per_cell": "Mean reads per cell (Total input reads / Passing cells)",
    "median_molecules_per_cell": "Median molecules per passing cell",
    "median_genes_per_cell": "Median genes per passing cell",
    "total_genes_detected": "Total genes detected",
}

# Samplesheet columns that locate files rather than describe sample metadata.
SAMPLESHEET_PATH_COLUMNS = {
    "sample_id", "molecule_info_h5", "scrna_metrics_csv", "filtered_barcodes_tsv",
    "features_tsv", "crispr_h5ad",
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--scrna-metrics-csvs",
        nargs="+",
        required=True,
        help="DRAGEN scRNA_metrics.csv files, one per sample, aligned with --sample-ids.",
    )
    parser.add_argument(
        "--sample-ids",
        required=True,
        help="Comma-separated sample IDs, aligned with --scrna-metrics-csvs.",
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="Batch samplesheet CSV, used only to pull per-sample metadata columns "
             "(path columns are ignored here).",
    )
    parser.add_argument(
        "--target-reads-per-cell",
        type=float,
        default=None,
        help="Target mean reads per cell for the batch. Defaults to the minimum "
             "'Mean reads per cell' observed across all samples.",
    )
    parser.add_argument(
        "--output",
        default="gex_downsample_summary.csv",
        help="Path to write the summary CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--target-depth-output",
        default="gex_target_depth.txt",
        help="Path to write the resolved integer target depth to (default: %(default)s)",
    )
    return parser.parse_args(argv)


def _to_number(value):
    """Cast a CSV string to int/float where possible, else leave as-is."""
    value = value.strip()
    try:
        if "." in value or "e" in value.lower():
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_scrna_metrics(csv_path):
    """Parse a DRAGEN *.scRNA_metrics.csv into {metric name: value}."""
    metrics = {}
    with open(csv_path, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) < 4:
                continue
            metric_name, value = row[2].strip(), row[3]
            metrics[metric_name] = _to_number(value)
    return metrics


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


def build_summary(scrna_metrics_csvs, sample_ids, samplesheet_path):
    if len(scrna_metrics_csvs) != len(sample_ids):
        raise ValueError(
            f"Got {len(scrna_metrics_csvs)} --scrna-metrics-csvs but "
            f"{len(sample_ids)} --sample-ids; these must be aligned 1:1."
        )
    metadata_by_sample = load_samplesheet_metadata(samplesheet_path)

    rows = []
    for sample_id, metrics_path in zip(sample_ids, scrna_metrics_csvs):
        metrics = parse_scrna_metrics(metrics_path)

        row = {"sample_id": sample_id}
        row.update(metadata_by_sample.get(sample_id, {}))
        row["metrics_file"] = metrics_path

        for out_field, metric_name in METRIC_FIELDS.items():
            if metric_name not in metrics:
                print(
                    f"WARNING: metric '{metric_name}' not found for sample "
                    f"'{sample_id}' in {metrics_path}",
                    file=sys.stderr,
                )
            row[out_field] = metrics.get(metric_name)
        rows.append(row)

    return pd.DataFrame(rows)


def add_downsampling_targets(df, target_reads_per_cell=None):
    """Add target_mean_reads_per_cell and downsample_fraction columns.

    target_reads_per_cell: optional override for the common target. If not
    given, defaults to the minimum 'mean_reads_per_cell' across all rows in
    df (i.e. normalize every sample down to the shallowest sample).
    """
    df = df.copy()
    if target_reads_per_cell is None:
        target_reads_per_cell = df["mean_reads_per_cell"].min()

    df["target_mean_reads_per_cell"] = target_reads_per_cell
    df["downsample_fraction"] = (
        target_reads_per_cell / df["mean_reads_per_cell"]
    ).clip(upper=1.0)
    return df, target_reads_per_cell


def main(argv=None):
    args = parse_args(argv)
    sample_ids = [s.strip() for s in args.sample_ids.split(",") if s.strip()]

    summary = build_summary(args.scrna_metrics_csvs, sample_ids, args.samplesheet)
    summary, target_reads_per_cell = add_downsampling_targets(
        summary, target_reads_per_cell=args.target_reads_per_cell
    )

    print(summary.to_string(index=False))
    summary.to_csv(args.output, index=False)
    print(f"\nWrote summary to {args.output}")

    target_depth = int(round(target_reads_per_cell))
    with open(args.target_depth_output, "w") as fh:
        fh.write(f"{target_depth}\n")
    print(f"Resolved GEX target depth: {target_depth} reads/cell -> {args.target_depth_output}")


if __name__ == "__main__":
    main()
