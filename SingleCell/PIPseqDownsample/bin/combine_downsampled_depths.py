#!/usr/bin/env python3
"""Combine one sample's downsampled GEX + gRNA (CRISPR Guide Capture) data
across multiple depths into a single AnnData, for within-sample depth-response
analysis (e.g. how per-cell metrics change across a depth series).

For each depth present in *both* the GEX matrix dir (written by
downsample_molecule_info.py under <gex-matrix-dir>/<depth>rpc/filtered_matrix/)
and the gRNA h5ad dir (written by downsample_crispr_anndata.py under
<crispr-h5ad-dir>/<depth>rpc.h5ad): normalize the trailing 10x barcode suffix
on both sides so barcodes line up, then merge the two modalities along the
feature axis (inner join on barcode) -- the same per-unit merge
combine_downsampled_batch.py does per sample, shared via
downsample_combine_common.py.

All depths are then concatenated along the cell axis into one combined
AnnData, with an integer `depth` column and a constant `sample_id` column
attached to .obs.

Usage:
    python combine_downsampled_depths.py \\
        --gex-matrix-dir sample1_gex/ \\
        --crispr-h5ad-dir sample1_crispr/ \\
        --sample-id sample1 \\
        --output-h5ad sample1.depths_combined.h5ad
"""

import argparse
import os
import re
import sys

import anndata as ad

from downsample_combine_common import load_gex_from_matrix_dir, load_downsampled_grna_from_h5ad, merge_gex_and_grna

GEX_DEPTH_DIR_RE = re.compile(r"^(\d+)rpc$")
CRISPR_DEPTH_FILE_RE = re.compile(r"^(\d+)rpc\.h5ad$")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--gex-matrix-dir",
        required=True,
        help="This sample's DOWNSAMPLE_MOLECULE_INFO output dir (one <depth>rpc/ subdir per requested depth).",
    )
    parser.add_argument(
        "--crispr-h5ad-dir",
        required=True,
        help="This sample's DOWNSAMPLE_CRISPR_ANNDATA output dir (one <depth>rpc.h5ad file per requested depth).",
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Sample ID, attached to the combined AnnData's .obs.",
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help="Path to write the combined AnnData to.",
    )
    return parser.parse_args(argv)


def discover_gex_depths(gex_matrix_dir):
    """Return {depth: filtered_matrix dir} for every <depth>rpc/filtered_matrix/
    subdir found directly under gex_matrix_dir."""
    depths = {}
    for name in sorted(os.listdir(gex_matrix_dir)):
        match = GEX_DEPTH_DIR_RE.match(name)
        if not match:
            continue
        filtered_dir = os.path.join(gex_matrix_dir, name, "filtered_matrix")
        if os.path.isdir(filtered_dir):
            depths[int(match.group(1))] = filtered_dir
    return depths


def discover_crispr_depths(crispr_h5ad_dir):
    """Return {depth: h5ad path} for every <depth>rpc.h5ad file found directly
    under crispr_h5ad_dir."""
    depths = {}
    for name in sorted(os.listdir(crispr_h5ad_dir)):
        match = CRISPR_DEPTH_FILE_RE.match(name)
        if match:
            depths[int(match.group(1))] = os.path.join(crispr_h5ad_dir, name)
    return depths


def main(argv=None):
    args = parse_args(argv)

    gex_depths = discover_gex_depths(args.gex_matrix_dir)
    crispr_depths = discover_crispr_depths(args.crispr_h5ad_dir)

    common_depths = sorted(set(gex_depths) & set(crispr_depths))
    gex_only = sorted(set(gex_depths) - set(crispr_depths))
    crispr_only = sorted(set(crispr_depths) - set(gex_depths))
    if gex_only or crispr_only:
        print(
            f"WARNING: depths present in only one modality are skipped -- "
            f"GEX-only: {gex_only}, gRNA-only: {crispr_only}",
            file=sys.stderr,
        )
    if not common_depths:
        raise FileNotFoundError(
            f"No depth present in both {args.gex_matrix_dir} and {args.crispr_h5ad_dir}."
        )

    per_depth = []
    for depth in common_depths:
        gex = load_gex_from_matrix_dir(gex_depths[depth])
        grna = load_downsampled_grna_from_h5ad(crispr_depths[depth])
        per_depth.append(merge_gex_and_grna(gex, grna, label=f"{args.sample_id} @ {depth}rpc"))

    combined = ad.concat(
        per_depth,
        axis=0,
        join="outer",
        label="depth",
        keys=[str(depth) for depth in common_depths],
        index_unique="_",
    )
    combined.obs["depth"] = combined.obs["depth"].astype(int)
    combined.obs["sample_id"] = args.sample_id

    print(f"\nFinal combined AnnData: {combined}")
    combined.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"\nWrote combined AnnData to {args.output_h5ad}")


if __name__ == "__main__":
    main()
