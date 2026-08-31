#!/usr/bin/env python3
"""Combine one sample's downsampled GEX + gRNA (CRISPR Guide Capture) data
across multiple depths into a single AnnData, for within-sample depth-response
analysis (e.g. how per-cell metrics change across a depth series).

Every depth written by either modality ends up in the output -- the GEX depths
(written by downsample_molecule_info.py under
<gex-matrix-dir>/<depth>rpc/filtered_matrix/) and the gRNA depths (written by
downsample_crispr_anndata.py under <crispr-h5ad-dir>/<depth>rpc.h5ad) are
independent params on different scales, so they are *not* required to line up.

Per depth, one cell-axis block is built from whatever that depth has:

  - both modalities  -> barcodes normalized on both sides, then merged along
                        the feature axis (inner join on barcode), the same
                        per-unit merge combine_downsampled_batch.py does per
                        sample (shared via downsample_combine_common.py);
                        .obs.modality == "gex+grna"
  - GEX only         -> that depth's GEX matrix alone; .obs.modality == "gex"
  - gRNA only        -> that depth's gRNA AnnData alone; .obs.modality == "grna"

Blocks are then concatenated along the cell axis (outer join on features,
missing modality zero-filled) into one AnnData, with integer `depth`,
`modality` and constant `sample_id` columns attached to .obs. Because a
single-modality block has no counts at all for the other modality's features,
`modality` is what distinguishes a structural zero there from a measured zero
-- filter on it before comparing feature sets across depths.

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

from downsample_combine_common import (
    load_gex_from_matrix_dir,
    load_downsampled_grna_from_h5ad,
    merge_gex_and_grna,
    normalize_barcodes,
)

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

    all_depths = sorted(set(gex_depths) | set(crispr_depths))
    if not all_depths:
        raise FileNotFoundError(
            f"No downsampled depth found in either {args.gex_matrix_dir} or "
            f"{args.crispr_h5ad_dir}."
        )
    print(
        f"GEX depths found: {sorted(gex_depths)}\n"
        f"gRNA depths found: {sorted(crispr_depths)}\n"
        f"Combining all {len(all_depths)} depth(s): {all_depths}"
    )

    per_depth = []
    modalities = []
    for depth in all_depths:
        label = f"{args.sample_id} @ {depth}rpc"
        gex_dir = gex_depths.get(depth)
        crispr_h5ad = crispr_depths.get(depth)
        if gex_dir is not None and crispr_h5ad is not None:
            gex = load_gex_from_matrix_dir(gex_dir)
            grna = load_downsampled_grna_from_h5ad(crispr_h5ad)
            block = merge_gex_and_grna(gex, grna, label=label)
            modality = "gex+grna"
        elif gex_dir is not None:
            block = normalize_barcodes(load_gex_from_matrix_dir(gex_dir))
            modality = "gex"
            print(f"{label}: GEX only {block.shape} (no gRNA h5ad at this depth)")
        else:
            block = normalize_barcodes(load_downsampled_grna_from_h5ad(crispr_h5ad))
            modality = "grna"
            print(f"{label}: gRNA only {block.shape} (no GEX matrix at this depth)")
        block.obs["modality"] = modality
        per_depth.append(block)
        modalities.append(modality)

    # join="outer" so single-modality depths keep their own features and get the other
    # modality's zero-filled; fill_value=0 keeps X sparse (the default would give NaN).
    combined = ad.concat(
        per_depth,
        axis=0,
        join="outer",
        label="depth",
        keys=[str(depth) for depth in all_depths],
        index_unique="_",
        fill_value=0,
    )
    combined.obs["depth"] = combined.obs["depth"].astype(int)
    combined.obs["sample_id"] = args.sample_id

    if set(modalities) != {"gex+grna"}:
        print(
            "WARNING: not every depth has both modalities -- "
            f"{modalities.count('gex')} GEX-only and {modalities.count('grna')} gRNA-only "
            "depth(s) are included with the other modality's features zero-filled. "
            "Filter on .obs['modality'] before treating those zeros as measured counts.",
            file=sys.stderr,
        )

    print(f"\nFinal combined AnnData: {combined}")
    combined.write_h5ad(args.output_h5ad, compression="gzip")
    print(f"\nWrote combined AnnData to {args.output_h5ad}")


if __name__ == "__main__":
    main()
