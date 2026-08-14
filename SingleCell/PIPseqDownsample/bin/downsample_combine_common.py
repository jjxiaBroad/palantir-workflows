"""Shared helpers for combining one unit's downsampled GEX + gRNA (CRISPR Guide
Capture) data into a single AnnData.

Used by combine_downsampled_batch.py (one AnnData per sample, concatenated
across samples at one common depth) and combine_downsampled_depths.py (one
AnnData per depth, concatenated across depths for one sample). Both scripts
merge the same two on-disk shapes -- a downsample_molecule_info.py GEX
filtered_matrix directory and a downsample_crispr_anndata.py <depth>rpc.h5ad
file -- so the barcode-normalization and feature-axis-merge logic lives here
once rather than being duplicated per script.
"""

import glob
import re
import sys

import anndata as ad
import scanpy as sc

BARCODE_SUFFIX_RE = re.compile(r"-\d+$")


def find_one(pattern, what):
    """Return the single file/dir matching a glob pattern, warning (not
    failing) if more than one match is found."""
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise FileNotFoundError(f"No {what} found for pattern: {pattern}")
    if len(matches) > 1:
        print(f"WARNING: multiple {what} matches for {pattern}; using first: {matches}", file=sys.stderr)
    return matches[0]


def load_gex_from_matrix_dir(filtered_matrix_dir):
    """Load a downsample_molecule_info.py filtered_matrix/ dir (10x-mtx format)."""
    adata = sc.read_10x_mtx(filtered_matrix_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    adata.var["feature_types"] = "Gene Expression"
    return adata


def load_downsampled_grna_from_h5ad(h5ad_path):
    """Load a downsample_crispr_anndata.py <depth>rpc.h5ad file."""
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names_make_unique()
    adata.var["feature_types"] = "CRISPR Guide Capture"
    return adata


def normalize_barcodes(adata):
    """Strip the trailing 10x barcode suffix (e.g. "-1") so GEX and gRNA
    barcodes -- which may have been suffixed independently -- line up."""
    adata.obs_names = adata.obs_names.str.replace(BARCODE_SUFFIX_RE, "", regex=True)
    return adata


def merge_gex_and_grna(gex, grna, label):
    """Normalize barcodes on both sides and merge GEX + gRNA along the
    feature axis (inner join on barcode). GEX and gRNA are expected to
    already share the same filtered-cell barcodes; a mismatch is reported
    (only shared barcodes are kept) rather than silently ignored.

    ``label`` is used only for the warning/progress messages (e.g. a sample
    ID, or a "<sample_id> @ <depth>rpc" string).
    """
    gex = normalize_barcodes(gex)
    grna = normalize_barcodes(grna)

    shared = gex.obs_names.intersection(grna.obs_names)
    if len(shared) != gex.n_obs or len(shared) != grna.n_obs:
        print(
            f"WARNING: {label}: GEX has {gex.n_obs} filtered cells, gRNA has "
            f"{grna.n_obs} filtered cells, only {len(shared)} barcodes are shared -- "
            "only shared barcodes will be kept.",
            file=sys.stderr,
        )

    combined = ad.concat([gex, grna], axis=1, join="inner", merge="unique")
    print(f"{label}: GEX {gex.shape} + downsampled gRNA {grna.shape} -> combined {combined.shape}")
    return combined
