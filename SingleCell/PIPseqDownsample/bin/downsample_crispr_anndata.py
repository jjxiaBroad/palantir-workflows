#!/usr/bin/env python3
"""Standalone script: downsample a cell x gRNA AnnData (.h5ad) via direct
binomial thinning of raw read counts.

This assumes ``.X`` (or ``--layer``) already holds a cells x guides matrix of
plain, undeduplicated read counts -- the same assumption used elsewhere in
this project: gRNA/CRISPR Guide Capture data has no meaningful binning index
and no molecule-level deduplication, so each read independently survives
downsampling with probability equal to the downsample factor. The
downsampled count per (cell, guide) is simply

    Binomial(original_count, downsample_factor)

applied directly to the sparse matrix's stored (nonzero) values -- there is
no unique-molecule tracking and no correction step, matching how the raw
counts were generated in the first place.

"Full depth" (the reads-per-cell baseline that --matrix-depths ratios are
computed against) is derived from the AnnData itself: total counts in
``.X`` (or ``--layer``) divided by the number of cells (``adata.n_obs``).
This intentionally does not require any external DRAGEN metrics file --
if you need to downsample against a *true* sequencing-library total instead
(e.g. because this AnnData only contains a post-QC cell subset), pass
--total-input-reads to override it.

For each requested reads-per-cell depth, this writes a full copy of the
input AnnData with the thinned counts in ``.X`` (or ``--layer``) to its own
.h5ad file, e.g. ``out_dir/5000rpc.h5ad``. All other AnnData content
(.obs, .var, .obsm, .uns, etc.) is carried through unchanged.

Usage:
    python downsample_crispr_anndata.py \\
        --input-h5ad /path/to/crispr_counts.h5ad \\
        --output-dir /path/to/out_crispr \\
        --matrix-depths 5000 10000 20000
"""

import argparse
import os

import numpy as np
from scipy import sparse
import anndata as ad


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            'Downsample a cell x gRNA AnnData (.h5ad) via direct binomial '
            'thinning of raw read counts. No IMI/IPM -- this modality has no '
            'meaningful binning index and no deduplication.'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--input-h5ad',
        required=True,
        help='Path to the cell x gRNA .h5ad file. .X (or --layer) must hold '
             'raw, undeduplicated read counts.',
    )
    parser.add_argument(
        '--layer',
        default=None,
        help='Use adata.layers[LAYER] as the raw-count source instead of .X. '
             'The downsampled output is written back to the same slot '
             '(.X or that layer) in each output file.',
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to write one .h5ad file per requested depth into. '
             'Created if it does not exist.',
    )
    parser.add_argument(
        '--matrix-depths',
        type=int,
        nargs='+',
        required=True,
        help='Reads-per-cell depths to downsample to. One <N>rpc.h5ad file '
             'is written per depth.',
    )
    parser.add_argument(
        '--total-input-reads',
        type=int,
        default=None,
        help='Override the full-depth read total used as the reads-per-cell '
             'denominator. Defaults to the sum of all counts in .X (or '
             '--layer). Use this if the AnnData is a post-QC cell subset and '
             'you want depths computed against the true sequencing-library '
             'total instead.',
    )
    parser.add_argument(
        '--min-reads-per-cell',
        type=int,
        default=1,
        help='Drop requested depths below this floor.',
    )
    parser.add_argument(
        '--random-seed',
        type=int,
        default=42,
        help='Seed for the binomial thinning draws, for reproducibility.',
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Silence per-step progress logging.',
    )
    return parser.parse_args(argv)


def get_counts_matrix(adata, layer=None):
    """Return the source counts matrix as CSR, whatever format it's stored
    in (dense, or any scipy sparse format)."""
    mat = adata.X if layer is None else adata.layers[layer]
    if sparse.issparse(mat):
        return mat.tocsr()
    return sparse.csr_matrix(mat)


def check_looks_like_raw_counts(mat, quiet=False):
    """Best-effort sanity check: warn (don't fail) if the matrix doesn't
    look like non-negative integer read counts -- binomial thinning is only
    meaningful on raw counts, not normalized/log-transformed values."""
    if mat.nnz == 0:
        return
    sample = mat.data[: min(len(mat.data), 100_000)]
    if np.any(sample < 0) or not np.allclose(sample, np.round(sample)):
        print(
            'WARNING: values in the input matrix do not look like raw, '
            'non-negative integer read counts (found negative and/or '
            'non-integer values). Binomial thinning assumes raw counts -- '
            'check --layer or that .X has not already been normalized.'
        )


def compute_downsample_targets(reads_per_cell_full, requested_rpcs, min_reads_per_cell):
    """Turn requested reads-per-cell values into (target_rpc, factor) pairs,
    deepest-first. Values below min_reads_per_cell, or above full depth
    (which would mean upsampling), are dropped. A requested depth exactly at
    full depth (factor == 1) is kept as a no-op "downsample" -- this is the
    normal case for whichever sample defines the batch's target depth."""
    target_rpcs = []
    factors = []
    if reads_per_cell_full > 0:
        for rpc in requested_rpcs:
            if rpc < min_reads_per_cell:
                continue
            factor = rpc / reads_per_cell_full
            if factor > 1:
                continue
            target_rpcs.append(rpc)
            factors.append(factor)
    pairs = sorted(zip(factors, target_rpcs), reverse=True)
    return [p[1] for p in pairs], [p[0] for p in pairs]


def downsample_binomial_sparse(mat, factor, rng):
    """Thin a sparse counts matrix's stored (nonzero) values by independent
    binomial draws: each read survives with probability = factor, so the
    downsampled count per entry is Binomial(original_count, factor).

    Implicit zeros in the sparse matrix are already zero and need no work --
    only ``.data`` (the nonzero entries) is touched. ``factor >= 1`` is a
    no-op copy (nothing to thin). Entries that thin down to exactly zero are
    dropped from the sparsity structure via ``eliminate_zeros()``.
    """
    mat = mat.tocsr()
    thinned = mat.copy()
    if factor < 1:
        thinned.data = rng.binomial(thinned.data.astype(np.int64), factor)
        thinned.eliminate_zeros()
    thinned.data = thinned.data.astype(np.int64)
    return thinned


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)

    if not args.quiet:
        print(f'Loading {args.input_h5ad} ...')
    adata = ad.read_h5ad(args.input_h5ad)
    counts = get_counts_matrix(adata, layer=args.layer)
    check_looks_like_raw_counts(counts, quiet=args.quiet)

    n_cells = adata.n_obs
    if args.total_input_reads is not None:
        total_reads = args.total_input_reads
        source = '--total-input-reads override'
    else:
        total_reads = int(counts.sum())
        source = f'sum of counts in {"X" if args.layer is None else args.layer!r}'
    reads_per_cell_full = total_reads / n_cells

    if not args.quiet:
        print(
            f'Full-depth reads/cell: {reads_per_cell_full:,.1f} '
            f'(total={total_reads:,} from {source}, n_cells={n_cells:,})'
        )

    requested = {int(rpc) for rpc in args.matrix_depths}
    target_rpcs, factors = compute_downsample_targets(
        reads_per_cell_full, requested, args.min_reads_per_cell
    )
    if not target_rpcs:
        raise SystemExit(
            'No requested depths remain after filtering -- check '
            '--min-reads-per-cell and the full-depth value printed above '
            '(all requested depths were either below the floor or at/above '
            'full depth).'
        )

    for target_rpc, factor in zip(target_rpcs, factors):
        depth_int = int(round(target_rpc))
        if not args.quiet:
            print(f'Downsampling to {depth_int:,} reads/cell (factor={factor:.4f}) ...')
        thinned = downsample_binomial_sparse(counts, factor, rng)

        out_adata = adata.copy()
        if args.layer is None:
            out_adata.X = thinned
        else:
            out_adata.layers[args.layer] = thinned

        out_path = os.path.join(args.output_dir, f'{depth_int}rpc.h5ad')
        out_adata.write_h5ad(out_path)
        if not args.quiet:
            print(f'  wrote {out_path}')

    print('Done.')


if __name__ == '__main__':
    main()
