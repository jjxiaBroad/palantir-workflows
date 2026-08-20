#!/usr/bin/env python3
"""Standalone command-line tool: downsample IMI counts from a DRAGEN scRNA
molecule-info h5 file.

Given a DRAGEN scRNA results directory (or a specific molecule-info h5), this
script can:

  1. Write out downsampled counts matrices at one or more reads-per-cell
     depths. By default both the raw (all barcodes) and cell-filtered layouts
     are produced, following DRAGEN's ``matrix.mtx.gz`` / ``barcodes.tsv.gz``
     / ``features.tsv.gz`` convention.
  2. Optionally run a saturation sweep over a built-in reads-per-cell ladder
     (plus any extra depths provided via ``--saturation-extra-depths``) and
     save a CSV describing median transcripts, median genes, and % sequencing
     saturation at each depth.

Only IMI counting is supported (DRAGEN 4.6+ with ``/binning_idx`` in the h5).
The script does not require an external features or barcodes file - both are
read directly from the h5.

Typical invocations:

    # Save a downsampled matrix at a single depth (raw + filtered layouts)
    python downsample_molecule_info.py \\
        --dragen-results-dir /path/to/results \\
        --output-dir /path/to/out \\
        --matrix-depths 20000

    # Add a saturation sweep and a few extra depths
    python downsample_molecule_info.py \\
        --dragen-results-dir /path/to/results \\
        --output-dir /path/to/out \\
        --matrix-depths 10000 20000 \\
        --saturation \\
        --saturation-extra-depths 40 80 160
"""

import argparse
import gzip
import glob
import multiprocessing
import os
import warnings

import numpy as np
import pandas as pd
import h5py
from scipy.sparse import coo_matrix, hstack, vstack


# =====================================================================
# I/O helpers: locate/read a DRAGEN scRNA results directory's well-known
# files, and parse "Total input reads" out of the metrics CSV.
# =====================================================================

DRAGEN_MOLECULE_INFO_SUFFIX = 'scRNA.moleculeInfo.h5'
DRAGEN_FILTERED_BARCODES_SUFFIX = 'scRNA.filtered.barcodes.tsv.gz'
DRAGEN_METRICS_SUFFIX = 'scRNA_metrics.csv'
DRAGEN_FEATURES_SUFFIX = 'scRNA.features.tsv.gz'


def smart_read_textfile(file_input_path):
    """Read a text file (optionally gzip-compressed) into a list of stripped lines.

    Trailing empty lines are dropped so the returned list matches ``wc -l`` on
    a well-formed newline-terminated file. Any ``.gz`` file is decompressed
    transparently based on its extension.
    """
    path = str(file_input_path)
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as fid:
        text = fid.read()
    lines = text.split('\n')
    # Drop the trailing empty entry produced by a final newline, if any.
    if lines and lines[-1] == '':
        lines.pop()
    return lines


def _files_with_suffix(dragen_results_dir, suffix):
    """Return the sorted list of files in ``dragen_results_dir`` whose basename
    ends with ``suffix``. Does not recurse; DRAGEN emits scRNA outputs directly
    in the results directory."""
    pattern = os.path.join(dragen_results_dir, f'*{suffix}')
    return sorted(glob.glob(pattern))


def detect_dragen_prefix(dragen_results_dir):
    """Detect the DRAGEN output-file-prefix used in a results directory.

    DRAGEN prepends an optional prefix (configured via ``output-file-prefix``)
    to every scRNA output file, followed by a ``.`` separator. For example,
    ``mysample.scRNA.moleculeInfo.h5``. If no prefix is configured, files are
    named ``scRNA.moleculeInfo.h5`` with an empty prefix.

    Returns the prefix as a string (including the trailing ``.`` when
    non-empty) so callers can concatenate it directly, e.g.
    ``prefix + DRAGEN_FILTERED_BARCODES_SUFFIX``. Returns an empty string if
    no molecule-info h5 file is found or the file has no prefix.
    """
    matches = _files_with_suffix(dragen_results_dir, DRAGEN_MOLECULE_INFO_SUFFIX)
    if not matches:
        return ''
    basename = os.path.basename(matches[0])
    prefix_part = basename[: -len(DRAGEN_MOLECULE_INFO_SUFFIX)]
    # prefix_part will end with '.' when non-empty (e.g. 'sample1.'), or be
    # empty for prefix-less runs. Either way, returning it verbatim is what
    # concatenating callers expect.
    return prefix_part


def find_molecule_info_h5(dragen_results_dir):
    """Return the path to the DRAGEN molecule-info h5 in ``dragen_results_dir``.

    Raises FileNotFoundError if there is no match, and RuntimeError if the
    directory contains multiple candidates (which would mean the user pointed
    at a folder containing more than one DRAGEN run).
    """
    matches = _files_with_suffix(dragen_results_dir, DRAGEN_MOLECULE_INFO_SUFFIX)
    if not matches:
        raise FileNotFoundError(
            f'No molecule-info h5 file (*{DRAGEN_MOLECULE_INFO_SUFFIX}) found in '
            f'{dragen_results_dir}. Pass --molecule-info-h5 explicitly, or check '
            f'that this is a DRAGEN scRNA results directory.'
        )
    if len(matches) > 1:
        raise RuntimeError(
            f'Multiple molecule-info h5 files found in {dragen_results_dir}: '
            f'{matches}. Pass --molecule-info-h5 explicitly to disambiguate.'
        )
    return matches[0]


def find_filtered_barcodes(dragen_results_dir):
    """Return the DRAGEN filtered-barcodes file in ``dragen_results_dir``, or
    ``None`` if the file is not present. The caller decides whether that's an
    error (the CLI treats absence as an error only when the filtered-barcodes
    list is actually needed)."""
    matches = _files_with_suffix(dragen_results_dir, DRAGEN_FILTERED_BARCODES_SUFFIX)
    if len(matches) > 1:
        raise RuntimeError(
            f'Multiple filtered-barcodes files found in {dragen_results_dir}: '
            f'{matches}. Pass --filtered-barcodes explicitly to disambiguate.'
        )
    return matches[0] if matches else None


def find_scrna_metrics_csv(dragen_results_dir):
    """Return the DRAGEN scRNA metrics CSV in ``dragen_results_dir``, or
    ``None`` if not present. See find_filtered_barcodes for the multi-match
    behavior."""
    matches = _files_with_suffix(dragen_results_dir, DRAGEN_METRICS_SUFFIX)
    if len(matches) > 1:
        raise RuntimeError(
            f'Multiple scRNA_metrics.csv files found in {dragen_results_dir}: '
            f'{matches}. Pass --scrna-metrics-csv explicitly to disambiguate.'
        )
    return matches[0] if matches else None


def find_features_tsv(dragen_results_dir):
    """Return the DRAGEN raw features.tsv.gz in ``dragen_results_dir``, or
    ``None`` if not present.

    DRAGEN emits both ``<prefix>.scRNA.features.tsv.gz`` (all barcodes) and
    ``<prefix>.scRNA.filtered.features.tsv.gz`` (filtered cells) - both end
    with ``scRNA.features.tsv.gz``. We deliberately pick the raw version by
    excluding paths whose basename contains ``.filtered.``, since the raw
    file lists every feature the molecule-info h5 references.
    """
    matches = _files_with_suffix(dragen_results_dir, DRAGEN_FEATURES_SUFFIX)
    matches = [m for m in matches if '.filtered.' not in os.path.basename(m)]
    if len(matches) > 1:
        raise RuntimeError(
            f'Multiple raw features.tsv.gz files found in {dragen_results_dir}: '
            f'{matches}. Pass --features-tsv explicitly to disambiguate.'
        )
    return matches[0] if matches else None


def parse_total_input_reads(scrna_metrics_csv):
    """Parse the "Total input reads" scalar from a DRAGEN scRNA_metrics.csv.

    The metrics CSV has no header row; each line is
    ``<metric_group>,<sample>,<metric>,<value>[,<pct>]``. We scan for the
    first row where ``metric == 'Total input reads'`` and return its value
    as an int.

    Returns ``None`` if the file has no such row (rather than raising) so the
    CLI can fall back to the h5-derived mapped-read total with a warning.
    """
    if not os.path.isfile(scrna_metrics_csv):
        raise FileNotFoundError(f'scRNA metrics CSV not found: {scrna_metrics_csv}')
    with open(scrna_metrics_csv, 'r', encoding='utf-8') as fid:
        for line in fid:
            parts = [p.strip() for p in line.rstrip('\n').split(',')]
            if len(parts) < 4:
                continue
            metric = parts[2]
            if metric == 'Total input reads':
                value = parts[3]
                try:
                    return int(value)
                except ValueError:
                    # DRAGEN occasionally writes floats/scientific notation.
                    return int(float(value))
    return None


# =====================================================================
# Matrix: minimal sparse counts matrix representation, enough to hold a
# downsampled result in memory and write it out in DRAGEN's
# matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz format.
# =====================================================================

# Standard MatrixMarket header used by DRAGEN's scRNA output.
_MTX_HEADER_LINE_1 = '%%MatrixMarket matrix coordinate integer general\n'
_MTX_HEADER_LINE_2 = '%\n'


class Matrix:
    """Base class for a sparse (barcodes x features) counts matrix.

    Subclasses (only ``MatrixData`` here) are responsible for populating
    ``self.matrix`` (CSR) plus the accompanying ``self.barcodes_seq`` and
    ``self.features`` metadata. The base class provides slice-to-filtered-
    cells and disk-serialization helpers.
    """

    def __init__(self, random_seed=42):
        self.random_seed = random_seed
        # Placeholders populated by subclasses; declared here to make the
        # expected attribute surface explicit.
        self.matrix = None
        self.barcodes_seq = None
        self.features = None
        self.filtered_cell_inds = None

    def get_filtered_cells(self, cell_barcodes_seq):
        """Restrict the matrix's ``filtered_cell_inds`` to the given barcodes.

        ``cell_barcodes_seq`` is an iterable of barcode strings (typically the
        DRAGEN filtered-barcodes list). Any barcode not present in
        ``self.barcodes_seq`` is silently ignored. The corresponding row
        indices are stored in ``self.filtered_cell_inds`` for later use by
        ``get_filtered_matrix``.
        """
        cell_inds_mask = np.isin(self.barcodes_seq, np.asarray(cell_barcodes_seq))
        self.filtered_cell_inds = np.where(cell_inds_mask)[0]
        return self.barcodes_seq[self.filtered_cell_inds]

    def get_filtered_matrix(self):
        """Return a new ``MatrixData`` restricted to ``self.filtered_cell_inds``.

        The returned matrix's barcodes list is exactly the filtered subset,
        so subsequent ``write_matrix`` calls emit the filtered layout DRAGEN
        would have produced with those cell calls.
        """
        if self.filtered_cell_inds is None:
            raise RuntimeError('Call get_filtered_cells(...) before get_filtered_matrix().')
        matrix_mask = self.matrix[self.filtered_cell_inds, :]
        matrix_barcodes_seq = self.barcodes_seq[self.filtered_cell_inds]
        matrix_mask.eliminate_zeros()
        barcodes_col, features_col = matrix_mask.nonzero()
        counts_col = matrix_mask.data
        mtx = MatrixData(
            barcodes_col=barcodes_col,
            features_col=features_col,
            counts_col=counts_col,
            barcodes_seq=matrix_barcodes_seq,
            features=self.features,
        )
        mtx.filtered_cell_inds = np.arange(len(self.filtered_cell_inds))
        return mtx

    def write_matrix(self, output_dir):
        """Write the matrix to ``output_dir`` in DRAGEN-compatible layout.

        Produces three gzip-compressed files in ``output_dir``:
        ``matrix.mtx.gz`` (MatrixMarket coordinate integer format, transposed
        from row-major CSR to the ``features x barcodes`` layout DRAGEN
        emits), ``barcodes.tsv.gz`` (one barcode per line), and
        ``features.tsv.gz`` (one feature per line, tab-separated columns).

        The output directory is created if it does not already exist.
        """
        os.makedirs(output_dir, exist_ok=True)

        # Barcodes: one per line, in the same order used by self.matrix rows.
        out_barcode_file = os.path.join(output_dir, 'barcodes.tsv.gz')
        with gzip.open(out_barcode_file, 'wt', compresslevel=1) as bc_fid:
            bc_fid.write('\n'.join(self.barcodes_seq))

        # Features: serialize the DataFrame with a tab separator, no header,
        # no index. This matches DRAGEN's features.tsv.gz format (three
        # columns: feature id, feature name, feature type).
        out_feature_file = os.path.join(output_dir, 'features.tsv.gz')
        self.features.to_csv(
            out_feature_file, sep='\t', compression='gzip', header=False, index=False
        )

        # Counts matrix in MatrixMarket coordinate integer format.
        out_counts_file = os.path.join(output_dir, 'matrix.mtx.gz')
        num_barcodes = len(self.barcodes_seq)
        num_features = len(self.features)
        self.matrix.eliminate_zeros()
        num_nonzero = self.matrix.nnz

        # DRAGEN's convention: rows are features, columns are barcodes.
        # Our CSR is (barcodes x features), so we iterate barcode-major and
        # emit "feature+1 barcode+1 count" per non-zero entry.
        line_3 = '%d %d %d\n' % (num_features, num_barcodes, num_nonzero)
        buffer_lines = [_MTX_HEADER_LINE_1, _MTX_HEADER_LINE_2, line_3]

        with gzip.open(out_counts_file, 'wt', compresslevel=1) as mat_fid:
            mat_fid.write(''.join(buffer_lines))
            buffer_lines = []
            block_size = 1_000_000
            buffer_len = 0
            num_indptr = len(self.matrix.indptr)
            for i in range(num_indptr - 1):
                start = self.matrix.indptr[i]
                end = self.matrix.indptr[i + 1]
                for idx in range(start, end):
                    feature = self.matrix.indices[idx]
                    count = self.matrix.data[idx]
                    buffer_lines.append('%d %d %d\n' % (feature + 1, i + 1, count))
                    buffer_len += 1
                if buffer_len >= block_size:
                    mat_fid.write(''.join(buffer_lines))
                    buffer_lines = []
                    buffer_len = 0
            if buffer_len > 0:
                mat_fid.write(''.join(buffer_lines))


class MatrixData(Matrix):
    """A concrete ``Matrix`` built from column-oriented COO arrays.

    ``barcodes_col``, ``features_col``, ``counts_col`` are 1-D arrays of the
    same length giving the row (barcode) index, column (feature) index, and
    count for each non-zero entry. ``barcodes_seq`` is the full barcode list
    and ``features`` is the full features DataFrame; both are needed so the
    matrix can be zero-padded to the expected shape before being converted
    to CSR.
    """

    def __init__(self, barcodes_col, features_col, counts_col, barcodes_seq, features):
        super().__init__()
        self.barcodes_seq = np.asarray(barcodes_seq)
        self.features = features

        if len(counts_col) > 0:
            sparse_mat = coo_matrix((counts_col, (barcodes_col, features_col)))

            # Downsampling may have zeroed out the last few features or
            # barcodes, in which case the derived COO shape is smaller than
            # the full metadata. Pad with zeros so the resulting matrix has
            # exactly the shape callers expect (num_barcodes x num_features).
            num_features = self.features.shape[0]
            num_barcodes = len(self.barcodes_seq)
            current_num_barcodes, current_num_features = sparse_mat.shape

            if current_num_features < num_features:
                padding = coo_matrix(
                    (sparse_mat.shape[0], num_features - current_num_features),
                    dtype=np.int64,
                )
                sparse_mat = hstack([sparse_mat, padding])
            if current_num_barcodes < num_barcodes:
                padding = coo_matrix(
                    (num_barcodes - current_num_barcodes, num_features),
                    dtype=np.int64,
                )
                sparse_mat = vstack([sparse_mat, padding])
        else:
            warnings.warn('Empty counts matrix')
            sparse_mat = coo_matrix((0, 0))

        self.matrix = sparse_mat.tocsr()


# =====================================================================
# MoleculeInfo: IMI-only molecule-info downsampler, ported from the parent
# pipeline's counting/minf.py.
#
# - Only DRAGEN's HDF5 molecule-info format (v4.6+) is supported. The CSV/text
#   parsers used for pre-4.6 outputs are not implemented here.
# - Feature and barcode metadata are read directly from the molecule-info h5,
#   so callers do not need to pass in external features/barcodes files
#   (except for the combined GEX+CRISPR case, see read_dragen_moleculeInfo_h5).
# - UMI counting is not supported. Every code path assumes IMI counting.
# - CRISPR support and get_reads_in_cells are not implemented (see
#   downsample_crispr_anndata.py for CRISPR/gRNA downsampling instead).
# - The parallel workers seed np.random deterministically per worker so
#   results are reproducible when a random_seed is provided.
# =====================================================================


class MoleculeInfo:
    """In-memory representation of a DRAGEN scRNA molecule-info h5, plus IMI
    downsampling routines.

    Instances hold the per-molecule (barcode, gene, bin, count) arrays
    extracted from the h5, along with the total mapped-read count and a
    pre-computed lookup table used by the IPM (Inferred Per Molecule)
    correction. Downsampling operates in place using these arrays.
    """

    def __init__(self, random_seed=42, parallel_pool=None, block_size=1_000_000, silent=False):
        self.random_seed = random_seed
        self.silent = silent
        self.parallel_pool = parallel_pool
        self.threads = parallel_pool._processes if parallel_pool else 1
        self.num_cores = self.threads
        self.block_size = block_size

        # Populated by read_dragen_moleculeInfo_h5.
        self.barcodes_seq = None
        self.features = None
        self.total_mi_reads = 0
        self.total_imi_reads = 0
        self.num_barcodes = 0

        self.barcodes = None
        self.genes = None
        self.bins = None
        self.counts = None
        self.ipm_bgcs = None

        # Expected number of molecules seen given N unique bins, used by the
        # IPM correction. Only valid for 5 <= n_bins <= 32; everything else
        # keeps the sentinel -1 and is filtered out downstream.
        self.expected_mol = {i: -1 for i in range(65)}
        expected = 0
        for n_bins in range(1, 33):
            expected += 64 / (64 - n_bins + 1)
            if n_bins >= 5:
                self.expected_mol[n_bins] = expected

    def read_dragen_moleculeInfo_h5(self, molecule_info_file, features_file=None):
        """Load a DRAGEN scRNA molecule-info h5 into memory (IMI mode only).

        Barcodes are always read directly from the h5. Features are read
        from the h5 when the file was produced by a pure-GEX run (which
        stores ``gene_ids`` / ``gene_names`` / ``gene_idx``). For combined
        GEX+CRISPR runs DRAGEN instead writes a single flat feature table
        under ``feature_ids`` / ``feature_names`` / ``feature_idx`` with no
        in-band feature-type column, and the caller must provide the
        DRAGEN ``features.tsv.gz`` file via ``features_file`` so we can
        filter down to just the "Gene Expression" rows.

        Only IMI-mode files (with ``binning_idx``) are supported. UMI-only
        outputs raise ``ValueError``.
        """
        if not os.path.isfile(molecule_info_file):
            raise FileNotFoundError(f'Molecule info file not found: {molecule_info_file}')

        with h5py.File(molecule_info_file, 'r') as h5:
            # Barcodes: DRAGEN stores the full whitelist as strings under
            # /barcodes, and per-molecule barcode indices under /barcode_idx.
            barcodes_seq = h5['barcodes'].asstr()[...]
            barcodes_idx = h5['barcode_idx'][...]

            # Binning must be present for IMI counting. If /binning_idx is
            # missing the file is UMI-mode, which we don't support.
            if 'binning_idx' not in h5:
                raise ValueError(
                    f'{molecule_info_file} has no /binning_idx dataset. This script only '
                    'supports IMI-counted DRAGEN molecule-info files.'
                )

            # Two possible feature-table layouts DRAGEN uses:
            #   1. Pure GEX runs: gene_ids / gene_names / gene_idx. All
            #      features are Gene Expression by construction, so no
            #      external features.tsv.gz is needed.
            #   2. Combined GEX+CRISPR (or CRISPR-only) runs: feature_ids /
            #      feature_names / feature_idx flatten every feature type
            #      into one table with no feature-type column. To pick out
            #      just the Gene Expression features we require the DRAGEN
            #      features.tsv.gz sidecar file.
            if 'gene_ids' in h5 and 'gene_idx' in h5:
                gene_ids = h5['gene_ids'].asstr()[...]
                gene_names = h5['gene_names'].asstr()[...] if 'gene_names' in h5 else gene_ids.copy()
                feature_idx_all = h5['gene_idx'][...]
                # Every entry in this h5 is GEX; no filtering needed.
                self.features = pd.DataFrame({
                    0: gene_ids,
                    1: gene_names,
                    2: np.full(len(gene_ids), 'Gene Expression'),
                })
                feature_idx_remapped = feature_idx_all
                molecule_mask = np.ones(len(feature_idx_all), dtype=bool)
            elif 'feature_ids' in h5 and 'feature_idx' in h5:
                if features_file is None:
                    raise ValueError(
                        f'{molecule_info_file} contains a combined feature table '
                        '(feature_ids / feature_idx) with no in-band feature-type '
                        'information. Pass features_file=... (the DRAGEN '
                        'scRNA.features.tsv.gz file) so GEX features can be '
                        'identified.'
                    )
                if not os.path.isfile(features_file):
                    raise FileNotFoundError(
                        f'Features TSV not found: {features_file}'
                    )

                feature_ids = h5['feature_ids'].asstr()[...]
                feature_names = h5['feature_names'].asstr()[...] if 'feature_names' in h5 else feature_ids.copy()
                feature_idx_all = h5['feature_idx'][...]

                # DRAGEN's features.tsv.gz has three tab-separated columns:
                # feature_id, feature_name, feature_type. We match by id
                # first, falling back to name if ids are not unique.
                tsv = pd.read_csv(features_file, header=None, sep='\t')
                tsv_ids = tsv.iloc[:, 0].astype(str).to_numpy()
                tsv_types = tsv.iloc[:, 2].astype(str).to_numpy()
                # Build a lookup from feature-id -> type. When there are
                # duplicate ids (rare but seen in some references) we fall
                # back to matching on name.
                if len(tsv_ids) == len(np.unique(tsv_ids)):
                    id_to_type = dict(zip(tsv_ids, tsv_types))
                    types_by_h5_row = np.array(
                        [id_to_type.get(fid, '') for fid in feature_ids], dtype=object
                    )
                else:
                    tsv_names = tsv.iloc[:, 1].astype(str).to_numpy()
                    name_to_type = dict(zip(tsv_names, tsv_types))
                    types_by_h5_row = np.array(
                        [name_to_type.get(fname, '') for fname in feature_names],
                        dtype=object,
                    )

                gex_mask = types_by_h5_row == 'Gene Expression'
                num_gex = int(np.sum(gex_mask))
                if num_gex == 0:
                    raise ValueError(
                        f'No "Gene Expression" features found in {features_file}. '
                        'Cannot downsample GEX-only from this run.'
                    )

                # Remap feature_idx: entries pointing to GEX features get
                # a new 0..num_gex-1 column position; entries pointing to
                # CRISPR (or other non-GEX) features are marked -1 and
                # dropped from the per-molecule arrays below.
                gex_positions = np.where(gex_mask)[0]
                # remap[i] = new column for original column i, or -1 if
                # original column i is not GEX.
                remap = np.full(len(feature_ids), -1, dtype=np.int64)
                remap[gex_positions] = np.arange(num_gex, dtype=np.int64)
                feature_idx_remapped = remap[feature_idx_all]
                molecule_mask = feature_idx_remapped != -1

                if not self.silent:
                    total = len(feature_ids)
                    print(
                        f'Feature table: {num_gex:,} GEX / {total - num_gex:,} non-GEX. '
                        f'Keeping {int(np.sum(molecule_mask)):,} of {len(feature_idx_all):,} '
                        'molecule entries (non-GEX dropped).'
                    )

                # self.features holds the GEX-only DataFrame so downstream
                # matrix code has the right column count.
                self.features = pd.DataFrame({
                    0: feature_ids[gex_positions],
                    1: feature_names[gex_positions],
                    2: np.full(num_gex, 'Gene Expression'),
                })
            else:
                raise ValueError(
                    f'{molecule_info_file} does not contain gene_ids/gene_idx nor '
                    'feature_ids/feature_idx datasets. This script only supports '
                    'DRAGEN scRNA molecule-info files.'
                )

            binning_idx_all = h5['binning_idx'][...]
            read_counts_all = h5['read_counts'][...]

        # Apply the (possibly all-True) molecule mask uniformly.
        self.barcodes_seq = np.asarray(barcodes_seq, dtype=str)
        self.barcodes = barcodes_idx[molecule_mask].astype(np.uint64)
        self.genes = feature_idx_remapped[molecule_mask].astype(np.uint32)
        self.bins = binning_idx_all[molecule_mask]
        self.counts = read_counts_all[molecule_mask]

        total_reads = int(np.sum(self.counts))
        self.total_mi_reads = total_reads
        self.total_imi_reads = total_reads
        self.num_barcodes = int(len(np.unique(self.barcodes)))

        # Precompute the per-(barcode, gene) mask of "eligible for IPM
        # calculation": at least 5 and at most 32 unique bins observed. The
        # pipeline uses this to average IPM values only over BGCs where the
        # bin coverage is dense enough to give a meaningful estimate.
        self.ipm_bgcs = _compute_ipm_bgcs(self.barcodes, self.genes, self.bins)

        if not self.silent:
            print(
                f'Loaded {len(self.barcodes):,} molecule entries across '
                f'{self.num_barcodes:,} barcodes and {len(self.features):,} features '
                f'({total_reads:,} mapped reads).'
            )

    def _resolve_cell_barcodes(self, cell_barcode_sequences):
        """Validate and convert the caller-provided cell-barcode iterable."""
        if cell_barcode_sequences is None:
            raise ValueError('Must provide filtered cell barcodes.')
        seqs = np.asarray(list(cell_barcode_sequences))
        if len(seqs) < 1:
            raise ValueError('Need at least one filtered cell barcode.')
        return seqs

    def _compute_downsample_factors(self, reads_per_cell_full, requested_rpcs, min_reads_per_cell):
        """Turn a set of requested reads-per-cell values into (target_rpc,
        downsample_factor) lists that are safe to feed into the workers.

        Values below ``min_reads_per_cell`` are dropped. Values that exceed
        the full-depth RPC would imply an upsample and are dropped as well.
        A requested depth exactly at full depth (ratio == 1) is kept as a
        no-op "downsample" -- this is the normal case for whichever sample
        defines the batch's target depth. The lists are returned in
        decreasing order of downsample factor (i.e. highest depth first) to
        match the parent pipeline's convention.
        """
        target_rpcs = []
        factors = []
        if reads_per_cell_full > 0:
            for rpc in requested_rpcs:
                if rpc < min_reads_per_cell:
                    continue
                ratio = rpc / reads_per_cell_full
                if ratio > 1:
                    continue
                target_rpcs.append(rpc)
                factors.append(ratio)
        # Preserve the deepest-first ordering the workers expect.
        pairs = sorted(zip(factors, target_rpcs), reverse=True)
        factors_sorted = [p[0] for p in pairs]
        rpcs_sorted = [p[1] for p in pairs]
        return factors_sorted, rpcs_sorted

    def run_downsample(
        self,
        cell_barcode_sequences,
        add_reads_per_cell=None,
        min_reads_per_cell=1,
        total_input_reads=None,
    ):
        """Compute IMI saturation metrics at a sweep of downsample depths.

        The sweep always includes the full-depth point plus a built-in ladder
        of reads-per-cell values (10 through 1,000,000). ``add_reads_per_cell``
        (a list of ints or a single int) can be used to add additional custom
        depths. Depths above the full-depth RPC are silently dropped.

        Returns a dict with per-depth arrays:
            - ``mean_reads_per_cell``: the target RPC (or full-depth value)
            - ``downsample_rate``: fraction of reads kept for that depth
            - ``pct_sequencing_saturation``: 100 * (1 - 1/dup_rate)
            - ``median_transcripts_in_cells``, ``total_transcripts_in_cells``
            - ``median_genes_in_cells``, ``pct_cells_with_features``
        """
        if not self.silent:
            print('Running downsample (IMI)')

        cell_barcode_sequences = self._resolve_cell_barcodes(cell_barcode_sequences)
        num_cells = len(cell_barcode_sequences)

        if total_input_reads is None:
            reads_per_cell_full = self.total_mi_reads / num_cells
        else:
            reads_per_cell_full = total_input_reads / num_cells

        # Built-in saturation depth ladder, plus optional caller-supplied extras.
        requested = {
            10, 100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 25000,
            30000, 40000, 50000, 75000, 100000, 125000, 150000, 200000,
            300000, 400000, 500000, 750000, 1000000,
        }
        if add_reads_per_cell is not None:
            if isinstance(add_reads_per_cell, int):
                requested.add(add_reads_per_cell)
            else:
                for rpc in add_reads_per_cell:
                    requested.add(int(rpc))

        downsample_factors, target_rpcs = self._compute_downsample_factors(
            reads_per_cell_full, requested, min_reads_per_cell
        )
        # Always emit the full-depth point at the top of the sweep.
        downsample_factors = [1.0] + downsample_factors
        target_rpcs = [reads_per_cell_full] + target_rpcs

        (
            median_transcripts,
            total_transcripts,
            median_genes,
            cells_with_features,
            _transcripts_df,
            _genes_df,
        ) = self.downsample_imi_counts(
            downsample_factors, cell_barcode_sequences=cell_barcode_sequences
        )

        # Sequencing saturation: 1 - 1/duplication_rate. duplication_rate is
        # (reads in cells at this depth) / (unique transcripts in cells at
        # this depth). Clamped at 0 so the CSV never shows a negative
        # saturation for very shallow depths.
        barcode_seq_inds = np.isin(self.barcodes_seq, cell_barcode_sequences)
        cell_barcode_seq_inds = np.where(barcode_seq_inds)[0]
        cell_barcode_inds_mask = np.isin(self.barcodes, cell_barcode_seq_inds)
        mapped_reads_in_cells = np.sum(self.counts[cell_barcode_inds_mask])
        downsampled_mapped_reads_in_cells = np.array(
            [mapped_reads_in_cells * d for d in downsample_factors]
        )
        dup_rate = downsampled_mapped_reads_in_cells / total_transcripts
        seq_sat = [max((1 - 1 / d) * 100, 0) for d in dup_rate]

        pct_cells_with_features = cells_with_features * 100 / num_cells

        return {
            'mean_reads_per_cell': target_rpcs,
            'downsample_rate': downsample_factors,
            'pct_sequencing_saturation': seq_sat,
            'median_transcripts_in_cells': median_transcripts,
            'total_transcripts_in_cells': total_transcripts,
            'median_genes_in_cells': median_genes,
            'pct_cells_with_features': pct_cells_with_features,
        }

    def run_downsample_matrix(
        self,
        reads_per_cell,
        cell_barcode_sequences,
        min_reads_per_cell=1,
        total_input_reads=None,
    ):
        """Return one ``MatrixData`` per requested downsampled depth.

        ``reads_per_cell`` is a list of target RPCs (a single int is also
        accepted). Depths above the full-depth RPC are silently dropped and
        the returned list is ordered deepest-first, matching the built-in
        depth ordering used elsewhere in the pipeline. Cells are needed only
        to compute the full-depth RPC baseline; the produced matrices are
        the *raw* downsampled matrices (all barcodes), not filtered.
        """
        if not self.silent:
            print('Running downsample matrix (IMI)')

        cell_barcode_sequences = self._resolve_cell_barcodes(cell_barcode_sequences)
        num_cells = len(cell_barcode_sequences)

        if isinstance(reads_per_cell, int):
            reads_per_cell = [reads_per_cell]
        requested = {int(rpc) for rpc in reads_per_cell}

        if total_input_reads is None:
            reads_per_cell_full = self.total_mi_reads / num_cells
        else:
            reads_per_cell_full = total_input_reads / num_cells

        downsample_factors, target_rpcs = self._compute_downsample_factors(
            reads_per_cell_full, requested, min_reads_per_cell
        )
        if not downsample_factors:
            if not self.silent:
                print('No requested matrix depths remain after filtering.')
            return [], []
        matrices = self.get_imi_matrix(downsample_factors=downsample_factors)
        return matrices, target_rpcs

    def get_ipms(self, downsample_factors, rand_block=100_000):
        """Estimate the global per-molecule IPM constant at each downsample factor.

        The IPM (Inferred Per Molecule) correction is applied when converting
        binning-index observations back into transcript counts: for each
        (barcode, gene) with enough bin coverage, we sum the observed count
        vs. an expected-count function of the number of unique bins seen.
        Averaging that ratio across all eligible BGCs gives the global IPM
        for a given depth. The returned array is one IPM value per entry in
        ``downsample_factors``, floored at 1 (the correction is a no-op when
        IPM <= 1).
        """
        ipm_totals = np.zeros(len(downsample_factors), dtype='uint64')
        ipm_counts = np.zeros(len(downsample_factors), dtype='int32')

        args = []
        start_pos = 0
        new_barcodes = np.where(
            np.append(np.append(True, np.diff(self.barcodes) != 0), True)
        )[0]
        num_barcodes = self.num_barcodes
        # Worker index used to seed each parallel process deterministically.
        worker_idx = 0

        for bc_idx in range(num_barcodes):
            end_pos = new_barcodes[bc_idx + 1]

            if (end_pos - start_pos >= self.block_size) or (bc_idx == num_barcodes - 1):
                # Reached the size of one thread block; extract the slice and
                # queue it for a worker.
                ipm_bgcs = self.ipm_bgcs[start_pos:end_pos]

                thread_bcs = self.barcodes[start_pos:end_pos][ipm_bgcs]
                thread_genes = self.genes[start_pos:end_pos][ipm_bgcs]
                thread_counts = self.counts[start_pos:end_pos][ipm_bgcs]
                thread_bis = self.bins[start_pos:end_pos][ipm_bgcs]
                start_pos = end_pos

                thread_args = {
                    'parallel': self.parallel_pool is not None,
                    'barcodes': thread_bcs,
                    'genes': thread_genes,
                    'counts': thread_counts,
                    'bin_indexes': thread_bis,
                    'rand_block': rand_block,
                    'expected_mol': self.expected_mol,
                    'downsample_factors': downsample_factors,
                    'random_seed': self.random_seed,
                    'worker_idx': worker_idx,
                }
                args.append(thread_args)
                worker_idx += 1

            if (len(args) == self.num_cores) or (bc_idx == num_barcodes - 1):
                if self.parallel_pool is None:
                    results = [get_ipm_counts_parallel(a) for a in args]
                else:
                    try:
                        results = self.parallel_pool.map(get_ipm_counts_parallel, args)
                    except Exception:
                        self.parallel_pool.close()
                        self.parallel_pool.join()
                        raise
                for thread_totals, thread_counts_arr in results:
                    ipm_totals = ipm_totals + thread_totals
                    ipm_counts = ipm_counts + thread_counts_arr
                args = []

        ipm_counts = np.where(ipm_counts == 0, -1, ipm_counts)
        ipms = ipm_totals / ipm_counts
        ipms = np.where(ipms < 1, 1, ipms)
        return ipms

    def downsample_imi_counts(self, downsample_factors, cell_barcode_sequences, rand_block=100_000):
        """Compute per-cell transcript / gene counts at each downsample factor.

        Returns a tuple of ``(median_transcripts, total_transcripts,
        median_genes, cells_with_features, transcripts_df, genes_df)`` where
        each ``*_df`` is a dict keyed by downsample factor with the per-cell
        arrays. The DataFrames are kept as dicts (rather than pandas frames)
        so the workers stay pandas-free.
        """
        ipms = self.get_ipms(downsample_factors)

        num_barcodes = self.num_barcodes
        barcode_seq_inds = np.isin(self.barcodes_seq, cell_barcode_sequences)
        cell_barcode_seq_inds = np.where(barcode_seq_inds)[0]
        num_cbs = len(cell_barcode_seq_inds)
        transcript_counts = np.tile(np.zeros(num_cbs), (len(downsample_factors), 1))
        gene_counts = np.tile(np.zeros(num_cbs), (len(downsample_factors), 1))

        args = []
        start_pos = 0
        new_barcodes = np.where(
            np.append(np.append(True, np.diff(self.barcodes) != 0), True)
        )[0]
        worker_idx = 0

        for bc_idx in range(num_barcodes):
            end_pos = new_barcodes[bc_idx + 1]

            if (end_pos - start_pos >= self.block_size) or (bc_idx == num_barcodes - 1):
                thread_bcs = self.barcodes[start_pos:end_pos]
                thread_barcodes_mask = np.isin(thread_bcs, cell_barcode_seq_inds)
                thread_bcs_mask = np.where(thread_barcodes_mask)[0]
                thread_bcs = thread_bcs[thread_bcs_mask]
                thread_genes = self.genes[start_pos:end_pos][thread_bcs_mask]
                thread_counts = self.counts[start_pos:end_pos][thread_bcs_mask]
                thread_bis = self.bins[start_pos:end_pos][thread_bcs_mask]
                start_pos = end_pos

                thread_args = {
                    'barcodes': thread_bcs,
                    'genes': thread_genes,
                    'counts': thread_counts,
                    'bin_indexes': thread_bis,
                    'rand_block': rand_block,
                    'ipms': ipms,
                    'downsample_factors': downsample_factors,
                    'random_seed': self.random_seed,
                    'worker_idx': worker_idx,
                }
                args.append(thread_args)
                worker_idx += 1

            if (len(args) == self.num_cores) or (bc_idx == num_barcodes - 1):
                if self.parallel_pool is None:
                    results = [get_imi_counts_parallel(a) for a in args]
                else:
                    try:
                        results = self.parallel_pool.map(get_imi_counts_parallel, args)
                    except Exception:
                        self.parallel_pool.close()
                        self.parallel_pool.join()
                        raise
                for thread_tc, thread_gc, thread_barcode_inds in results:
                    ind_mask = np.isin(cell_barcode_seq_inds, thread_barcode_inds)
                    array_inds = np.where(ind_mask)[0]
                    transcript_counts[:, array_inds] += thread_tc
                    gene_counts[:, array_inds] += thread_gc
                args = []

        transcripts_df = {
            downsample_factors[d]: transcript_counts[d, :]
            for d in range(len(downsample_factors))
        }
        transcripts_df['Barcode'] = self.barcodes_seq[cell_barcode_seq_inds]
        genes_df = {
            downsample_factors[d]: gene_counts[d, :]
            for d in range(len(downsample_factors))
        }
        genes_df['Barcode'] = self.barcodes_seq[cell_barcode_seq_inds]

        median_transcripts = np.median(transcript_counts, axis=1)
        total_transcripts = np.sum(transcript_counts, axis=1)
        median_genes = np.median(gene_counts, axis=1)
        cells_with_features = np.sum(gene_counts > 0, axis=1)

        return (
            median_transcripts,
            total_transcripts,
            median_genes,
            cells_with_features,
            transcripts_df,
            genes_df,
        )

    def get_imi_matrix(self, downsample_factors=(1,), rand_block=100_000):
        """Produce one ``MatrixData`` per downsample factor.

        For each entry in ``downsample_factors``, this returns a raw
        (unfiltered) barcodes x features counts matrix representing the
        downsampled data. The IPM correction is applied per-BGC to convert
        binning observations back into transcript counts. The returned list
        is in the same order as ``downsample_factors``.
        """
        downsample_factors = list(downsample_factors)
        args = []
        start_pos = 0

        max_nonzero = len(self.counts)
        barcodes_col = np.zeros(max_nonzero, dtype='uint64')
        genes_col = np.zeros(max_nonzero, dtype='uint32')
        unique_bins = np.zeros((len(downsample_factors), max_nonzero), dtype='uint8')
        counts = np.zeros((len(downsample_factors), max_nonzero), dtype='uint32')
        thread_start_idx = 0

        new_barcodes = np.where(
            np.append(np.append(True, np.diff(self.barcodes) != 0), True)
        )[0]
        num_barcodes = self.num_barcodes
        worker_idx = 0
        thread_end_idx = 0

        for bc_idx in range(num_barcodes):
            end_pos = new_barcodes[bc_idx + 1]

            if (end_pos - start_pos >= self.block_size) or (bc_idx == num_barcodes - 1):
                thread_bcs = self.barcodes[start_pos:end_pos]
                thread_genes = self.genes[start_pos:end_pos]
                thread_counts = self.counts[start_pos:end_pos]
                thread_bis = self.bins[start_pos:end_pos]
                start_pos = end_pos

                thread_args = {
                    'barcodes': thread_bcs,
                    'genes': thread_genes,
                    'counts': thread_counts,
                    'bin_indexes': thread_bis,
                    'rand_block': rand_block,
                    'downsample_factors': downsample_factors,
                    'random_seed': self.random_seed,
                    'worker_idx': worker_idx,
                }
                args.append(thread_args)
                worker_idx += 1

            if (len(args) == self.num_cores) or (bc_idx == num_barcodes - 1):
                if self.parallel_pool is None:
                    results = [get_imi_matrix_parallel(a) for a in args]
                else:
                    try:
                        results = self.parallel_pool.map(get_imi_matrix_parallel, args)
                    except Exception:
                        self.parallel_pool.close()
                        self.parallel_pool.join()
                        raise
                for result_barcodes, result_genes, result_read_counts, result_unique_bins in results:
                    thread_end_idx = thread_start_idx + len(result_barcodes)
                    barcodes_col[thread_start_idx:thread_end_idx] = result_barcodes
                    genes_col[thread_start_idx:thread_end_idx] = result_genes
                    counts[:, thread_start_idx:thread_end_idx] = result_read_counts
                    unique_bins[:, thread_start_idx:thread_end_idx] = result_unique_bins
                    thread_start_idx = thread_end_idx
                args = []

        barcodes_col = barcodes_col[:thread_end_idx]
        genes_col = genes_col[:thread_end_idx]
        counts = counts[:, :thread_end_idx]
        unique_bins = unique_bins[:, :thread_end_idx]

        # Global IPM per downsample factor.
        ipms = np.ones(len(downsample_factors))
        for i in range(len(downsample_factors)):
            ipm_idxs = np.where(
                np.logical_and(unique_bins[i, :] >= 5, unique_bins[i, :] <= 32)
            )[0]
            expected = np.array([self.expected_mol[b] for b in unique_bins[i, ipm_idxs]])
            bgc_ipms = counts[i, ipm_idxs] / expected
            if len(bgc_ipms) > 0:
                ipms[i] = max(1, sum(bgc_ipms) / len(bgc_ipms))
            else:
                ipms[i] = 1

        # Apply the IPM correction. For BGCs with fewer than 10 unique bins
        # we fall back to using the raw unique-bin count directly, since the
        # IPM correction is unreliable at low coverage. Corrected values are
        # also floored at the unique-bin count so a single low-coverage BGC
        # can never produce a count of 0.
        for i in range(len(downsample_factors)):
            corrected = np.floor_divide(counts[i, :], ipms[i])
            n_bins_row = unique_bins[i, :]
            corrected = np.where(n_bins_row > 10, corrected, n_bins_row)
            corrected = np.maximum(corrected, n_bins_row)
            counts[i, :] = corrected

        matrices = []
        for i in range(len(downsample_factors)):
            matrices.append(
                MatrixData(
                    barcodes_col=barcodes_col,
                    features_col=genes_col,
                    counts_col=counts[i],
                    barcodes_seq=self.barcodes_seq,
                    features=self.features,
                )
            )
        return matrices


def _compute_ipm_bgcs(barcodes, genes, bins):
    """Vectorized replacement for the per-(barcode, gene) IPM-eligibility mask.

    Returns a boolean array with the same length as ``barcodes`` / ``genes`` /
    ``bins`` where each row is ``True`` iff its BGC (contiguous run of same
    ``barcode`` and same ``gene``) observed between 5 and 32 unique
    ``bins`` values inclusive.

    Assumes ``barcodes`` is non-decreasing and, within each barcode, ``genes``
    is non-decreasing (which is DRAGEN's storage order). ``bins`` must be an
    integer array with values in ``[0, 63]``.

    The implementation is O(n) in time and O(num_bgcs) in extra memory: a
    per-BGC 64-bit bitmask records which of the 64 possible bin values each
    BGC has seen, and the final unique-bin count is a popcount over that
    bitmask. This replaces a Python while-loop that scaled poorly for the
    hundreds-of-millions-of-molecules h5 files DRAGEN emits on large runs.
    """
    n = len(barcodes)
    if n == 0:
        return np.zeros(0, dtype=bool)

    # Identify BGC boundaries: a row starts a new BGC iff its (barcode, gene)
    # pair differs from the previous row.
    new_bgc_flags = np.empty(n, dtype=bool)
    new_bgc_flags[0] = True
    np.logical_or(
        np.diff(barcodes) != 0,
        np.diff(genes) != 0,
        out=new_bgc_flags[1:],
    )
    # bgc_id[i] gives the 0-based BGC index of row i.
    bgc_id = np.cumsum(new_bgc_flags, dtype=np.int64) - 1
    num_bgcs = int(bgc_id[-1]) + 1

    # OR together (1 << bin) into a per-BGC uint64 bitmask so we can popcount
    # unique bin values without touching Python.
    bit_flags = (np.uint64(1) << bins.astype(np.uint64))
    bgc_bitmask = np.zeros(num_bgcs, dtype=np.uint64)
    np.bitwise_or.at(bgc_bitmask, bgc_id, bit_flags)

    if hasattr(np, 'bitwise_count'):
        # numpy 2.0+ - straight popcount over the uint64 mask.
        n_unique_per_bgc = np.bitwise_count(bgc_bitmask).astype(np.int64)
    else:
        # Fallback for older numpy: unpackbits interprets each byte as 8 bits
        # and sums per row.
        byte_view = bgc_bitmask.view(np.uint8).reshape(num_bgcs, 8)
        n_unique_per_bgc = np.unpackbits(byte_view, axis=1).sum(axis=1).astype(np.int64)

    # Eligibility mask per BGC, then scatter back to per-row.
    bgc_eligible = (n_unique_per_bgc >= 5) & (n_unique_per_bgc <= 32)
    return bgc_eligible[bgc_id]


def _seed_worker(random_seed, worker_idx):
    """Seed ``np.random`` deterministically for a given worker.

    Called at the top of every parallel worker so that repeated runs of the
    same script with the same ``--random-seed`` produce identical
    downsampled counts. Without this, workers inherit whatever np.random
    state the child process was spawned with, which is not deterministic
    across processes.
    """
    np.random.seed(int(random_seed) + int(worker_idx))


def get_ipm_counts_parallel(args):
    """Worker: accumulate IPM totals/counts for one block of BGC data.

    Runs the downsample simulation for every requested downsample factor
    over the block's (barcode, gene, count, bin) arrays, and returns per-
    factor totals of (a) sum of BGC IPM values and (b) number of BGCs
    that contributed to that sum. The caller averages across blocks.
    """
    _seed_worker(args['random_seed'], args['worker_idx'])
    barcodes = args['barcodes']
    genes = args['genes']
    counts = args['counts']
    bin_indexes = args['bin_indexes']
    expected_mol = args['expected_mol']
    downsample_factors = args['downsample_factors']

    ipm_totals = np.zeros(len(downsample_factors), dtype='uint64')
    ipm_counts = np.zeros(len(downsample_factors), dtype='int32')

    num_elements = len(counts)
    if num_elements == 0:
        return ipm_totals, ipm_counts

    same_barcode_and_gene = (
        [False]
        + np.logical_and(np.diff(barcodes) == 0, np.diff(genes) == 0).tolist()
        + [False]
    )
    idx = 0
    rand_cursor = 0
    rand_selection = np.zeros((1, 0))
    num_results = len(downsample_factors)

    while idx < num_elements:
        same_idx = idx + 1
        while same_barcode_and_gene[same_idx]:
            same_idx += 1

        counts_slice = counts[idx:same_idx]
        bin_index_slice = bin_indexes[idx:same_idx]
        rand_nums_needed = int(sum(counts_slice))

        if rand_cursor + rand_nums_needed >= rand_selection.shape[1] - 1:
            rand_nums = np.random.rand(max(args['rand_block'], rand_nums_needed) + 1)
            rand_selection = []
            for ds in range(num_results):
                rand_selection.append(
                    np.cumsum(rand_nums < downsample_factors[ds], dtype='uint32')
                )
                rand_selection[-1][-1] = 0
            rand_selection = np.array(rand_selection)
            rand_cursor = -1

        seen_bins = np.zeros((num_results, 64))
        imi_counts = np.zeros(num_results)
        for i in range(len(counts_slice)):
            curr_count = int(counts_slice[i])
            mol_counts = (
                rand_selection[:, rand_cursor + curr_count]
                - rand_selection[:, rand_cursor]
            )
            curr_imi_counts = mol_counts > 0
            imi_counts += curr_imi_counts
            seen_bins[curr_imi_counts > 0, bin_index_slice[i]] = 1
            rand_cursor += curr_count

        num_unique_bins = np.sum(seen_bins, axis=1)
        bgc_ipms = imi_counts / np.array([expected_mol[n] for n in num_unique_bins])
        bgc_ipms = np.where(bgc_ipms > 0, bgc_ipms, 0)
        bgc_counts = np.where(bgc_ipms > 0, 1, 0)
        ipm_totals = ipm_totals + bgc_ipms
        ipm_counts = ipm_counts + bgc_counts

        idx = same_idx
    return ipm_totals, ipm_counts


def get_imi_counts_parallel(args):
    """Worker: compute per-cell transcript / gene counts for one block.

    Applies the IPM correction using the globally averaged IPM values
    passed in via ``args['ipms']``. Returns per-cell arrays with one row
    per downsample factor and one column per unique cell barcode seen in
    the block, plus a parallel array giving those cell barcode indices so
    the caller can scatter results into the global cell-barcode layout.
    """
    _seed_worker(args['random_seed'], args['worker_idx'])
    barcodes = args['barcodes']
    genes = args['genes']
    counts = args['counts']
    bin_indexes = args['bin_indexes']
    ipms = args['ipms']
    downsample_factors = args['downsample_factors']
    rand_block = args['rand_block']

    num_barcodes = len(np.unique(barcodes))
    transcript_counts = np.tile(np.zeros(num_barcodes), (len(downsample_factors), 1))
    gene_counts = np.tile(np.zeros(num_barcodes), (len(downsample_factors), 1))
    barcode_inds = np.zeros(num_barcodes, dtype=barcodes.dtype)

    num_elements = len(counts)
    if num_elements == 0:
        return transcript_counts, gene_counts, barcode_inds

    same_barcode_and_gene = (
        [False]
        + np.logical_and(np.diff(barcodes) == 0, np.diff(genes) == 0).tolist()
        + [False]
    )
    idx = 0
    curr_barcode = barcodes[0]
    curr_barcode_idx = 0

    rand_cursor = 0
    rand_selection = np.zeros((1, 0))
    num_results = len(downsample_factors)

    while idx < num_elements:
        same_idx = idx + 1
        while same_barcode_and_gene[same_idx]:
            same_idx += 1

        counts_slice = counts[idx:same_idx]
        bin_index_slice = bin_indexes[idx:same_idx]
        rand_nums_needed = int(sum(counts_slice))

        if rand_cursor + rand_nums_needed >= rand_selection.shape[1] - 1:
            rand_nums = np.random.rand(max(rand_block, rand_nums_needed) + 1)
            rand_selection = []
            for ds in range(num_results):
                rand_selection.append(
                    np.cumsum(rand_nums < downsample_factors[ds], dtype='uint32')
                )
                rand_selection[-1][-1] = 0
            rand_selection = np.array(rand_selection)
            rand_cursor = -1

        if (same_idx - idx) == 1:
            curr_count = int(counts_slice[0])
            mol_counts = (
                rand_selection[:, rand_cursor + curr_count]
                - rand_selection[:, rand_cursor]
            )
            corrected_imi_counts = mol_counts > 0
            transcript_counts[:, curr_barcode_idx] += corrected_imi_counts
            gene_counts[:, curr_barcode_idx] += corrected_imi_counts
        else:
            seen_bins = np.zeros((num_results, 64))
            imi_counts = np.zeros(num_results)
            for i in range(len(counts_slice)):
                curr_count = int(counts_slice[i])
                mol_counts = (
                    rand_selection[:, rand_cursor + curr_count]
                    - rand_selection[:, rand_cursor]
                )
                curr_imi_counts = mol_counts > 0
                imi_counts += curr_imi_counts
                seen_bins[curr_imi_counts, bin_index_slice[i]] = 1
                rand_cursor += curr_count

            num_unique_bins = np.sum(seen_bins, axis=1)
            corrected_imi_counts = np.floor_divide(imi_counts, ipms)
            corrected_imi_counts = np.where(
                num_unique_bins > 10, corrected_imi_counts, num_unique_bins
            )
            corrected_imi_counts = np.maximum(corrected_imi_counts, num_unique_bins)
            gene_present = corrected_imi_counts > 0

            transcript_counts[:, curr_barcode_idx] += corrected_imi_counts
            gene_counts[:, curr_barcode_idx] += gene_present

        if barcodes[idx] != curr_barcode:
            barcode_inds[curr_barcode_idx] = curr_barcode
            curr_barcode = barcodes[idx]
            curr_barcode_idx += 1
        idx = same_idx
    barcode_inds[curr_barcode_idx] = curr_barcode
    return transcript_counts, gene_counts, barcode_inds


def get_imi_matrix_parallel(args):
    """Worker: build the downsampled matrix (as COO-style arrays) for a block.

    Emits, for each BGC in the block, one entry per downsample factor with
    the corrected count and the number of unique bins seen. The caller
    concatenates blocks and applies the global IPM correction.
    """
    _seed_worker(args['random_seed'], args['worker_idx'])
    barcodes = args['barcodes']
    genes = args['genes']
    counts = args['counts']
    bin_indexes = args['bin_indexes']
    downsample_factors = args['downsample_factors']
    rand_block = args['rand_block']

    max_nonzero = len(counts)
    new_counts = np.tile(np.zeros(max_nonzero), (len(downsample_factors), 1))
    unique_bins = np.tile(np.zeros(max_nonzero), (len(downsample_factors), 1))
    new_barcodes = np.zeros(max_nonzero)
    new_genes = np.zeros(max_nonzero)
    new_count_idx = 0

    same_barcode_and_gene = (
        [False]
        + np.logical_and(np.diff(barcodes) == 0, np.diff(genes) == 0).tolist()
        + [False]
    )
    idx = 0
    num_elements = len(counts)

    rand_cursor = 0
    rand_selection = np.zeros((1, 0))
    num_results = len(downsample_factors)

    while idx < num_elements:
        same_idx = idx + 1
        while same_barcode_and_gene[same_idx]:
            same_idx += 1

        barcode = barcodes[idx]
        gene = genes[idx]
        counts_slice = counts[idx:same_idx]
        bin_index_slice = bin_indexes[idx:same_idx]
        rand_nums_needed = int(sum(counts_slice))

        if rand_cursor + rand_nums_needed >= rand_selection.shape[1] - 1:
            rand_nums = np.random.rand(max(rand_block, rand_nums_needed) + 1)
            rand_selection = []
            for ds in range(num_results):
                rand_selection.append(
                    np.cumsum(rand_nums < downsample_factors[ds], dtype='uint32')
                )
                rand_selection[-1][-1] = 0
            rand_selection = np.array(rand_selection)
            rand_cursor = -1

        if (same_idx - idx) == 1:
            curr_count = int(counts_slice[0])
            mol_counts = (
                rand_selection[:, rand_cursor + curr_count]
                - rand_selection[:, rand_cursor]
            )
            n_bins_row = mol_counts > 0
            new_counts[:, new_count_idx] = mol_counts
            unique_bins[:, new_count_idx] = n_bins_row
            new_barcodes[new_count_idx] = barcode
            new_genes[new_count_idx] = gene
        else:
            seen_bins = np.zeros((num_results, 64))
            imi_counts = np.zeros(num_results)
            for i in range(len(counts_slice)):
                curr_count = int(counts_slice[i])
                mol_counts = (
                    rand_selection[:, rand_cursor + curr_count]
                    - rand_selection[:, rand_cursor]
                )
                curr_imi_counts = mol_counts > 0
                imi_counts += curr_imi_counts
                seen_bins[curr_imi_counts > 0, bin_index_slice[i]] = 1
                rand_cursor += curr_count

            n_bins_row = np.sum(seen_bins, axis=1)
            new_counts[:, new_count_idx] = imi_counts
            unique_bins[:, new_count_idx] = n_bins_row
            new_barcodes[new_count_idx] = barcode
            new_genes[new_count_idx] = gene

        idx = same_idx
        new_count_idx += 1

    return (
        new_barcodes[:new_count_idx],
        new_genes[:new_count_idx],
        new_counts[:, :new_count_idx],
        unique_bins[:, :new_count_idx],
    )


# =====================================================================
# CLI
# =====================================================================


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            'Downsample IMI counts from a DRAGEN scRNA molecule-info h5. '
            'Saves downsampled matrices and/or a saturation sweep CSV.'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input selection: users typically point at a full DRAGEN results
    # directory and let the script auto-detect the h5, metrics CSV, and
    # filtered-barcodes file. Direct-file overrides exist for edge cases.
    parser.add_argument(
        '--dragen-results-dir',
        default=None,
        help='Path to a DRAGEN scRNA results directory. Used to auto-detect the '
             'molecule-info h5, filtered-barcodes list, and scRNA_metrics.csv.',
    )
    parser.add_argument(
        '--molecule-info-h5',
        default=None,
        help='Path to a DRAGEN scRNA moleculeInfo.h5. Overrides auto-detection '
             'from --dragen-results-dir.',
    )
    parser.add_argument(
        '--filtered-barcodes',
        default=None,
        help='Path to a filtered-barcodes list (one barcode per line, plain text '
             'or .gz). Defaults to the DRAGEN filtered-barcodes file in '
             '--dragen-results-dir.',
    )
    parser.add_argument(
        '--scrna-metrics-csv',
        default=None,
        help='Path to a DRAGEN scRNA_metrics.csv file used to look up '
             '"Total input reads". Defaults to the auto-detected metrics file '
             'in --dragen-results-dir.',
    )
    parser.add_argument(
        '--features-tsv',
        default=None,
        help='Path to the DRAGEN raw scRNA.features.tsv.gz file. Required for '
             'combined GEX+CRISPR runs so GEX features can be identified; '
             'ignored for pure-GEX runs (feature metadata is read from the h5 '
             'directly). Auto-detected from --dragen-results-dir when possible.',
    )
    parser.add_argument(
        '--total-input-reads',
        type=int,
        default=None,
        help='Total input reads to use when computing reads-per-cell. Overrides '
             'the value parsed from --scrna-metrics-csv. If neither is provided, '
             'the mapped-read total from the h5 is used as a lower-bound fallback.',
    )

    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to write outputs into. Created if it does not exist.',
    )

    # Depth selection.
    parser.add_argument(
        '--matrix-depths',
        type=int,
        nargs='+',
        default=None,
        help='Reads-per-cell depths at which to save downsampled matrices. '
             'One <N>rpc/ folder is created per depth. If omitted, no matrices '
             'are written and --saturation must be provided.',
    )
    parser.add_argument(
        '--matrix-format',
        choices=('raw', 'filtered', 'both'),
        default='both',
        help='Which layout(s) to write per matrix depth. "raw" writes the full '
             'unfiltered matrix, "filtered" writes only the cells in the '
             'filtered-barcodes list, "both" writes both under raw_matrix/ and '
             'filtered_matrix/ subdirs.',
    )
    parser.add_argument(
        '--saturation',
        action='store_true',
        help='Also run a saturation sweep and write saturation.csv with '
             'per-depth median transcripts, median genes, and % sequencing '
             'saturation.',
    )
    parser.add_argument(
        '--saturation-extra-depths',
        type=int,
        nargs='+',
        default=None,
        help='Extra reads-per-cell depths (in addition to the built-in ladder) '
             'to include in the saturation sweep. Ignored when --saturation is '
             'not set.',
    )

    parser.add_argument(
        '--min-reads-per-cell',
        type=int,
        default=1,
        help='Drop requested depths below this floor. Matches the parent '
             'pipeline default.',
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of worker processes for the downsampling loops. Defaults '
             'to 1 (single-threaded).',
    )
    parser.add_argument(
        '--random-seed',
        type=int,
        default=42,
        help='Base random seed. Each worker seeds np.random with '
             'random_seed + worker_index so results are reproducible across runs '
             'that use the same --threads value.',
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Silence per-step progress logging.',
    )

    args = parser.parse_args(argv)

    if not args.matrix_depths and not args.saturation:
        parser.error(
            'Nothing to do. Pass --matrix-depths and/or --saturation.'
        )
    return args


def resolve_inputs(args):
    """Turn the CLI args into concrete file paths, using --dragen-results-dir
    for auto-detection when explicit paths were not provided."""
    results_dir = args.dragen_results_dir
    molecule_info_h5 = args.molecule_info_h5
    filtered_barcodes = args.filtered_barcodes
    scrna_metrics_csv = args.scrna_metrics_csv
    features_tsv = args.features_tsv

    if molecule_info_h5 is None:
        if results_dir is None:
            raise SystemExit(
                'Must provide either --molecule-info-h5 or --dragen-results-dir.'
            )
        molecule_info_h5 = find_molecule_info_h5(results_dir)
        print(f'Auto-detected molecule-info h5: {molecule_info_h5}')

    # For any auto-detected companion files we only fall back when the user
    # gave us a results directory to search. Missing auto-detected files are
    # not fatal unless we actually need them later.
    if filtered_barcodes is None and results_dir is not None:
        auto = find_filtered_barcodes(results_dir)
        if auto is not None:
            filtered_barcodes = auto
            print(f'Auto-detected filtered-barcodes file: {filtered_barcodes}')

    if scrna_metrics_csv is None and results_dir is not None:
        auto = find_scrna_metrics_csv(results_dir)
        if auto is not None:
            scrna_metrics_csv = auto
            print(f'Auto-detected scRNA metrics CSV: {scrna_metrics_csv}')

    if features_tsv is None and results_dir is not None:
        auto = find_features_tsv(results_dir)
        if auto is not None:
            features_tsv = auto
            print(f'Auto-detected features TSV: {features_tsv}')

    return molecule_info_h5, filtered_barcodes, scrna_metrics_csv, features_tsv


def resolve_total_input_reads(args, scrna_metrics_csv, h5_total_mi_reads):
    """Determine the "Total input reads" value used to compute reads-per-cell.

    Precedence:
      1. --total-input-reads (explicit flag)
      2. Parsed from scRNA_metrics.csv when available
      3. h5 total_mi_reads (with a warning; this is mapped-reads-only and
         under-reports the true input by the unmapped-reads fraction)
    """
    if args.total_input_reads is not None:
        print(f'Using --total-input-reads override: {args.total_input_reads:,}')
        return args.total_input_reads
    if scrna_metrics_csv is not None:
        value = parse_total_input_reads(scrna_metrics_csv)
        if value is not None:
            print(
                f'Using Total input reads from {os.path.basename(scrna_metrics_csv)}: '
                f'{value:,}'
            )
            return value
        print(
            f'WARNING: {scrna_metrics_csv} does not contain a "Total input reads" row. '
            f'Falling back to h5 mapped-read total.'
        )
    else:
        print(
            'WARNING: no --total-input-reads or scRNA_metrics.csv available. '
            'Falling back to h5 mapped-read total (under-reports true RPC).'
        )
    print(f'Using h5 mapped-read total: {h5_total_mi_reads:,}')
    return h5_total_mi_reads


def resolve_filtered_barcodes(filtered_barcodes_path):
    """Read the filtered-barcodes file into a numpy string array."""
    if filtered_barcodes_path is None:
        raise SystemExit(
            'No filtered-barcodes file available. Pass --filtered-barcodes '
            'explicitly, or point --dragen-results-dir at a directory that '
            'contains scRNA.filtered.barcodes.tsv.gz.'
        )
    seqs = np.asarray(smart_read_textfile(filtered_barcodes_path))
    print(f'Loaded {len(seqs):,} filtered-cell barcodes from {filtered_barcodes_path}')
    return seqs


def write_saturation_csv(results, output_dir):
    """Write the saturation sweep CSV. Matches the columns used by the parent
    pipeline's ``data.csv`` output so downstream tooling remains compatible."""
    df = pd.DataFrame({
        'reads_per_cell': np.round(results['mean_reads_per_cell']).astype(int),
        'median_transcripts_per_cell': np.round(
            results['median_transcripts_in_cells']
        ).astype(int),
        'median_genes_per_cell': np.round(results['median_genes_in_cells']).astype(int),
        'sequencing_saturation': np.round(results['pct_sequencing_saturation'], 2),
    })
    out_path = os.path.join(output_dir, 'saturation.csv')
    df.to_csv(out_path, index=False)
    print(f'Saturation sweep written to {out_path}')


def write_matrices(minf, target_rpcs, matrices, filtered_barcodes_seq, output_dir, matrix_format):
    """Write raw and/or filtered matrix layouts for each downsample depth."""
    if not matrices:
        print('No matrix depths remained after filtering; nothing to write.')
        return

    base_dir = os.path.join(output_dir, 'downsampled_matrix')
    for depth_rpc, mtx in zip(target_rpcs, matrices):
        depth_int = int(round(depth_rpc))
        depth_dir = os.path.join(base_dir, f'{depth_int}rpc')
        if matrix_format in ('raw', 'both'):
            raw_dir = os.path.join(depth_dir, 'raw_matrix')
            mtx.write_matrix(output_dir=raw_dir)
            print(f'  wrote raw matrix ({depth_int} rpc) -> {raw_dir}')
        if matrix_format in ('filtered', 'both'):
            mtx.get_filtered_cells(cell_barcodes_seq=filtered_barcodes_seq)
            filtered_mtx = mtx.get_filtered_matrix()
            filt_dir = os.path.join(depth_dir, 'filtered_matrix')
            filtered_mtx.write_matrix(output_dir=filt_dir)
            print(f'  wrote filtered matrix ({depth_int} rpc) -> {filt_dir}')


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)

    molecule_info_h5, filtered_barcodes_path, scrna_metrics_csv, features_tsv = resolve_inputs(args)

    # Optional confirmation that the DRAGEN prefix is what we think it is.
    if args.dragen_results_dir is not None:
        prefix = detect_dragen_prefix(args.dragen_results_dir)
        if prefix:
            print(f'DRAGEN output-file-prefix: {prefix.rstrip(".")}')

    # Set up an optional multiprocessing pool. Using 'spawn' matches the
    # parent pipeline; 'fork' would leak the whole importer state.
    parallel_pool = None
    if args.threads > 1:
        parallel_pool = multiprocessing.get_context('spawn').Pool(args.threads)

    try:
        minf = MoleculeInfo(
            random_seed=args.random_seed,
            parallel_pool=parallel_pool,
            silent=args.quiet,
        )
        minf.read_dragen_moleculeInfo_h5(molecule_info_h5, features_file=features_tsv)

        total_input_reads = resolve_total_input_reads(
            args, scrna_metrics_csv, minf.total_mi_reads
        )

        # The filtered-barcodes list is required whenever we compute a
        # reads-per-cell (both saturation and matrix paths depend on
        # num_cells). Load it once and reuse.
        filtered_barcodes_seq = resolve_filtered_barcodes(filtered_barcodes_path)

        if args.saturation:
            print('Running saturation sweep...')
            sat_results = minf.run_downsample(
                cell_barcode_sequences=filtered_barcodes_seq,
                add_reads_per_cell=args.saturation_extra_depths,
                min_reads_per_cell=args.min_reads_per_cell,
                total_input_reads=total_input_reads,
            )
            write_saturation_csv(sat_results, args.output_dir)

        if args.matrix_depths:
            print('Building downsampled matrices...')
            matrices, target_rpcs = minf.run_downsample_matrix(
                reads_per_cell=list(args.matrix_depths),
                cell_barcode_sequences=filtered_barcodes_seq,
                min_reads_per_cell=args.min_reads_per_cell,
                total_input_reads=total_input_reads,
            )
            write_matrices(
                minf,
                target_rpcs,
                matrices,
                filtered_barcodes_seq,
                args.output_dir,
                args.matrix_format,
            )
    finally:
        if parallel_pool is not None:
            parallel_pool.close()
            parallel_pool.join()

    print('Done.')


if __name__ == '__main__':
    main()
