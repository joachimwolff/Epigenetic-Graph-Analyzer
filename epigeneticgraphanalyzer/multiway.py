"""Ligation-free multi-way chromatin contact data (SPRITE, GAM).

Both assays are represented by a sparse binary incidence matrix X with one
row per genomic bin and one column per group: a SPRITE cluster or a GAM tube
(nuclear profile). X[i, g] = 1 if bin i was detected in group g. Each group
carries a weight; SPRITE clusters with n bins are down-weighted by 2/n as in
Quinodoz et al. (2018), GAM tubes have weight 1.

Co-occurrence of a set of bins is the (weighted) number of groups that contain
all of them; for GAM pairs the normalised pointwise mutual information (NPMI)
of Beagrie et al. (2017) is available as well.
"""

import logging
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import parallel
from .genome import BinIndex, IntervalIndex, normalize_chrom, open_text

log = logging.getLogger(__name__)


class MultiwayData:

    def __init__(self, bins, incidence, group_weight, kind):
        self.bins = bins
        self.incidence = incidence.tocsr().astype(np.float32)
        self.incidence.sort_indices()
        self.group_weight = np.asarray(group_weight, dtype=np.float64)
        self.kind = kind
        self.n_groups = self.incidence.shape[1]
        self.coverage_raw = np.asarray(self.incidence.sum(axis=1)).ravel()
        self.coverage_weighted = self.incidence @ self.group_weight
        self._expected = {}
        self._band = None
        self._band_offset = 0
        self._ones = None
        self.threads = 1
        self.batch_rows = None   # rows per sparse query batch; None: from the matrix density (see _batches)

    def set_threads(self, threads):
        """Worker processes for pair and set queries (results do not depend on the number)."""
        self.threads = parallel.resolve_threads(threads)
        return self

    def unit_weights(self):
        if self._ones is None:
            self._ones = np.ones(self.n_groups)
        return self._ones

    # ------------------------------------------------------------------ io
    def save(self, path):
        extra = {}
        if not self.bins.fixed:
            extra = dict(interval_chroms=np.array([self.bins.chrom_names[i] for i in self.bins._chrom_of_bin]),
                         interval_starts=self.bins.starts, interval_ends=self.bins.ends)
        np.savez_compressed(path, indptr=self.incidence.indptr, indices=self.incidence.indices,
                            shape=np.array(self.incidence.shape), group_weight=self.group_weight,
                            kind=np.array(self.kind), bin_size=np.array(self.bins.bin_size),
                            chroms=np.array(list(self.bins.chrom_sizes.keys())),
                            sizes=np.array(list(self.bins.chrom_sizes.values())), **extra)

    @classmethod
    def load(cls, path):
        data = np.load(path, allow_pickle=False)
        shape = tuple(int(x) for x in data['shape'])
        incidence = sp.csr_matrix((np.ones(len(data['indices']), dtype=np.float32), data['indices'],
                                   data['indptr']), shape=shape)
        if 'interval_starts' in data:
            bins = IntervalIndex([str(c) for c in data['interval_chroms']], data['interval_starts'],
                                 data['interval_ends'])
        else:
            bins = BinIndex(dict(zip([str(c) for c in data['chroms']], [int(s) for s in data['sizes']])),
                            int(data['bin_size']))
        return cls(bins, incidence, data['group_weight'], str(data['kind']))

    # ------------------------------------------------------------ queries
    def _batches(self, n, *bin_arrays):
        mean_nnz = max(1.0, self.incidence.nnz / max(1, self.incidence.shape[0]))
        # about 3e5 stored entries per batch: small batches stay in cache and spread over the workers
        batch = self.batch_rows or int(max(256, min(200000, 3e5 / mean_nnz)))
        for start in range(0, n, batch):
            yield start, [array[start:start + batch] for array in bin_arrays]

    def build_band(self, max_offset):
        """Pairwise co-occurrence of all bins up to `max_offset` bins apart on the same chromosome.

        band[i, k] holds the weighted (index 0) and unweighted (index 1) number of
        groups containing bins i and i + k. One sparse element-wise product per
        offset and chromosome; pair queries within the band become lookups.
        """
        if self._band is not None and self._band_offset >= max_offset:
            return
        if not self.bins.fixed or self.incidence.shape[0] * (max_offset + 1) > 4e8:
            # variable intervals or a band too large for memory: pair queries use the sparse path
            return
        band = np.full((2, self.incidence.shape[0], max_offset + 1), np.nan, dtype=np.float64)
        ones = self.unit_weights()
        for chrom in self.bins.chrom_names:
            start = self.bins.offsets[chrom]
            n_chrom = self.bins.n_chrom_bins(chrom)
            block = self.incidence[start:start + n_chrom]
            for offset in range(0, min(max_offset, n_chrom - 1) + 1):
                if offset == 0:
                    product = block
                else:
                    product = block[:-offset].multiply(block[offset:]).tocsr()
                band[0, start:start + n_chrom - offset, offset] = product @ self.group_weight
                band[1, start:start + n_chrom - offset, offset] = product @ ones
        self._band = band
        self._band_offset = max_offset
        log.info('co-occurrence band up to %d bins computed', max_offset)

    def cooccurrence(self, *bin_arrays, weighted=True):
        """Number of groups containing all given bins, element-wise over equally long bin arrays."""
        arrays = [np.asarray(a, dtype=np.int64) for a in bin_arrays]
        n = len(arrays[0])
        result = np.zeros(n, dtype=np.float64)
        weights = self.group_weight if weighted else self.unit_weights()
        if len(arrays) == 2 and self._band is not None and n:
            low = np.minimum(arrays[0], arrays[1])
            offset = np.abs(arrays[1] - arrays[0])
            in_band = (offset <= self._band_offset)
            in_band &= self.bins._chrom_of_bin[low] == self.bins._chrom_of_bin[np.maximum(arrays[0], arrays[1])]
            result[in_band] = self._band[0 if weighted else 1, low[in_band], offset[in_band]]
            rest = np.nonzero(~in_band)[0]
            if len(rest):
                result[rest] = self._sparse_cooccurrence([arrays[0][rest], arrays[1][rest]], weights)
            return result
        return self._sparse_cooccurrence(arrays, weights)

    def _sparse_cooccurrence(self, arrays, weights):
        n = len(arrays[0])
        result = np.zeros(n, dtype=np.float64)
        batches = list(self._batches(n, *arrays))
        if self.threads > 1 and len(batches) > 1:
            # each batch is computed independently (same batches as the serial path), so the values are identical
            parallel.share(multiway_data=self, multiway_weights=weights)
            values = parallel.map_chunks(_cooccurrence_batch, [chunk for _, chunk in batches], self.threads)
            for (start, chunk), value in zip(batches, values):
                result[start:start + len(chunk[0])] = value
            return result
        for start, chunk in batches:
            result[start:start + len(chunk[0])] = _cooccurrence_of(self.incidence, chunk, weights)
        return result

    def _sampled_expected(self, max_offset, metric):
        return _sampled_expected_impl(self, max_offset, metric)

    def coverage(self, bins, weighted=True):
        bins = np.asarray(bins, dtype=np.int64)
        return (self.coverage_weighted if weighted else self.coverage_raw)[bins]

    def npmi(self, bin_a, bin_b):
        """Normalised pointwise mutual information of co-segregation (unweighted groups)."""
        n_ab = self.cooccurrence(bin_a, bin_b, weighted=False)
        n_a = self.coverage(bin_a, weighted=False)
        n_b = self.coverage(bin_b, weighted=False)
        return npmi_from_counts(n_ab, n_a, n_b, self.n_groups)

    def set_cooccurrence(self, bin_set, min_members, weighted=False):
        """Number of groups containing at least `min_members` of the bins in `bin_set`.

        Unweighted by default: a group either holds the multi-way contact or not,
        the 2/n weighting of SPRITE is meant for pairwise contact frequencies.
        """
        bin_set = np.unique(np.asarray(bin_set, dtype=np.int64))
        groups, hits = np.unique(self.incidence[bin_set].indices, return_counts=True)
        selected = groups[hits >= min_members]
        return float(self.group_weight[selected].sum()) if weighted else float(len(selected))

    def expected_by_offset(self, max_offset, metric='cooc'):
        """Mean metric of all covered bin pairs at each bin offset 0..max_offset (index = offset)."""
        key = (max_offset, metric)
        if key in self._expected:
            return self._expected[key]
        self.build_band(max_offset)
        if self._band is None:
            clock = time.time()
            self._expected[key] = self._sampled_expected(max_offset, metric)
            log.info('expected %s by offset (0..%d, sampled pairs): %.1f s', metric, max_offset, time.time() - clock)
            return self._expected[key]
        covered = self.coverage_raw > 0
        expected = np.full(max_offset + 1, np.nan)
        for offset in range(1, max_offset + 1):
            left = np.nonzero(covered[:-offset] & covered[offset:])[0]
            left = left[self.bins._chrom_of_bin[left] == self.bins._chrom_of_bin[left + offset]]
            if len(left) == 0:
                continue
            if metric == 'cooc':
                values = self._band[0, left, offset]
            elif metric == 'npmi':
                values = npmi_from_counts(self._band[1, left, offset], self.coverage_raw[left],
                                          self.coverage_raw[left + offset], self.n_groups)
            else:
                raise ValueError(metric)
            expected[offset] = np.nanmean(values)
        self._expected[key] = expected
        return expected


def _cooccurrence_of(incidence, chunk, weights):
    product = incidence[chunk[0]]
    for other in chunk[1:]:
        product = product.multiply(incidence[other])
    return sp.csr_matrix(product) @ weights


def _cooccurrence_batch(chunk):
    return _cooccurrence_of(parallel.shared('multiway_data').incidence, chunk, parallel.shared('multiway_weights'))


def _sampled_expected_impl(data, max_offset, metric, samples_per_offset=20000, seed=0):
    """Mean metric of covered pairs at each offset, from a random sample of pairs per offset.

    The samples of all offsets are drawn first (one random stream, offset by offset) and queried together, so that
    the pair queries run in parallel batches."""
    rng = np.random.default_rng(seed)
    covered = np.nonzero(data.coverage_raw > 0)[0]
    chrom_of = data.bins._chrom_of_bin
    expected = np.full(max_offset + 1, np.nan)
    covered_mask = data.coverage_raw > 0
    offsets, lefts = [], []
    for offset in range(1, max_offset + 1):
        left = covered[covered + offset < len(chrom_of)]
        left = left[covered_mask[left + offset] & (chrom_of[left] == chrom_of[left + offset])]
        if len(left) == 0:
            continue
        if len(left) > samples_per_offset:
            left = rng.choice(left, samples_per_offset, replace=False)
        offsets.append(offset)
        lefts.append(left)
    if not offsets:
        return expected
    left = np.concatenate(lefts)
    right = left + np.repeat(np.array(offsets, dtype=np.int64), [len(x) for x in lefts])
    values = data.cooccurrence(left, right) if metric == 'cooc' else data.npmi(left, right)
    bounds = np.cumsum([0] + [len(x) for x in lefts])
    for k, offset in enumerate(offsets):
        expected[offset] = np.nanmean(values[bounds[k]:bounds[k + 1]])
    return expected


def npmi_from_counts(n_ab, n_a, n_b, n_total):
    """NPMI = log(p_ab / (p_a p_b)) / -log(p_ab); -1 if the pair never co-occurs, NaN if a bin is empty."""
    n_ab = np.asarray(n_ab, dtype=float)
    n_a = np.asarray(n_a, dtype=float)
    n_b = np.asarray(n_b, dtype=float)
    result = np.full(n_ab.shape, np.nan)
    defined = (n_a > 0) & (n_b > 0)
    never = defined & (n_ab == 0)
    result[never] = -1.0
    both = defined & (n_ab > 0)
    p_ab = n_ab[both] / n_total
    p_a = n_a[both] / n_total
    p_b = n_b[both] / n_total
    with np.errstate(divide='ignore', invalid='ignore'):
        pmi = np.log(p_ab / (p_a * p_b))
        denominator = -np.log(p_ab)
        values = np.where(denominator > 0, pmi / denominator, 1.0)
    result[both] = values
    return result


def read_sprite_clusters(path, chrom_sizes, bin_size, min_size=2, max_size=1000, chunk_reads=5000000, index=None):
    """Parse a SPRITE cluster file (barcode TAB chrom:pos TAB ...) into MultiwayData.

    Clusters are filtered by their number of reads (`min_size`..`max_size`),
    reads are binned and de-duplicated per cluster. The group weight is 2/n with
    n the number of distinct bins of the cluster. With `index` (an IntervalIndex,
    e.g. restriction fragments) reads are assigned to the interval containing them.
    """
    bins = index if index is not None else BinIndex(chrom_sizes, bin_size)
    known = set(bins.chrom_sizes)
    row_parts, col_parts = [], []
    chroms, positions, groups = [], [], []
    n_clusters = 0

    def flush():
        if not positions:
            return
        bin_ids = bins.bins_of(np.array(chroms), np.array(positions, dtype=np.int64))
        group_ids = np.array(groups, dtype=np.int64)
        keep = bin_ids >= 0
        row_parts.append(bin_ids[keep])
        col_parts.append(group_ids[keep])
        chroms.clear()
        positions.clear()
        groups.clear()

    with open_text(path) as handle:
        for line in handle:
            tokens = line.rstrip('\n').split('\t')
            n_reads = len(tokens) - 1
            if n_reads < min_size or n_reads > max_size:
                continue
            for token in tokens[1:]:
                chrom, _, position = token.partition(':')
                if chrom not in known:
                    continue
                chroms.append(chrom)
                positions.append(int(position))
                groups.append(n_clusters)
            n_clusters += 1
            if len(positions) >= chunk_reads:
                flush()
    flush()
    rows = np.concatenate(row_parts) if row_parts else np.zeros(0, dtype=np.int64)
    cols = np.concatenate(col_parts) if col_parts else np.zeros(0, dtype=np.int64)
    incidence = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                              shape=(bins.n_bins, n_clusters))
    incidence.data[:] = 1.0  # duplicates of a bin within a cluster were summed
    # drop clusters that lost all but one bin after binning
    bins_per_cluster = np.asarray((incidence > 0).sum(axis=0)).ravel()
    keep = bins_per_cluster >= 2
    incidence = incidence[:, np.nonzero(keep)[0]]
    weight = 2.0 / bins_per_cluster[keep]
    log.info('SPRITE: %d clusters kept (%d with >= 2 %s)', n_clusters, int(keep.sum()),
             'bins of {} bp'.format(bin_size) if bins.fixed else 'fragments')
    return MultiwayData(bins, incidence, weight, 'sprite')


def read_gam_segregation(paths, chrom_sizes):
    """Read one or more GAM segregation tables (chrom, start, stop, one 0/1 column per tube).

    Tables of replicates with identical windows are concatenated column-wise.
    """
    blocks = []
    bin_size = None
    windows = None
    for path in paths:
        table = pd.read_csv(path, sep='\t', compression='infer')
        table = table.rename(columns={table.columns[0]: 'chrom', table.columns[1]: 'start', table.columns[2]: 'stop'})
        table['chrom'] = table['chrom'].map(normalize_chrom)
        size = int(np.median(table['stop'] - table['start']))
        if bin_size is None:
            bin_size = size
            windows = table[['chrom', 'start']]
        elif size != bin_size or not windows.equals(table[['chrom', 'start']]):
            raise ValueError('GAM tables use different windows: {}'.format(path))
        values = table.iloc[:, 3:].to_numpy(dtype=np.float32)
        values = (values > 0).astype(np.float32)
        blocks.append(values)
    bins = BinIndex(chrom_sizes, bin_size)
    bin_ids = bins.bins_of(windows['chrom'].to_numpy(), windows['start'].to_numpy())
    matrix = np.hstack(blocks)
    keep = bin_ids >= 0
    matrix = sp.csr_matrix(matrix[keep])
    rows = bin_ids[keep]
    incidence = sp.csr_matrix((bins.n_bins, matrix.shape[1]), dtype=np.float32)
    coo = matrix.tocoo()
    incidence = sp.csr_matrix((coo.data, (rows[coo.row], coo.col)), shape=(bins.n_bins, matrix.shape[1]))
    # tubes without any detected window carry no information
    used = np.asarray(incidence.sum(axis=0)).ravel() > 0
    incidence = incidence[:, np.nonzero(used)[0]]
    log.info('GAM: %d tubes, %d windows of %d bp', incidence.shape[1], int(keep.sum()), bin_size)
    return MultiwayData(bins, incidence, np.ones(incidence.shape[1]), 'gam')
