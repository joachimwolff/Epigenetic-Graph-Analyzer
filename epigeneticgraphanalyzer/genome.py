"""Genomic coordinate helpers: node identifiers, chromosome sizes, binning, liftover."""

import gzip
import logging
import re
from collections import OrderedDict

import numpy as np

log = logging.getLogger(__name__)

_NODE_RE = re.compile(r'^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$')


def normalize_chrom(chrom):
    """Return a UCSC style chromosome name ('1' -> 'chr1', 'MT' -> 'chrM')."""
    chrom = str(chrom)
    if chrom.startswith('chr'):
        return chrom
    if chrom in ('MT', 'M'):
        return 'chrM'
    return 'chr' + chrom


def node_id(chrom, start, end):
    """Canonical node identifier of a genomic fragment."""
    return '{}:{}-{}'.format(chrom, int(start), int(end))


def parse_node_id(identifier):
    """Split 'chr:start-end' into (chrom, start, end)."""
    match = _NODE_RE.match(str(identifier))
    if match is None:
        raise ValueError('Not a genomic node identifier: {}'.format(identifier))
    return match.group('chrom'), int(match.group('start')), int(match.group('end'))


def midpoint(start, end):
    return (int(start) + int(end)) // 2


def open_text(path):
    """Open plain or gzip/bzip2 compressed text files transparently."""
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    if path.endswith('.bz2'):
        import bz2
        return bz2.open(path, 'rt')
    return open(path, 'r')


def read_chrom_sizes(path, keep_regex=r'^chr([0-9]+|X|Y)$'):
    """Read a two column chromosome size file. Unplaced contigs are dropped by default."""
    sizes = OrderedDict()
    pattern = re.compile(keep_regex) if keep_regex else None
    with open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            chrom = normalize_chrom(fields[0])
            if pattern is not None and not pattern.match(chrom):
                continue
            sizes[chrom] = int(fields[1])
    return sizes


class BinIndex:
    """Fixed size genomic bins with a global integer id per bin."""

    def __init__(self, chrom_sizes, bin_size):
        self.bin_size = int(bin_size)
        self.chrom_sizes = OrderedDict(chrom_sizes)
        self.offsets = OrderedDict()
        current = 0
        for chrom, size in self.chrom_sizes.items():
            self.offsets[chrom] = current
            current += (size + self.bin_size - 1) // self.bin_size
        self.n_bins = current
        self._chrom_of_bin = np.empty(self.n_bins, dtype=np.int32)
        self.chrom_names = list(self.chrom_sizes.keys())
        for i, chrom in enumerate(self.chrom_names):
            start = self.offsets[chrom]
            self._chrom_of_bin[start:start + self.n_chrom_bins(chrom)] = i

    def n_chrom_bins(self, chrom):
        return (self.chrom_sizes[chrom] + self.bin_size - 1) // self.bin_size

    def bin_of(self, chrom, position):
        """Global bin id of a position or -1 if the position is outside the index."""
        offset = self.offsets.get(chrom)
        if offset is None or position < 0 or position >= self.chrom_sizes[chrom]:
            return -1
        return offset + int(position) // self.bin_size

    def bins_of(self, chroms, positions):
        """Vectorised bin_of for arrays of chromosome names and positions."""
        chroms = np.asarray(chroms)
        positions = np.asarray(positions, dtype=np.int64)
        result = np.full(len(positions), -1, dtype=np.int64)
        for chrom in np.unique(chroms):
            if chrom not in self.offsets:
                continue
            mask = chroms == chrom
            pos = positions[mask]
            valid = (pos >= 0) & (pos < self.chrom_sizes[chrom])
            values = np.full(len(pos), -1, dtype=np.int64)
            values[valid] = self.offsets[chrom] + pos[valid] // self.bin_size
            result[mask] = values
        return result

    def chrom_index_of_bins(self, bins):
        return self._chrom_of_bin[np.asarray(bins, dtype=np.int64)]

    fixed = True

    def bin_start(self, bins):
        bins = np.asarray(bins, dtype=np.int64)
        chrom_idx = self._chrom_of_bin[bins]
        offsets = np.array([self.offsets[c] for c in self.chrom_names], dtype=np.int64)
        return (bins - offsets[chrom_idx]) * self.bin_size

    def bin_end(self, bins):
        bins = np.asarray(bins, dtype=np.int64)
        sizes = np.array([self.chrom_sizes[c] for c in self.chrom_names], dtype=np.int64)
        return np.minimum(self.bin_start(bins) + self.bin_size, sizes[self._chrom_of_bin[bins]])


class IntervalIndex:
    """Genomic intervals of variable length (e.g. restriction fragments) with a global integer id each.

    Offers the interface of BinIndex; ids are ordered by chromosome and start, so the difference of two ids
    on one chromosome is the number of intervals between them. Positions outside every interval map to -1.
    """

    fixed = False

    def __init__(self, chroms, starts, ends):
        order = np.lexsort((np.asarray(starts), np.asarray(chroms)))
        chroms = np.asarray(chroms)[order]
        self.starts = np.asarray(starts, dtype=np.int64)[order]
        self.ends = np.asarray(ends, dtype=np.int64)[order]
        self.chrom_names = []
        self.offsets = OrderedDict()
        self._counts = {}
        for chrom in pd_unique(chroms):
            self.chrom_names.append(str(chrom))
        position = 0
        self._chrom_of_bin = np.empty(len(chroms), dtype=np.int32)
        for i, chrom in enumerate(self.chrom_names):
            count = int(np.sum(chroms == chrom))
            self.offsets[chrom] = position
            self._counts[chrom] = count
            self._chrom_of_bin[position:position + count] = i
            position += count
        self.n_bins = position
        self.chrom_sizes = OrderedDict((c, int(self.ends[self.offsets[c] + self._counts[c] - 1]))
                                       for c in self.chrom_names)
        self.bin_size = int(np.median(self.ends - self.starts))

    @classmethod
    def from_bed(cls, path):
        chroms, starts, ends = [], [], []
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split()
                chroms.append(normalize_chrom(fields[0]))
                starts.append(int(fields[1]))
                ends.append(int(fields[2]))
        return cls(chroms, starts, ends)

    def n_chrom_bins(self, chrom):
        return self._counts[chrom]

    def bin_of(self, chrom, position):
        return int(self.bins_of(np.array([chrom]), np.array([position]))[0])

    def bins_of(self, chroms, positions):
        chroms = np.asarray(chroms)
        positions = np.asarray(positions, dtype=np.int64)
        result = np.full(len(positions), -1, dtype=np.int64)
        for chrom in np.unique(chroms):
            if chrom not in self.offsets:
                continue
            first = self.offsets[chrom]
            starts = self.starts[first:first + self._counts[chrom]]
            ends = self.ends[first:first + self._counts[chrom]]
            mask = chroms == chrom
            index = np.searchsorted(starts, positions[mask], side='right') - 1
            valid = (index >= 0) & (positions[mask] < ends[np.clip(index, 0, len(ends) - 1)])
            values = np.full(len(index), -1, dtype=np.int64)
            values[valid] = first + index[valid]
            result[mask] = values
        return result

    def chrom_index_of_bins(self, bins):
        return self._chrom_of_bin[np.asarray(bins, dtype=np.int64)]

    def bin_start(self, bins):
        return self.starts[np.asarray(bins, dtype=np.int64)]

    def bin_end(self, bins):
        return self.ends[np.asarray(bins, dtype=np.int64)]


def pd_unique(values):
    """Unique values in order of first occurrence."""
    seen = OrderedDict()
    for value in values:
        seen.setdefault(value, None)
    return list(seen)


class LiftOver:
    """Midpoint liftover of genomic positions with a UCSC chain file."""

    def __init__(self, chain_file):
        from pyliftover import LiftOver as _LiftOver
        self._lift = _LiftOver(chain_file)

    def position(self, chrom, position):
        hits = self._lift.convert_coordinate(chrom, int(position))
        if not hits:
            return None
        new_chrom, new_pos = hits[0][0], hits[0][1]
        return new_chrom, int(new_pos)

    def hits(self, chrom, position):
        """All lifted positions (chrom, position) of a position; position() uses the first."""
        return [(h[0], int(h[1])) for h in (self._lift.convert_coordinate(chrom, int(position)) or [])]
