"""TAD / domain boundaries crossed by the edges of a capture graph.

Readers (`read_tads`, auto-detected or with an explicit format):

* 'hicexplorer': hicFindTADs output, given as the prefix or any of <prefix>_boundaries.bed, _domains.bed,
  _boundaries.gff, _tad_score.bm, _score.bedgraph. Boundary score = TAD-separation score (lower: stronger
  insulation), strength = -score; the gff adds delta and pvalue. With only a score file, boundaries are the local
  minima of the score (`boundaries_from_track`).
* 'arrowhead': Juicer Arrowhead domain lists and other 2D domain lists (4DN domain calls): chr1 x1 x2 chr2 y1 y2
  [color corner_score ...]. Domains (chr1, x1, x2); boundaries are the domain ends, strength = the largest corner
  score of the domains ending there (else the number of domains).
* 'cooltools': insulation table (TSV with log2_insulation_score_<w>, boundary_strength_<w>, is_boundary_<w>); one
  window (`window`, default the largest); strength = boundary_strength.
* 'bed' (auto: boundaries or domains from the interval widths), 'bed-boundaries', 'bed-domains': generic / 4DN BED,
  optional score column (5th) = strength.
* 'ega': egaBoundaries output (<prefix>_bins.tsv, _boundaries.bed, _domains.bed, _capture_insulation.bedGraph);
  strength = rule score.

Chromosome names are harmonised to those of the graph ('6' -> 'chr6' and back, `chrom_mapper`); an optional chain
file lifts the boundaries (midpoints) and the domain ends to the assembly of the graph.

Crossing (`edge_crossings`): an edge between the midpoints a < b of its two fragments crosses a boundary whose
midpoint p lies strictly between them (a < p < b). Edges are 'same domain' when one domain holds both ends,
'cross domain' when the ends lie in different domains and 'outside domains' when an end lies in no domain (gaps of
the domain calls); without domain calls the domains are the segments between consecutive boundaries. Boundary
strength is also given as its percentile among all boundaries of the input ('strong': top decile).

Per gene (`gene_crossings`): gained / lost / shared contacts of each promoter fragment that cross a boundary, the
strongest boundary crossed by at least half of them, the region of the partners beyond it, their element classes
and mark changes, a sentence. Genome-wide (`promoter_ranking`, `boundary_ranking`): per promoter, the enrichment of
gained (or lost) edges among its crossing edges against its non-crossing edges (one-sided Fisher / hypergeometric,
Benjamini-Hochberg); per boundary, gained (lost) crossing edges against the expectation from all union edges of the
same length (20 length quantiles; exact tail of the sum of binomials).
"""

import logging
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

from .genome import LiftOver, normalize_chrom, open_text, parse_node_id

log = logging.getLogger(__name__)

FORMATS = ('auto', 'hicexplorer', 'arrowhead', 'cooltools', 'bed', 'bed-boundaries', 'bed-domains', 'ega')
HICEXPLORER_SUFFIXES = ('_boundaries.bed', '_domains.bed', '_boundaries.gff', '_tad_score.bm', '_score.bedgraph',
                        '_zscore_matrix.h5')
EGA_SUFFIXES = ('_bins.tsv', '_boundaries.bed', '_domains.bed', '_capture_insulation.bedGraph')
STATUSES = ('gained', 'lost', 'shared')
STRONG_PCT = 0.9


# ----------------------------------------------------------------------------------------------------------------
# container

class Tads:
    """Boundaries (chrom, start, end, pos, score, strength, name, strength_pct, id), optional domains (chrom, start,
    end, name) and an optional insulation / TAD-separation track (chrom, start, end, value)."""

    def __init__(self, boundaries, domains=None, track=None, fmt='', source=''):
        b = boundaries.copy()
        if 'pos' not in b:
            b['pos'] = ((b['start'] + b['end']) // 2).astype(np.int64)
        for column in ('score', 'strength'):
            if column not in b:
                b[column] = np.nan
        if 'name' not in b:
            b['name'] = ''
        b = b.sort_values(['chrom', 'pos'], kind='stable').drop_duplicates(['chrom', 'pos']).reset_index(drop=True)
        strength = b['strength'].astype(float)
        b['strength_pct'] = strength.rank(pct=True, method='average') if strength.notna().any() else np.nan
        b['strength_rank'] = strength.rank(ascending=False, method='min') if strength.notna().any() else np.nan
        b['id'] = ['{}:{}'.format(c, int(p)) for c, p in zip(b['chrom'], b['pos'])]
        self.boundaries = b
        self.domains = None if domains is None or not len(domains) else \
            domains.sort_values(['chrom', 'start']).reset_index(drop=True)
        self.track = track
        self.fmt = fmt
        self.source = source
        self._chroms = {}
        for chrom, g in b.groupby('chrom', sort=False):
            self._chroms[chrom] = (g.index.values, g['pos'].values.astype(np.int64),
                                   g['strength'].values.astype(float))

    def __len__(self):
        return len(self.boundaries)

    @property
    def chroms(self):
        return list(self._chroms)

    def map_chroms(self, mapper):
        """Rename chromosomes with `mapper(name) -> name or None` (None drops)."""
        def apply(df):
            if df is None:
                return None
            names = df['chrom'].map(lambda c: mapper(c))
            return df.assign(chrom=names)[names.notna()]
        return Tads(apply(self.boundaries).drop(columns=['id', 'strength_pct', 'strength_rank']), apply(self.domains),
                    apply(self.track), self.fmt, self.source)

    def region(self, chrom, start, end):
        b = self.boundaries
        return b[(b.chrom == chrom) & (b.pos > start) & (b.pos < end)]


# ----------------------------------------------------------------------------------------------------------------
# readers

def _read_bed(path, min_columns=3):
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < min_columns:
                fields = line.split()
            if len(fields) < min_columns:
                continue
            try:
                int(float(fields[1]))
                int(float(fields[2]))
            except ValueError:
                continue  # header
            rows.append(fields)
    return rows


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def boundaries_from_domains(domains, scores=None):
    """Domain ends as boundaries; strength = largest score of the domains ending there (else their number)."""
    if domains is None or not len(domains):
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'pos', 'score', 'strength', 'name'])
    ends = pd.concat([domains[['chrom', 'start']].rename(columns={'start': 'pos'}).assign(s=scores),
                      domains[['chrom', 'end']].rename(columns={'end': 'pos'}).assign(s=scores)], ignore_index=True)
    grouped = ends.groupby(['chrom', 'pos'], sort=False)
    out = grouped.size().rename('n').reset_index()
    best = grouped['s'].max().values if scores is not None else np.full(len(out), np.nan)
    out['score'] = best
    out['strength'] = np.where(np.isfinite(best.astype(float)), best, out['n'].values)
    out['start'] = out['pos']
    out['end'] = out['pos'] + 1
    out['name'] = ''
    return out.drop(columns='n')


def domains_from_boundaries(boundaries, chrom_sizes=None):
    """Segments between consecutive boundaries of a chromosome (from 0 and to the chromosome end when sizes are
    given)."""
    rows = []
    for chrom, g in boundaries.groupby('chrom', sort=False):
        pos = np.sort(g['pos'].values.astype(np.int64))
        edges = list(pos)
        if chrom_sizes and chrom in chrom_sizes:
            edges = [0] + edges + [int(chrom_sizes[chrom])]
        for s, e in zip(edges[:-1], edges[1:]):
            if e > s:
                rows.append((chrom, int(s), int(e)))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


def boundaries_from_track(track, prominence=0.3):
    """Local minima of an insulation-like track (lower: stronger insulation) with a minimum prominence."""
    from scipy import signal
    rows = []
    for chrom, g in track.groupby('chrom', sort=False):
        v = g['value'].values.astype(float)
        filled = np.where(np.isfinite(v), v, np.nanmax(v) if np.isfinite(v).any() else 0.0)
        peaks, props = signal.find_peaks(-filled, prominence=prominence)
        for k, prom in zip(peaks, props['prominences']):
            rows.append((chrom, int(g['start'].values[k]), int(g['end'].values[k]), float(v[k]), float(prom)))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'score', 'strength'])


def _prefix(path, suffixes):
    for suffix in suffixes:
        if str(path).endswith(suffix):
            return str(path)[:-len(suffix)]
    return str(path)


def read_hicexplorer(path):
    prefix = _prefix(path, HICEXPLORER_SUFFIXES)
    exists = {s: os.path.exists(prefix + s) for s in HICEXPLORER_SUFFIXES}
    track = None
    if exists['_score.bedgraph']:
        rows = _read_bed(prefix + '_score.bedgraph', 4)
        track = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                              'end': [int(r[2]) for r in rows], 'value': [_float(r[3]) for r in rows]})
    elif exists['_tad_score.bm']:
        rows = _read_bed(prefix + '_tad_score.bm', 4)
        track = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                              'end': [int(r[2]) for r in rows],
                              'value': [np.nanmean([_float(x) for x in r[3:]]) for r in rows]})
    domains = None
    if exists['_domains.bed']:
        rows = _read_bed(prefix + '_domains.bed')
        domains = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                                'end': [int(r[2]) for r in rows],
                                'name': [r[3] if len(r) > 3 else '' for r in rows]})
    if exists['_boundaries.bed']:
        rows = _read_bed(prefix + '_boundaries.bed')
        b = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                          'end': [int(r[2]) for r in rows], 'name': [r[3] if len(r) > 3 else '' for r in rows],
                          'score': [_float(r[4]) if len(r) > 4 else np.nan for r in rows]})
        b['strength'] = -b['score']
    elif domains is not None:
        b = boundaries_from_domains(domains)
        b['strength'] = np.nan
    elif track is not None:
        b = boundaries_from_track(track)
    else:
        raise FileNotFoundError('no hicFindTADs output found for prefix {}'.format(prefix))
    if exists['_boundaries.gff']:
        extra = {}
        with open_text(prefix + '_boundaries.gff') as handle:
            for line in handle:
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or line.startswith('#'):
                    continue
                attrs = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
                extra[(f[0], int(f[3]), int(f[4]))] = (_float(attrs.get('delta')), _float(attrs.get('pvalue')))
        keys = list(zip(b['chrom'], b['start'], b['end']))
        b['delta'] = [extra.get(k, (np.nan, np.nan))[0] for k in keys]
        b['pvalue'] = [extra.get(k, (np.nan, np.nan))[1] for k in keys]
    return Tads(b, domains, track, 'hicexplorer', prefix)


def read_arrowhead(path):
    rows = _read_bed(path, 6)
    rows = [r for r in rows if len(r) >= 6]
    corner = []
    for r in rows:
        corner.append(_float(r[7]) if len(r) > 7 else np.nan)
    domains = pd.DataFrame({'chrom': [r[0] for r in rows],
                            'start': [int(float(r[1])) for r in rows],
                            'end': [int(float(r[5])) if r[5] else int(float(r[2])) for r in rows]})
    domains['end'] = np.maximum(domains['end'].values, [int(float(r[2])) for r in rows])
    domains['name'] = ['D{}'.format(k + 1) for k in range(len(domains))]
    scores = np.asarray(corner, dtype=float)
    b = boundaries_from_domains(domains, scores if np.isfinite(scores).any() else None)
    return Tads(b, domains, None, 'arrowhead', str(path))


def read_cooltools(path, window=None):
    table = pd.read_csv(path, sep='\t')
    windows = sorted({int(c.rsplit('_', 1)[1]) for c in table.columns if c.startswith('is_boundary_')})
    if not windows:
        raise ValueError('{}: no is_boundary_<window> column'.format(path))
    w = int(window) if window else windows[-1]
    if w not in windows:
        raise ValueError('{}: window {} not among {}'.format(path, w, windows))
    flag = table['is_boundary_{}'.format(w)].astype(str).str.lower().isin(('true', '1'))
    b = table.loc[flag, ['chrom', 'start', 'end']].copy()
    b['score'] = table.loc[flag, 'log2_insulation_score_{}'.format(w)].values \
        if 'log2_insulation_score_{}'.format(w) in table else np.nan
    b['strength'] = table.loc[flag, 'boundary_strength_{}'.format(w)].values \
        if 'boundary_strength_{}'.format(w) in table else np.nan
    track = None
    if 'log2_insulation_score_{}'.format(w) in table:
        track = table[['chrom', 'start', 'end']].assign(value=table['log2_insulation_score_{}'.format(w)].values)
    tads = Tads(b, None, track, 'cooltools', str(path))
    tads.window = w
    return tads


def read_bed(path, kind='auto', min_domain=100000):
    rows = _read_bed(path)
    df = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(float(r[1])) for r in rows],
                       'end': [int(float(r[2])) for r in rows],
                       'name': [r[3] if len(r) > 3 else '' for r in rows],
                       'score': [_float(r[4]) if len(r) > 4 else np.nan for r in rows]})
    if kind == 'auto':
        width = (df['end'] - df['start']).median() if len(df) else 0
        kind = 'domains' if width >= min_domain else 'boundaries'
    if kind == 'domains':
        b = boundaries_from_domains(df[['chrom', 'start', 'end']], df['score'].values
                                    if df['score'].notna().any() else None)
        return Tads(b, df[['chrom', 'start', 'end', 'name']], None, 'bed-domains', str(path))
    df['strength'] = df['score']
    return Tads(df, None, None, 'bed-boundaries', str(path))


def read_ega(path):
    prefix = _prefix(path, EGA_SUFFIXES)
    track = domains = None
    if os.path.exists(prefix + '_bins.tsv'):
        bins = pd.read_csv(prefix + '_bins.tsv', sep='\t')
        flag = bins['boundary'].astype(str).str.lower().isin(('true', '1'))
        b = bins.loc[flag, ['chrom', 'start', 'end']].copy()
        b['score'] = bins.loc[flag, 'score'].values if 'score' in bins else np.nan
        if 'capture_insulation' in bins:
            track = bins[['chrom', 'start', 'end']].assign(value=bins['capture_insulation'].values)
    elif os.path.exists(prefix + '_boundaries.bed'):
        rows = _read_bed(prefix + '_boundaries.bed')
        b = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                          'end': [int(r[2]) for r in rows],
                          'score': [_float(r[4]) if len(r) > 4 else np.nan for r in rows]})
    else:
        raise FileNotFoundError('no egaBoundaries output for prefix {}'.format(prefix))
    b['strength'] = b['score']
    b['name'] = 'boundary'
    if os.path.exists(prefix + '_domains.bed'):
        rows = _read_bed(prefix + '_domains.bed')
        domains = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                                'end': [int(r[2]) for r in rows]})
    if track is None and os.path.exists(prefix + '_capture_insulation.bedGraph'):
        rows = _read_bed(prefix + '_capture_insulation.bedGraph', 4)
        track = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                              'end': [int(r[2]) for r in rows], 'value': [_float(r[3]) for r in rows]})
    return Tads(b, domains, track, 'ega', prefix)


def detect_format(path):
    path = str(path)
    if not os.path.exists(path):
        if any(os.path.exists(path + s) for s in ('_bins.tsv', '_capture_insulation.bedGraph')):
            return 'ega'
        if any(os.path.exists(path + s) for s in HICEXPLORER_SUFFIXES):
            return 'hicexplorer'
        raise FileNotFoundError(path)
    if path.endswith(('_tad_score.bm', '_score.bedgraph', '_boundaries.gff', '_zscore_matrix.h5')):
        return 'hicexplorer'
    if path.endswith(('_bins.tsv', '_capture_insulation.bedGraph')):
        return 'ega'
    with open_text(path) as handle:
        lines = [line for _, line in zip(range(50), handle)]
    header = lines[0] if lines else ''
    if 'is_boundary_' in header:
        return 'cooltools'
    if re.match(r'^#?\s*chr1\s+x1\s+x2\s+chr2\s+y1\s+y2', header):
        return 'arrowhead'
    data = [line.rstrip('\n').split() for line in lines if line.strip() and not line.startswith(('#', 'track'))]
    for fields in data:
        if len(fields) >= 3 and not fields[1].replace('.', '').isdigit():
            continue  # header
        if len(fields) >= 6 and fields[3] == fields[0] and all(f.replace('.', '').isdigit() for f in fields[4:6]):
            return 'arrowhead'
        break
    for suffix in ('_boundaries.bed', '_domains.bed'):
        if path.endswith(suffix):
            prefix = path[:-len(suffix)]
            if any(os.path.exists(prefix + s) for s in ('_tad_score.bm', '_boundaries.gff', '_score.bedgraph')):
                return 'hicexplorer'
            if any(os.path.exists(prefix + s) for s in ('_bins.tsv', '_capture_insulation.bedGraph')):
                return 'ega'
            names = [f[3] for f in data if len(f) > 3]
            if names and all(re.match(r'^(B\d+|ID_)', n) for n in names[:20]):
                return 'hicexplorer'
            if names and all(n == 'boundary' for n in names[:20]):
                return 'ega'
    return 'bed'


def chrom_mapper(targets=None):
    """Map chromosome names to the style of `targets` ('6' <-> 'chr6', 'MT' <-> 'chrM'); without targets to UCSC
    style. Names without a counterpart map to None."""
    if targets is None:
        return normalize_chrom
    targets = set(targets)

    def mapper(chrom):
        chrom = str(chrom)
        if chrom in targets:
            return chrom
        ucsc = normalize_chrom(chrom)
        if ucsc in targets:
            return ucsc
        bare = chrom[3:] if chrom.startswith('chr') else chrom
        for candidate in (bare, 'MT' if bare == 'M' else bare):
            if candidate in targets:
                return candidate
        return None
    return mapper


def liftover_tads(tads, chain):
    """Lift boundaries (midpoints) and domain ends with a chain file (UCSC chromosome names)."""
    lift = LiftOver(chain) if isinstance(chain, str) else chain
    b = tads.boundaries
    rows = []
    for r in b.itertuples(index=False):
        hit = lift.position(normalize_chrom(r.chrom), int(r.pos))
        if hit is None:
            continue
        half = (int(r.end) - int(r.start)) // 2
        d = r._asdict()
        d.update(chrom=hit[0], pos=hit[1], start=hit[1] - half, end=hit[1] - half + (int(r.end) - int(r.start)))
        rows.append(d)
    nb = pd.DataFrame(rows).drop(columns=['id', 'strength_pct', 'strength_rank'], errors='ignore')
    nd = None
    if tads.domains is not None:
        drows = []
        for r in tads.domains.itertuples(index=False):
            s = lift.position(normalize_chrom(r.chrom), int(r.start))
            e = lift.position(normalize_chrom(r.chrom), int(r.end) - 1)
            if s is None or e is None or s[0] != e[0] or e[1] + 1 <= s[1]:
                continue
            d = r._asdict()
            d.update(chrom=s[0], start=s[1], end=e[1] + 1)
            drows.append(d)
        nd = pd.DataFrame(drows)
    nt = None
    if tads.track is not None:
        trows = []
        for r in tads.track.itertuples(index=False):
            s = lift.position(normalize_chrom(r.chrom), (int(r.start) + int(r.end)) // 2)
            if s is not None:
                half = (int(r.end) - int(r.start)) // 2
                trows.append((s[0], s[1] - half, s[1] + half, r.value))
        nt = pd.DataFrame(trows, columns=['chrom', 'start', 'end', 'value']).sort_values(['chrom', 'start'])
    log.info('liftover: %d of %d boundaries, %s of %s domains', len(nb), len(b),
             len(nd) if nd is not None else 0, len(tads.domains) if tads.domains is not None else 0)
    return Tads(nb, nd, nt, tads.fmt, tads.source)


def read_tads(path, fmt='auto', targets=None, chain=None, window=None):
    """Read TAD / boundary results in any supported format, lift them over with `chain` and harmonise the
    chromosome names to `targets` (the chromosomes of the graph)."""
    fmt = fmt or 'auto'
    if fmt == 'auto':
        fmt = detect_format(path)
    if fmt == 'hicexplorer':
        tads = read_hicexplorer(path)
    elif fmt == 'arrowhead':
        tads = read_arrowhead(path)
    elif fmt == 'cooltools':
        tads = read_cooltools(path, window)
    elif fmt in ('bed', 'bed-boundaries', 'bed-domains'):
        tads = read_bed(path, {'bed': 'auto', 'bed-boundaries': 'boundaries', 'bed-domains': 'domains'}[fmt])
    elif fmt == 'ega':
        tads = read_ega(path)
    else:
        raise ValueError('unknown TAD format {}'.format(fmt))
    if chain:
        tads = liftover_tads(tads.map_chroms(normalize_chrom), chain)
    if targets is not None:
        tads = tads.map_chroms(chrom_mapper(targets))
    log.info('%s (%s): %d boundaries on %d chromosomes%s', path, tads.fmt, len(tads), len(tads.chroms),
             ', {} domains'.format(len(tads.domains)) if tads.domains is not None else '')
    return tads


# ----------------------------------------------------------------------------------------------------------------
# crossing

def _sparse_argmax(values):
    """Sparse table for range-argmax queries over `values` (NaN as -inf)."""
    v = np.where(np.isfinite(values), values, -np.inf)
    table = [np.arange(len(v))]
    k = 1
    while (1 << k) <= len(v):
        prev = table[-1]
        half = 1 << (k - 1)
        left, right = prev[:len(v) - (1 << k) + 1], prev[half:half + len(v) - (1 << k) + 1]
        table.append(np.where(v[left] >= v[right], left, right))
        k += 1
    return v, table


def _range_argmax(v, table, lo, hi):
    """Argmax of v[lo:hi] (hi > lo) for arrays of ranges."""
    length = hi - lo
    k = np.floor(np.log2(np.maximum(length, 1))).astype(int)
    out = np.empty(len(lo), dtype=np.int64)
    for level in np.unique(k):
        sel = k == level
        a = table[level][lo[sel]]
        b = table[level][hi[sel] - (1 << level)]
        out[sel] = np.where(v[a] >= v[b], a, b)
    return out


def _domain_lookup(domains):
    """Per chromosome: starts sorted, prefix max of ends (nested or tiling domains)."""
    out = {}
    if domains is None:
        return out
    for chrom, g in domains.groupby('chrom', sort=False):
        g = g.sort_values('start')
        out[chrom] = (g['start'].values.astype(np.int64), np.maximum.accumulate(g['end'].values.astype(np.int64)))
    return out


def _contained(lookup, chrom, lo, hi):
    """Whether a domain holds [lo, hi]: among domains starting at or before lo, the largest end reaches hi."""
    if chrom not in lookup:
        return np.zeros(len(lo), dtype=bool)
    starts, maxend = lookup[chrom]
    k = np.searchsorted(starts, lo, side='right') - 1
    ok = k >= 0
    out = np.zeros(len(lo), dtype=bool)
    out[ok] = maxend[k[ok]] >= hi[ok]
    return out


def edge_positions(edges):
    """Midpoints of the two ends of union edges (node ids chrom:start-end): chrom, a, b (a <= b), trans edges
    dropped (chrom = '')."""
    def mids(ids):
        parsed = [parse_node_id(x) for x in ids]
        return [p[0] for p in parsed], np.array([(p[1] + p[2]) // 2 for p in parsed], dtype=np.int64)
    ca, ma = mids(edges['node_a'].values)
    cb, mb = mids(edges['node_b'].values)
    chrom = np.where(np.asarray(ca) == np.asarray(cb), ca, '')
    return chrom, np.minimum(ma, mb), np.maximum(ma, mb)


def edge_crossings(edges, tads, min_strength_pct=0.0, max_names=50):
    """Per edge: boundaries crossed (a < p < b), their number, the strongest one, the domain relation."""
    chrom, lo, hi = edge_positions(edges)
    n = len(edges)
    out = pd.DataFrame({'chrom': chrom, 'a_mid': lo, 'b_mid': hi}, index=edges.index)
    n_cross = np.zeros(n, dtype=np.int64)
    first = np.full(n, -1, dtype=np.int64)
    strongest = np.full(n, -1, dtype=np.int64)
    same = np.zeros(n, dtype=bool)
    in_a = np.zeros(n, dtype=bool)
    in_b = np.zeros(n, dtype=bool)
    b = tads.boundaries
    keep = (b['strength_pct'].fillna(1.0) >= min_strength_pct).values if min_strength_pct > 0 else \
        np.ones(len(b), dtype=bool)
    lookup = _domain_lookup(tads.domains)
    for c in pd.unique(chrom):
        if not c:
            continue
        rows = np.flatnonzero(chrom == c)
        g = b[(b.chrom == c).values & keep]
        idx, pos, strength = g.index.values, g['pos'].values.astype(np.int64), g['strength'].values.astype(float)
        if len(pos):
            i = np.searchsorted(pos, lo[rows], side='right')
            j = np.searchsorted(pos, hi[rows], side='left')
            k = np.maximum(j - i, 0)
            n_cross[rows] = k
            has = k > 0
            first[rows[has]] = idx[i[has]]
            if has.any():
                v, table = _sparse_argmax(strength)
                strongest[rows[has]] = idx[_range_argmax(v, table, i[has], j[has])]
        if lookup:
            same[rows] = _contained(lookup, c, lo[rows], hi[rows])
            in_a[rows] = _contained(lookup, c, lo[rows], lo[rows])
            in_b[rows] = _contained(lookup, c, hi[rows], hi[rows])
    if not lookup:
        same = n_cross == 0
        in_a[:] = in_b[:] = True
    relation = np.where(chrom == '', 'trans', np.where(same, 'same domain', np.where(in_a & in_b, 'cross domain',
                                                                                         'outside domains')))
    out['n_boundaries'] = n_cross
    out['crosses'] = n_cross > 0
    out['domain_relation'] = relation
    sb = strongest >= 0
    out['strongest_boundary'] = ''
    out.loc[sb, 'strongest_boundary'] = b['id'].values[strongest[sb]]
    out['strongest_strength'] = np.nan
    out.loc[sb, 'strongest_strength'] = b['strength'].values[strongest[sb]]
    out['strongest_pct'] = np.nan
    out.loc[sb, 'strongest_pct'] = b['strength_pct'].values[strongest[sb]]
    # boundary ids of each crossing edge (capped), for the per-gene and per-boundary tables
    ids = b['id'].values
    pos_all = {c: (g.index.values, g['pos'].values) for c, g in b[keep].groupby('chrom', sort=False)}
    names = np.full(n, '', dtype=object)
    for r in np.flatnonzero(n_cross > 0):
        gi, gp = pos_all[chrom[r]]
        i = np.searchsorted(gp, lo[r], side='right')
        names[r] = ','.join(ids[gi[i:i + min(n_cross[r], max_names)]])
    out['boundaries'] = names
    return out


def condition_names(edges):
    names = [c[len('score_'):] for c in edges.columns if c.startswith('score_')]
    return names[:2] if len(names) >= 2 else ['A', 'B']


def presence(edges, names):
    """Whether each union edge is called in condition A and in B."""
    a, b = names
    if 'significant_{}'.format(a) in edges and 'significant_{}'.format(b) in edges:
        def flag(col):
            return edges[col].astype(str).str.lower().isin(('true', '1')).values
        return flag('significant_{}'.format(a)), flag('significant_{}'.format(b))
    status = edges['status'].values
    return status != 'gained', status != 'lost'


def read_union_edges(prefix_or_path):
    """Union edges of egaDifferential: <prefix>_edges.tsv.gz or a union graphml."""
    path = str(prefix_or_path)
    if path.endswith('.graphml'):
        from .graph import read_graph
        g = read_graph(path)
        rows = [dict(node_a=u, node_b=v, **d) for u, v, d in g.edges(data=True)]
        return pd.DataFrame(rows)
    if not os.path.exists(path):
        for suffix in ('_edges.tsv.gz', '_edges.tsv'):
            if os.path.exists(path + suffix):
                path = path + suffix
                break
    return pd.read_csv(path, sep='\t', low_memory=False)


def annotate_edges(edges, tads, names=None, min_strength_pct=0.0):
    """Union edges with the crossing columns and presence per condition."""
    names = names or condition_names(edges)
    cross = edge_crossings(edges, tads, min_strength_pct)
    out = pd.concat([edges[['node_a', 'node_b', 'status']].reset_index(drop=True),
                     cross.reset_index(drop=True)], axis=1)
    out['distance'] = out['b_mid'] - out['a_mid']
    pa, pb = presence(edges, names)
    out['in_{}'.format(names[0])] = pa
    out['in_{}'.format(names[1])] = pb
    for column in ('log2FC', 'score_{}'.format(names[0]), 'score_{}'.format(names[1])):
        if column in edges:
            out[column] = edges[column].values
    return out


def condition_summary(annotated, names):
    """Fraction of boundary-crossing edges per condition and per status."""
    rows = []
    cis = annotated[annotated['chrom'] != '']
    for name in names:
        part = cis[cis['in_{}'.format(name)]]
        rows.append(dict(set='condition', group=name, edges=len(part), crossing=int(part['crosses'].sum()),
                         fraction_crossing=part['crosses'].mean() if len(part) else np.nan,
                         same_domain=int((part['domain_relation'] == 'same domain').sum()),
                         cross_domain=int((part['domain_relation'] == 'cross domain').sum()),
                         outside_domains=int((part['domain_relation'] == 'outside domains').sum()),
                         median_distance=float(part['distance'].median()) if len(part) else np.nan))
    for status in STATUSES:
        part = cis[cis['status'] == status]
        rows.append(dict(set='status', group=status, edges=len(part), crossing=int(part['crosses'].sum()),
                         fraction_crossing=part['crosses'].mean() if len(part) else np.nan,
                         same_domain=int((part['domain_relation'] == 'same domain').sum()),
                         cross_domain=int((part['domain_relation'] == 'cross domain').sum()),
                         outside_domains=int((part['domain_relation'] == 'outside domains').sum()),
                         median_distance=float(part['distance'].median()) if len(part) else np.nan))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------------------
# statistics

def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * ok.sum() / np.arange(1, ok.sum() + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(q)
    res[order] = np.minimum(ranked, 1.0)
    out[ok] = res
    return out


def fisher_greater(x, row, col, total):
    """One-sided Fisher P (enrichment) of x successes in the crossing row: hypergeometric upper tail."""
    x, row, col, total = (np.asarray(v, dtype=np.int64) for v in (x, row, col, total))
    return stats.hypergeom.sf(x - 1, total, col, row)


def sum_binomial_sf(observed, counts, probs):
    """P(X >= observed) for X = sum of Binomial(counts[k], probs[k])."""
    pmf = np.array([1.0])
    for n, p in zip(counts, probs):
        n = int(n)
        if n == 0:
            continue
        pmf = np.convolve(pmf, stats.binom.pmf(np.arange(n + 1), n, p))
    if observed <= 0:
        return 1.0
    return float(min(1.0, pmf[int(observed):].sum())) if observed < len(pmf) else 0.0


def promoter_ranking(annotated, bait_nodes, genes=None):
    """Per promoter fragment: crossing / non-crossing edges by status, Fisher enrichment of gained and of lost
    edges among the crossing edges, BH over the promoters with a test."""
    cis = annotated[annotated['chrom'] != '']
    rows = []
    long = pd.concat([cis[['node_a', 'status', 'crosses', 'distance']].rename(columns={'node_a': 'node'}),
                      cis[['node_b', 'status', 'crosses', 'distance']].rename(columns={'node_b': 'node'})])
    long = long[long['node'].isin(bait_nodes)]
    counts = long.groupby(['node', 'crosses', 'status']).size().unstack('status', fill_value=0)
    for s in STATUSES:
        if s not in counts:
            counts[s] = 0
    counts = counts.unstack('crosses', fill_value=0)
    table = pd.DataFrame(index=counts.index)
    for s in STATUSES:
        table['{}_crossing'.format(s)] = counts[(s, True)] if (s, True) in counts else 0
        table['{}_within'.format(s)] = counts[(s, False)] if (s, False) in counts else 0
    dist = long.groupby(['node', 'crosses'])['distance'].median().unstack('crosses')
    table['median_distance_crossing'] = dist[True] if True in dist else np.nan
    table['median_distance_within'] = dist[False] if False in dist else np.nan
    crossing = sum(table['{}_crossing'.format(s)] for s in STATUSES)
    within = sum(table['{}_within'.format(s)] for s in STATUSES)
    total = crossing + within
    table['edges'] = total
    table['crossing'] = crossing
    for s in ('gained', 'lost'):
        x = table['{}_crossing'.format(s)]
        col = x + table['{}_within'.format(s)]
        testable = (crossing > 0) & (within > 0) & (col > 0)
        p = np.full(len(table), np.nan)
        p[testable.values] = fisher_greater(x[testable], crossing[testable], col[testable], total[testable])
        a, b = x + 0.5, crossing - x + 0.5
        c, d = table['{}_within'.format(s)] + 0.5, within - table['{}_within'.format(s)] + 0.5
        table['{}_odds_ratio'.format(s)] = (a * d) / (b * c)
        table['{}_p'.format(s)] = p
        table['{}_q'.format(s)] = benjamini_hochberg(p)
        table['{}_rank'.format(s)] = pd.Series(p, index=table.index).rank(method='min')
    table = table.reset_index().rename(columns={'index': 'node'})
    if genes is not None:
        table.insert(1, 'genes', table['node'].map(genes).fillna(''))
    return table.sort_values('gained_p', na_position='last').reset_index(drop=True)


def length_bins(distance, n_bins=20):
    edges = np.unique(np.quantile(np.log10(np.maximum(distance, 1)), np.linspace(0, 1, n_bins + 1)))
    return edges


def boundary_ranking(annotated, tads, n_bins=20, min_strength_pct=0.0):
    """Per boundary (at or above `min_strength_pct`): union edges crossing it, gained and lost among them against
    the expectation from edges of the same length (genome-wide status fraction per length quantile)."""
    cis = annotated[annotated['chrom'] != ''].reset_index(drop=True)
    bins = length_bins(cis['distance'].values, n_bins)
    k = np.clip(np.searchsorted(bins, np.log10(np.maximum(cis['distance'].values, 1)), side='right') - 1, 0,
                max(len(bins) - 2, 0))
    rates = {s: np.array([(cis['status'].values[k == j] == s).mean() if (k == j).any() else 0.0
                          for j in range(max(len(bins) - 1, 1))]) for s in STATUSES}
    b = tads.boundaries
    if min_strength_pct > 0:
        b = b[b['strength_pct'].fillna(1.0) >= min_strength_pct]
    rows = []
    for chrom, g in b.groupby('chrom', sort=False):
        sel = np.flatnonzero(cis['chrom'].values == chrom)
        if not len(sel):
            continue
        lo, hi = cis['a_mid'].values[sel], cis['b_mid'].values[sel]
        st, kk = cis['status'].values[sel], k[sel]
        for r in g.itertuples():
            p = r.pos
            m = (lo < p) & (hi > p)
            n = int(m.sum())
            row = dict(id=r.id, chrom=chrom, pos=int(p), start=int(r.start), end=int(r.end), score=r.score,
                       strength=r.strength, strength_pct=r.strength_pct, strength_rank=r.strength_rank, crossing=n)
            if n:
                per_bin = np.bincount(kk[m], minlength=len(rates['gained']))
                for s in STATUSES:
                    obs = int((st[m] == s).sum())
                    exp = float((per_bin * rates[s]).sum())
                    row['{}_crossing'.format(s)] = obs
                    row['{}_expected'.format(s)] = exp
                    if s != 'shared':
                        row['{}_enrichment'.format(s)] = (obs + 0.5) / (exp + 0.5)
                        row['{}_p'.format(s)] = sum_binomial_sf(obs, per_bin, rates[s])
            rows.append(row)
    out = pd.DataFrame(rows)
    for s in ('gained', 'lost'):
        if '{}_p'.format(s) not in out:
            out['{}_p'.format(s)] = np.nan
        out['{}_q'.format(s)] = benjamini_hochberg(out['{}_p'.format(s)].values)
        out['{}_rank'.format(s)] = out['{}_p'.format(s)].rank(method='min')
    return out.sort_values('gained_p', na_position='last').reset_index(drop=True)


# ----------------------------------------------------------------------------------------------------------------
# genes

def read_regions(path):
    """Named regions (BED: chrom start end name) used to name partner regions."""
    if not path:
        return None
    rows = _read_bed(path, 3)
    return pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [int(r[1]) for r in rows],
                         'end': [int(r[2]) for r in rows],
                         'name': [r[3] if len(r) > 3 else '{}:{}-{}'.format(*r[:3]) for r in rows]})


def read_elements(prefix):
    """Element classes per node and condition from <prefix>_elements.tsv (declared roles, fragment mean)."""
    path = prefix + '_elements.tsv'
    if not os.path.exists(path):
        return {}
    table = pd.read_csv(path, sep='\t', usecols=lambda c: c in ('node', 'condition', 'class_declared_mean'))
    if 'class_declared_mean' not in table:
        return {}
    return {(n, c): k for n, c, k in table[['node', 'condition', 'class_declared_mean']].itertuples(index=False)}


def fmt_mb(pos, digits=3):
    return '{:.{}f} Mb'.format(pos / 1e6, digits)


def strength_word(pct):
    if pct is None or not np.isfinite(pct):
        return ''
    return 'strong' if pct >= STRONG_PCT else 'intermediate' if pct >= 0.5 else 'weak'


def _plural(n, word, plural=None):
    return '{} {}'.format(n, word if n == 1 else (plural or word + 's'))


def _class_phrase(classes, total):
    counts = pd.Series([c for c in classes if isinstance(c, str) and c]).value_counts()
    counts = counts[~counts.index.isin(['unmarked element', 'unmarked promoter'])]
    if not len(counts):
        return ''
    items = ['{} {}'.format(n, cls + ('s' if n != 1 else '')) for cls, n in counts.head(4).items()]
    text = ', '.join(items[:-1]) + ' and ' + items[-1] if len(items) > 1 else items[0]
    return 'the {} partners include {}'.format(total, text)


class GeneCrossing:
    """Crossing information of every promoter fragment of a union graph."""

    def __init__(self, annotated, tads, names, graph=None, elements=None, regions=None, marks=()):
        self.annotated = annotated
        self.tads = tads
        self.names = list(names)
        self.graph = graph
        self.elements = elements or {}
        self.regions = regions
        self.marks = list(marks)
        cis = annotated[annotated['chrom'] != '']
        self.by_node = {}
        for col, other in (('node_a', 'node_b'), ('node_b', 'node_a')):
            part = cis[[col, other, 'status', 'crosses', 'n_boundaries', 'boundaries', 'strongest_boundary',
                        'domain_relation', 'distance']].rename(columns={col: 'node', other: 'partner'})
            for node, g in part.groupby('node', sort=False):
                self.by_node.setdefault(node, []).append(g)
        self.by_node = {n: pd.concat(v, ignore_index=True) for n, v in self.by_node.items()}
        self.bnd = tads.boundaries.set_index('id')

    def node_mid(self, node):
        c, s, e = parse_node_id(node)
        return c, (s + e) // 2

    def region_names(self, chrom, start, end):
        if self.regions is None:
            return []
        r = self.regions
        hit = r[(r.chrom == chrom) & (r.start < end) & (r.end > start)]
        return list(hit['name'])

    def mark_change(self, partners, mark):
        if self.graph is None:
            return np.nan
        a, b = self.names
        values = []
        for n in partners:
            d = self.graph.nodes.get(n, {}) if n in self.graph else {}
            va, vb = d.get('{}_{}'.format(mark, a)), d.get('{}_{}'.format(mark, b))
            if va is None or vb is None:
                continue
            values.append(np.log2((float(vb) + 1.0) / (float(va) + 1.0)))
        return float(np.mean(values)) if values else np.nan

    def status_info(self, node, status):
        g = self.by_node.get(node)
        if g is None:
            return None
        part = g[g['status'] == status]
        cross = part[part['crosses']]
        info = dict(n=len(part), crossing=len(cross), within=len(part) - len(cross),
                    same_domain=int((part['domain_relation'] == 'same domain').sum()))
        if not len(cross):
            return info
        counts = pd.Series(','.join(cross['boundaries']).split(',')).value_counts()
        frequent = counts[counts >= max(1, int(np.ceil(len(cross) / 2.0)))]
        cand = self.bnd.loc[[i for i in frequent.index if i in self.bnd.index]]
        if len(cand):
            best = cand['strength'].astype(float).fillna(-np.inf).idxmax()
            info.update(boundary=best, boundary_pos=int(self.bnd.at[best, 'pos']),
                        boundary_strength=self.bnd.at[best, 'strength'],
                        boundary_pct=self.bnd.at[best, 'strength_pct'],
                        boundary_rank=self.bnd.at[best, 'strength_rank'], boundary_edges=int(counts[best]))
        info['boundaries_crossed'] = len(counts)
        chrom, mid = self.node_mid(node)
        if 'boundary' in info:
            # the partner region: partners of the edges that cross the named boundary
            p = info['boundary_pos']
            pm = np.array([(s + e) // 2 for _, s, e in (parse_node_id(x) for x in cross['partner'])])
            over = (np.minimum(pm, mid) < p) & (np.maximum(pm, mid) > p)
            cross = cross[over]
            info['boundary_edges'] = int(over.sum())
        spans = [parse_node_id(p) for p in cross['partner']]
        info['partner_chrom'] = chrom
        info['partner_start'] = min(s for _, s, _ in spans)
        info['partner_end'] = max(e for _, _, e in spans)
        info['partner_side'] = 'upstream' if np.median([(s + e) / 2 for _, s, e in spans]) < mid else 'downstream'
        names = self.region_names(chrom, info['partner_start'], info['partner_end'])
        in_named = {}
        if self.regions is not None:
            for _, s, e in spans:
                for nm in self.region_names(chrom, s, e):
                    in_named[nm] = in_named.get(nm, 0) + 1
        info['partner_regions'] = names
        info['partner_region_counts'] = in_named
        condition = self.names[1] if status == 'gained' else self.names[0]
        classes = [self.elements.get((p, condition), '') for p in cross['partner']]
        info['partner_classes'] = pd.Series([c for c in classes if c]).value_counts().to_dict()
        info['class_phrase'] = _class_phrase(classes, len(cross))
        info['marks'] = {m: self.mark_change(list(cross['partner']), m) for m in self.marks}
        return info

    def gene_row(self, node, gene=''):
        a, b = self.names
        row = dict(node=node, gene=gene)
        infos = {s: self.status_info(node, s) for s in STATUSES}
        for s, info in infos.items():
            info = info or dict(n=0, crossing=0, within=0)
            row['{}'.format(s)] = info['n']
            row['{}_crossing'.format(s)] = info['crossing']
            row['{}_within'.format(s)] = info['within']
            row['{}_boundary'.format(s)] = info.get('boundary', '')
            row['{}_boundary_pct'.format(s)] = info.get('boundary_pct', np.nan)
            row['{}_boundaries_crossed'.format(s)] = info.get('boundaries_crossed', 0)
            if info.get('partner_start') is not None:
                row['{}_partner_region'.format(s)] = '{}:{}-{}'.format(info['partner_chrom'], info['partner_start'],
                                                                       info['partner_end'])
                row['{}_partner_named'.format(s)] = ','.join(info.get('partner_regions', []))
                row['{}_partner_classes'.format(s)] = ';'.join('{}={}'.format(k, v) for k, v in
                                                               info.get('partner_classes', {}).items())
                for m, v in info.get('marks', {}).items():
                    row['{}_partner_{}_log2'.format(s, m)] = v
        row['sentence'] = self.sentence(gene or node, infos)
        row['_infos'] = infos
        return row

    def boundary_phrase(self, info):
        pct = info.get('boundary_pct', np.nan)
        word = strength_word(pct)
        chrom = info.get('partner_chrom', '')
        text = 'the {}boundary at {}:{}'.format(word + ' ' if word else '', chrom, fmt_mb(info['boundary_pos']))
        rank = info.get('boundary_rank', np.nan)
        if rank is not None and np.isfinite(rank):
            text += ' (strength rank {} of {})'.format(int(rank), len(self.tads))
        return text

    def partner_phrase(self, info):
        text = 'towards {}:{:.3f}-{:.3f} Mb'.format(info['partner_chrom'], info['partner_start'] / 1e6,
                                                   info['partner_end'] / 1e6)
        named = info.get('partner_region_counts', {})
        if named:
            text += ' ({})'.format(', '.join('{}: {} partner{}'.format(k, v, '' if v == 1 else 's')
                                             for k, v in named.items()))
        extras = []
        if info.get('class_phrase'):
            extras.append(info['class_phrase'])
        changes = [(m, v) for m, v in info.get('marks', {}).items() if np.isfinite(v) and abs(v) >= 0.3]
        if changes:
            a, b = self.names
            extras.append('mean log2 change {} / {} at these partners {}'.format(
                b, a, ', '.join('{} {:+.1f}'.format(m, v) for m, v in changes)))
        if extras:
            text += '; ' + '; '.join(extras)
        return text

    def status_phrase(self, status, cond, info):
        n = info['n']
        if not info['crossing']:
            return 'its {} {} contact{} ({} only) stay within its own domain'.format(
                n, status, '' if n == 1 else 's', cond)
        text = 'in {}, {} of its {} {} contacts cross a boundary'.format(cond, info['crossing'], n, status)
        if 'boundary' in info:
            text += '; {} of them cross {}, {}'.format(info['boundary_edges'], self.boundary_phrase(info),
                                                      self.partner_phrase(info))
        if info['within']:
            text += '; {} stay{} within its own domain'.format(info['within'], 's' if info['within'] == 1 else '')
        return text

    def sentence(self, gene, infos):
        a, b = self.names
        parts = []
        for status, cond in (('gained', b), ('lost', a)):
            info = infos.get(status)
            if info and info['n']:
                parts.append(self.status_phrase(status, cond, info))
        shared = infos.get('shared')
        if shared and shared['n']:
            parts.append('{} of its {} shared contacts cross a boundary'.format(shared['crossing'], shared['n']))
        if not parts:
            return ''
        return '{}: {}.'.format(gene, '; '.join(parts))

    def table(self, nodes_genes):
        rows = [self.gene_row(n, g) for n, g in nodes_genes]
        return pd.DataFrame([{k: v for k, v in r.items() if k != '_infos'} for r in rows])


def bait_genes(graph):
    out = {}
    for n, d in graph.nodes(data=True):
        if d.get('bait'):
            genes = d.get('genes', '')
            out[n] = '' if genes is None or (isinstance(genes, float) and np.isnan(genes)) else str(genes)
    return out


def graph_marks(graph, names):
    marks = set()
    for _, d in list(graph.nodes(data=True))[:200]:
        for k in d:
            for name in names:
                if k.endswith('_' + name):
                    marks.add(k[:-len(name) - 1])
    return sorted(marks - {'score', 'significant', 'enrichment'})


class Analysis:
    """Everything of one comparison: the annotated union edges, per-gene crossing, the rankings."""

    def __init__(self, prefix, tads_path, fmt='auto', chain=None, window=None, regions=None, min_strength_pct=0.0,
                 names=None, graph=None):
        from .graph import read_graph
        self.prefix = prefix
        graph_path = prefix + '_union.graphml'
        self.graph = graph if graph is not None else (read_graph(graph_path) if os.path.exists(graph_path) else None)
        self.edges = read_union_edges(prefix)
        self.names = list(names) if names else condition_names(self.edges)
        targets = set(self.graph.nodes[n]['chrom'] for n in self.graph.nodes) if self.graph is not None else \
            set(edge_positions(self.edges)[0]) - {''}
        self.tads = read_tads(tads_path, fmt, targets=targets, chain=chain, window=window)
        self.min_strength_pct = min_strength_pct
        self.annotated = annotate_edges(self.edges, self.tads, self.names, min_strength_pct)
        self.genes = bait_genes(self.graph) if self.graph is not None else {}
        marks = graph_marks(self.graph, self.names) if self.graph is not None else []
        self.crossing = GeneCrossing(self.annotated, self.tads, self.names, self.graph, read_elements(prefix),
                                     read_regions(regions) if isinstance(regions, str) else regions, marks)

    def gene_nodes(self, gene):
        return [n for n, g in self.genes.items() if gene in str(g).split(',')]

    def gene_table(self):
        return self.crossing.table(sorted(self.genes.items()))

    def promoters(self):
        return promoter_ranking(self.annotated, set(self.genes), self.genes)

    def promoter_table(self):
        """promoters(), computed once."""
        if getattr(self, '_promoters', None) is None:
            self._promoters = self.promoters()
        return self._promoters

    def boundaries(self):
        return boundary_ranking(self.annotated, self.tads, min_strength_pct=self.min_strength_pct)

    def summary(self):
        return condition_summary(self.annotated, self.names)

    def write(self, out_prefix, genes_table=True):
        from ._cli import write_table
        write_table(self.annotated, out_prefix + '_edges.tsv.gz')
        write_table(self.summary(), out_prefix + '_summary.tsv')
        promoters = self.promoter_table()
        write_table(promoters, out_prefix + '_promoters.tsv')
        bounds = self.boundaries()
        write_table(bounds, out_prefix + '_boundaries.tsv')
        if genes_table:
            write_table(self.gene_table(), out_prefix + '_genes.tsv')
        return promoters, bounds


def run_for_prefix(prefix, tads_path, out_prefix=None, **kwargs):
    """Used by egaDifferential --tads: the crossing tables of a finished run (<prefix>_tad_crossing_*)."""
    analysis = Analysis(prefix, tads_path, **kwargs)
    analysis.write(out_prefix or prefix + '_tad_crossing')
    return analysis
