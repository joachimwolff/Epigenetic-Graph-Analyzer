"""Domain boundaries from capture contacts, CTCF / cohesin binding and CTCF motif orientation.

Three evidence sources per fixed genomic bin:

* capture insulation (`capture_insulation`): the number (and summed weight) of graph edges whose span covers the
  bin centre, against a distance-matched expectation that keeps every edge at its bait and keeps its length but lets
  it point to either side (mirror null: expected = (observed + coverage of the mirrored edges) / 2). Bait density,
  bait degree and the length distribution of every bait are therefore part of the expectation; only the side to which
  the edges of a bait point is tested. A second expectation spreads every bait's edges with the length distribution
  of all edges of the graph (`global` null). Insulation = log2((observed + c) / (expected + c)); a bin is covered when
  its expected count reaches `min_expected`. Boundaries: local minima with a prominence threshold (`local_minima`).
* CTCF / cohesin sites (`call_summits`, `scan_motif`, `ctcf_sites`): summits of a CTCF ChIP signal (100 bp windows,
  enrichment over the genome-wide mean and over a 20 kb local mean), cohesin (RAD21 / SMC1A / SMC3) enrichment at the
  summit, best CTCF motif hit (JASPAR PFM, log-odds, P-value threshold from the exact score distribution) within
  +-150 bp on either strand: '+' (forward, JASPAR MA0139 orientation), '-' (reverse) or 'both'.
* loops: graph edges whose left end carries a forward and whose right end a reverse cohesin-bound CTCF site
  (convergent pair, `classify_loops`).

`bin_features` turns the three into per-bin features (site strength and counts, the orientation pattern of a
boundary: reverse sites at and left of the bin, forward sites at and right of it, capture insulation and its local
dip, convergent loop anchors and loops spanning the bin). `LogisticModel` is a plain L2-penalised logistic regression
on z-scored features; `peak_calls` turns a per-bin score into boundary calls (local maxima above a threshold chosen on
the training data, `best_threshold`). Evaluation: `match_boundaries` (precision, recall, F1 at a tolerance in bins),
`average_precision`, `variation_of_information` and `local_rand` between domain partitions, and the diamond
insulation of a Hi-C matrix (`diamond_insulation`) as reference.
"""

import gzip
import logging

import numpy as np
import pandas as pd
from scipy import ndimage, optimize, signal

log = logging.getLogger(__name__)

BASES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


# ----------------------------------------------------------------------------------------------------------------
# bins

def bin_index(bins):
    """Per chromosome: (first row, number of bins, bin size) of a bins table (chrom, start, end) in order."""
    out = {}
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        size = int(bins.end.values[idx[0]] - bins.start.values[idx[0]])
        out[chrom] = (int(idx[0]), len(idx), size)
    return out


def chrom_codes(bins):
    return pd.factorize(bins.chrom.values)[0]


# ----------------------------------------------------------------------------------------------------------------
# capture insulation

def _add_spans(diff, lo, hi, weight, res):
    """Add weight to every bin whose centre lies strictly inside (lo, hi) (difference array of length n + 1)."""
    n = len(diff) - 1
    first = np.floor((lo - res / 2.0) / res).astype(np.int64) + 1
    last = np.ceil((hi - res / 2.0) / res).astype(np.int64) - 1
    first = np.clip(first, 0, n)
    last = np.clip(last, -1, n - 1)
    ok = last >= first
    np.add.at(diff, first[ok], weight[ok])
    np.add.at(diff, last[ok] + 1, -weight[ok])


def edge_table(nodes, edges, position=None):
    """Edges as (chrom, a, b, anchor_a, anchor_b, weight) from egaCreateGraph node and edge tables; a and b are node
    midpoints (or `position(node_ids) -> (chrom, mid)` for lifted coordinates); anchors are baits. Trans edges and
    unmapped nodes are dropped."""
    nd = nodes.set_index('node')
    if position is None:
        chrom = nd.chrom
        mid = (nd.start + nd.end) // 2
    else:
        chrom, mid = position(nd.index.values)
        chrom = pd.Series(chrom, index=nd.index)
        mid = pd.Series(mid, index=nd.index)
    e = edges[edges.node_a.isin(nd.index) & edges.node_b.isin(nd.index)]
    out = pd.DataFrame({'chrom': chrom.loc[e.node_a.values].values,
                        'chrom_b': chrom.loc[e.node_b.values].values,
                        'a': mid.loc[e.node_a.values].values.astype(np.int64),
                        'b': mid.loc[e.node_b.values].values.astype(np.int64),
                        'anchor_a': nd.bait.loc[e.node_a.values].values.astype(bool),
                        'anchor_b': nd.bait.loc[e.node_b.values].values.astype(bool),
                        'weight': e.weight.values.astype(float),
                        'node_a': e.node_a.values, 'node_b': e.node_b.values})
    keep = (out.chrom == out.chrom_b) & (out.a >= 0) & (out.b >= 0) & (out.a != out.b)
    return out[keep].drop(columns='chrom_b').reset_index(drop=True)


def capture_insulation(bins, edges, min_expected=5.0, pseudocount=1.0):
    """Capture insulation per bin (see the module docstring).

    bins: chrom, start, end (fixed size per chromosome, in order). edges: chrom, a, b, anchor_a, anchor_b, weight.
    Returns a DataFrame aligned to bins: observed, expected (mirror null), expected_global, observed_w, expected_w,
    baits (number of bait positions in the bin), insulation, insulation_w, insulation_global, covered.
    """
    n_all = len(bins)
    cols = {k: np.zeros(n_all) for k in ('observed', 'expected', 'expected_global', 'observed_w', 'expected_w',
                                         'baits')}
    weights = edges.weight.values.astype(float)
    pseudo_w = float(np.median(weights)) * pseudocount if len(weights) else pseudocount
    lengths_all = np.abs(edges.b.values - edges.a.values)
    for chrom, (first, n, res) in bin_index(bins).items():
        e = edges[edges.chrom == chrom]
        if e.empty:
            continue
        a, b = e.a.values.astype(float), e.b.values.astype(float)
        w = e.weight.values.astype(float)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        anc_a, anc_b = e.anchor_a.values, e.anchor_b.values
        none = ~anc_a & ~anc_b
        anc_a, anc_b = anc_a | none, anc_b | none
        share = 1.0 / (anc_a.astype(float) + anc_b.astype(float))
        for kind, wt in (('', np.ones(len(e))), ('_w', w)):
            obs = np.zeros(n + 1)
            _add_spans(obs, lo, hi, wt, res)
            flipped = np.zeros(n + 1)
            for anchor, other, use in ((a, b, anc_a), (b, a, anc_b)):
                mirror = 2 * anchor - other
                _add_spans(flipped, np.minimum(anchor, mirror)[use], np.maximum(anchor, mirror)[use],
                           (wt * share)[use], res)
            o = np.cumsum(obs)[:n]
            f = np.cumsum(flipped)[:n]
            key = 'observed' + kind
            cols[key][first:first + n] = o
            cols['expected' + kind][first:first + n] = 0.5 * (o + f)
        # global null: every anchor's edges spread with the length distribution of all edges (both sides)
        density = np.zeros(n)
        for anchor, use in ((a, anc_a), (b, anc_b)):
            k = np.clip((anchor[use] // res).astype(np.int64), 0, n - 1)
            np.add.at(density, k, share[use])
        uniq = np.unique(np.concatenate([a[e.anchor_a.values], b[e.anchor_b.values]]))
        cols['baits'][first:first + n] = np.bincount(np.clip((uniq // res).astype(np.int64), 0, n - 1), minlength=n)
        max_d = int(min(n, np.ceil(lengths_all.max() / res) + 1))
        offsets = np.arange(-max_d, max_d + 1)
        dist = np.maximum(np.abs(offsets) * res, 0)
        srt = np.sort(lengths_all)
        survival = 1.0 - np.searchsorted(srt, dist, side='right') / len(srt)
        kernel = 0.5 * survival
        kernel[max_d] = 0.5
        conv = signal.fftconvolve(density, kernel, mode='full')[max_d:max_d + n]
        cols['expected_global'][first:first + n] = np.clip(conv, 0, None)
    out = pd.DataFrame(cols)
    out['insulation'] = np.log2((out.observed + pseudocount) / (out.expected + pseudocount))
    out['insulation_w'] = np.log2((out.observed_w + pseudo_w) / (out.expected_w + pseudo_w))
    out['insulation_global'] = np.log2((out.observed + pseudocount) / (out.expected_global + pseudocount))
    out['covered'] = out.expected >= min_expected
    for col in ('insulation', 'insulation_w', 'insulation_global'):
        out.loc[~out.covered, col] = np.nan
    return out


def local_minima(values, codes, prominence=0.3, distance=1):
    """Indices of local minima of a per-bin track with at least `prominence` (same units as the track), searched
    within each chromosome and within runs of finite values."""
    values = np.asarray(values, dtype=float)
    out = []
    for c in np.unique(codes):
        idx = np.flatnonzero(codes == c)
        v = values[idx]
        finite = np.isfinite(v)
        if not finite.any():
            continue
        # runs of finite values
        edges = np.flatnonzero(np.diff(np.concatenate([[0], finite.astype(int), [0]])))
        for s, e in zip(edges[::2], edges[1::2]):
            if e - s < 3:
                continue
            peaks, _ = signal.find_peaks(-v[s:e], prominence=prominence, distance=distance)
            out.append(idx[s + peaks])
    return np.sort(np.concatenate(out)) if out else np.array([], dtype=np.int64)


# ----------------------------------------------------------------------------------------------------------------
# Hi-C insulation (reference)

def band_from_dense(matrix, max_d):
    """band[i, d] = M[i, i + d] for d = 0..max_d (NaN outside the matrix)."""
    n = matrix.shape[0]
    band = np.full((n, max_d + 1), np.nan)
    for d in range(max_d + 1):
        if d < n:
            band[:n - d, d] = np.diagonal(matrix, d)
    return band


def band_from_pixels(n, i, j, v, max_d, valid=None):
    """Band from sparse upper-triangle pixels (i <= j); missing pixels are 0 between valid bins, NaN otherwise."""
    band = np.zeros((n, max_d + 1))
    d = j - i
    k = (d >= 0) & (d <= max_d)
    band[i[k], d[k]] = v[k]
    rows = np.arange(n)[:, None]
    cols = rows + np.arange(max_d + 1)[None, :]
    band[cols >= n] = np.nan
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        ok = valid[:, None] & valid[np.clip(cols, 0, n - 1)]
        band[~ok] = np.nan
    return band


def diamond_insulation(band, window, min_valid=0.5):
    """Diamond insulation per bin i: mean of M[a, b] for a in [i - window, i - 1], b in [i + 1, i + window];
    log2 of the ratio to the chromosome mean; NaN where fewer than `min_valid` of the pixels are valid."""
    n, width = band.shape
    if width < 2 * window + 1:
        raise ValueError('band too narrow for the window')
    total = np.zeros(n)
    count = np.zeros(n)
    i = np.arange(n)
    for d in range(2, 2 * window + 1):
        col = band[:, d]
        valid = np.isfinite(col)
        cs = np.concatenate([[0.0], np.cumsum(np.where(valid, col, 0.0))])
        cc = np.concatenate([[0.0], np.cumsum(valid)])
        lo = np.maximum(i - window, i + 1 - d)
        hi = np.minimum(i - 1, i + window - d)
        lo = np.clip(lo, 0, n)
        hi = np.clip(hi + 1, 0, n)
        ok = hi > lo
        total[ok] += cs[hi[ok]] - cs[lo[ok]]
        count[ok] += cc[hi[ok]] - cc[lo[ok]]
    n_pix = window * window
    with np.errstate(invalid='ignore', divide='ignore'):
        ins = total / count
    ins[(count < min_valid * n_pix) | (i < window) | (i >= n - window)] = np.nan
    mean = np.nanmean(ins) if np.isfinite(ins).any() else np.nan
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.log2(ins / mean)
    out[~np.isfinite(out)] = np.nan
    return out


# ----------------------------------------------------------------------------------------------------------------
# CTCF / cohesin sites and motif orientation

def bigwig_windows(handle, chrom, size, window=100):
    """Mean signal per window of a pyBigWig handle (0 where no data)."""
    n = int(size // window)
    if n == 0:
        return np.zeros(0)
    out = np.zeros(n)
    step = 2000 * window * 100        # chunks of 20 Mb (100 bp windows)
    for lo in range(0, n * window, step):
        hi = min(lo + step, n * window)
        vals = np.nan_to_num(handle.values(chrom, lo, hi, numpy=True)).astype(float)
        out[lo // window:hi // window] = vals.reshape(-1, window).mean(axis=1)
    return out


def call_summits(values, genome_mean, window=100, min_enrichment=5.0, min_local=3.0, local_bp=10000,
                 summit_bp=500, smooth=3):
    """Summits of a windowed ChIP signal: local maxima (+-summit_bp) of the `smooth`-window running mean that are at
    least `min_enrichment` times the genome-wide mean and `min_local` times the mean over +-local_bp.
    Returns (window indices, enrichment over the genome-wide mean)."""
    v = np.asarray(values, dtype=float)
    if len(v) == 0:
        return np.array([], dtype=np.int64), np.array([])
    s = ndimage.uniform_filter1d(v, smooth, mode='constant')
    local = ndimage.uniform_filter1d(v, 2 * (local_bp // window) + 1, mode='nearest')
    mx = ndimage.maximum_filter1d(s, 2 * (summit_bp // window) + 1, mode='constant')
    cand = (s >= mx) & (s > 0) & (s >= min_enrichment * genome_mean) & (s >= min_local * local)
    idx = np.flatnonzero(cand)
    if len(idx) > 1:   # plateaus: keep the first of each run of candidates closer than the summit window
        keep = np.concatenate([[True], np.diff(idx) > summit_bp // window])
        idx = idx[keep]
    return idx, s[idx] / genome_mean


def signal_at(values, idx, genome_mean, flank=2, local_bp=10000, window=100):
    """Maximum of a windowed signal (300 bp running mean) within +-flank windows of each index, over the genome-wide
    mean and over the mean within +-local_bp. Returns (genome ratio, local ratio)."""
    v = np.asarray(values, dtype=float)
    if len(v) == 0:
        return np.zeros(len(idx)), np.zeros(len(idx))
    k = np.clip(idx, 0, len(v) - 1)
    mx = ndimage.maximum_filter1d(ndimage.uniform_filter1d(v, 3, mode='constant'), 2 * flank + 1, mode='constant')[k]
    local = ndimage.uniform_filter1d(v, 2 * (local_bp // window) + 1, mode='nearest')[k]
    with np.errstate(invalid='ignore', divide='ignore'):
        return mx / genome_mean, np.where(local > 0, mx / np.maximum(local, 1e-12), 0.0)


def read_jaspar(path, motif_id):
    """Counts matrix (L x 4, columns A C G T) of one motif of a JASPAR-format file (prefix match on the ID)."""
    rows, found = [], False
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if found:
                    break
                found = line[1:].split()[0].startswith(motif_id)
                continue
            if found and line.strip():
                nums = line.split('[')[1].split(']')[0].split()
                rows.append([float(x) for x in nums])
    if not found or len(rows) != 4:
        raise ValueError('motif {} not found in {}'.format(motif_id, path))
    return np.array(rows).T


def log_odds(counts, pseudocount=0.8, background=(0.25, 0.25, 0.25, 0.25)):
    counts = np.asarray(counts, dtype=float)
    freq = (counts + pseudocount * np.asarray(background)) / (counts.sum(axis=1, keepdims=True) + pseudocount)
    return np.log2(freq / np.asarray(background))


def score_threshold(pwm, pvalue=1e-4, background=(0.25, 0.25, 0.25, 0.25), scale=100):
    """Smallest log-odds score whose upper tail probability under the background is <= pvalue (exact distribution of
    the integer-rounded scores)."""
    q = np.round(pwm * scale).astype(np.int64)
    lo = q.min(axis=1)
    dist = np.array([1.0])
    offset = 0
    for k in range(q.shape[0]):
        width = int(q[k].max() - lo[k]) + 1
        step = np.zeros(width)
        for b in range(4):
            step[q[k, b] - lo[k]] += background[b]
        dist = np.convolve(dist, step)
        offset += lo[k]
    tail = np.cumsum(dist[::-1])[::-1]
    ok = np.flatnonzero((tail <= pvalue) & (dist > 0))
    if len(ok) == 0:
        return float(q.max(axis=1).sum()) / scale + 1
    # the smallest attainable (rounded) score with tail <= pvalue; minus the largest rounding error of a sum so that
    # unrounded scores of that atom pass
    return (ok[0] + offset) / scale - 0.5 * q.shape[0] / scale


def encode(seq):
    table = np.full(256, 4, dtype=np.uint8)
    for base, code in BASES.items():
        table[ord(base)] = code
        table[ord(base.lower())] = code
    return table[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]


def scan_motif(codes, pwm):
    """Best forward and reverse-strand scores of a PWM in each row of a code array (n_seq x width; 4 = N).
    Returns (best_forward, best_reverse) arrays."""
    codes = np.atleast_2d(codes)
    n, width = codes.shape
    length = pwm.shape[0]
    if width < length:
        return np.full(n, -np.inf), np.full(n, -np.inf)
    ext = np.vstack([pwm.T, np.full(length, -1e9)])            # 5 x L
    rc = pwm[::-1, ::-1]                                        # reverse complement (A<->T, C<->G)
    ext_rc = np.vstack([rc.T, np.full(length, -1e9)])
    m = width - length + 1
    fwd = np.zeros((n, m))
    rev = np.zeros((n, m))
    for k in range(length):
        sub = codes[:, k:k + m]
        fwd += ext[sub, k]
        rev += ext_rc[sub, k]
    return fwd.max(axis=1), rev.max(axis=1)


def orientation(best_fwd, best_rev, threshold):
    f = np.asarray(best_fwd) >= threshold
    r = np.asarray(best_rev) >= threshold
    out = np.full(len(f), '', dtype=object)
    out[f & ~r] = '+'
    out[r & ~f] = '-'
    out[f & r] = 'both'
    return out


def iter_fasta(path, wanted):
    """(chrom, sequence) for the wanted chromosomes of a (gzipped) FASTA, streamed."""
    opener = gzip.open if str(path).endswith('.gz') else open
    chrom, chunks = None, []
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if chrom in wanted:
                    yield chrom, ''.join(chunks)
                chrom, chunks = line[1:].split()[0], []
            elif chrom in wanted:
                chunks.append(line.rstrip())
        if chrom in wanted:
            yield chrom, ''.join(chunks)


def annotate_motifs(sites, fasta, pwm, threshold, flank=150):
    """Adds motif_fwd, motif_rev (best scores) and orientation ('+', '-', 'both', '') to a sites table (chrom,
    summit)."""
    sites = sites.copy()
    sites['motif_fwd'] = -np.inf
    sites['motif_rev'] = -np.inf
    wanted = set(sites.chrom)
    for chrom, seq in iter_fasta(fasta, wanted):
        idx = np.flatnonzero(sites.chrom.values == chrom)
        if len(idx) == 0:
            continue
        codes = encode(seq)
        starts = np.clip(sites.summit.values[idx] - flank, 0, max(len(codes) - 2 * flank - 1, 0))
        win = codes[starts[:, None] + np.arange(2 * flank + 1)[None, :]]
        f, r = scan_motif(win, pwm)
        sites.loc[sites.index[idx], 'motif_fwd'] = f
        sites.loc[sites.index[idx], 'motif_rev'] = r
    sites['orientation'] = orientation(sites.motif_fwd.values, sites.motif_rev.values, threshold)
    return sites


def read_windows(path, chrom_sizes, window=100):
    """Windowed signal of a bigWig on the given chromosomes and its genome-wide mean (missing data count as 0)."""
    import pyBigWig
    h = pyBigWig.open(path)
    chroms = h.chroms()
    vals = {c: bigwig_windows(h, c, min(s, chroms[c]), window) if c in chroms else np.zeros(int(s // window))
            for c, s in chrom_sizes.items()}
    h.close()
    total = float(sum(len(v) for v in vals.values())) * window
    return vals, sum(v.sum() for v in vals.values()) * window / total


def ctcf_sites(ctcf_paths, cohesin_paths, chrom_sizes, window=100, min_enrichment=10.0, min_local=4.0,
               min_cohesin=3.0, min_cohesin_local=2.0, windows=None):
    """CTCF summits (signals of several files are averaged after division by their genome-wide means), CTCF
    enrichment and cohesin enrichment (mean over the cohesin files of the maximum near the summit / genome mean).
    `windows` (path -> (values per chromosome, genome mean)) reuses signals read with `read_windows`.
    Cohesin-bound: cohesin >= min_cohesin x its genome mean and >= min_cohesin_local x its mean over +-10 kb.
    Returns chrom, summit (bp), ctcf, cohesin, cohesin_local, cohesin_bound."""
    windows = dict(windows or {})
    for p in list(ctcf_paths) + list(cohesin_paths):
        if p not in windows:
            windows[p] = read_windows(p, chrom_sizes, window)
    rows = []
    for chrom in chrom_sizes:
        c_sig = np.mean([windows[p][0][chrom] / windows[p][1] for p in ctcf_paths], axis=0)
        idx, enr = call_summits(c_sig, 1.0, window, min_enrichment, min_local)
        if len(cohesin_paths):
            at = [signal_at(windows[p][0][chrom], idx, windows[p][1]) for p in cohesin_paths]
            coh = np.mean([a[0] for a in at], axis=0)
            coh_local = np.mean([a[1] for a in at], axis=0)
        else:
            coh = coh_local = np.full(len(idx), np.nan)
        rows.append(pd.DataFrame({'chrom': chrom, 'summit': idx * window + window // 2, 'ctcf': enr,
                                  'cohesin': coh, 'cohesin_local': coh_local}))
    out = pd.concat(rows, ignore_index=True)
    out['cohesin_bound'] = (out.cohesin >= min_cohesin) & (out.cohesin_local >= min_cohesin_local)
    return out


# ----------------------------------------------------------------------------------------------------------------
# loops

def classify_loops(ends, sites, pad=1000, cohesin_only=True, return_ambiguous=False):
    """Orientation class of each edge from the CTCF sites at its two ends.

    ends: chrom, left_start, left_end, right_start, right_end (left = smaller coordinate). A site belongs to an end
    when its summit lies within the end interval +- pad. Classes: 'convergent' (forward site at the left end, reverse
    at the right), 'divergent', 'tandem_forward', 'tandem_reverse' (an end with 'both' counts as either), 'one_end',
    'none'. Returns the class and the summed CTCF enrichment at the two ends (and, with return_ambiguous, whether an
    end carries both orientations, which makes the class depend on the precedence convergent > divergent > tandem).
    """
    s = sites[sites.cohesin_bound] if cohesin_only else sites
    s = s[s.orientation != '']
    by_chrom = {c: g.sort_values('summit') for c, g in s.groupby('chrom')}

    def end_info(chrom, start, end):
        fwd = np.zeros(len(chrom), dtype=bool)
        rev = np.zeros(len(chrom), dtype=bool)
        strength = np.zeros(len(chrom))
        for c in pd.unique(chrom):
            k = np.flatnonzero(chrom == c)
            if c not in by_chrom:
                continue
            g = by_chrom[c]
            pos = g.summit.values
            cs_f = np.concatenate([[0], np.cumsum(g.orientation.isin(['+', 'both']).values)])
            cs_r = np.concatenate([[0], np.cumsum(g.orientation.isin(['-', 'both']).values)])
            cs_s = np.concatenate([[0], np.cumsum(g.ctcf.values)])
            lo = np.searchsorted(pos, start[k] - pad, side='left')
            hi = np.searchsorted(pos, end[k] + pad, side='right')
            fwd[k] = cs_f[hi] - cs_f[lo] > 0
            rev[k] = cs_r[hi] - cs_r[lo] > 0
            strength[k] = cs_s[hi] - cs_s[lo]
        return fwd, rev, strength

    chrom = ends.chrom.values
    lf, lr, ls = end_info(chrom, ends.left_start.values, ends.left_end.values)
    rf, rr, rs = end_info(chrom, ends.right_start.values, ends.right_end.values)
    cls = np.full(len(ends), 'none', dtype=object)
    left_any, right_any = lf | lr, rf | rr
    cls[left_any ^ right_any] = 'one_end'
    both = left_any & right_any
    cls[both & lr & rr] = 'tandem_reverse'
    cls[both & lf & rf] = 'tandem_forward'
    cls[both & lr & rf] = 'divergent'
    cls[both & lf & rr] = 'convergent'      # takes precedence when an end carries both orientations
    if return_ambiguous:
        return cls, ls + rs, (lf & lr) | (rf & rr)
    return cls, ls + rs


# ----------------------------------------------------------------------------------------------------------------
# features

def _count_in_bins(n, positions, res, weights=None):
    k = np.clip((np.asarray(positions) // res).astype(np.int64), 0, n - 1)
    out = np.zeros(n)
    np.add.at(out, k, 1.0 if weights is None else weights)
    return out


def _max_in_bins(n, positions, res, values):
    k = np.clip((np.asarray(positions) // res).astype(np.int64), 0, n - 1)
    out = np.zeros(n)
    np.maximum.at(out, k, values)
    return out


def _window_sum(x, left, right):
    """sum of x[i - left .. i + right] (clipped at the ends)."""
    cs = np.concatenate([[0.0], np.cumsum(x)])
    i = np.arange(len(x))
    lo = np.clip(i - left, 0, len(x))
    hi = np.clip(i + right + 1, 0, len(x))
    return cs[hi] - cs[lo]


FEATURE_SETS = {
    # bait / edge density only (no information on the side to which edges point)
    'density': ['covered', 'log_expected', 'log_expected_global', 'bait_n'],
    'insulation': ['insulation', 'insulation_w', 'insulation_global', 'insulation_dip'],
    'ctcf': ['ctcf_max', 'ctcf_n', 'cb_n', 'cohesin_max', 'cb_window'],
    'orientation': ['fwd_n', 'rev_n', 'rev_left', 'fwd_right', 'flip', 'antiflip'],
    'loops': ['loop_left_anchor', 'loop_right_anchor', 'loop_cross', 'loop_strength'],
}
METHODS = {
    'density': FEATURE_SETS['density'],
    'capture': FEATURE_SETS['density'] + FEATURE_SETS['insulation'],
    'ctcf_cohesin': FEATURE_SETS['ctcf'],
    'ctcf_orientation': FEATURE_SETS['ctcf'] + FEATURE_SETS['orientation'],
    'ctcf_orientation_density': FEATURE_SETS['ctcf'] + FEATURE_SETS['orientation'] + FEATURE_SETS['density'],
    'combined': FEATURE_SETS['ctcf'] + FEATURE_SETS['orientation'] + FEATURE_SETS['density'] +
    FEATURE_SETS['insulation'] + FEATURE_SETS['loops'],
}


def bin_features(bins, sites=None, insulation=None, loops=None, context_bp=50000):
    """Per-bin features (see FEATURE_SETS). sites: chrom, summit, ctcf, cohesin, cohesin_bound, orientation.
    insulation: output of capture_insulation. loops: chrom, left_mid, right_mid, cls, strength (convergent loops are
    used). context_bp sets the window (in bins: max(1, context_bp / bin size)) of the orientation pattern, the site
    sum and the insulation dip."""
    out = pd.DataFrame(index=bins.index)
    cols = {c: np.zeros(len(bins)) for s in FEATURE_SETS.values() for c in s}
    for chrom, (first, n, res) in bin_index(bins).items():
        sl = slice(first, first + n)
        w = max(1, int(round(context_bp / res)))
        if sites is not None:
            g = sites[sites.chrom == chrom]
            cb = g[g.cohesin_bound]
            cols['ctcf_max'][sl] = np.log2(1 + _max_in_bins(n, g.summit.values, res, g.ctcf.values))
            cols['ctcf_n'][sl] = np.log1p(_count_in_bins(n, g.summit.values, res))
            cols['cb_n'][sl] = np.log1p(_count_in_bins(n, cb.summit.values, res))
            coh = np.nan_to_num(g.cohesin.values) if len(g) else np.zeros(0)
            cols['cohesin_max'][sl] = np.log2(1 + _max_in_bins(n, g.summit.values, res, coh))
            cols['cb_window'][sl] = np.log1p(_window_sum(_count_in_bins(n, cb.summit.values, res, cb.ctcf.values),
                                                         w, w))
            fwd = _count_in_bins(n, cb.summit.values, res, cb.orientation.isin(['+', 'both']).values.astype(float))
            rev = _count_in_bins(n, cb.summit.values, res, cb.orientation.isin(['-', 'both']).values.astype(float))
            cols['fwd_n'][sl] = np.log1p(fwd)
            cols['rev_n'][sl] = np.log1p(rev)
            rl, fr = _window_sum(rev, w, 0), _window_sum(fwd, 0, w)
            fl, rr = _window_sum(fwd, w, 0), _window_sum(rev, 0, w)
            cols['rev_left'][sl] = np.log1p(rl)
            cols['fwd_right'][sl] = np.log1p(fr)
            cols['flip'][sl] = np.log1p(rl) * np.log1p(fr)
            cols['antiflip'][sl] = np.log1p(fl) * np.log1p(rr)
        if insulation is not None:
            ins = insulation.iloc[sl]
            v = ins.insulation.values
            cov = ins.covered.values.astype(bool)
            filled = np.where(cov, v, 0.0)
            cols['insulation'][sl] = filled
            cols['insulation_w'][sl] = np.where(cov, ins.insulation_w.values, 0.0)
            cols['insulation_global'][sl] = np.where(cov, np.nan_to_num(ins.insulation_global.values), 0.0)
            cols['covered'][sl] = cov.astype(float)
            cols['log_expected'][sl] = np.log1p(ins.expected.values)
            cols['log_expected_global'][sl] = np.log1p(ins.expected_global.values)
            if 'baits' in ins:
                cols['bait_n'][sl] = np.log1p(ins.baits.values)
            # dip: insulation minus the mean of the flanks (bins i-2w..i-1 and i+1..i+2w)
            flank_sum = _window_sum(filled, 2 * w, 2 * w) - filled
            flank_n = _window_sum(cov.astype(float), 2 * w, 2 * w) - cov
            with np.errstate(invalid='ignore', divide='ignore'):
                dip = np.where(cov & (flank_n > 0), filled - flank_sum / np.maximum(flank_n, 1), 0.0)
            cols['insulation_dip'][sl] = dip
        if loops is not None:
            g = loops[(loops.chrom == chrom) & (loops.cls == 'convergent')]
            cols['loop_left_anchor'][sl] = np.log1p(_count_in_bins(n, g.left_mid.values, res))
            cols['loop_right_anchor'][sl] = np.log1p(_count_in_bins(n, g.right_mid.values, res))
            st = _count_in_bins(n, g.left_mid.values, res, g.strength.values) + \
                _count_in_bins(n, g.right_mid.values, res, g.strength.values)
            cols['loop_strength'][sl] = np.log1p(st)
            diff = np.zeros(n + 1)
            if len(g):
                _add_spans(diff, g.left_mid.values.astype(float) + res, g.right_mid.values.astype(float) - res,
                           np.ones(len(g)), res)
            cols['loop_cross'][sl] = np.log1p(np.cumsum(diff)[:n])
    for c, v in cols.items():
        out[c] = v
    return out


# ----------------------------------------------------------------------------------------------------------------
# model and calls

class LogisticModel:
    """L2-penalised logistic regression on z-scored features (z-scores from the data the model is applied to, so a
    model transfers between cell types and assays of different scale)."""

    def __init__(self, features, l2=1.0):
        self.features = list(features)
        self.l2 = l2
        self.coef = None
        self.intercept = None

    @staticmethod
    def standardise(x):
        x = np.asarray(x, dtype=float)
        mu = np.nanmean(x, axis=0)
        sd = np.nanstd(x, axis=0)
        sd[~(sd > 0)] = 1.0
        return np.nan_to_num((x - mu) / sd)

    def _matrix(self, table):
        return self.standardise(table[self.features].values)

    def fit(self, table, y):
        x = self._matrix(table)
        y = np.asarray(y, dtype=float)
        n, p = x.shape

        def loss(beta):
            z = beta[0] + x @ beta[1:]
            nll = np.sum(np.logaddexp(0, z) - y * z)
            pr = 1.0 / (1.0 + np.exp(-z))
            grad = np.concatenate([[np.sum(pr - y)], x.T @ (pr - y)])
            reg = 0.5 * self.l2 * np.sum(beta[1:] ** 2)
            grad[1:] += self.l2 * beta[1:]
            return nll + reg, grad

        p0 = np.clip(y.mean(), 1e-6, 1 - 1e-6)
        beta0 = np.zeros(p + 1)
        beta0[0] = np.log(p0 / (1 - p0))
        res = optimize.minimize(loss, beta0, jac=True, method='L-BFGS-B')
        self.intercept, self.coef = float(res.x[0]), res.x[1:]
        self.loglik = -float(np.sum(np.logaddexp(0, res.x[0] + x @ res.x[1:]) - y * (res.x[0] + x @ res.x[1:])))
        return self

    def decision(self, table):
        return self.intercept + self._matrix(table) @ self.coef

    def predict(self, table):
        return 1.0 / (1.0 + np.exp(-self.decision(table)))

    def to_dict(self):
        return {'features': self.features, 'l2': self.l2, 'intercept': self.intercept,
                'coef': [float(c) for c in self.coef]}

    @classmethod
    def from_dict(cls, d):
        m = cls(d['features'], d.get('l2', 1.0))
        m.intercept = d['intercept']
        m.coef = np.asarray(d['coef'], dtype=float)
        return m


def peak_calls(score, codes, threshold, radius=1):
    """Bins whose score is >= threshold and the maximum within +-radius bins of the same chromosome."""
    score = np.asarray(score, dtype=float)
    out = []
    for c in np.unique(codes):
        idx = np.flatnonzero(codes == c)
        v = np.nan_to_num(score[idx], nan=-np.inf)
        mx = ndimage.maximum_filter1d(v, 2 * radius + 1, mode='constant', cval=-np.inf)
        k = np.flatnonzero((v >= mx) & (v >= threshold))
        if len(k) > 1:
            k = k[np.concatenate([[True], np.diff(k) > radius])]
        out.append(idx[k])
    return np.sort(np.concatenate(out)) if out else np.array([], dtype=np.int64)


def match_boundaries(called, reference, codes, tolerance=1, mask=None):
    """Precision (called within `tolerance` bins of a reference boundary), recall (reference boundaries within
    tolerance of a call) and F1; with a bin mask both sets are restricted to masked bins."""
    called = np.asarray(called, dtype=np.int64)
    reference = np.asarray(reference, dtype=np.int64)
    if mask is not None:
        called = called[mask[called]]
        reference = reference[mask[reference]]
    # chromosome-aware distances: offset every chromosome far apart
    offs = np.asarray(codes, dtype=np.int64) * (10 ** 7)
    c = np.sort(called + offs[called])
    r = np.sort(reference + offs[reference])

    def hits(q, t):
        if len(t) == 0 or len(q) == 0:
            return np.zeros(len(q), dtype=bool)
        k = np.searchsorted(t, q)
        d1 = np.abs(q - t[np.clip(k, 0, len(t) - 1)])
        d0 = np.abs(q - t[np.clip(k - 1, 0, len(t) - 1)])
        return np.minimum(d0, d1) <= tolerance

    tp_c = hits(c, r).sum()
    tp_r = hits(r, c).sum()
    precision = tp_c / len(c) if len(c) else np.nan
    recall = tp_r / len(r) if len(r) else np.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
    return {'n_called': int(len(c)), 'n_reference': int(len(r)), 'precision': float(precision),
            'recall': float(recall), 'f1': float(f1)}


def best_threshold(score, reference, codes, tolerance=1, radius=1, n_grid=60):
    """Threshold of `peak_calls` that maximises F1 against the reference (training data)."""
    score = np.asarray(score, dtype=float)
    finite = score[np.isfinite(score)]
    grid = np.unique(np.quantile(finite, np.linspace(0.5, 0.999, n_grid)))
    best = (-1.0, grid[0])
    for t in grid:
        f1 = match_boundaries(peak_calls(score, codes, t, radius), reference, codes, tolerance)['f1']
        if f1 > best[0]:
            best = (f1, t)
    return best[1], best[0]


def average_precision(score, label):
    """Area under the precision-recall curve (step-wise average precision, ties broken conservatively)."""
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=bool)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n_pos = label.sum()
    if n_pos == 0:
        return np.nan
    order = np.argsort(-score, kind='mergesort')
    s, lab = score[order], label[order]
    tp = np.cumsum(lab)
    # evaluate at the last index of each tie group
    last = np.concatenate([np.flatnonzero(np.diff(s) != 0), [len(s) - 1]])
    tp_t = tp[last]
    prec = tp_t / (last + 1)
    rec = tp_t / n_pos
    rec_prev = np.concatenate([[0.0], rec[:-1]])
    return float(np.sum((rec - rec_prev) * prec))


def bin_labels(n_bins, boundary_bins, tolerance=0, codes=None):
    y = np.zeros(n_bins, dtype=bool)
    b = np.asarray(boundary_bins, dtype=np.int64)
    for d in range(-tolerance, tolerance + 1):
        k = b + d
        ok = (k >= 0) & (k < n_bins)
        if codes is not None:
            ok &= codes[np.clip(k, 0, n_bins - 1)] == codes[b]
        y[k[ok]] = True
    return y


def positions_to_bins(bins, positions):
    """Bin row index of positions: dict chrom -> array of bp (positions outside the bins are dropped)."""
    idx = bin_index(bins)
    out = []
    for chrom, pos in positions.items():
        if chrom not in idx:
            continue
        first, n, res = idx[chrom]
        k = (np.asarray(pos) // res).astype(np.int64)
        k = k[(k >= 0) & (k < n)]
        out.append(first + k)
    return np.unique(np.concatenate(out)) if out else np.array([], dtype=np.int64)


# ----------------------------------------------------------------------------------------------------------------
# domain agreement

def segment_labels(codes, boundary_bins):
    """Domain label per bin: segments between consecutive boundaries within each chromosome."""
    codes = np.asarray(codes)
    marks = np.zeros(len(codes), dtype=np.int64)
    marks[np.asarray(boundary_bins, dtype=np.int64)] = 1
    new_chrom = np.concatenate([[1], (np.diff(codes) != 0).astype(np.int64)])
    return np.cumsum(np.maximum(marks, new_chrom))


def variation_of_information(a, b, normalise=True):
    a = pd.factorize(np.asarray(a))[0]
    b = pd.factorize(np.asarray(b))[0]
    n = len(a)
    joint = pd.Series(1, index=pd.MultiIndex.from_arrays([a, b])).groupby(level=[0, 1]).sum().values / n
    pa = np.bincount(a) / n
    pb = np.bincount(b) / n
    h_ab = -np.sum(joint * np.log(joint))
    h_a = -np.sum(pa[pa > 0] * np.log(pa[pa > 0]))
    h_b = -np.sum(pb[pb > 0] * np.log(pb[pb > 0]))
    vi = 2 * h_ab - h_a - h_b
    return vi / np.log(n) if normalise else vi


def local_rand(a, b, codes, max_offset):
    """Fraction of bin pairs up to `max_offset` bins apart (same chromosome) on which two partitions agree about
    being in the same domain; also returns the fraction of same-domain pairs of each partition."""
    a, b, codes = np.asarray(a), np.asarray(b), np.asarray(codes)
    agree = total = same_a = same_b = 0
    for d in range(1, max_offset + 1):
        ok = codes[d:] == codes[:-d]
        sa = (a[d:] == a[:-d])[ok]
        sb = (b[d:] == b[:-d])[ok]
        agree += np.sum(sa == sb)
        total += ok.sum()
        same_a += sa.sum()
        same_b += sb.sum()
    return agree / total, same_a / total, same_b / total


def random_boundaries(codes, counts, rng):
    """Random boundary bins with the given number per chromosome code."""
    out = []
    for c, k in counts.items():
        idx = np.flatnonzero(codes == c)
        if k > 0 and len(idx):
            out.append(np.sort(rng.choice(idx, size=min(k, len(idx)), replace=False)))
    return np.concatenate(out) if out else np.array([], dtype=np.int64)
