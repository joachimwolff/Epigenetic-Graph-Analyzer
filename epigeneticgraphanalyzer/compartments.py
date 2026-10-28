"""A/B compartment estimates from chromatin tracks, the Hi-C compartment eigenvector and domain statistics of the
interaction graph.

Compartment estimates per fixed genomic bin (default 100 kb):

* chromatin, unsupervised: first principal component of the binned log enrichment of all tracks, sign oriented
  so that the tracks with an active role load positively (`chromatin_pc1`);
* chromatin, declared roles: mean z-score of the active tracks minus mean z-score of the repressive tracks
  (`role_score`);
* Hi-C: per chromosome observed / expected matrix, Pearson correlation of its rows, first eigenvector
  (`contact_eigenvector`; the reference eigenvector of a Hi-C or Micro-C matrix), bins without coverage NaN;
* profiles (Fortin and Hansen 2015 style): bin-bin correlation of binned signal profiles over several tracks or
  samples, first eigenvector (`profile_eigenvector`).

Every estimate is compared with a reference eigenvector by `compare_to_reference` (Pearson, Spearman, A/B sign
agreement and the AUC of the estimate for the reference A label, per chromosome and genome-wide).

Domain statistics: community boundaries along the chromosome, their enrichment at TAD boundaries against circular
shifts and against a node-matched null, the fraction of spans inside one domain and the fraction of edges crossing
a domain boundary against shifted copies.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)

ACTIVE_ROLES = ('active', 'enhancer', 'transcription')
REPRESSIVE_ROLES = ('repressive', 'heterochromatin')


# ----------------------------------------------------------------------------------------------------------------
# binning

def make_bins(chrom_sizes, bin_size, chroms=None):
    """DataFrame of fixed bins (chrom, start, end) in chromosome order."""
    rows = []
    for chrom, size in chrom_sizes.items():
        if chroms is not None and chrom not in chroms:
            continue
        starts = np.arange(0, size, bin_size, dtype=np.int64)
        rows.append(pd.DataFrame({'chrom': chrom, 'start': starts,
                                  'end': np.minimum(starts + bin_size, size)}))
    return pd.concat(rows, ignore_index=True)


def bin_bigwig(path, bins):
    """Mean signal of a bigWig in every bin (NaN where the file has no data or no such chromosome)."""
    import pyBigWig
    handle = pyBigWig.open(path)
    sizes = handle.chroms()
    out = np.full(len(bins), np.nan)
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        if chrom not in sizes:
            continue
        for i in idx:
            start, end = int(bins.start.iat[i]), min(int(bins.end.iat[i]), sizes[chrom])
            if end <= start:
                continue
            value = handle.stats(chrom, start, end, type='mean', exact=False)[0]
            out[i] = np.nan if value is None else value
    handle.close()
    return out


def interval_density(bins, chroms, starts, ends=None):
    """Number of intervals (by midpoint) per bin, e.g. gene density."""
    mids = np.asarray(starts, dtype=np.int64) if ends is None else \
        (np.asarray(starts, dtype=np.int64) + np.asarray(ends, dtype=np.int64)) // 2
    chroms = np.asarray(chroms)
    out = np.zeros(len(bins))
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        starts_c = bins.start.values[idx]
        pos = np.sort(mids[chroms == chrom])
        if len(pos) == 0:
            continue
        ends_c = bins.end.values[idx]
        out[idx] = np.searchsorted(pos, ends_c, side='left') - np.searchsorted(pos, starts_c, side='left')
    return out


def log_enrichment(values, pseudocount_fraction=0.1):
    """log2 of the signal over its genome-wide mean, with a pseudocount of a fraction of that mean."""
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values)
    if not np.isfinite(mean) or mean <= 0:
        return np.full_like(values, np.nan)
    pseudo = pseudocount_fraction * mean
    return np.log2((np.clip(values, 0, None) + pseudo) / (mean + pseudo))


def smooth_per_chrom(values, chroms, window):
    """Running mean over 2 * window + 1 bins within each chromosome (NaN ignored); window 0 returns a copy."""
    values = np.asarray(values, dtype=float)
    if window <= 0:
        return values.copy()
    out = np.full_like(values, np.nan)
    chroms = np.asarray(chroms)
    kernel = np.ones(2 * window + 1)
    for chrom in pd.unique(chroms):
        idx = np.flatnonzero(chroms == chrom)
        x = values[idx]
        finite = np.isfinite(x)
        num = np.convolve(np.where(finite, x, 0.0), kernel, mode='same')
        den = np.convolve(finite.astype(float), kernel, mode='same')
        with np.errstate(invalid='ignore', divide='ignore'):
            y = num / den
        y[~finite] = np.nan
        out[idx] = y
    return out


# ----------------------------------------------------------------------------------------------------------------
# chromatin scores

def _zscore_columns(features):
    x = np.asarray(features, dtype=float)
    mean = np.nanmean(x, axis=0)
    std = np.nanstd(x, axis=0)
    std[std == 0] = 1.0
    return (x - mean) / std


def chromatin_pc1(features, active=None, orient_by=None):
    """First principal component of z-scored features (bins x tracks).

    Rows with a missing value get NaN. Sign: the mean loading of the `active` columns (boolean mask or indices) is
    positive; without active columns the score correlates positively with `orient_by` (per-bin vector).
    Returns (score, loadings, explained variance fraction).
    """
    z = _zscore_columns(features)
    valid = np.all(np.isfinite(z), axis=1)
    score = np.full(z.shape[0], np.nan)
    if valid.sum() < 3:
        return score, np.full(z.shape[1], np.nan), np.nan
    zv = z[valid] - z[valid].mean(axis=0)
    _, singular, vt = np.linalg.svd(zv, full_matrices=False)
    loadings = vt[0]
    pc = zv @ loadings
    sign = 1.0
    if active is not None and np.any(np.asarray(active)):
        sign = np.sign(np.mean(loadings[np.asarray(active)])) or 1.0
    elif orient_by is not None:
        o = np.asarray(orient_by, dtype=float)[valid]
        ok = np.isfinite(o)
        if ok.sum() > 2:
            sign = np.sign(np.corrcoef(pc[ok], o[ok])[0, 1]) or 1.0
    score[valid] = sign * pc
    explained = singular[0] ** 2 / np.sum(singular ** 2)
    return score, sign * loadings, explained


def role_score(features, roles):
    """Mean z-score of the active tracks minus mean z-score of the repressive tracks.

    roles: one role string per column; 'active', 'enhancer', 'transcription' count as active, 'repressive' and
    'heterochromatin' as repressive, anything else is ignored. Without repressive tracks the score is the active mean.
    """
    z = _zscore_columns(features)
    roles = [str(r).lower() for r in roles]
    act = [i for i, r in enumerate(roles) if r in ACTIVE_ROLES]
    rep = [i for i, r in enumerate(roles) if r in REPRESSIVE_ROLES]
    if not act:
        raise ValueError('role_score needs at least one active track')
    score = np.mean(z[:, act], axis=1)
    if rep:
        score = score - np.mean(z[:, rep], axis=1)
    return score


def loco_regression(features, target, chroms, ridge=1.0):
    """Leave-one-chromosome-out ridge regression of the reference eigenvector on the features.

    A supervised upper bound for what the tracks can predict (it needs Hi-C for training); returns the
    out-of-chromosome predictions.
    """
    z = _zscore_columns(features)
    target = np.asarray(target, dtype=float)
    chroms = np.asarray(chroms)
    valid = np.all(np.isfinite(z), axis=1) & np.isfinite(target) & (target != 0)
    pred = np.full(len(target), np.nan)
    design = np.column_stack([np.ones(len(z)), np.nan_to_num(z)])
    penalty = ridge * np.eye(design.shape[1])
    penalty[0, 0] = 0
    for chrom in pd.unique(chroms):
        test = chroms == chrom
        train = valid & ~test
        if train.sum() < design.shape[1] + 2:
            continue
        a = design[train]
        coef = np.linalg.solve(a.T @ a + penalty, a.T @ target[train])
        rows = test & np.all(np.isfinite(z), axis=1)
        pred[rows] = design[rows] @ coef
    return pred


# ----------------------------------------------------------------------------------------------------------------
# contact matrices

def observed_expected(matrix, valid=None):
    """Distance-normalised matrix: each diagonal divided by its mean over pairs of valid bins.

    matrix: symmetric dense array of one chromosome; valid: boolean mask of bins to use (default: row sum > 0).
    Entries of invalid bins are NaN.
    """
    m = np.asarray(matrix, dtype=float)
    n = m.shape[0]
    if valid is None:
        valid = np.nansum(m, axis=1) > 0
    valid = np.asarray(valid, dtype=bool)
    oe = np.full_like(m, np.nan)
    for k in range(n):
        i = np.arange(n - k)
        j = i + k
        ok = valid[i] & valid[j]
        if not ok.any():
            continue
        diag = m[i[ok], j[ok]]
        mean = np.nanmean(diag)
        if not np.isfinite(mean) or mean <= 0:
            continue
        oe[i[ok], j[ok]] = diag / mean
        oe[j[ok], i[ok]] = diag / mean
    return oe


def _leading_eigenvectors_of_rows(x, n=1):
    """Leading n eigenvectors of the row-row Pearson correlation of x (rows = bins); rows without variance NaN.

    Returns (bins x n array, eigenvalue fractions)."""
    x = np.where(np.isfinite(x), x, 0.0)
    keep = np.nanstd(x, axis=1) > 0
    vecs = np.full((x.shape[0], n), np.nan)
    if keep.sum() < max(3, n + 1):
        return vecs, np.full(n, np.nan)
    corr = np.nan_to_num(np.corrcoef(x[keep]))
    values, vectors = np.linalg.eigh(corr)
    order = np.argsort(values)[::-1][:n]
    vecs[keep] = vectors[:, order]
    return vecs, values[order] / np.sum(np.abs(values))


def _first_eigenvector_of_rows(x):
    vecs, frac = _leading_eigenvectors_of_rows(x, 1)
    return vecs[:, 0], frac[0]


def _orient(vec, orient_by):
    if orient_by is None:
        return vec
    o = np.asarray(orient_by, dtype=float)
    ok = np.isfinite(vec) & np.isfinite(o)
    if ok.sum() > 2 and np.std(vec[ok]) > 0 and np.std(o[ok]) > 0:
        if np.corrcoef(vec[ok], o[ok])[0, 1] < 0:
            return -vec
    return vec


def _select(vecs, phasing, margin=1.5):
    """Index of the compartment eigenvector: the first one, unless another of the leading ones correlates at least
    `margin` times more strongly (absolute Pearson r) with the phasing track. A plain maximum is too greedy:
    on shallow matrices a second eigenvector often correlates slightly better with GC than the first one although
    the first one is the compartment pattern."""
    if phasing is None or vecs.shape[1] == 1:
        return vecs[:, 0], 0
    o = np.asarray(phasing, dtype=float)
    r = np.zeros(vecs.shape[1])
    for k in range(vecs.shape[1]):
        v = vecs[:, k]
        ok = np.isfinite(v) & np.isfinite(o)
        if ok.sum() > 2 and np.std(v[ok]) > 0 and np.std(o[ok]) > 0:
            r[k] = abs(np.corrcoef(v[ok], o[ok])[0, 1])
    best = int(np.argmax(r))
    if best != 0 and r[best] >= margin * r[0]:
        return vecs[:, best], best
    return vecs[:, 0], 0


def contact_eigenvector(matrix, valid=None, orient_by=None, min_contacts=1, n_candidates=1):
    """Compartment eigenvector of the correlation of the observed / expected matrix of one chromosome.

    Bins with fewer than `min_contacts` non-zero entries (or outside `valid`) are NaN. With n_candidates > 1 another
    of the leading eigenvectors replaces the first one when it correlates at least 1.5 times more strongly with
    `orient_by` (_select), so that a chromosome-arm or other leading pattern does not replace the compartments. The sign
    is chosen so that the vector correlates positively with `orient_by` (gene density, GC, chromatin score).
    Returns (vector, index of the chosen eigenvector).
    """
    m = np.asarray(matrix, dtype=float)
    nonzero = np.sum(np.nan_to_num(m) > 0, axis=1)
    covered = nonzero >= min_contacts
    if valid is not None:
        covered &= np.asarray(valid, dtype=bool)
    vec = np.full(m.shape[0], np.nan)
    if covered.sum() < max(3, n_candidates + 1):
        return vec, -1
    oe = observed_expected(m, covered)
    sub = oe[np.ix_(covered, covered)]
    vecs, _ = _leading_eigenvectors_of_rows(sub, n_candidates)
    phasing = None if orient_by is None else np.asarray(orient_by, dtype=float)[covered]
    chosen, index = _select(vecs, phasing)
    vec[covered] = chosen
    return _orient(vec, orient_by), index


def profile_eigenvector(profiles, orient_by=None):
    """Fortin and Hansen style estimate: bins x samples (or tracks) of one chromosome; first eigenvector of the
    bin-bin correlation of the column-centred profiles. Needs at least three columns."""
    x = np.asarray(profiles, dtype=float)
    if x.ndim != 2 or x.shape[1] < 3:
        raise ValueError('profile_eigenvector needs at least three samples or tracks')
    valid = np.all(np.isfinite(x), axis=1)
    vec = np.full(x.shape[0], np.nan)
    xv = x[valid] - x[valid].mean(axis=0)
    v, _ = _first_eigenvector_of_rows(xv)
    vec[valid] = v
    return _orient(vec, orient_by)


def per_chrom(function, values_by_bin, bins, **kwargs):
    """Apply a per-chromosome eigenvector function to a bins x k array; returns the concatenated vector."""
    out = np.full(len(bins), np.nan)
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        sub = {k: (None if v is None else np.asarray(v)[idx]) for k, v in kwargs.items()}
        result = function(np.asarray(values_by_bin)[idx], **sub)
        out[idx] = result[0] if isinstance(result, tuple) else result
    return out


# ----------------------------------------------------------------------------------------------------------------
# evaluation against a reference eigenvector

def _auc(score, label):
    label = np.asarray(label, dtype=bool)
    n_pos, n_neg = label.sum(), (~label).sum()
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = stats.rankdata(score)
    return (ranks[label].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def _metrics(pred, ref, threshold):
    ok = np.isfinite(pred) & np.isfinite(ref) & (ref != 0)
    p, r = pred[ok], ref[ok]
    row = {'n_bins': int(ok.sum())}
    if ok.sum() < 5 or np.std(p) == 0:
        row.update(pearson=np.nan, spearman=np.nan, sign_agreement=np.nan, auc=np.nan, frac_A_reference=np.nan)
        return row
    cut = np.median(p) if threshold == 'median' else 0.0
    row.update(pearson=stats.pearsonr(p, r)[0], spearman=stats.spearmanr(p, r)[0],
               sign_agreement=np.mean((p > cut) == (r > 0)), auc=_auc(p, r > 0), frac_A_reference=np.mean(r > 0))
    return row


def compare_to_reference(pred, ref, chroms, threshold='zero', min_bins=20):
    """Per chromosome and genome-wide agreement of an estimate with a reference eigenvector.

    Bins where either value is missing or the reference is exactly 0 (masked) are skipped. threshold: 'zero' calls
    A where the estimate is > 0, 'median' where it exceeds the per-chromosome median of the estimate. The genome-wide
    row pools all bins (estimates centred per chromosome when threshold is 'median'); 'median_of_chroms' is the median
    of the per-chromosome values.
    """
    pred = np.asarray(pred, dtype=float)
    ref = np.asarray(ref, dtype=float)
    chroms = np.asarray(chroms)
    rows = []
    pooled = pred.copy()
    for chrom in pd.unique(chroms):
        idx = chroms == chrom
        row = _metrics(pred[idx], ref[idx], threshold)
        if threshold == 'median':
            finite = idx & np.isfinite(pred)
            if finite.any():
                pooled[finite] = pred[finite] - np.median(pred[finite])
        if row['n_bins'] >= min_bins:
            row['chrom'] = chrom
            rows.append(row)
    table = pd.DataFrame(rows)
    genome = _metrics(pooled, ref, 'zero')
    genome['chrom'] = 'genome'
    summary = {'chrom': 'median_of_chroms', 'n_bins': int(table.n_bins.sum()) if len(table) else 0}
    for key in ('pearson', 'spearman', 'sign_agreement', 'auc', 'frac_A_reference'):
        summary[key] = float(np.nanmedian(table[key])) if len(table) else np.nan
    return pd.concat([table, pd.DataFrame([genome, summary])], ignore_index=True)


# ----------------------------------------------------------------------------------------------------------------
# domains

def read_domains(path):
    """Domains as a DataFrame (chrom, start, end), sorted; accepts BED and the Rao et al. domain list."""
    rows = []
    with open(path) if not str(path).endswith('.gz') else __import__('gzip').open(path, 'rt') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3 or line.startswith(('#', 'track', 'browser')):
                continue
            try:
                start, end = int(float(fields[1])), int(float(fields[2]))
            except ValueError:
                continue
            chrom = fields[0] if fields[0].startswith('chr') else 'chr' + fields[0]
            rows.append((chrom, start, end))
    df = pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
    return df.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)


def flatten_domains(domains):
    """Nested or overlapping domains as non-overlapping segments between consecutive boundaries (kept where at
    least one domain covers the segment). Two positions share a segment exactly when no domain boundary lies
    between them and both are inside some domain."""
    rows = []
    for chrom, sub in domains.groupby('chrom', sort=False):
        cuts = np.unique(np.concatenate([sub.start.values, sub.end.values]))
        mids = (cuts[:-1] + cuts[1:]) // 2
        starts, ends = np.sort(sub.start.values), np.sort(sub.end.values)
        covered = (np.searchsorted(starts, mids, side='right') - np.searchsorted(ends, mids, side='right')) > 0
        rows.append(pd.DataFrame({'chrom': chrom, 'start': cuts[:-1][covered], 'end': cuts[1:][covered]}))
    if not rows:
        return domains.iloc[:0][['chrom', 'start', 'end']]
    return pd.concat(rows, ignore_index=True)


def domain_boundaries(domains):
    """Unique boundary positions (domain starts and ends) per chromosome: dict chrom -> sorted array."""
    out = {}
    for chrom, sub in domains.groupby('chrom'):
        out[chrom] = np.unique(np.concatenate([sub.start.values, sub.end.values]))
    return out


def domain_id(domains, chrom, positions):
    """Index of the domain containing each position (the innermost if nested, -1 outside every domain).

    For non-nested domains (hicFindTADs, insulation) the domain is unique.
    """
    positions = np.asarray(positions, dtype=np.int64)
    out = np.full(len(positions), -1, dtype=np.int64)
    sub = domains[domains.chrom == chrom]
    if sub.empty:
        return out
    starts, ends, ids = sub.start.values, sub.end.values, sub.index.values
    order = np.argsort(starts, kind='stable')
    starts, ends, ids = starts[order], ends[order], ids[order]
    if np.all(ends[:-1] <= starts[1:]):    # non-overlapping: binary search
        k = np.searchsorted(starts, positions, side='right') - 1
        hit = (k >= 0) & (positions < ends[np.clip(k, 0, None)])
        out[hit] = ids[k[hit]]
        return out
    lengths = ends - starts
    order = np.argsort(-lengths)          # larger first, so smaller (inner) domains overwrite
    for k in order:
        inside = (positions >= starts[k]) & (positions < ends[k])
        out[inside] = ids[k]
    return out


def community_boundaries(nodes):
    """Positions where consecutive nodes along a chromosome switch community.

    nodes: DataFrame with chrom, start, end, community. Nodes without a community are skipped. Returns a DataFrame
    (chrom, position, left, right) with the midpoint between the two nodes as position.
    """
    rows = []
    df = nodes.dropna(subset=['community']).copy()
    df['mid'] = (df.start.values + df.end.values) // 2
    for chrom, sub in df.sort_values(['chrom', 'mid']).groupby('chrom', sort=False):
        com = sub.community.values
        mid = sub.mid.values
        switch = np.flatnonzero(com[1:] != com[:-1])
        for i in switch:
            rows.append((chrom, (sub.end.values[i] + sub.start.values[i + 1]) // 2, mid[i], mid[i + 1]))
    return pd.DataFrame(rows, columns=['chrom', 'position', 'left', 'right'])


def _near(positions, boundaries, window):
    """Boolean: position within `window` bp of any boundary (boundaries sorted)."""
    if len(boundaries) == 0:
        return np.zeros(len(positions), dtype=bool)
    idx = np.searchsorted(boundaries, positions)
    left = np.abs(positions - boundaries[np.clip(idx - 1, 0, len(boundaries) - 1)])
    right = np.abs(boundaries[np.clip(idx, 0, len(boundaries) - 1)] - positions)
    return np.minimum(left, right) <= window


def boundary_enrichment(query, reference, spans, window=50000, n_shifts=1000, seed=0, null_positions=None):
    """Fraction of query positions within `window` of a reference boundary, against two nulls.

    query, reference: dict chrom -> positions; spans: dict chrom -> (lo, hi), the covered range in which positions
    are circularly shifted (the same random offset for all query positions of a chromosome per shift).
    null_positions: optional dict chrom -> candidate positions (e.g. all gaps between consecutive nodes); then a
    second null draws the same number of query positions from the candidates per chromosome.
    Returns a dict with observed, shift mean, enrichment, empirical P (one-sided, +1 corrected) and the same for the
    matched null.
    """
    rng = np.random.default_rng(seed)
    chroms = [c for c in query if c in reference and c in spans and len(query[c])]
    n = sum(len(query[c]) for c in chroms)
    observed = sum(_near(np.asarray(query[c]), reference[c], window).sum() for c in chroms)
    shifted = np.zeros(n_shifts)
    matched = np.zeros(n_shifts) if null_positions is not None else None
    for s in range(n_shifts):
        hits = 0
        hits_m = 0
        for c in chroms:
            lo, hi = spans[c]
            length = max(hi - lo, 1)
            q = (np.asarray(query[c]) - lo + rng.integers(0, length)) % length + lo
            hits += _near(q, reference[c], window).sum()
            if matched is not None:
                cand = np.asarray(null_positions[c])
                pick = rng.choice(cand, size=min(len(query[c]), len(cand)), replace=False)
                hits_m += _near(pick, reference[c], window).sum()
        shifted[s] = hits
        if matched is not None:
            matched[s] = hits_m
    out = {'n_query': int(n), 'observed_fraction': observed / n if n else np.nan,
           'shift_fraction': shifted.mean() / n if n else np.nan,
           'shift_enrichment': observed / shifted.mean() if shifted.mean() > 0 else np.nan,
           'shift_p': (1 + np.sum(shifted >= observed)) / (1 + n_shifts)}
    if matched is not None:
        out.update(matched_fraction=matched.mean() / n if n else np.nan,
                   matched_enrichment=observed / matched.mean() if matched.mean() > 0 else np.nan,
                   matched_p=(1 + np.sum(matched >= observed)) / (1 + n_shifts))
    return out


def _pair_domains(domains, chrom, a, b):
    """(same domain, both ends inside a domain) for position pairs."""
    ia, ib = domain_id(domains, chrom, a), domain_id(domains, chrom, b)
    both = (ia >= 0) & (ib >= 0)
    return (ia == ib) & both, both


def _shift_test(groups, domains, chrom_spans, n_shifts, seed, first, second):
    """Observed and circularly shifted fractions of pairs in one domain, among pairs with both ends in a domain.

    first / second: functions of a group returning the two positions of every pair. Each shift moves all pairs of a
    chromosome by one random offset within chrom_spans (distances kept); pairs leaving the span are dropped.
    """
    rng = np.random.default_rng(seed)
    same = inside = 0
    for c, g in groups:
        s, b = _pair_domains(domains, c, first(g), second(g))
        same += s.sum()
        inside += b.sum()
    shifted = np.zeros(n_shifts)
    for k in range(n_shifts):
        s_tot = b_tot = 0
        for c, g in groups:
            lo, hi = chrom_spans[c]
            length = max(hi - lo, 1)
            offset = rng.integers(0, length)
            x = (first(g) - lo + offset) % length + lo
            y = x + (second(g) - first(g))
            ok = (y >= lo) & (y < hi)
            s, b = _pair_domains(domains, c, x[ok], y[ok])
            s_tot += s.sum()
            b_tot += b.sum()
        shifted[k] = s_tot / max(b_tot, 1)
    return same / max(inside, 1), int(inside), shifted


def spans_within_domain(spans, domains, chrom_spans, n_shifts=1000, seed=0):
    """Fraction of spans (chrom, start, end) without a domain boundary inside, observed and under circular shifts.

    A span is inside one domain when its start and end - 1 fall in the same domain; only spans whose two ends lie
    in some domain count (for observed and shifted spans alike), so gaps between domains (centromeres, unmappable
    regions) do not enter. Shifts move every span of a chromosome by one random offset within chrom_spans[chrom].
    """
    spans = spans[spans.chrom.isin(chrom_spans.keys())]
    groups = list(spans.groupby('chrom'))
    obs, n_inside, shifted = _shift_test(groups, domains, chrom_spans, n_shifts, seed,
                                         lambda g: g.start.values, lambda g: g.end.values - 1)
    return {'n_spans': int(len(spans)), 'n_spans_in_domains': n_inside, 'observed_fraction': obs,
            'shift_fraction': shifted.mean(), 'ratio': obs / shifted.mean() if shifted.mean() > 0 else np.nan,
            'shift_p': (1 + np.sum(shifted >= obs)) / (1 + n_shifts)}


def edges_crossing(edges, domains, chrom_spans, n_shifts=1000, seed=0):
    """Fraction of edges whose two ends (midpoints mid_a, mid_b) lie in different domains, among edges with both
    ends in a domain, observed and for the same edges shifted circularly along the chromosome (distances kept, so the
    expectation is distance matched). Returns the observed fraction, the shift mean, ratio and the empirical P of a
    deficit of crossings (one-sided)."""
    edges = edges[edges.chrom.isin(chrom_spans.keys())]
    groups = list(edges.groupby('chrom'))
    same, n_inside, shifted = _shift_test(groups, domains, chrom_spans, n_shifts, seed,
                                          lambda g: g.mid_a.values, lambda g: g.mid_b.values)
    obs = 1 - same
    crossing = 1 - shifted
    return {'n_edges': int(len(edges)), 'n_edges_in_domains': n_inside, 'observed_fraction': obs,
            'shift_fraction': crossing.mean(), 'ratio': obs / crossing.mean() if crossing.mean() > 0 else np.nan,
            'p_deficit': (1 + np.sum(crossing <= obs)) / (1 + n_shifts)}


# ----------------------------------------------------------------------------------------------------------------
# condition-specific compartment switches

def compartment_switches(score_a, score_b, threshold=0.5):
    """Switch call per bin from two compartment scores on comparable scales (e.g. per-chromosome standardised).

    Returns an array of 'B->A', 'A->B', 'A', 'B' or 'NA': a switch needs a sign change and |b - a| >= threshold.
    """
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    out = np.full(len(a), 'NA', dtype=object)
    ok = np.isfinite(a) & np.isfinite(b)
    delta = b - a
    big = ok & (np.abs(delta) >= threshold)
    out[ok & (a > 0)] = 'A'
    out[ok & (a <= 0)] = 'B'
    out[big & (a <= 0) & (b > 0)] = 'B->A'
    out[big & (a > 0) & (b <= 0)] = 'A->B'
    return out


def quantile_normalise_pair(a, b, chroms):
    """Per chromosome, map both scores onto their mean quantile distribution (symmetric quantile normalisation).

    Removes differences in location and spread between the two conditions (e.g. ChIP efficiency), so that a switch
    reflects a change of rank along the chromosome. Bins missing in either condition stay NaN in both.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out_a = np.full_like(a, np.nan)
    out_b = np.full_like(b, np.nan)
    chroms = np.asarray(chroms)
    for chrom in pd.unique(chroms):
        idx = np.flatnonzero((chroms == chrom) & np.isfinite(a) & np.isfinite(b))
        if len(idx) < 3:
            continue
        target = (np.sort(a[idx]) + np.sort(b[idx])) / 2
        out_a[idx] = target[stats.rankdata(a[idx], method='ordinal').astype(int) - 1]
        out_b[idx] = target[stats.rankdata(b[idx], method='ordinal').astype(int) - 1]
    return out_a, out_b


def standardise_per_chrom(values, chroms):
    """Centre on the per-chromosome median and divide by the per-chromosome standard deviation."""
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan)
    chroms = np.asarray(chroms)
    for chrom in pd.unique(chroms):
        idx = chroms == chrom
        x = values[idx]
        ok = np.isfinite(x)
        if ok.sum() < 3:
            continue
        sd = np.std(x[ok])
        out[idx] = (x - np.median(x[ok])) / (sd if sd > 0 else 1.0)
    return out
