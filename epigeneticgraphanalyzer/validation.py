"""Validation of graph-derived contacts with ligation-free multi-way data.

Three levels are tested, each against controls that share the genomic
distance of the tested loci:

1. Pairs (direct edges, indirect pairs): co-occurrence of the two bins is
   compared with control pairs at the identical bin offset whose ends are drawn
   from the same pools of graph loci (bait or other-end) but that are not
   connected within two hops. This controls for distance decay and for the
   general enrichment of captured loci in active, co-localising chromatin.
2. Triplets u - v - w: the number of groups containing all three bins is
   compared with control triplets that keep the mediator v and replace u and w
   by graph loci at nearly the same offsets that are not adjacent to v. The
   per-triplet table also lists the ratio C(u,v,w) N(v) / (C(u,v) C(v,w)) as a
   descriptive quantity.
3. Communities: the number of groups containing at least k distinct bins of a
   community is compared with the same bin set circularly shifted along the
   chromosome, which preserves all pairwise distances.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from .genome import parse_node_id

log = logging.getLogger(__name__)

DISTANCE_STRATA = (0, 100000, 250000, 500000, 1000000, 2000001)


# --------------------------------------------------------------------- mapping
def map_nodes_to_bins(graph, data, liftover=None):
    """Bin id of every graph node (midpoint, optionally lifted to the assembly of `data`)."""
    mapping = {}
    for node, attrs in graph.nodes(data=True):
        chrom = attrs.get('chrom')
        if chrom is None:
            chrom, start, end = parse_node_id(node)
        else:
            start, end = int(attrs['start']), int(attrs['end'])
        position = (start + end) // 2
        if liftover is not None:
            lifted = liftover.position(chrom, position)
            if lifted is None:
                continue
            chrom, position = lifted
        bin_id = data.bins.bin_of(chrom, position)
        if bin_id >= 0 and data.coverage_raw[bin_id] > 0:
            mapping[node] = bin_id
    return mapping


def pairs_to_bins(pairs, node_bins):
    """Add bin_a/bin_b to a pair table; drops unmapped pairs and pairs inside one bin."""
    a = pairs['node_a'].map(node_bins)
    b = pairs['node_b'].map(node_bins)
    out = pairs.assign(bin_a=a, bin_b=b).dropna(subset=['bin_a', 'bin_b'])
    out['bin_a'] = out['bin_a'].astype(np.int64)
    out['bin_b'] = out['bin_b'].astype(np.int64)
    out = out[out['bin_a'] != out['bin_b']].copy()
    swap = out['bin_a'] > out['bin_b']
    out.loc[swap, ['bin_a', 'bin_b']] = out.loc[swap, ['bin_b', 'bin_a']].to_numpy()
    out['offset'] = out['bin_b'] - out['bin_a']
    return out.reset_index(drop=True)


# --------------------------------------------------------------------- metrics
def pair_metrics(data, bin_a, bin_b, max_offset):
    """Co-occurrence, distance-normalised co-occurrence and (GAM) NPMI of bin pairs."""
    bin_a = np.asarray(bin_a, dtype=np.int64)
    bin_b = np.asarray(bin_b, dtype=np.int64)
    offset = np.abs(bin_b - bin_a)
    data.build_band(max_offset)
    cooc = data.cooccurrence(bin_a, bin_b)
    expected = data.expected_by_offset(max_offset, 'cooc')
    exp_values = expected[np.clip(offset, 0, max_offset)]
    pseudo = np.nanmin(expected[1:][expected[1:] > 0]) * 0.1 if np.any(expected[1:] > 0) else 1e-9
    metrics = {
        'cooc': cooc,
        'cooc_oe': np.log2((cooc + pseudo) / (exp_values + pseudo)),
    }
    # NPMI normalises co-occurrence by how often each locus is detected at all, which removes the
    # dependence on locus length and mappability (computed for SPRITE clusters and GAM tubes alike)
    npmi = data.npmi(bin_a, bin_b)
    npmi_expected = data.expected_by_offset(max_offset, 'npmi')
    metrics['npmi'] = npmi
    metrics['npmi_delta'] = npmi - npmi_expected[np.clip(offset, 0, max_offset)]
    return metrics


# ------------------------------------------------------------- control pools
class LocusPools:
    """Bins hosting bait and other-end nodes, per chromosome, for control sampling."""

    def __init__(self, graph, node_bins, data):
        self.data = data
        self.bait = set()
        self.other = set()
        for node, bin_id in node_bins.items():
            if graph.nodes[node].get('bait'):
                self.bait.add(bin_id)
            else:
                self.other.add(bin_id)
        self.bait_array = np.array(sorted(self.bait), dtype=np.int64)
        self.other_array = np.array(sorted(self.other), dtype=np.int64)
        self.chrom_index = data.bins._chrom_of_bin

    def pool(self, kind):
        return self.bait_array if kind == 'bait' else self.other_array


def _class_pools(pair_class):
    if pair_class == 'bait-bait':
        return 'bait', 'bait'
    if pair_class == 'oe-oe':
        return 'other', 'other'
    return 'bait', 'other'


def coverage_classes(pools, n_classes=5):
    """Coverage class (quantile of the multi-way coverage over all graph loci) of every bin."""
    coverage = pools.data.coverage_raw
    loci = np.union1d(pools.bait_array, pools.other_array)
    edges = np.quantile(coverage[loci], np.linspace(0, 1, n_classes + 1)[1:-1])
    return np.searchsorted(edges, coverage, side='right')


def sample_matched_controls(pairs, pools, excluded, n_controls=10, seed=0, coverage_class=None, tad_id=None):
    """Control pairs at the same bin offset with ends from the same locus pools.

    `excluded` is a set of (bin_a, bin_b) tuples (bin_a < bin_b) that must not
    be drawn (direct and indirect contacts). With `coverage_class` (per bin, see
    coverage_classes) a control must also match the unordered pair of coverage classes
    of the tested pair; pairs without such a candidate fall back to offset matching. `coverage_class` can
    combine several per-bin classes (see domain_classes). With `tad_id` (per bin, -1 outside domains) a control
    must also match whether the two ends lie in the same topological domain.
    Returns a DataFrame with test_index, bin_a, bin_b.
    """
    rng = np.random.default_rng(seed)
    chrom_of = pools.chrom_index
    n_fallback = 0
    rows_index, rows_a, rows_b = [], [], []
    for (pair_class, offset), group in pairs.groupby(['pair_class', 'offset']):
        left_kind, right_kind = _class_pools(pair_class)
        left = pools.pool(left_kind)
        right = pools.pool(right_kind)
        # orientation 1: left pool upstream, right pool downstream
        starts = left[np.isin(left + offset, right)]
        starts = starts[chrom_of[starts] == chrom_of[np.minimum(starts + offset, len(chrom_of) - 1)]]
        candidates = [(s, s + offset) for s in starts]
        if left_kind != right_kind:
            starts2 = right[np.isin(right + offset, left)]
            starts2 = starts2[chrom_of[starts2] == chrom_of[np.minimum(starts2 + offset, len(chrom_of) - 1)]]
            candidates.extend((s, s + offset) for s in starts2)
        candidates = [c for c in candidates if c not in excluded]
        if not candidates:
            continue
        candidates = np.array(candidates, dtype=np.int64)
        by_class = None
        if coverage_class is not None or tad_id is not None:
            codes = _pair_code(coverage_class, tad_id, candidates[:, 0], candidates[:, 1])
            by_class = {code: candidates[codes == code] for code in np.unique(codes)}
            test_codes = _pair_code(coverage_class, tad_id, group['bin_a'].to_numpy(), group['bin_b'].to_numpy())
        for i, test_index in enumerate(group.index):
            pool = candidates
            if by_class is not None:
                pool = by_class.get(test_codes[i])
                if pool is None:
                    pool = candidates
                    n_fallback += 1
            choice = rng.integers(0, len(pool), size=n_controls)
            rows_index.extend([test_index] * n_controls)
            rows_a.extend(pool[choice, 0])
            rows_b.extend(pool[choice, 1])
    if coverage_class is not None or tad_id is not None:
        log.info('class-matched controls: %d of %d pairs fell back to offset matching', n_fallback, len(pairs))
    return pd.DataFrame({'test_index': np.array(rows_index, dtype=np.int64),
                         'bin_a': np.array(rows_a, dtype=np.int64),
                         'bin_b': np.array(rows_b, dtype=np.int64)})


def _pair_code(bin_class, tad_id, bin_a, bin_b):
    """Code of a pair: unordered pair of bin classes and, with tad_id, whether both ends share a domain."""
    code = _class_code(bin_class, bin_a, bin_b) if bin_class is not None else np.zeros(len(bin_a), dtype=np.int64)
    if tad_id is not None:
        same = (tad_id[bin_a] == tad_id[bin_b]) & (tad_id[bin_a] >= 0)
        code = code * 2 + same.astype(np.int64)
    return code


def _class_code(coverage_class, bin_a, bin_b):
    a, b = coverage_class[bin_a], coverage_class[bin_b]
    return np.minimum(a, b) * 1000 + np.maximum(a, b)


def _read_intervals(path):
    """Intervals of a BED-like file (first three columns; chromosome names with or without 'chr') with the name
    column if present. Rao et al. Arrowhead domain lists (chr1 x1 x2 chr2 y1 y2 ...) are read the same way."""
    rows = []
    opener = __import__('gzip').open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser', 'chr1\tx1')):
                continue
            fields = line.rstrip('\n').split('\t')
            try:
                start, end = int(float(fields[1])), int(float(fields[2]))
            except (ValueError, IndexError):
                continue
            chrom = fields[0] if fields[0].startswith('chr') else 'chr' + fields[0]
            rows.append((chrom, start, end, fields[3] if len(fields) > 3 else ''))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'name'])


def domain_annotation(bins, compartments=None, domains=None):
    """Per bin of a multi-way data set: A/B compartment (1 A, 0 B, -1 unknown) and topological domain id (-1 outside
    any domain), from the midpoint of the bin.

    `compartments`: BED whose name column starts with A or B (e.g. subcompartments A1, B2) or a bedGraph whose fourth
    column is the eigenvector (positive A). `domains`: BED of domains; for nested domains the smallest domain that
    contains the midpoint is used, so that two loci share a domain if they lie in the same innermost domain."""
    n = bins.n_bins
    chrom_names = np.array(bins.chrom_names)[bins.chrom_index_of_bins(np.arange(n))]
    middle = (np.asarray(bins.bin_start(np.arange(n))) + np.asarray(bins.bin_end(np.arange(n)))) // 2
    compartment = np.full(n, -1, dtype=np.int64)
    tad_id = np.full(n, -1, dtype=np.int64)
    if compartments is not None:
        table = _read_intervals(compartments)
        for chrom, sub in table.groupby('chrom'):
            on = np.nonzero(chrom_names == chrom)[0]
            if len(on) == 0:
                continue
            sub = sub.sort_values('start')
            index = np.searchsorted(sub['start'].to_numpy(), middle[on], side='right') - 1
            inside = (index >= 0) & (middle[on] < sub['end'].to_numpy()[np.clip(index, 0, None)])
            names = sub['name'].astype(str).to_numpy()
            label = np.full(len(on), -1, dtype=np.int64)
            for k in np.nonzero(inside)[0]:
                name = names[index[k]]
                if name[:1] in ('A', 'B'):
                    label[k] = 1 if name[0] == 'A' else 0
                else:
                    try:
                        value = float(name)
                    except ValueError:
                        continue
                    if np.isfinite(value) and value != 0:
                        label[k] = 1 if value > 0 else 0
            compartment[on] = label
    if domains is not None:
        table = _read_intervals(domains)
        table = table.assign(size=table['end'] - table['start']).sort_values('size', ascending=False)
        for k, (chrom, start, end) in enumerate(zip(table['chrom'], table['start'], table['end'])):
            on = np.nonzero((chrom_names == chrom) & (middle >= start) & (middle < end))[0]
            tad_id[on] = k        # smaller domains, written later, overwrite larger ones
    return compartment, tad_id


def domain_classes(coverage_class, compartment):
    """Per-bin class combining the coverage class and the compartment (A, B or unknown)."""
    base = coverage_class if coverage_class is not None else np.zeros(len(compartment), dtype=np.int64)
    return base * 3 + (compartment + 1)


def pair_concordance(test_values, control_values, control_index, n_tests):
    """Per tested item, the fraction of its matched controls with a lower value (ties one half)."""
    test_values = np.asarray(test_values, dtype=float)
    control_values = np.asarray(control_values, dtype=float)
    control_index = np.asarray(control_index, dtype=np.int64)
    keep = ~np.isnan(control_values) & ~np.isnan(test_values[control_index])
    control_values, control_index = control_values[keep], control_index[keep]
    reference = test_values[control_index]
    wins = (reference > control_values) + 0.5 * (reference == control_values)
    counts = np.bincount(control_index, minlength=n_tests)
    score = np.full(n_tests, np.nan)
    has = counts > 0
    score[has] = np.bincount(control_index, weights=wins, minlength=n_tests)[has] / counts[has]
    return score


def block_bootstrap(values, blocks, n_boot=1000, seed=0, statistic=np.nanmean):
    """Percentile 95% interval of a statistic over items, resampling whole blocks (e.g. chromosomes) with replacement,
    so that items sharing a block, and their dependence, stay together."""
    values = np.asarray(values, dtype=float)
    blocks = np.asarray(blocks)
    labels, index = np.unique(blocks, return_inverse=True)
    members = [np.nonzero(index == k)[0] for k in range(len(labels))]
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_boot)
    for b in range(n_boot):
        chosen = rng.integers(0, len(labels), len(labels))
        estimates[b] = statistic(values[np.concatenate([members[k] for k in chosen])])
    return float(np.nanpercentile(estimates, 2.5)), float(np.nanpercentile(estimates, 97.5))


def distance_stratum(distance):
    labels = ['{}-{}kb'.format(DISTANCE_STRATA[i] // 1000, DISTANCE_STRATA[i + 1] // 1000)
              for i in range(len(DISTANCE_STRATA) - 1)]
    index = np.clip(np.searchsorted(DISTANCE_STRATA, distance, side='right') - 1, 0, len(labels) - 1)
    return np.array(labels)[index]


def summarize_against_controls(test_values, control_values, control_index, n_tests):
    """Compare every tested item with its own matched controls.

    The concordance of an item is the fraction of its controls with a lower value, counting ties as one half;
    its mean over items, `fraction_above_control`, is the probability that a tested item exceeds a matched
    control and is 0.5 under the null hypothesis also for sparse, discrete values such as co-occurrence
    counts. The Wilcoxon signed-rank test compares the item concordances with 0.5. Mean values of items and
    controls are reported alongside.
    """
    test_values = np.asarray(test_values, dtype=float)
    control_values = np.asarray(control_values, dtype=float)
    control_index = np.asarray(control_index, dtype=np.int64)
    keep = ~np.isnan(control_values) & ~np.isnan(test_values[control_index])
    control_values, control_index = control_values[keep], control_index[keep]
    reference = test_values[control_index]
    wins = (reference > control_values) + 0.5 * (reference == control_values)
    counts = np.bincount(control_index, minlength=n_tests)
    score = np.full(n_tests, np.nan)
    has = counts > 0
    score[has] = np.bincount(control_index, weights=wins, minlength=n_tests)[has] / counts[has]
    control_mean = np.full(n_tests, np.nan)
    control_mean[has] = np.bincount(control_index, weights=control_values, minlength=n_tests)[has] / counts[has]
    valid = has & ~np.isnan(test_values)
    t, c, sc = test_values[valid], control_mean[valid], score[valid]
    result = {
        'n': int(valid.sum()),
        'test_mean': float(np.mean(t)) if len(t) else np.nan,
        'test_median': float(np.median(t)) if len(t) else np.nan,
        'control_mean': float(np.mean(c)) if len(c) else np.nan,
        'control_median': float(np.median(c)) if len(c) else np.nan,
        'mean_difference': float(np.mean(t - c)) if len(t) else np.nan,
        'fraction_above_control': float(np.mean(sc)) if len(sc) else np.nan,
    }
    if len(sc) >= 10 and np.any(sc != 0.5):
        result['wilcoxon_p'] = float(stats.wilcoxon(sc - 0.5, zero_method='zsplit').pvalue)
    else:
        result['wilcoxon_p'] = np.nan
    return result, control_mean


def validate_pairs(pairs, data, pools, excluded, label, n_controls=10, max_offset=None, seed=0,
                   coverage_matched=0, compartment=None, tad_id=None):
    """Metrics of test pairs and their matched controls; returns (per-pair table, summary table).

    coverage_matched > 0 additionally matches controls on that many coverage classes of their ends; `compartment`
    and `tad_id` (per bin, see domain_annotation) on the compartments of the ends and on whether both ends share a
    topological domain. The per-pair table holds the concordance of every pair (<metric>_concordance).
    """
    if len(pairs) == 0:
        return pairs, pd.DataFrame()
    if max_offset is None:
        max_offset = int(pairs['offset'].max())
    test = pair_metrics(data, pairs['bin_a'], pairs['bin_b'], max_offset)
    classes = coverage_classes(pools, coverage_matched) if coverage_matched else None
    if compartment is not None:
        classes = domain_classes(classes, compartment)
    controls = sample_matched_controls(pairs, pools, excluded, n_controls=n_controls, seed=seed,
                                       coverage_class=classes, tad_id=tad_id)
    control = pair_metrics(data, controls['bin_a'], controls['bin_b'], max_offset)
    table = pairs.copy()
    rows = []
    strata = distance_stratum(table['distance'].to_numpy())
    table['stratum'] = strata
    for metric, values in test.items():
        table[metric] = values
        summary, control_mean = summarize_against_controls(values, control[metric],
                                                           controls['test_index'].to_numpy(), len(table))
        table[metric + '_control'] = control_mean
        table[metric + '_concordance'] = pair_concordance(values, control[metric], controls['test_index'].to_numpy(),
                                                          len(table))
        rows.append(dict(set=label, pair_class='all', stratum='all', metric=metric, **summary))
        for pair_class in sorted(table['pair_class'].unique()):
            mask = (table['pair_class'] == pair_class).to_numpy()
            sub_index = np.nonzero(mask)[0]
            remap = -np.ones(len(table), dtype=np.int64)
            remap[sub_index] = np.arange(len(sub_index))
            keep = remap[controls['test_index'].to_numpy()] >= 0
            summary, _ = summarize_against_controls(values[mask], control[metric][keep],
                                                    remap[controls['test_index'].to_numpy()[keep]], len(sub_index))
            rows.append(dict(set=label, pair_class=pair_class, stratum='all', metric=metric, **summary))
            for stratum in sorted(np.unique(strata[mask])):
                mask_s = mask & (strata == stratum)
                sub_index = np.nonzero(mask_s)[0]
                remap = -np.ones(len(table), dtype=np.int64)
                remap[sub_index] = np.arange(len(sub_index))
                keep = remap[controls['test_index'].to_numpy()] >= 0
                summary, _ = summarize_against_controls(values[mask_s], control[metric][keep],
                                                        remap[controls['test_index'].to_numpy()[keep]],
                                                        len(sub_index))
                rows.append(dict(set=label, pair_class=pair_class, stratum=stratum, metric=metric, **summary))
    return table, pd.DataFrame(rows)


def _offset_bucket(offset, exact=16, ratio=1.2):
    """Offsets up to `exact` stay separate, larger ones are grouped in geometric buckets."""
    offset = np.asarray(offset, dtype=float)
    bucket = np.floor(np.log(np.maximum(offset, exact) / exact) / np.log(ratio)).astype(np.int64)
    return np.where(offset <= exact, offset.astype(np.int64), exact + 1 + bucket)


def confirm_pairs(pairs, data, pools, excluded, label, max_offset, n_background=20, n_classes=5, alpha=0.05,
                  fdr=0.1, seed=0):
    """Test every pair individually against background pairs of the same distance and coverage.

    The statistic is the co-occurrence of the two bins divided by the mean co-occurrence at their offset. The
    background of a pair are the matched controls (same pair class, coverage classes of both ends, offset) of
    all pairs in its stratum (pair class, offset bucket, coverage classes); its empirical P-value is
    (1 + background values >= statistic) / (1 + background size), so a pair that never co-occurs has P = 1.
    One further held-out control per pair passes the same test and estimates the rate of false confirmation.
    Returns (per-pair table, summary table).
    """
    if len(pairs) == 0:
        return pairs, pd.DataFrame()
    classes = coverage_classes(pools, n_classes)
    expected = data.expected_by_offset(max_offset, 'cooc')
    floor = np.nanmin(expected[1:][expected[1:] > 0]) if np.any(expected[1:] > 0) else 1.0
    expected = np.where(np.isfinite(expected) & (expected > 0), expected, floor)

    def statistic(bin_a, bin_b):
        bin_a, bin_b = np.asarray(bin_a, dtype=np.int64), np.asarray(bin_b, dtype=np.int64)
        offset = np.clip(np.abs(bin_b - bin_a), 0, max_offset)
        return data.cooccurrence(bin_a, bin_b, weighted=False), data.cooccurrence(bin_a, bin_b) / expected[offset]

    table = pairs.copy().reset_index(drop=True)
    n_groups, test_stat = statistic(table['bin_a'], table['bin_b'])
    background = sample_matched_controls(table, pools, excluded, n_controls=n_background, seed=seed,
                                         coverage_class=classes)
    held_out = sample_matched_controls(table, pools, excluded, n_controls=1, seed=seed + 7919,
                                       coverage_class=classes)
    _, background_stat = statistic(background['bin_a'], background['bin_b'])
    held_groups, held_stat = statistic(held_out['bin_a'], held_out['bin_b'])

    class_code = _class_code(classes, table['bin_a'].to_numpy(), table['bin_b'].to_numpy())
    kind = pd.factorize(table['pair_class'])[0]
    key = (kind * 1000 + _offset_bucket(table['offset'].to_numpy())) * 1000000 + class_code
    p_test = np.full(len(table), np.nan)
    p_held = np.full(len(table), np.nan)
    held_index = held_out['test_index'].to_numpy()
    background_key = key[background['test_index'].to_numpy()]
    order = np.argsort(background_key, kind='stable')
    sorted_keys = background_key[order]
    for k in np.unique(key):
        lo, hi = np.searchsorted(sorted_keys, k, 'left'), np.searchsorted(sorted_keys, k, 'right')
        if hi == lo:
            continue
        values = np.sort(background_stat[order[lo:hi]])
        n = len(values)
        members = np.nonzero(key == k)[0]
        p_test[members] = (1 + n - np.searchsorted(values, test_stat[members], 'left')) / (n + 1)
        held = np.nonzero(key[held_index] == k)[0]
        p_held[held_index[held]] = (1 + n - np.searchsorted(values, held_stat[held], 'left')) / (n + 1)
    p_test[n_groups == 0] = 1.0
    held_zero = np.zeros(len(table), dtype=bool)
    held_zero[held_index[held_groups == 0]] = True
    p_held[held_zero] = 1.0

    table['groups_both'] = n_groups
    table['confirm_p'] = p_test
    table['confirm_q'] = _bh(p_test)
    table['control_groups_both'] = np.nan
    table.loc[held_index, 'control_groups_both'] = held_groups
    table['control_confirm_p'] = p_held
    table['control_confirm_q'] = _bh(p_held)
    table['stratum'] = distance_stratum(table['distance'].to_numpy())

    def row(mask, pair_class, stratum):
        t = table[mask]
        ok = t['confirm_p'].notna() & t['control_confirm_p'].notna()
        t = t[ok]
        return dict(set=label, pair_class=pair_class, stratum=stratum, n=len(t),
                    observed=float(np.mean(t['groups_both'] > 0)) if len(t) else np.nan,
                    control_observed=float(np.mean(t['control_groups_both'] > 0)) if len(t) else np.nan,
                    confirmed_p=float(np.mean(t['confirm_p'] < alpha)) if len(t) else np.nan,
                    control_confirmed_p=float(np.mean(t['control_confirm_p'] < alpha)) if len(t) else np.nan,
                    confirmed_fdr=float(np.mean(t['confirm_q'] < fdr)) if len(t) else np.nan,
                    control_confirmed_fdr=float(np.mean(t['control_confirm_q'] < fdr)) if len(t) else np.nan)

    rows = [row(np.ones(len(table), dtype=bool), 'all', 'all')]
    for pair_class in sorted(table['pair_class'].unique()):
        mask = (table['pair_class'] == pair_class).to_numpy()
        rows.append(row(mask, pair_class, 'all'))
        for stratum in sorted(table.loc[mask, 'stratum'].unique()):
            rows.append(row(mask & (table['stratum'] == stratum).to_numpy(), pair_class, stratum))
    if 'weight' in table and table['weight'].nunique() > 4:
        quartile = pd.qcut(table['weight'].rank(method='first'), 4, labels=False).to_numpy()
        for q in range(4):
            rows.append(row(quartile == q, 'all', 'weight_quartile_{}'.format(q + 1)))
    log.info('%s: %.3f of pairs confirmed at P < %g (held-out controls %.3f)', label, rows[0]['confirmed_p'],
             alpha, rows[0]['control_confirmed_p'])
    return table, pd.DataFrame(rows)


def _bh(p):
    """Benjamini-Hochberg adjusted P-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        values = p[ok]
        order = np.argsort(values)
        ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(len(values))
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q


# --------------------------------------------------------------------- triplets
def triplet_cooperativity(data, bin_u, bin_v, bin_w):
    """Cooperativity ratio R = C(uvw) N(v) / (C(uv) C(vw)) on unweighted group counts (NaN if undefined)."""
    c_uvw = data.cooccurrence(bin_u, bin_v, bin_w, weighted=False)
    c_uv = data.cooccurrence(bin_u, bin_v, weighted=False)
    c_vw = data.cooccurrence(bin_v, bin_w, weighted=False)
    n_v = data.coverage(bin_v, weighted=False)
    denominator = c_uv * c_vw
    ratio = np.full(len(c_uvw), np.nan)
    ok = denominator > 0
    ratio[ok] = c_uvw[ok] * n_v[ok] / denominator[ok]
    return c_uvw, c_uv, c_vw, n_v, ratio


def validate_triplets(triplets, node_bins, data, pools, graph, n_controls=5, max_triplets=500000,
                      max_shift=20, seed=0):
    """Cooperativity of graph triplets versus control triplets with the same mediator."""
    rng = np.random.default_rng(seed)
    t = triplets.copy()
    t['bin_u'] = t['node_a'].map(node_bins)
    t['bin_v'] = t['mediator'].map(node_bins)
    t['bin_w'] = t['node_b'].map(node_bins)
    t = t.dropna(subset=['bin_u', 'bin_v', 'bin_w'])
    for column in ('bin_u', 'bin_v', 'bin_w'):
        t[column] = t[column].astype(np.int64)
    t = t[(t['bin_u'] != t['bin_v']) & (t['bin_v'] != t['bin_w']) & (t['bin_u'] != t['bin_w'])].reset_index(drop=True)
    if len(t) > max_triplets:
        t = t.sample(n=max_triplets, random_state=seed).reset_index(drop=True)
    if len(t) == 0:
        return t, pd.DataFrame()
    c_uvw, c_uv, c_vw, n_v, ratio = triplet_cooperativity(data, t['bin_u'], t['bin_v'], t['bin_w'])
    t['c_uvw'] = c_uvw
    t['c_uv'] = c_uv
    t['c_vw'] = c_vw
    t['n_v'] = n_v
    with np.errstate(divide='ignore'):
        t['log2_cooperativity'] = np.log2(ratio)

    # controls: keep the mediator bin, shift u and w jointly by a random number of bins so that both
    # offsets to v stay close to the test offsets; u and w must be graph loci that are not neighbours of v
    n_bins_total = len(pools.chrom_index)
    neighbour_codes = []
    for u_node, w_node in graph.edges:
        if u_node in node_bins and w_node in node_bins:
            a, b = node_bins[u_node], node_bins[w_node]
            neighbour_codes.append(a * n_bins_total + b)
            neighbour_codes.append(b * n_bins_total + a)
    neighbour_codes = np.unique(np.array(neighbour_codes, dtype=np.int64))
    all_loci = np.union1d(pools.bait_array, pools.other_array)
    chrom_of = pools.chrom_index
    tries = n_controls * 8
    index = np.repeat(np.arange(len(t)), tries)
    v = np.repeat(t['bin_v'].to_numpy(), tries)
    shift = rng.integers(-max_shift, max_shift + 1, size=len(index))
    u = np.repeat(t['bin_u'].to_numpy(), tries) + shift
    w = np.repeat(t['bin_w'].to_numpy(), tries) + shift
    ok = (shift != 0) & (u >= 0) & (w >= 0) & (u < n_bins_total) & (w < n_bins_total)
    ok &= (u != v) & (w != v) & (u != w)
    u_c, w_c = np.clip(u, 0, n_bins_total - 1), np.clip(w, 0, n_bins_total - 1)
    ok &= (chrom_of[u_c] == chrom_of[v]) & (chrom_of[w_c] == chrom_of[v])
    ok &= np.isin(u_c, all_loci) & np.isin(w_c, all_loci)
    ok &= ~np.isin(v * n_bins_total + u_c, neighbour_codes) & ~np.isin(v * n_bins_total + w_c, neighbour_codes)
    candidates = pd.DataFrame({'test_index': index[ok], 'bin_u': u_c[ok], 'bin_v': v[ok], 'bin_w': w_c[ok]})
    controls = candidates.groupby('test_index').head(n_controls).reset_index(drop=True)
    rows = []
    if len(controls):
        ctrl_c_uvw = data.cooccurrence(controls['bin_u'], controls['bin_v'], controls['bin_w'], weighted=False)
        for metric, test_values, control_values in (('c_uvw', t['c_uvw'].to_numpy(), ctrl_c_uvw),):
            finite_t = np.where(np.isfinite(test_values), test_values, np.nan)
            finite_c = np.where(np.isfinite(control_values), control_values, np.nan)
            ok_c = ~np.isnan(finite_c)
            summary, control_mean = summarize_against_controls(
                finite_t, finite_c[ok_c], controls['test_index'].to_numpy()[ok_c], len(t))
            t[metric + '_control'] = control_mean
            rows.append(dict(set='triplets', pair_class='all', stratum='all', metric=metric, **summary))
            for mediator_is_bait in (True, False):
                mask = (t['mediator_is_bait'] == mediator_is_bait).to_numpy()
                sub = np.nonzero(mask)[0]
                remap = -np.ones(len(t), dtype=np.int64)
                remap[sub] = np.arange(len(sub))
                index = controls['test_index'].to_numpy()[ok_c]
                keep = remap[index] >= 0
                summary, _ = summarize_against_controls(finite_t[mask], finite_c[ok_c][keep], remap[index[keep]], len(sub))
                rows.append(dict(set='triplets', pair_class='mediator_bait' if mediator_is_bait else 'mediator_oe',
                                 stratum='all', metric=metric, **summary))
    return t, pd.DataFrame(rows)


# ------------------------------------------------------------------ communities
def validate_communities(communities, node_bins, data, min_bins=3, min_members=3, n_shifts=20,
                         max_bins=200, seed=0, min_locus_fraction=0.5, bait_nodes=None):
    """Multi-way co-occurrence of community bin sets versus circularly shifted sets.

    Two nulls are drawn per community: 'shift' accepts any shift whose bins are
    covered by the multi-way data (>= 80 %), 'locus' additionally requires that
    at least `min_locus_fraction` of the shifted bins host graph loci, so the
    shifted set lies in captured chromatin like the community itself.
    """
    rng = np.random.default_rng(seed)
    covered = data.coverage_raw > 0
    locus = np.zeros(len(covered), dtype=bool)
    locus[np.fromiter(node_bins.values(), dtype=np.int64)] = True
    rows = []
    for community_id, members in enumerate(communities):
        bins = np.unique([node_bins[n] for n in members if n in node_bins])
        if len(bins) < min_bins or len(bins) > max_bins:
            continue
        chrom_index = data.bins.chrom_index_of_bins(bins[:1])[0]
        chrom = data.bins.chrom_names[chrom_index]
        if np.any(data.bins.chrom_index_of_bins(bins) != chrom_index):
            continue
        first = data.bins.offsets[chrom]
        n_chrom = data.bins.n_chrom_bins(chrom)
        observed = data.set_cooccurrence(bins, min_members)
        local = bins - first
        null = {'shift': [], 'locus': []}
        attempts = 0
        while (len(null['shift']) < n_shifts or len(null['locus']) < n_shifts) and attempts < n_shifts * 200:
            attempts += 1
            shift = int(rng.integers(1, n_chrom))
            shifted = (local + shift) % n_chrom + first
            if covered[shifted].mean() < 0.8:
                continue
            value = None
            if len(null['shift']) < n_shifts:
                value = data.set_cooccurrence(shifted, min_members)
                null['shift'].append(value)
            if len(null['locus']) < n_shifts and locus[shifted].mean() >= min_locus_fraction:
                null['locus'].append(value if value is not None else data.set_cooccurrence(shifted, min_members))
        row = {'community': community_id, 'chrom': chrom, 'n_nodes': len(members),
               'n_baits': int(sum(1 for n in members if bait_nodes is not None and n in bait_nodes)),
               'n_bins': len(bins), 'span_bins': int(bins.max() - bins.min()), 'observed': observed}
        for name, values in null.items():
            values = np.array(values, dtype=float)
            row['null_{}_n'.format(name)] = len(values)
            row['null_{}_mean'.format(name)] = float(values.mean()) if len(values) else np.nan
            row['log2_enrichment_{}'.format(name)] = float(np.log2((observed + 1e-3) / (values.mean() + 1e-3))) \
                if len(values) else np.nan
            row['empirical_p_{}'.format(name)] = float((np.sum(values >= observed) + 1) / (len(values) + 1)) \
                if len(values) else np.nan
            row['concordance_{}'.format(name)] = float(np.mean((observed > values) + 0.5 * (observed == values))) \
                if len(values) else np.nan
        rows.append(row)
    table = pd.DataFrame(rows)
    summary = []
    for name in ('shift', 'locus'):
        if not len(table):
            break
        sub = table.dropna(subset=['null_{}_mean'.format(name)])
        sub = sub[sub['null_{}_n'.format(name)] >= 5]
        if not len(sub):
            continue
        diff = np.log2(sub['observed'] + 1e-3) - np.log2(sub['null_{}_mean'.format(name)] + 1e-3)
        concordance = sub['concordance_{}'.format(name)].to_numpy()
        summary.append({'set': 'communities', 'pair_class': 'null_' + name, 'stratum': 'all',
                        'metric': 'groups_with_{}_members'.format(min_members), 'n': len(sub),
                        'test_mean': float(sub['observed'].mean()), 'test_median': float(sub['observed'].median()),
                        'control_mean': float(sub['null_{}_mean'.format(name)].mean()),
                        'control_median': float(sub['null_{}_mean'.format(name)].median()),
                        'mean_difference': float(diff.mean()),
                        'fraction_above_control': float(np.mean(concordance)),
                        'wilcoxon_p': float(stats.wilcoxon(concordance - 0.5, zero_method='zsplit').pvalue)
                        if len(sub) >= 10 and np.any(concordance != 0.5) else np.nan})
    return table, pd.DataFrame(summary)
