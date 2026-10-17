"""Time course of the explanation cascade over ordered conditions (e.g. differentiation days).

egaDifferential compares two conditions. With three or more ordered conditions (days d_0 < d_1 < ... ) the pairwise
runs of consecutive days and of the first and last day are combined here:

* Trajectories. All days share one graph, the union of the called contacts of every day (edge strength r_d per day,
  0 where the pair is not called at day d). Every graph step is defined on this union graph, so a node is at the same
  step of a gene at every day. Per gene node and day: expression (log2 of the mean normalised count + 1), promoter
  mark signal, summed contact strength of the promoter, and the contact-weighted mark activity at step k
  (H_d^(k) @ signal_d, the operator of egaDifferential), each as log2(value + 1). The pairwise runs add whether each
  level explains the gene in each pair.
* Ordering of events. Units are a gene with one of its explaining nodes (union over the pairwise runs; the promoter
  fragment is the unit of step 0). Per unit three layers: contact (path weight of the element at its step, summed
  promoter contact strength at step 0), the mark at the element and the expression of the gene. A layer changes in
  interval i (d_i -> d_i+1) when |log2((v_i+1 + 1) / (v_i + 1))| reaches its threshold (expression: |log2FC| at least
  the RNA threshold and padj <= alpha in the contrast of the interval); its event time is the first such interval.
  Units with an event in every layer of a layer set are classified (three layers contact, mark, expression: forward
  = contact <= mark <= expression and not all equal; simultaneous; reversed = contact >= mark >= expression and not
  all equal; other) and compared with a permutation null that shuffles the event intervals among the layers within
  each gene (the expression slot is shared by the gene's units), which keeps how early or late a gene changes and
  tests only the order of the layers. A second null shuffles each layer's intervals across units and keeps how
  often each layer changes in each interval: an order caused by one layer changing early everywhere (e.g. through
  a technical difference of one sample) is part of it. With two intervals the order is resolved only as 'same interval' versus
  'earlier interval'; changes within an interval are simultaneous by construction.
* Knockdown test. For a factor F knocked down: H3K27ac change (knockdown versus control, log2 with pseudo count 1 on
  the genome-mean scale) at the explaining nodes of the genes whose candidate trigger (egaTriggers) is F, against the
  explaining nodes of the other explained genes; Mann-Whitney U, difference of medians and Cliff's delta with
  intervals from a bootstrap over genes.
* UMI-4C. Viewpoints (bins with the most read pairs) that overlap a gene node with explaining nodes: contact
  frequency of the viewpoint with the explaining nodes per condition.
"""

import logging
import re

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)

ORDER3 = ('forward', 'simultaneous', 'reversed', 'other')
ORDER2 = ('before', 'same', 'after')
LAYER_SETS = {
    'contact_mark_expression': ('contact', 'mark', 'expression'),
    'contact_expression': ('contact', 'expression'),
    'mark_expression': ('mark', 'expression'),
    'contact_mark': ('contact', 'mark'),
}


def intervals_of(days):
    """Consecutive pairs of ordered days."""
    return list(zip(days[:-1], days[1:]))


def log2_change(before, after, pseudo_count=1.0):
    return np.log2((np.asarray(after, dtype=float) + pseudo_count) / (np.asarray(before, dtype=float) + pseudo_count))


# --------------------------------------------------------------------------------------------------------------------
# events and orderings
# --------------------------------------------------------------------------------------------------------------------

def first_interval(changes, threshold, significant=None):
    """Event time per row: index of the first interval (column) with |change| >= threshold (and `significant` True
    where given), -1 when no interval qualifies. `changes` is (n_rows, n_intervals)."""
    changes = np.atleast_2d(np.asarray(changes, dtype=float))
    hit = np.abs(np.nan_to_num(changes, nan=0.0)) >= threshold
    if significant is not None:
        hit &= np.atleast_2d(np.asarray(significant, dtype=bool))
    first = np.argmax(hit, axis=1)
    return np.where(hit.any(axis=1), first, -1)


def classify_order(times):
    """Order category codes of event times (n_units, n_layers), every time >= 0.

    Two layers (x, y): 0 before (x earlier than y), 1 same, 2 after. Three layers (contact, mark, expression):
    0 forward (contact <= mark <= expression, not all equal), 1 simultaneous, 2 reversed (contact >= mark >=
    expression, not all equal), 3 other."""
    t = np.asarray(times)
    if t.ndim != 2 or t.shape[1] not in (2, 3):
        raise ValueError('times need 2 or 3 columns')
    if t.shape[1] == 2:
        return np.where(t[:, 0] < t[:, 1], 0, np.where(t[:, 0] == t[:, 1], 1, 2))
    c, m, e = t[:, 0], t[:, 1], t[:, 2]
    equal = (c == m) & (m == e)
    forward = (c <= m) & (m <= e) & ~equal
    reverse = (c >= m) & (m >= e) & ~equal
    return np.where(equal, 1, np.where(forward, 0, np.where(reverse, 2, 3)))


def category_names(n_layers):
    return ORDER2 if n_layers == 2 else ORDER3


def order_slots(units, layers, gene_column='gene'):
    """Slot layout of the within-gene permutation for units with an event in every layer.

    Returns (labels, groups, slot_index, genes): `labels` the observed event times of every slot, `groups` the gene
    code of every slot (sorted), `slot_index` (n_units, n_layers) the slot of each unit and layer (the expression
    slot is shared by all units of a gene), `genes` the gene code of every unit."""
    units = units.reset_index(drop=True)
    gene_codes, gene_names = pd.factorize(units[gene_column], sort=True)
    order = np.argsort(gene_codes, kind='stable')
    units = units.iloc[order].reset_index(drop=True)
    gene_codes = gene_codes[order]
    labels, groups = [], []
    slot_index = np.zeros((len(units), len(layers)), dtype=np.int64)
    expression_slot = {}
    for i, (code, row) in enumerate(zip(gene_codes, units.itertuples(index=False))):
        row = row._asdict()
        for j, layer in enumerate(layers):
            if layer == 'expression':
                if code not in expression_slot:
                    expression_slot[code] = len(labels)
                    labels.append(int(row['t_expression']))
                    groups.append(code)
                slot_index[i, j] = expression_slot[code]
            else:
                slot_index[i, j] = len(labels)
                labels.append(int(row['t_' + layer]))
                groups.append(code)
    labels = np.array(labels, dtype=np.int64)
    groups = np.array(groups, dtype=np.int64)
    # slots of one gene are contiguous already (units sorted by gene, slots appended in unit order)
    return labels, groups, slot_index, gene_codes, units


def shuffle_within_groups(labels, groups, rng):
    """Permutation of `labels` within each group (groups sorted ascending, contiguous)."""
    order = np.lexsort((rng.random(len(labels)), groups))
    return labels[order]


def order_test(units, layers, n_permutations=10000, n_bootstrap=2000, seed=0, gene_column='gene'):
    """Counts of the order categories of `units` (columns t_<layer>, gene) against two permutation nulls.

    Units need an event (t >= 0) in every layer. Within-gene null (columns null_*, ratio, p): the event intervals of a
    gene are shuffled among its layer slots, which keeps how early a gene changes and makes the layers exchangeable
    within the gene. Layer-marginal null (columns marginal_*): the event intervals of each layer are shuffled across
    units (expression across genes), which keeps how often each layer changes in each interval and removes only the
    coupling of the layers within a unit; an order produced by one layer changing earlier everywhere (e.g. a
    technical difference of one day) is part of this null. Returns one row per category: observed count and
    fraction, bootstrap interval over genes (2.5 / 97.5 percentiles), null means and 2.5 / 97.5 percentiles, ratios
    observed / null mean and two-sided permutation p (deviation from the null mean)."""
    columns = ['t_' + layer for layer in layers]
    names = category_names(len(layers))
    eligible = units[(units[columns] >= 0).all(axis=1)]
    empty = pd.DataFrame({'category': names, 'units': len(eligible), 'genes': 0, 'observed': 0,
                          'fraction': np.nan, 'boot_low': np.nan, 'boot_high': np.nan, 'null_mean': np.nan,
                          'null_low': np.nan, 'null_high': np.nan, 'ratio': np.nan, 'p': np.nan,
                          'marginal_mean': np.nan, 'marginal_low': np.nan, 'marginal_high': np.nan,
                          'marginal_ratio': np.nan, 'marginal_p': np.nan})
    if len(eligible) == 0:
        return empty
    labels, groups, slot_index, gene_codes, ordered = order_slots(eligible, layers, gene_column)
    observed_codes = classify_order(labels[slot_index])
    n_cat = len(names)
    observed = np.bincount(observed_codes, minlength=n_cat)
    rng = np.random.default_rng(seed)
    null = np.zeros((n_permutations, n_cat))
    marginal = np.zeros((n_permutations, n_cat))
    layer_slots = [np.unique(slot_index[:, j]) for j in range(len(layers))]
    for b in range(n_permutations):
        permuted = shuffle_within_groups(labels, groups, rng)
        null[b] = np.bincount(classify_order(permuted[slot_index]), minlength=n_cat)
        permuted = labels.copy()
        for idx in layer_slots:
            permuted[idx] = labels[idx][rng.permutation(len(idx))]
        marginal[b] = np.bincount(classify_order(permuted[slot_index]), minlength=n_cat)
    n_genes = int(gene_codes.max()) + 1
    per_gene = np.zeros((n_genes, n_cat))
    np.add.at(per_gene, (gene_codes, observed_codes), 1)
    boot = np.zeros((n_bootstrap, n_cat))
    for b in range(n_bootstrap):
        boot[b] = per_gene[rng.integers(0, n_genes, n_genes)].sum(axis=0)
    total = observed.sum()

    def summary(values):
        mean = values.mean(axis=0)
        deviation = np.abs(values - mean)
        p = (1 + (deviation >= np.abs(observed - mean) - 1e-9).sum(axis=0)) / (1 + n_permutations)
        ratio = np.where(mean > 0, observed / np.where(mean > 0, mean, 1), np.nan)
        return mean, np.percentile(values, 2.5, axis=0), np.percentile(values, 97.5, axis=0), ratio, p

    null_mean, null_low, null_high, ratio, p = summary(null)
    m_mean, m_low, m_high, m_ratio, m_p = summary(marginal)
    return pd.DataFrame({
        'category': names, 'units': len(eligible), 'genes': n_genes, 'observed': observed,
        'fraction': observed / total, 'boot_low': np.percentile(boot, 2.5, axis=0),
        'boot_high': np.percentile(boot, 97.5, axis=0), 'null_mean': null_mean,
        'null_low': null_low, 'null_high': null_high, 'ratio': ratio, 'p': p, 'marginal_mean': m_mean,
        'marginal_low': m_low, 'marginal_high': m_high, 'marginal_ratio': m_ratio, 'marginal_p': m_p})


def event_patterns(units, layers):
    """Counts of the joint event-time patterns (e.g. contact 0, mark 0, expression 1) of units with all events."""
    columns = ['t_' + layer for layer in layers]
    eligible = units[(units[columns] >= 0).all(axis=1)]
    if not len(eligible):
        return pd.DataFrame(columns=columns + ['units'])
    return eligible.groupby(columns).size().rename('units').reset_index()


# --------------------------------------------------------------------------------------------------------------------
# graph over all days
# --------------------------------------------------------------------------------------------------------------------

def day_tables(paths_by_day, file_format='ibed', baits=None, mode='all', threshold=None):
    """{day: consensus table of the day's replicate files} (differential.consensus_scores)."""
    from .differential import consensus_scores
    from .io import read_interactions
    result = {}
    for day, paths in paths_by_day.items():
        tables = [read_interactions([p], file_format, baits=baits) for p in paths]
        if len(tables) == 1 and threshold is not None:
            result[day] = consensus_scores([tables[0][tables[0]['score'] >= threshold]], 'union')
        else:
            result[day] = consensus_scores(tables, mode, threshold)
    return result


def timecourse_graph(tables, reference_score=5.0):
    """Union graph of the called pairs of every day; edge attribute r_<day> = score / reference_score (0 where not
    called), node attributes chrom, start, end."""
    graph = nx.Graph()
    days = list(tables)
    for day, table in tables.items():
        for row in table.itertuples(index=False):
            a = '{}:{}-{}'.format(row.chrom, int(row.bait_start), int(row.bait_end))
            b = '{}:{}-{}'.format(row.chrom, int(row.other_start), int(row.other_end))
            for node, start, end in ((a, row.bait_start, row.bait_end), (b, row.other_start, row.other_end)):
                if node not in graph:
                    graph.add_node(node, chrom=row.chrom, start=int(start), end=int(end))
            if not graph.has_edge(a, b):
                graph.add_edge(a, b, **{'r_' + d: 0.0 for d in days})
            graph.edges[a, b]['r_' + day] = float(row.score) / reference_score
    return graph


def day_operators(graph, sources, days, max_hops=3):
    """Per day the k-step path weight operators {day: {k: H}} on the union graph (differential.ActivityOperators, two
    days per construction). Returns (operators of the last construction for the node / source order, {day: hops})."""
    from .differential import ActivityOperators
    hops = {}
    ops = None
    pending = list(days)
    while pending:
        first = pending.pop(0)
        second = pending.pop(0) if pending else first
        for _, _, data in graph.edges(data=True):
            data['r_a'] = data['r_' + first]
            data['r_b'] = data['r_' + second]
        ops = ActivityOperators(graph, sources, max_hops=max_hops)
        hops[first] = ops.hop_a
        hops[second] = ops.hop_b
    for _, _, data in graph.edges(data=True):
        data.pop('r_a', None)
        data.pop('r_b', None)
    return ops, hops


def node_signals(graph_nodes, node_table, tracks, threads=1):
    """{mark: {day: array over graph_nodes}}: mean signal over every node fragment, averaged over the files of a day.
    `tracks` is {mark: {day: [paths]}}; `node_table` holds chrom, start, end indexed by node."""
    from .differential import _interval_signals
    chroms = node_table.loc[graph_nodes, 'chrom'].tolist()
    starts = node_table.loc[graph_nodes, 'start'].to_numpy(dtype=np.int64)
    ends = node_table.loc[graph_nodes, 'end'].to_numpy(dtype=np.int64)
    paths = list(dict.fromkeys(p for by_day in tracks.values() for files in by_day.values() for p in files))
    read = _interval_signals(paths, chroms, starts, ends, threads)
    return {mark: {day: np.mean([read[p] for p in files], axis=0) for day, files in by_day.items()}
            for mark, by_day in tracks.items()}


# --------------------------------------------------------------------------------------------------------------------
# expression
# --------------------------------------------------------------------------------------------------------------------

def node_gene_choice(node_genes, tables, rna_threshold=1.0, alpha=0.05):
    """Per gene node one gene: among the node's genes present in the expression tables, the one with the largest
    |log2FC| over all contrasts (significant genes first). `tables` {contrast: DataFrame indexed by gene name with
    log2FoldChange, padj}."""
    lfc = pd.concat({c: t['log2FoldChange'] for c, t in tables.items()}, axis=1)
    padj = pd.concat({c: t['padj'] for c, t in tables.items()}, axis=1)
    significant = ((padj <= alpha) & (lfc.abs() >= rna_threshold)).any(axis=1)
    magnitude = lfc.abs().max(axis=1)
    choice = {}
    for node, names in node_genes.items():
        candidates = [g for g in str(names).split(',') if g and g in lfc.index]
        if not candidates:
            continue
        choice[node] = max(candidates, key=lambda g: (bool(significant.get(g, False)), float(magnitude.get(g, 0.0))))
    return choice


def expression_pattern(lfc, padj, rna_threshold=1.0, alpha=0.05):
    """Per gene and interval '+', '-' or '0' (significant change of at least the threshold), joined."""
    lfc = np.atleast_2d(lfc)
    padj = np.atleast_2d(padj)
    up = (lfc >= rna_threshold) & (padj <= alpha)
    down = (lfc <= -rna_threshold) & (padj <= alpha)
    symbols = np.where(up, '+', np.where(down, '-', '0'))
    return [''.join(row) for row in symbols]


# --------------------------------------------------------------------------------------------------------------------
# knockdown
# --------------------------------------------------------------------------------------------------------------------

def cliffs_delta(x, y):
    """P(x > y) - P(x < y) from the Mann-Whitney U of x against y."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not len(x) or not len(y):
        return np.nan
    u = stats.mannwhitneyu(x, y, alternative='two-sided').statistic
    return 2.0 * u / (len(x) * len(y)) - 1.0


def knockdown_test(values, node_genes_candidate, node_genes_control, n_bootstrap=2000, seed=0):
    """Compare a per-node knockdown effect between the explaining nodes of candidate genes and of control genes.

    `values` {node: effect}; `node_genes_*` DataFrames with columns gene_node, node (explaining nodes per gene). Nodes
    of candidate genes are removed from the control group. Returns a dict with sizes, medians, the difference of
    medians (candidate - control), Cliff's delta, both with 95% intervals from a bootstrap over genes (each group
    resampled by gene, nodes of a resampled gene taken with it), and the Mann-Whitney p (two-sided and candidate
    lower)."""
    cand = node_genes_candidate[node_genes_candidate['node'].isin(values)]
    cand_nodes = set(cand['node'])
    ctrl = node_genes_control[node_genes_control['node'].isin(values) & ~node_genes_control['node'].isin(cand_nodes)]
    x = np.array([values[n] for n in sorted(cand_nodes)], dtype=float)
    y = np.array([values[n] for n in sorted(set(ctrl['node']))], dtype=float)
    out = {'genes_candidate': cand['gene_node'].nunique(), 'genes_control': ctrl['gene_node'].nunique(),
           'nodes_candidate': len(x), 'nodes_control': len(y)}
    if len(x) < 3 or len(y) < 3:
        return dict(out, median_candidate=np.nan, median_control=np.nan, median_difference=np.nan, md_low=np.nan,
                    md_high=np.nan, cliffs_delta=np.nan, delta_low=np.nan, delta_high=np.nan, p_two_sided=np.nan,
                    p_less=np.nan)
    p_two = stats.mannwhitneyu(x, y, alternative='two-sided').pvalue
    p_less = stats.mannwhitneyu(x, y, alternative='less').pvalue
    rng = np.random.default_rng(seed)

    def grouped(table):
        nodes_of = table.groupby('gene_node')['node'].apply(lambda s: np.array([values[n] for n in s], dtype=float))
        return list(nodes_of)

    # every node in one cluster only (its first gene), so a resample holds each node at most as often as its gene
    cand_groups, ctrl_groups = grouped(cand.drop_duplicates('node')), grouped(ctrl.drop_duplicates('node'))
    md, delta = np.zeros(n_bootstrap), np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        bx = np.concatenate([cand_groups[i] for i in rng.integers(0, len(cand_groups), len(cand_groups))])
        by = np.concatenate([ctrl_groups[i] for i in rng.integers(0, len(ctrl_groups), len(ctrl_groups))])
        md[b] = np.median(bx) - np.median(by)
        delta[b] = cliffs_delta(bx, by)
    return dict(out, median_candidate=float(np.median(x)), median_control=float(np.median(y)),
                median_difference=float(np.median(x) - np.median(y)), md_low=float(np.percentile(md, 2.5)),
                md_high=float(np.percentile(md, 97.5)), cliffs_delta=float(cliffs_delta(x, y)),
                delta_low=float(np.percentile(delta, 2.5)), delta_high=float(np.percentile(delta, 97.5)),
                p_two_sided=float(p_two), p_less=float(p_less))


def bigwig_scale(path):
    """1 / genome-wide mean of a bigWig (header sumData over the summed chromosome lengths)."""
    import pyBigWig
    handle = pyBigWig.open(path)
    try:
        return sum(handle.chroms().values()) / float(handle.header()['sumData'])
    finally:
        handle.close()


# --------------------------------------------------------------------------------------------------------------------
# UMI-4C
# --------------------------------------------------------------------------------------------------------------------

UMI4C_NAME = re.compile(r'UMI4C_KD_(?P<condition>[A-Za-z0-9]+?)i_Day(?P<day>\d+)_Rep(?P<rep>\d+)')


def read_umi4c(path):
    """UMI-4C pair counts (BED6-like: chrom1 start1 end1 chrom2 start2 end2 count) between fixed bins."""
    table = pd.read_csv(path, sep='\t', header=None, usecols=range(7), compression='infer',
                        names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count'])
    return table


def bin_totals(table):
    """Read pairs per bin (both ends)."""
    a = table.groupby(['chrom1', 'start1', 'end1'])['count'].sum()
    a.index.names = ['chrom', 'start', 'end']
    b = table.groupby(['chrom2', 'start2', 'end2'])['count'].sum()
    b.index.names = ['chrom', 'start', 'end']
    return a.add(b, fill_value=0)


def find_viewpoints(totals, min_fraction=0.002, flank=100000):
    """Viewpoint bins: at least `min_fraction` of all pair ends and the largest bin within +-`flank`."""
    totals = totals.sort_values(ascending=False)
    limit = min_fraction * totals.sum()
    chosen = []
    for (chrom, start, end), value in totals.items():
        if value < limit:
            break
        if any(c == chrom and abs(s - start) <= flank for c, s, _ in chosen):
            continue
        chosen.append((chrom, int(start), int(end)))
    return chosen


def viewpoint_profile(table, viewpoint, window=1000000):
    """Reads of the viewpoint bin with every other cis bin within +-window: Series indexed by (start, end)."""
    chrom, start, end = viewpoint
    first = (table['chrom1'] == chrom) & (table['start1'] == start)
    second = (table['chrom2'] == chrom) & (table['start2'] == start)
    a = table[first & (table['chrom2'] == chrom)][['start2', 'end2', 'count']].set_axis(['start', 'end', 'count'],
                                                                                        axis=1)
    b = table[second & (table['chrom1'] == chrom)][['start1', 'end1', 'count']].set_axis(['start', 'end', 'count'],
                                                                                         axis=1)
    both = pd.concat([a, b])
    both = both[(both['start'] != start) & ((both['start'] - start).abs() <= window)]
    return both.groupby(['start', 'end'])['count'].sum()


def profile_at(profile, start, end):
    """Sum of profile reads of bins overlapping [start, end)."""
    return profile_over(profile, [(start, end)])


def profile_over(profile, intervals):
    """Sum of profile reads of the bins overlapping any of `intervals` [(start, end)], every bin counted once."""
    if not len(profile) or not len(intervals):
        return 0.0
    s = profile.index.get_level_values(0).to_numpy()
    e = profile.index.get_level_values(1).to_numpy()
    hit = np.zeros(len(profile), dtype=bool)
    for start, end in intervals:
        hit |= (s < end) & (e > start)
    return float(profile.to_numpy()[hit].sum())
