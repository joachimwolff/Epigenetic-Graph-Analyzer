"""Differential analysis of two capture Hi-C graphs with chromatin and expression changes.

Edges
    The union graph holds the interactions called (e.g. by CHiCAGO) in either
    condition. The strength of an edge in a condition is its score divided by
    the calling threshold (1 at the threshold, 0 where the pair was not
    called). A Chicdiff results table, when given, adds the tested change of
    each interaction (log2 fold change and adjusted p-value).

Genes
    For every bait gene the chromatin change around it is summarised on the
    promoter and on k = 1 .. n hops of the union graph (n = max_hops, default
    2), each weighting a node's mark signal by the edge strength r in the
    respective condition:

    * promoter:          the mark at the bait fragment itself,
    * direct (1 hop):    sum over neighbours n of r(bait, n) * signal(n),
    * indirect (2 hops): sum over nodes m at exactly two hops of
                         sum over mediators n of r(bait, n) * r(n, m) * signal(m),
    * hop<k> (k hops):   sum over nodes m at exactly k hops (shortest-path
                         distance k) of the sum over the shortest paths from the
                         bait to m of the product of r along the path, times signal(m).

    The feature is log2(activity_B / activity_A). Expression changes are
    related to these features by nested linear models (one level added at a
    time) and by a hierarchical classification (promoter, 1 hop, 2 hops, ...,
    n hops, unexplained). The contribution of the levels at k >= 2 hops is
    assessed by permuting node signals.
"""

import logging
import multiprocessing
import os
import re
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)

PAIR_KEY = ['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end']


def _key(df):
    return (df['chrom'] + ':' + df['bait_start'].astype(str) + '-' + df['bait_end'].astype(str) + '|' +
            df['chrom'] + ':' + df['other_start'].astype(str) + '-' + df['other_end'].astype(str))


REPLICATE_MODES = ('all', 'intersect', 'union')


def consensus_scores(tables, mode='all', threshold=None):
    """Combine replicate score tables of one condition into one table with the mean score per pair.

    'all' keeps pairs listed in every replicate whose score reaches `threshold` in at least one of them (the other
    replicates may hold the pair below the threshold, e.g. CHiCAGO exports with a lower cutoff); 'intersect' keeps
    pairs reaching `threshold` in every replicate, 'union' pairs reaching it in any. Without `threshold` every
    listed pair counts as called, so 'all' and 'intersect' coincide (published files hold called pairs only).
    """
    if mode not in REPLICATE_MODES:
        raise ValueError('Unknown replicate mode {}'.format(mode))
    keyed = []
    for table in tables:
        t = table.assign(key=_key(table)).sort_values('score', ascending=False).drop_duplicates('key')
        keyed.append(t.set_index('key'))
    listed = keyed[0].index
    for t in keyed[1:]:
        listed = listed.union(t.index)
    scores = pd.concat([t['score'].reindex(listed) for t in keyed], axis=1)
    present = scores.notna()
    called = scores >= threshold if threshold is not None else present
    if mode == 'all':
        keep = present.all(axis=1) & called.any(axis=1)
    elif mode == 'intersect':
        keep = called.all(axis=1)
    else:
        keep = called.any(axis=1)
    keys = listed[keep.to_numpy()]
    base = pd.concat([t.drop(columns=['score']) for t in keyed])
    base = base[~base.index.duplicated()].reindex(keys)
    base['score'] = scores.loc[keys].mean(axis=1, skipna=True).to_numpy()
    base['n_replicates_called'] = called.loc[keys].sum(axis=1).to_numpy()
    base['significant'] = True
    return base.reset_index(drop=True)


def differential_edges_from_scores(table_a, table_b, name_a='A', name_b='B', reference_score=5.0):
    """Union of significant pairs of two conditions described by interaction scores only (e.g. CHiCAGO).

    The edge strength of a condition is score / reference_score (1 at the
    calling threshold) and 0 where the pair was not called. No test is
    possible without counts; the status records gained, lost and shared pairs.
    """
    a = table_a.assign(key=_key(table_a)).drop_duplicates('key').set_index('key')
    b = table_b.assign(key=_key(table_b)).drop_duplicates('key').set_index('key')
    keys = a.index.union(b.index)
    base = pd.concat([a, b])
    base = base[~base.index.duplicated()].reindex(keys)
    score_a = a['score'].reindex(keys).fillna(0).to_numpy(dtype=float)
    score_b = b['score'].reindex(keys).fillna(0).to_numpy(dtype=float)
    sig_a = score_a > 0
    sig_b = score_b > 0
    out = pd.DataFrame({
        'chrom': base['chrom'].to_numpy(),
        'node_a': (base['chrom'] + ':' + base['bait_start'].astype(np.int64).astype(str) + '-' +
                   base['bait_end'].astype(np.int64).astype(str)).to_numpy(),
        'node_b': (base['chrom'] + ':' + base['other_start'].astype(np.int64).astype(str) + '-' +
                   base['other_end'].astype(np.int64).astype(str)).to_numpy(),
        'other_is_bait': base['other_is_bait'].astype(bool).to_numpy(),
        'distance': base['distance'].to_numpy(),
        'score_' + name_a: score_a, 'score_' + name_b: score_b,
        'enrichment_' + name_a: score_a / reference_score, 'enrichment_' + name_b: score_b / reference_score,
        'log2FC': np.log2((score_b + 1) / (score_a + 1)),
        'significant_' + name_a: sig_a, 'significant_' + name_b: sig_b,
        'status': np.where(sig_a & sig_b, 'shared', np.where(sig_a, 'lost', 'gained')),
        'pvalue': np.nan, 'qvalue': np.nan,
    })
    return out


def annotate_chicdiff(edges, chicdiff, padj_threshold=0.05):
    """Add the Chicdiff test of each union edge (bait identical, other end overlapping a Chicdiff region).

    Edges overlapping several regions take the one with the smallest adjusted p-value. Columns log2FC,
    pvalue and qvalue are replaced where a test exists; chicdiff_significant marks padj < `padj_threshold`.
    """
    edges = edges.copy()
    bait = edges['node_a'].str.split(r'[:-]', expand=True)
    other = edges['node_b'].str.split(r'[:-]', expand=True)
    frame = pd.DataFrame({'chrom': bait[0], 'bait_start': bait[1].astype(np.int64), 'bait_end': bait[2].astype(np.int64),
                          'other_start': other[1].astype(np.int64), 'other_end': other[2].astype(np.int64),
                          'edge': np.arange(len(edges))})
    merged = frame.merge(chicdiff, on=['chrom', 'bait_start', 'bait_end'], suffixes=('', '_region'))
    merged = merged[(merged['other_start'] <= merged['other_end_region']) &
                    (merged['other_end'] >= merged['other_start_region'])]
    merged = merged.sort_values('padj').drop_duplicates('edge')
    tested = merged.set_index('edge')
    edges['chicdiff_tested'] = False
    edges.loc[tested.index, 'chicdiff_tested'] = True
    edges.loc[tested.index, 'log2FC'] = tested['log2FC'].to_numpy()
    edges.loc[tested.index, 'pvalue'] = tested['pvalue'].to_numpy()
    edges.loc[tested.index, 'qvalue'] = tested['padj'].to_numpy()
    edges['chicdiff_significant'] = edges['chicdiff_tested'] & (edges['qvalue'] < padj_threshold)
    log.info('Chicdiff: %d of %d union edges tested, %d with padj < %g', len(tested), len(edges),
             int(edges['chicdiff_significant'].sum()), padj_threshold)
    return edges


def union_graph(edges, name_a='A', name_b='B', bait_nodes=None):
    """Graph of the union of significant edges, carrying both enrichments."""
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        row = row._asdict()
        graph.add_edge(row['node_a'], row['node_b'],
                       r_a=float(row['enrichment_' + name_a]), r_b=float(row['enrichment_' + name_b]),
                       log2FC=float(row['log2FC']), qvalue=float(row['qvalue']), status=row['status'])
    baits = set(edges['node_a']) | set(edges.loc[edges['other_is_bait'], 'node_b'])
    if bait_nodes is not None:
        baits |= set(bait_nodes)
    for node in graph.nodes:
        chrom, rest = node.split(':')
        start, end = rest.split('-')
        graph.nodes[node].update(chrom=chrom, start=int(start), end=int(end), bait=node in baits)
    return graph


def level_name(hops):
    """Name of the level at `hops` graph steps: 1 'direct', 2 'indirect', k >= 3 'hop<k>'."""
    if hops == 1:
        return 'direct'
    if hops == 2:
        return 'indirect'
    return 'hop{}'.format(hops)


def level_hops(level):
    """Inverse of level_name; the promoter level has 0 hops."""
    if level == 'promoter':
        return 0
    if level == 'direct':
        return 1
    if level == 'indirect':
        return 2
    return int(level[3:])


def hop_levels(features, marks):
    """Graph levels present in a feature table: direct, indirect and hop<k> in increasing order of k."""
    found = set()
    for column in features.columns:
        for mark in marks:
            match = re.fullmatch(re.escape(mark) + r'_hop(\d+)', column)
            if match:
                found.add(int(match.group(1)))
    return ['direct', 'indirect'] + ['hop{}'.format(k) for k in sorted(found) if k >= 3]


class ActivityOperators:
    """Sparse operators mapping a node signal vector to the k-hop activities of source nodes, k = 1 .. max_hops.

    The activity of source s at k hops in condition c is H_c^(k) @ signal_c with
    H_c^(k)[s, m] = sum over the shortest paths s = n_0, n_1, ..., n_k = m of prod_i r_c(n_{i-1}, n_i),
    restricted to nodes m at exactly k steps (shortest-path distance k in the union graph). It is computed as
    H^(1) = W[sources] and H^(k) = (H^(k-1) @ W) masked to the nodes at distance k: a walk of k steps that ends at
    distance k is a shortest path, and its first k-1 steps are a shortest path to a node at distance k-1, so the
    recursion sums exactly the path products. direct_c = H_c^(1) and indirect_c = H_c^(2) (sum over mediators n of
    r_c(s, n) r_c(n, m), m at exactly two hops) as before; hop_a[k], hop_b[k] and n_hop[k] hold every level.
    """

    def __init__(self, graph, sources, max_hops=2):
        import scipy.sparse as sp
        if max_hops < 2:
            raise ValueError('max_hops must be at least 2')
        self.max_hops = int(max_hops)
        self.levels = [level_name(k) for k in range(1, self.max_hops + 1)]
        self.nodes = list(graph.nodes)
        self.index = {n: i for i, n in enumerate(self.nodes)}
        self.sources = [s for s in sources if s in self.index]
        n = len(self.nodes)
        rows, cols, w_a, w_b = [], [], [], []
        for u, v, data in graph.edges(data=True):
            i, j = self.index[u], self.index[v]
            rows += [i, j]
            cols += [j, i]
            w_a += [data['r_a'], data['r_a']]
            w_b += [data['r_b'], data['r_b']]
        W_a = sp.csr_matrix((w_a, (rows, cols)), shape=(n, n))
        W_b = sp.csr_matrix((w_b, (rows, cols)), shape=(n, n))
        A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
        src = np.array([self.index[s] for s in self.sources], dtype=np.int64)
        A_s = A[src]
        self_mask = sp.csr_matrix((np.ones(len(src)), (np.arange(len(src)), src)), shape=(len(src), n))
        # nodes at distance < k of every source (the source itself and its direct neighbours to start with)
        visited = A_s + self_mask
        visited.data[:] = 1.0
        frontier = A_s
        self.hop_a = {1: W_a[src]}
        self.hop_b = {1: W_b[src]}
        self.n_hop = {}
        for k in range(2, self.max_hops + 1):
            reach = (frontier @ A).tocsr()
            reach.data[:] = 1.0
            # nodes at exactly k steps: reachable in k steps and not at a smaller distance
            mask = reach - reach.multiply(visited)
            mask.eliminate_zeros()
            self.hop_a[k] = (self.hop_a[k - 1] @ W_a).multiply(mask).tocsr()
            self.hop_b[k] = (self.hop_b[k - 1] @ W_b).multiply(mask).tocsr()
            self.n_hop[k] = np.asarray((mask > 0).sum(axis=1)).ravel()
            if k < self.max_hops:
                visited = visited + mask
                visited.data[:] = 1.0
                frontier = mask
        self.direct_a = self.hop_a[1]
        self.direct_b = self.hop_b[1]
        self.indirect_a = self.hop_a[2]
        self.indirect_b = self.hop_b[2]
        self.n_direct = np.asarray((self.direct_a > 0).sum(axis=1)).ravel()
        self.n_hop[1] = self.n_direct
        self.n_indirect = self.n_hop[2]
        self.source_index = src

    def level(self, level):
        """(H_a, H_b, number of nodes) of a level name (direct, indirect, hop<k>)."""
        k = level_hops(level)
        return self.hop_a[k], self.hop_b[k], self.n_hop[k]

    def vector(self, signal):
        return np.array([signal.get(n, 0.0) for n in self.nodes], dtype=float)


def gene_features(graph, marks, expression_change, pseudo_count=1.0, operators=None, mark_vectors=None, max_hops=2):
    """Per bait node: expression change and promoter and k-hop log2 activity changes per mark, k = 1 .. max_hops.

    `marks` maps mark name -> (signal_a dict, signal_b dict);
    `expression_change` maps bait node -> RNA log2 fold change. Precomputed
    `operators` (which then set the number of hops) and `mark_vectors` (name -> (vector_a, vector_b)) are used
    when given. Columns per level (direct, indirect, hop3, ...): n_<level> (nodes at that distance),
    contact_<level> and <mark>_<level>.
    """
    if operators is None:
        operators = ActivityOperators(graph, list(expression_change), max_hops=max_hops)
    if mark_vectors is None:
        mark_vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    sources = operators.sources
    levels = operators.levels
    columns = {
        'node': sources,
        'genes': [graph.nodes[n].get('genes', '') for n in sources],
        'rna_log2FC': [expression_change[n] for n in sources],
    }
    for level in levels:
        columns['n_' + level] = operators.level(level)[2]
    features = pd.DataFrame(columns)
    ones = np.ones(len(operators.nodes))
    for level in levels:
        H_a, H_b, _ = operators.level(level)
        features['contact_' + level] = np.log2((H_b @ ones + pseudo_count) / (H_a @ ones + pseudo_count))
    src = operators.source_index
    for mark, (vector_a, vector_b) in mark_vectors.items():
        features[mark + '_promoter'] = np.log2((vector_b[src] + pseudo_count) / (vector_a[src] + pseudo_count))
        for level in levels:
            H_a, H_b, _ = operators.level(level)
            features[mark + '_' + level] = np.log2((H_b @ vector_b + pseudo_count) / (H_a @ vector_a + pseudo_count))
    return features


def default_threads(cap=32):
    """Cores available to this process (os.sched_getaffinity), at most `cap`."""
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:
        available = os.cpu_count() or 1
    return max(1, min(cap, available))


_SHARED = {}


def _pool_map(func, items, threads, shared):
    """map(func, items) over `threads` forked worker processes. `shared` is visible to the workers as the module
    global _SHARED (set before the pool starts, so it is inherited read-only and not pickled); results come back in
    the order of `items`."""
    global _SHARED
    items = list(items)
    _SHARED = shared
    try:
        if threads <= 1 or len(items) <= 1:
            return [func(item) for item in items]
        context = multiprocessing.get_context('fork')
        with context.Pool(min(threads, len(items))) as pool:
            return pool.map(func, items, chunksize=1)
    finally:
        _SHARED = {}


def _bigwig_chunk(task):
    from .annotation import bigwig_interval_signal
    path, start, end = task
    s = _SHARED
    return bigwig_interval_signal(path, s['chroms'][start:end], s['starts'][start:end], s['ends'][start:end])


def _interval_signals(paths, chroms, starts, ends, threads, chunk=20000):
    """{bigWig path: mean signal of every interval}, read in parallel over files and chunks of intervals (the values
    are those of annotation.bigwig_interval_signal on all intervals at once)."""
    tasks = [(p, i, min(i + chunk, len(starts))) for p in paths for i in range(0, len(starts), chunk)]
    parts = _pool_map(_bigwig_chunk, tasks, threads, dict(chroms=chroms, starts=starts, ends=ends))
    signals = {p: [] for p in paths}
    for (p, _, _), values in zip(tasks, parts):
        signals[p].append(values)
    return {p: np.concatenate(v) if v else np.zeros(0) for p, v in signals.items()}


def mirrored_features(graph, operators, tracks_a, tracks_b, pseudo_count=1.0, threads=1):
    """Linear-proximity control of the k-hop features (direct, indirect, hop3, ...).

    Every node m at k hops of a source s is replaced by the interval at the same
    distance on the opposite side of s (mirrored around the midpoint of s),
    keeping the contact (path-product) weights. The mark signal of the mirrored intervals is
    read from the bigWig files (in parallel over `threads` processes), so the control has the genomic distance
    profile of the real neighbourhood without its 3D topology.
    """
    src_nodes = operators.sources
    node_start = np.array([graph.nodes[n]['start'] for n in operators.nodes], dtype=np.int64)
    node_end = np.array([graph.nodes[n]['end'] for n in operators.nodes], dtype=np.int64)
    source_chrom = np.array([graph.nodes[n]['chrom'] for n in src_nodes], dtype=object)
    source_centre = np.array([graph.nodes[n]['start'] + graph.nodes[n]['end'] for n in src_nodes], dtype=np.int64)
    columns = {}
    for level in operators.levels:
        W_a, W_b, _ = operators.level(level)
        # paths of either condition: gained neighbours exist in B only, lost ones in A only
        pattern = ((abs(W_a) > 0).astype(np.int8) + (abs(W_b) > 0).astype(np.int8)).tocoo()
        rows, cols = pattern.row, pattern.col
        weight_a = np.asarray(W_a.tocsr()[rows, cols]).ravel()
        weight_b = np.asarray(W_b.tocsr()[rows, cols]).ravel()
        chroms = list(source_chrom[rows])
        starts = source_centre[rows] - node_end[cols]
        ends = source_centre[rows] - node_start[cols]
        paths = list(dict.fromkeys(p for tracks in (tracks_a, tracks_b) for mark in tracks for p in tracks[mark]))
        read = _interval_signals(paths, chroms, starts, ends, threads) if len(rows) else {}
        for mark in tracks_a:
            signal_a = np.mean([read[p] for p in tracks_a[mark]], axis=0) if len(rows) else np.zeros(0)
            signal_b = np.mean([read[p] for p in tracks_b[mark]], axis=0) if len(rows) else np.zeros(0)
            activity_a = np.bincount(rows, weights=weight_a * signal_a, minlength=len(src_nodes))
            activity_b = np.bincount(rows, weights=weight_b * signal_b, minlength=len(src_nodes))
            columns['{}_{}_mirror'.format(mark, level)] = np.log2((activity_b + pseudo_count) /
                                                                   (activity_a + pseudo_count))
    return pd.DataFrame(columns, index=src_nodes)


def _ols_r2(y, X):
    X = np.column_stack([np.ones(len(y)), X]) if X.size else np.ones((len(y), 1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    residual = y - X @ beta
    rss = float(residual @ residual)
    tss = float(((y - y.mean()) ** 2).sum())
    return 1 - rss / tss if tss > 0 else 0.0, rss, X.shape[1]


def nested_models(features, marks, levels=None):
    """R^2 of RNA change explained by promoter, + 1 hop (direct), + 2 hops (indirect), + 3 hops ... features;
    partial F tests of every added level."""
    levels = ['promoter'] + (levels if levels is not None else hop_levels(features, marks))
    columns = []
    data = features.dropna(subset=['rna_log2FC'])
    y = data['rna_log2FC'].to_numpy(dtype=float)
    rows = []
    previous = None
    for level in levels:
        columns = columns + [m + '_' + level for m in marks]
        X = np.nan_to_num(data[columns].to_numpy(dtype=float))
        r2, rss, k = _ols_r2(y, X)
        row = {'model': '+'.join(levels[:levels.index(level) + 1]), 'n_genes': len(y), 'r2': r2,
               'n_parameters': k}
        if previous is not None:
            prev_r2, prev_rss, prev_k = previous
            df1 = k - prev_k
            df2 = len(y) - k
            f_stat = ((prev_rss - rss) / df1) / (rss / df2) if rss > 0 and df2 > 0 else np.nan
            row['delta_r2'] = r2 - prev_r2
            row['f_statistic'] = f_stat
            row['f_pvalue'] = float(stats.f.sf(f_stat, df1, df2)) if not np.isnan(f_stat) else np.nan
        rows.append(row)
        previous = (r2, rss, k)
    return pd.DataFrame(rows)


def model_r2(features, columns):
    data = features.dropna(subset=['rna_log2FC'])
    y = data['rna_log2FC'].to_numpy(dtype=float)
    X = np.nan_to_num(data[columns].to_numpy(dtype=float)) if columns else np.zeros((len(y), 0))
    return _ols_r2(y, X)[0]


def mirror_comparison(features, marks, n_bootstrap=1000, seed=0, levels=None):
    """Variance in RNA change gained by graph neighbours versus their mirrored loci, per level.

    Returns the R^2 of the nested models (promoter, + direct, + indirect, + hop3 ...) and the gain of every level
    over the model of the closer levels, for the graph features and for the mirror features of that level, with
    paired bootstrap intervals (over genes) of graph minus mirror.
    """
    levels = levels if levels is not None else hop_levels(features, marks)
    promoter = [m + '_promoter' for m in marks]

    def gains(df):
        r2 = {'r2_promoter': model_r2(df, promoter)}
        out_gains = {}
        columns, name, previous = list(promoter), 'r2_promoter', r2['r2_promoter']
        for level in levels:
            mirror = model_r2(df, columns + [m + '_' + level + '_mirror' for m in marks])
            columns = columns + [m + '_' + level for m in marks]
            name = name + '_' + level
            r2[name] = model_r2(df, columns)
            out_gains['gain_{}_graph'.format(level)] = r2[name] - previous
            out_gains['gain_{}_mirror'.format(level)] = mirror - previous
            previous = r2[name]
        r2.update(out_gains)
        return r2

    observed = gains(features)
    rng = np.random.default_rng(seed)
    boot = {level: [] for level in levels}
    for _ in range(n_bootstrap):
        sample = features.iloc[rng.integers(0, len(features), len(features))]
        g = gains(sample)
        for level in levels:
            boot[level].append(g['gain_{}_graph'.format(level)] - g['gain_{}_mirror'.format(level)])
    observed['n_genes'] = len(features)
    for level in levels:
        values = np.array(boot[level])
        observed.update({
            '{}_graph_minus_mirror'.format(level): observed['gain_{}_graph'.format(level)] -
            observed['gain_{}_mirror'.format(level)],
            '{}_ci_low'.format(level): float(np.quantile(values, 0.025)),
            '{}_ci_high'.format(level): float(np.quantile(values, 0.975)),
            '{}_fraction_bootstrap_above_0'.format(level): float(np.mean(values > 0)),
        })
    observed['n_bootstrap'] = n_bootstrap
    return observed


def feature_correlations(features, marks):
    rows = []
    graph_levels = hop_levels(features, marks)
    levels = ['promoter'] + graph_levels
    if any(m + '_direct_mirror' in features for m in marks):
        levels += [level + '_mirror' for level in graph_levels]
    for level in levels:
        for mark in marks:
            column = mark + '_' + level
            if column not in features:
                continue
            data = features[['rna_log2FC', column]].dropna()
            if len(data) < 3:
                continue
            rho, p = stats.spearmanr(data['rna_log2FC'], data[column])
            rows.append({'mark': mark, 'level': level, 'n': len(data), 'spearman_rho': rho, 'pvalue': p})
    for level in graph_levels:
        column = 'contact_' + level
        if column in features:
            data = features[['rna_log2FC', column]].dropna()
            if len(data) >= 3:
                rho, p = stats.spearmanr(data['rna_log2FC'], data[column])
                rows.append({'mark': 'contact', 'level': level, 'n': len(data), 'spearman_rho': rho, 'pvalue': p})
    return pd.DataFrame(rows)


def read_proximity_support(path):
    """Step curve of analysis/hop_distance.py (columns steps, concordance, ci_low, ci_high) as
    {steps: (concordance, ci_low, ci_high)}: the probability that graph neighbours at that many steps are closer in
    the nucleus than distance-matched pairs (0.5 without relation)."""
    table = pd.read_csv(path, sep='\t')
    missing = {'steps', 'concordance', 'ci_low', 'ci_high'} - set(table.columns)
    if missing:
        raise ValueError('{} lacks the columns {}'.format(path, sorted(missing)))
    return {int(r.steps): (float(r.concordance), float(r.ci_low), float(r.ci_high)) for r in table.itertuples()}


def proximity_columns(hops, support):
    """proximity_support and proximity_support_ci of explaining levels (0 hops = promoter: 1.0 and 'NA'; levels
    outside the step curve and unexplained genes: NaN and '')."""
    values, intervals = [], []
    for h in hops:
        if pd.isna(h):
            values.append(np.nan)
            intervals.append('')
        elif int(h) == 0:
            values.append(1.0)
            intervals.append('NA')
        elif int(h) in support:
            value, low, high = support[int(h)]
            values.append(value)
            intervals.append('{:.2f}-{:.2f}'.format(low, high))
        else:
            values.append(np.nan)
            intervals.append('')
    return values, intervals


def _level_size(row, level):
    """Number of nodes tested at a level (the bait itself for the promoter)."""
    if level == 'promoter':
        return 1
    value = row.get('n_' + level, np.nan)
    return np.nan if pd.isna(value) else int(value)


def classify_genes(features, marks, mark_signs, rna_threshold=1.0, chromatin_threshold=1.0, proximity=None):
    """Hierarchical explanation of expression changes: promoter > 1 hop (direct) > 2 hops (indirect) > 3 hops (hop3)
    > ... > unexplained.

    A level explains a gene if at least one mark changes by at least
    `chromatin_threshold` (log2) on that level in the direction expected from
    the mark sign (+1 activating, -1 repressive) and the RNA change; the gene takes the first such level.
    explaining_hops is its distance (0 at the promoter) and n_hop_nodes_tested the number of nodes at that level:
    more hops reach more nodes and so more chances of a concordant change by chance, which the per-level
    permutation test (gene_level_hop_test) accounts for. With `proximity` (read_proximity_support) every explained
    gene gets the 3D support of its level (proximity_support, proximity_support_ci).
    """
    levels = ['promoter'] + hop_levels(features, marks)
    de = features[features['rna_log2FC'].abs() >= rna_threshold].copy()
    direction = np.sign(de['rna_log2FC'].to_numpy())
    explanation = np.array(['unexplained'] * len(de), dtype=object)
    explaining_marks = [''] * len(de)
    for level in reversed(levels):
        hit = np.zeros(len(de), dtype=bool)
        level_marks = [[] for _ in range(len(de))]
        for mark in marks:
            column = mark + '_' + level
            values = np.nan_to_num(de[column].to_numpy(dtype=float))
            concordant = (np.sign(values) * mark_signs.get(mark, 1) == direction) & (np.abs(values) >= chromatin_threshold)
            hit |= concordant
            for i in np.nonzero(concordant)[0]:
                level_marks[i].append(mark)
        explanation[hit] = level
        for i in np.nonzero(hit)[0]:
            explaining_marks[i] = ','.join(level_marks[i])
    de['explanation'] = explanation
    de['explaining_marks'] = explaining_marks
    de['explaining_hops'] = pd.array([None if e == 'unexplained' else level_hops(e) for e in explanation],
                                     dtype='Int64')
    sizes = np.full(len(de), np.nan)
    sizes[explanation == 'promoter'] = 1
    for level in levels[1:]:
        if 'n_' + level in de:
            hit = explanation == level
            sizes[hit] = de['n_' + level].to_numpy(dtype=float)[hit]
    de['n_hop_nodes_tested'] = pd.array([None if np.isnan(v) else int(v) for v in sizes], dtype='Int64')
    if proximity is not None:
        de['proximity_support'], de['proximity_support_ci'] = proximity_columns(de['explaining_hops'], proximity)
    return de


def cascade_label(hops):
    """Cascade level of a number of steps: 0 'gene', k 'step<k>'."""
    return 'gene' if hops == 0 else 'step{}'.format(hops)


def explanation_cascade(classes, marks, mark_signs, chromatin_threshold=1.0, contact_threshold=None):
    """Where the expression change of every differentially expressed gene finds its first explanation.

    `classes` is the output of classify_genes (differential genes with their features and mirror features). Level
    'gene': a concordant change of a mark at the promoter fragment (the promoter rule of classify_genes) and/or of
    the promoter's own interactions, i.e. |contact_direct| (log2 B/A of the summed strength of its direct
    interactions) >= `contact_threshold` (default `chromatin_threshold`) with the sign of the expression change.
    Genes not explained at the gene go on to step 1, 2, ..., n (n = the levels in `classes`): the first step at which
    a mark changes concordantly (the rule of classify_genes); the rest stay unexplained. An explanation is the
    first concordant change given the measured layers, not a proof of regulation.

    Returns a DataFrame (indexed by node) with explained_at ('gene', 'step1', ..., 'unexplained'), explained_hops
    (0 at the gene), explained_by (at the gene: 'marks', 'interactions' or 'both'; at a step: the concordant marks),
    beyond_mirror (at a step: not every explaining mark changes concordantly at the mirrored loci of that step; NA
    without mirror features) and nodes_tested (1 at the gene, else the number of nodes at that step).
    """
    if contact_threshold is None:
        contact_threshold = chromatin_threshold
    levels = hop_levels(classes, marks)
    direction = np.sign(classes['rna_log2FC'].to_numpy(dtype=float))

    def concordant(column, sign=1):
        if column not in classes:
            return np.zeros(len(classes), dtype=bool)
        values = np.nan_to_num(classes[column].to_numpy(dtype=float))
        return (np.sign(values) * sign == direction) & (np.abs(values) >= chromatin_threshold)

    by_marks = np.zeros(len(classes), dtype=bool)
    for mark in marks:
        by_marks |= concordant(mark + '_promoter', mark_signs.get(mark, 1))
    contact = np.nan_to_num(classes['contact_direct'].to_numpy(dtype=float)) if 'contact_direct' in classes \
        else np.zeros(len(classes))
    by_contact = (np.sign(contact) == direction) & (np.abs(contact) >= contact_threshold)
    at = np.array(['unexplained'] * len(classes), dtype=object)
    hops = np.full(len(classes), np.nan)
    explained_by = np.array([''] * len(classes), dtype=object)
    beyond = np.full(len(classes), np.nan)
    nodes = np.full(len(classes), np.nan)
    gene = by_marks | by_contact
    at[gene] = 'gene'
    hops[gene] = 0
    nodes[gene] = 1
    explained_by[by_marks & by_contact] = 'both'
    explained_by[by_marks & ~by_contact] = 'marks'
    explained_by[~by_marks & by_contact] = 'interactions'
    for level in levels:
        k = level_hops(level)
        open_ = at == 'unexplained'
        hit_marks = {m: concordant('{}_{}'.format(m, level), mark_signs.get(m, 1)) & open_ for m in marks}
        hit = np.zeros(len(classes), dtype=bool)
        for values in hit_marks.values():
            hit |= values
        at[hit] = cascade_label(k)
        hops[hit] = k
        if 'n_' + level in classes:
            nodes[hit] = classes['n_' + level].to_numpy(dtype=float)[hit]
        have_mirror = all('{}_{}_mirror'.format(m, level) in classes for m in marks)
        mirrored = {m: hit_marks[m] & concordant('{}_{}_mirror'.format(m, level), mark_signs.get(m, 1))
                    for m in marks} if have_mirror else {}
        for i in np.nonzero(hit)[0]:
            explaining = [m for m in marks if hit_marks[m][i]]
            explained_by[i] = ','.join(explaining)
            if have_mirror:
                beyond[i] = float(sum(mirrored[m][i] for m in explaining) < len(explaining))
    return pd.DataFrame({'explained_at': at, 'explained_hops': pd.array([None if np.isnan(h) else int(h) for h in hops],
                                                                        dtype='Int64'),
                         'explained_by': explained_by,
                         'beyond_mirror': pd.array([None if np.isnan(b) else bool(b) for b in beyond], dtype='boolean'),
                         'nodes_tested': pd.array([None if np.isnan(v) else int(v) for v in nodes], dtype='Int64')},
                        index=classes['node'].to_numpy() if 'node' in classes else classes.index)


def cascade_table(cascade, genes=None, pvalues=None, proximity=None, max_hops=None):
    """Counts of the explanation cascade (explanation_cascade), one row per level.

    Rows: gene (with gene_marks_only, gene_interactions_only, gene_both), step1 .. step<max_hops>, unexplained.
    genes_reaching: genes not explained at an earlier level (all differential genes for the gene level);
    fraction_of_remaining = genes / genes_reaching. Per step: beyond_mirror, and with `pvalues` (node -> per-gene
    permutation P of the step at which the gene is first explained) p_lt_0.05, beyond_mirror_p_lt_0.05 and the
    numbers at FDR < 0.1 and < 0.05 (Benjamini-Hochberg across the genes first explained at that step, names of
    `genes` (node -> gene names) at FDR < 0.05 listed); proximity_support and proximity_ci of the step from a
    read_proximity_support curve; median_nodes_tested.
    """
    total = len(cascade)
    if max_hops is None:
        max_hops = int(cascade['explained_hops'].max()) if cascade['explained_hops'].notna().any() else 0
    rows = []
    reaching = total

    def row(level, mask, reaching_count):
        n = int(mask.sum())
        return {'level': level, 'genes': n, 'genes_reaching': reaching_count,
                'fraction_of_differential': n / total if total else np.nan,
                'fraction_of_remaining': n / reaching_count if reaching_count else np.nan}

    at = cascade['explained_at'].to_numpy()
    gene = at == 'gene'
    r = row('gene', gene, reaching)
    if proximity is not None:
        r.update(proximity_support=1.0, proximity_ci='NA')
    r['median_nodes_tested'] = 1.0
    rows.append(r)
    for name, value in (('gene_marks_only', 'marks'), ('gene_interactions_only', 'interactions'), ('gene_both', 'both')):
        r = row(name, gene & (cascade['explained_by'].to_numpy() == value), reaching)
        rows.append(r)
    reaching -= int(gene.sum())
    for k in range(1, max_hops + 1):
        level = cascade_label(k)
        mask = at == level
        part = cascade[mask]
        r = row(level, mask, reaching)
        beyond = part['beyond_mirror'].fillna(False).to_numpy(dtype=bool)
        r['beyond_mirror'] = int(beyond.sum()) if part['beyond_mirror'].notna().any() or not len(part) else np.nan
        if pvalues is not None:
            p = np.array([pvalues.get(n, np.nan) for n in part.index], dtype=float)
            q = _benjamini_hochberg(p) if len(p) else p
            r['p_lt_0.05'] = int((p < 0.05).sum())
            r['beyond_mirror_p_lt_0.05'] = int(((p < 0.05) & beyond).sum())
            r['q_lt_0.1'] = int((q < 0.1).sum())
            r['q_lt_0.05'] = int((q < 0.05).sum())
            names = [str(genes.get(n, n)) if genes is not None else n for n in part.index[q < 0.05]]
            r['genes_q_lt_0.05'] = ';'.join(names)
        if proximity is not None:
            values, intervals = proximity_columns([k], proximity)
            r['proximity_support'], r['proximity_ci'] = values[0], intervals[0]
        r['median_nodes_tested'] = float(part['nodes_tested'].astype(float).median()) if len(part) else np.nan
        rows.append(r)
        reaching -= int(mask.sum())
    rows.append(row('unexplained', at == 'unexplained', reaching))
    columns = ['level', 'genes', 'genes_reaching', 'fraction_of_differential', 'fraction_of_remaining',
               'beyond_mirror', 'p_lt_0.05', 'beyond_mirror_p_lt_0.05', 'q_lt_0.1', 'q_lt_0.05', 'proximity_support',
               'proximity_ci', 'median_nodes_tested', 'genes_q_lt_0.05']
    table = pd.DataFrame(rows)
    return table[[c for c in columns if c in table]]


def _concordant(values, mark_sign, direction, threshold):
    values = np.nan_to_num(np.asarray(values, dtype=float))
    return (np.sign(values) * mark_sign == direction) & (np.abs(values) >= threshold)


def _path_mediators(graph, source, target, hops, max_paths=3):
    """Intermediate nodes of up to `max_paths` shortest paths, labelled by their genes, as 'a>b' chains."""
    chains = []
    for path in nx.all_shortest_paths(graph, source, target):
        if len(path) != hops + 1:
            break
        chains.append('>'.join(graph.nodes[m].get('genes', '') or m for m in path[1:-1]))
        if len(chains) >= max_paths:
            break
    return '+'.join(sorted(set(chains)))


def indirect_only_genes(graph, operators, mark_vectors, features, marks, mark_signs, rna_threshold=1.0,
                        chromatin_threshold=1.0, genes=None, top=3, proximity=None):
    """Differential genes first explained at k >= 2 hops (indirect neighbours or further).

    A gene qualifies if its expression changes at least `rna_threshold`
    (log2), no mark changes concordantly at its promoter or at any level closer than k, and at least one mark
    changes concordantly at k hops; explaining_hops holds k and n_hop_nodes_tested the number of nodes at k hops.
    `beyond_mirror` marks genes whose explaining marks show no concordant change at the mirrored loci of that level.
    For each explaining mark the k-hop nodes with the largest concordant contribution are listed with their
    annotated genes and the mediators that connect them to the gene (for k = 2 the shared neighbours, for k >= 3
    up to three shortest paths, intermediate nodes joined by '>').
    """
    levels = ['promoter'] + hop_levels(features, marks)
    extra_levels = [level for level in levels if level_hops(level) >= 3]
    data = features if genes is None else features[features['node'].isin(genes)]
    data = data[data['rna_log2FC'].abs() >= rna_threshold]
    row_of = {node: i for i, node in enumerate(operators.sources)}
    rows = []
    for _, gene in data.iterrows():
        direction = np.sign(gene['rna_log2FC'])
        first, explaining = None, []
        for level in levels:
            explaining = [m for m in marks
                          if _concordant(gene[m + '_' + level], mark_signs.get(m, 1), direction, chromatin_threshold)]
            if explaining:
                first = level
                break
        if first is None or level_hops(first) < 2:
            continue
        hops = level_hops(first)
        mirrored = [m for m in explaining if m + '_' + first + '_mirror' in gene and
                    _concordant(gene[m + '_' + first + '_mirror'], mark_signs.get(m, 1), direction,
                                chromatin_threshold)]
        i = row_of[gene['node']]
        H_a, H_b, _ = operators.level(first)
        drivers = []
        for mark in explaining:
            vector_a, vector_b = mark_vectors[mark]
            row_a = H_a.getrow(i)
            row_b = H_b.getrow(i)
            columns = np.union1d(row_a.indices, row_b.indices)
            weight_a = np.asarray(row_a[:, columns].todense()).ravel()
            weight_b = np.asarray(row_b[:, columns].todense()).ravel()
            contribution = (weight_b * vector_b[columns] - weight_a * vector_a[columns]) * mark_signs.get(mark, 1) * direction
            for k in np.argsort(contribution)[::-1][:top]:
                if contribution[k] <= 0:
                    break
                node = operators.nodes[columns[k]]
                if hops == 2:
                    mediators = set(graph.neighbors(gene['node'])) & set(graph.neighbors(node))
                    via = '+'.join(sorted(graph.nodes[m].get('genes', '') or m for m in mediators))
                else:
                    via = _path_mediators(graph, gene['node'], node, hops)
                drivers.append('{}:{}[{}] via {}'.format(
                    mark, node, graph.nodes[node].get('genes', '') or ('promoter' if graph.nodes[node].get('bait') else 'OE'),
                    via))
        row = {'node': gene['node'], 'genes': gene['genes'], 'rna_log2FC': gene['rna_log2FC'],
               'n_direct': gene['n_direct'], 'n_indirect': gene['n_indirect'],
               'explaining_marks': ','.join(explaining), 'mirror_concordant_marks': ','.join(mirrored),
               'beyond_mirror': len(mirrored) < len(explaining), 'drivers': '; '.join(drivers)}
        for mark in marks:
            for level in ['promoter', 'direct', 'indirect', 'indirect_mirror'] + \
                    [name for extra in extra_levels for name in (extra, extra + '_mirror')]:
                column = '{}_{}'.format(mark, level)
                if column in gene:
                    row[column] = gene[column]
        row['explaining_hops'] = hops
        row['n_hop_nodes_tested'] = _level_size(gene, first)
        for level in extra_levels:
            row['n_' + level] = gene['n_' + level]
        if proximity is not None:
            support, interval = proximity_columns([hops], proximity)
            row['proximity_support'], row['proximity_support_ci'] = support[0], interval[0]
        rows.append(row)
    result = pd.DataFrame(rows)
    if len(result):
        result = result.sort_values(['beyond_mirror', 'rna_log2FC'], ascending=[False, False],
                                    key=lambda c: c.abs() if c.name == 'rna_log2FC' else c)
    return result


def permutation_test(graph, marks, expression_change, mark_signs, n_permutations=100, seed=0,
                     rna_threshold=1.0, chromatin_threshold=1.0, max_hops=2):
    """Null distribution of the contribution of every level at k >= 2 hops.

    Mark signals of all non-source nodes are permuted jointly (A and B values of
    a node move together, so the distribution of changes is preserved) while the
    graph and the signals at the source promoters stay fixed. The indirect (2-hop) level keeps its keys
    (observed_indirect_delta_r2, observed_fraction_indirect, ...); levels at k >= 3 hops add the same statistics as
    observed_hop<k>_delta_r2, observed_fraction_hop<k>, ...
    """
    rng = np.random.default_rng(seed)
    operators = ActivityOperators(graph, list(expression_change), max_hops=max_hops)
    levels = operators.levels
    tested = levels[1:]   # indirect, hop3, ...
    position = {level: levels.index(level) + 1 for level in tested}   # row of the level in nested_models
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    observed = gene_features(graph, marks, expression_change, operators=operators, mark_vectors=vectors)
    observed_models = nested_models(observed, list(marks))
    observed_class = classify_genes(observed, list(marks), mark_signs, rna_threshold, chromatin_threshold)
    obs_delta = {level: float(observed_models['delta_r2'].iloc[position[level]]) for level in tested}
    obs_fraction = {level: float((observed_class['explanation'] == level).mean()) if len(observed_class) else np.nan
                    for level in tested}
    is_source = np.zeros(len(operators.nodes), dtype=bool)
    is_source[operators.source_index] = True
    others = np.nonzero(~is_source)[0]
    null_delta = {level: [] for level in tested}
    null_fraction = {level: [] for level in tested}
    for _ in range(n_permutations):
        donors = others[rng.permutation(len(others))]
        permuted = {}
        for mark, (vector_a, vector_b) in vectors.items():
            new_a = vector_a.copy()
            new_b = vector_b.copy()
            new_a[others] = vector_a[donors]
            new_b[others] = vector_b[donors]
            permuted[mark] = (new_a, new_b)
        features = gene_features(graph, marks, expression_change, operators=operators, mark_vectors=permuted)
        models = nested_models(features, list(marks))
        classes = classify_genes(features, list(marks), mark_signs, rna_threshold, chromatin_threshold)
        for level in tested:
            null_delta[level].append(float(models['delta_r2'].iloc[position[level]]))
            null_fraction[level].append(float((classes['explanation'] == level).mean()) if len(classes) else np.nan)
    result = {'n_genes': int(len(observed)), 'n_differential_genes': int(len(observed_class))}
    for level in tested:
        delta = np.array(null_delta[level])
        fraction = np.array(null_fraction[level])
        summary = {
            'observed_indirect_delta_r2': obs_delta[level],
            'null_indirect_delta_r2_mean': float(np.mean(delta)),
            'null_indirect_delta_r2_sd': float(np.std(delta)),
            'empirical_p_delta_r2': float((np.sum(delta >= obs_delta[level]) + 1) / (len(delta) + 1)),
            'observed_fraction_indirect': obs_fraction[level],
            'null_fraction_indirect_mean': float(np.nanmean(fraction)),
            'null_fraction_indirect_sd': float(np.nanstd(fraction)),
            'empirical_p_fraction': float((np.sum(fraction >= obs_fraction[level]) + 1) / (len(fraction) + 1)),
        }
        if level == 'indirect':
            result.update(summary)
            result['n_permutations'] = int(n_permutations)
        else:
            for key, value in summary.items():
                key = key.replace('indirect', level)
                if key.startswith('empirical_p_'):
                    key = key + '_' + level
                result[key] = value
    return result


PERMUTATION_BLOCK = 1000
BLOCK_MEMORY = 256e6   # bytes of permuted mark values (and their node indices) per worker


def _hop_statistics(state, index):
    """Statistic of every tested gene and level for the permuted signals given by `index` (used nodes x
    permutations: the node whose signal each used node carries). One sparse-dense product per mark and level."""
    best = {level: np.full((len(test['rows']), index.shape[1]), -np.inf) for level, test in state['tests'].items()}
    for mark in state['marks']:
        vector_a, vector_b = state['vectors'][mark]
        values_a = vector_a[index]
        values_b = vector_b[index]
        sign = state['signs'].get(mark, 1)
        for level, test in state['tests'].items():
            change = np.log2((test['H_b'] @ values_b + state['pseudo_count']) /
                             (test['H_a'] @ values_a + state['pseudo_count']))
            np.maximum(best[level], change * (sign * test['direction'])[:, None], out=best[level])
    return best


def _hop_block(block):
    """Exceedance counts of one block of permutations, drawn from its own generator (seed, block)."""
    state = _SHARED
    start = block * PERMUTATION_BLOCK
    size = min(PERMUTATION_BLOCK, state['n_permutations'] - start)
    rng = np.random.default_rng(np.random.SeedSequence([state['seed'], block]))
    others, used, position = state['others'], state['used'], state['position']
    in_others = position >= 0
    exceed = {level: np.zeros(len(test['rows'])) for level, test in state['tests'].items()}
    chunk = state['chunk']
    for offset in range(0, size, chunk):
        width = min(chunk, size - offset)
        index = np.repeat(used[:, None], width, axis=1)
        for j in range(width):
            donors = others[rng.permutation(len(others))]
            index[in_others, j] = donors[position[in_others]]
        best = _hop_statistics(state, index)
        for level in exceed:
            exceed[level] += (best[level] >= state['observed'][level][:, None]).sum(axis=1)
    return exceed


def gene_level_hop_test(operators, mark_vectors, features, marks, mark_signs, rna_threshold=1.0,
                        n_permutations=1000, seed=0, pseudo_count=1.0, genes=None, threads=1):
    """Per-gene permutation test of the levels at k >= 2 hops for differential genes.

    `genes` maps a level name (indirect, hop3, ...) to the nodes tested at that level (None: every differential
    gene); without it the indirect level is tested for every differential gene. The statistic of a gene at level
    k is the largest concordant k-hop change over the marks, max_m sign_m * sign(RNA) * <mark>_<level>
    (log2 B/A). Its null distribution comes from permuting the mark signals of all non-source nodes jointly (as in
    permutation_test), which keeps the graph and the distribution of changes but breaks which nodes carry them; all
    levels share the same permutations. The test conditions on the nodes of the gene at that level: more hops reach
    more nodes and so give more chances of a large concordant change by chance, and the permuted statistic is
    computed over the same nodes and path weights, so the P-value accounts for the number of nodes tested. BH
    q-values are computed across the genes tested at each level. Returns {level: DataFrame with node,
    <level>_statistic, <level>_p, <level>_q}.

    The permutations come in blocks of PERMUTATION_BLOCK; block b draws from
    np.random.default_rng(np.random.SeedSequence([seed, b])), so the result does not depend on `threads` (worker
    processes over blocks). Within a block the permuted signals form a matrix (nodes x permutations) and every mark
    and level takes one sparse-dense product.
    """
    if genes is None:
        genes = {'indirect': None}
    de = features[features['rna_log2FC'].abs() >= rna_threshold]
    row_of = {node: i for i, node in enumerate(operators.sources)}
    tests = {}
    for level, nodes in genes.items():
        subset = de if nodes is None else de[de['node'].isin(set(nodes))]
        rows = np.array([row_of[n] for n in subset['node']], dtype=np.int64)
        H_a, H_b, n_level = operators.level(level)
        tests[level] = dict(nodes=subset['node'].to_numpy(), rows=rows, H_a=H_a[rows].tocsr(),
                            H_b=H_b[rows].tocsr(), direction=np.sign(subset['rna_log2FC'].to_numpy()),
                            n=n_level[rows])
    # nodes with a weight for any tested gene; the operators are restricted to them
    used = np.unique(np.concatenate([np.concatenate([t['H_a'].indices, t['H_b'].indices]) for t in tests.values()]
                                    + [np.zeros(0, dtype=np.int64)])).astype(np.int64)
    for test in tests.values():
        test['H_a'] = test['H_a'][:, used].tocsr()
        test['H_b'] = test['H_b'][:, used].tocsr()
    is_source = np.zeros(len(operators.nodes), dtype=bool)
    is_source[operators.source_index] = True
    others = np.nonzero(~is_source)[0]
    position = np.full(len(operators.nodes), -1, dtype=np.int64)
    position[others] = np.arange(len(others))
    state = dict(tests=tests, marks=list(marks), signs=mark_signs, vectors=mark_vectors, pseudo_count=pseudo_count,
                 used=used, position=position[used], others=others, seed=seed, n_permutations=n_permutations)
    observed = _hop_statistics(state, used[:, None])
    state['observed'] = {level: values[:, 0] for level, values in observed.items()}
    # index plus two value matrices of used nodes x permutations per chunk
    state['chunk'] = int(max(1, min(PERMUTATION_BLOCK, BLOCK_MEMORY // (24 * max(1, len(used))))))
    blocks = range((n_permutations + PERMUTATION_BLOCK - 1) // PERMUTATION_BLOCK)
    exceed = {level: np.zeros(len(test['rows'])) for level, test in tests.items()}
    if any(len(test['rows']) for test in tests.values()):
        for counts in _pool_map(_hop_block, blocks, threads, state):
            for level in exceed:
                exceed[level] += counts[level]
    result = {}
    for level, test in tests.items():
        pvalues = (exceed[level] + 1) / (n_permutations + 1)
        # genes without nodes at this distance carry no information on the level
        pvalues[test['n'] == 0] = 1.0
        result[level] = pd.DataFrame({'node': test['nodes'], level + '_statistic': state['observed'][level],
                                      level + '_p': pvalues, level + '_q': _benjamini_hochberg(pvalues)})
    return result


def gene_level_indirect_test(operators, mark_vectors, features, marks, mark_signs, rna_threshold=1.0,
                             n_permutations=1000, seed=0, pseudo_count=1.0, threads=1):
    """Per-gene test of the indirect (2-hop) level for every differential gene (gene_level_hop_test).

    The statistic of a gene is the largest concordant indirect change over the marks,
    max_m sign_m * sign(RNA) * indirect_m (log2 B/A). Its null distribution comes from permuting the mark signals
    of all non-source nodes jointly (as in permutation_test), which keeps the graph and the distribution of
    changes but breaks which nodes carry them. Returns node, statistic, empirical p-value and BH q-value.
    """
    return gene_level_hop_test(operators, mark_vectors, features, marks, mark_signs, rna_threshold, n_permutations,
                               seed, pseudo_count, genes={'indirect': None}, threads=threads)['indirect']


def _benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(len(p))
    q[order] = np.minimum(ranked, 1.0)
    return q


def community_transitions(communities_a, communities_b, nodes_a, nodes_b, stable_jaccard=0.5, part_fraction=0.2,
                          present_fraction=0.5):
    """How every community of condition A continues in condition B, and which B communities are new.

    Communities are partitions of the two condition graphs computed independently with the same method. For a
    community of A, the fragments still present in graph B are distributed over the B communities. Classes, in this
    order: 'dissolved' (less than `present_fraction` of its fragments are in graph B), 'split' (at least two B
    communities receive >= `part_fraction` of its present fragments), 'merged' (its main B community draws
    >= `part_fraction` of its fragments from at least two A communities), 'stable' (reciprocal best match with node
    Jaccard >= `stable_jaccard`) and 'remodelled' otherwise. A B community with less than `present_fraction` of its fragments in
    graph A is 'formed'. Returns (table, dict A index -> main B index).
    """
    of_b = {n: i for i, members in enumerate(communities_b) for n in members}
    of_a = {n: i for i, members in enumerate(communities_a) for n in members}
    best_a_for_b = {}
    for index_b, members in enumerate(communities_b):
        counts = defaultdict(int)
        for n in members:
            if n in of_a:
                counts[of_a[n]] += 1
        if counts:
            # ties broken by the lower community index, so that the result does not depend on set order
            best_a_for_b[index_b] = max(sorted(counts), key=lambda k: counts[k] / len(members | communities_a[k]))
    rows, main = [], {}
    for index_a, members in enumerate(communities_a):
        present = [n for n in members if n in nodes_b]
        counts = defaultdict(int)
        for n in present:
            counts[of_b[n]] += 1
        row = {'side': 'A', 'community_a': index_a, 'size_a': len(members),
               'fraction_present_in_b': len(present) / len(members), 'community_b': -1, 'size_b': 0,
               'jaccard': 0.0, 'n_b_parts': 0}
        if counts:
            best = max(sorted(counts), key=lambda k: counts[k] / len(members | communities_b[k]))
            main[index_a] = best
            parts = [k for k, c in counts.items() if c >= part_fraction * len(present)]
            sources = defaultdict(int)
            for n in communities_b[best]:
                if n in of_a:
                    sources[of_a[n]] += 1
            merged_sources = [k for k, c in sources.items() if c >= part_fraction * len(communities_b[best])]
            jaccard = counts[best] / len(members | communities_b[best])
            row.update(community_b=best, size_b=len(communities_b[best]), jaccard=jaccard, n_b_parts=len(parts))
            if row['fraction_present_in_b'] < present_fraction:
                row['transition'] = 'dissolved'
            elif len(parts) >= 2:
                row['transition'] = 'split'
            elif len(merged_sources) >= 2:
                row['transition'] = 'merged'
            elif best_a_for_b.get(best) == index_a and jaccard >= stable_jaccard:
                row['transition'] = 'stable'
            else:
                row['transition'] = 'remodelled'
        else:
            row['transition'] = 'dissolved'
        rows.append(row)
    for index_b, members in enumerate(communities_b):
        present = sum(1 for n in members if n in nodes_a)
        if present < present_fraction * len(members):
            rows.append({'side': 'B', 'community_a': -1, 'size_a': 0, 'fraction_present_in_b': np.nan,
                         'community_b': index_b, 'size_b': len(members), 'jaccard': 0.0, 'n_b_parts': 0,
                         'transition': 'formed', 'fraction_present_in_a': present / len(members)})
    return pd.DataFrame(rows), main


def gene_community_changes(communities_a, communities_b, transitions, node_genes, bait_nodes, expression=None,
                           differential=None, switch_jaccard=0.3):
    """Community membership of every bait (promoter) fragment in both conditions.

    A bait 'switches' when the Jaccard index of its community co-members in A and in B is below `switch_jaccard`;
    baits present in only one condition graph are 'lost' or 'gained'. Joined and left co-members are listed by their
    annotated genes. With `expression` (node -> RNA log2FC) and `differential` (set of DE nodes), the enrichment of
    DE genes among switching baits is tested (Fisher exact test, baits present in both graphs).
    """
    of_a = {n: i for i, members in enumerate(communities_a) for n in members}
    of_b = {n: i for i, members in enumerate(communities_b) for n in members}
    transition_of = dict(zip(transitions.loc[transitions['side'] == 'A', 'community_a'],
                             transitions.loc[transitions['side'] == 'A', 'transition']))

    def genes_of(nodes):
        names = sorted({g for n in nodes for g in str(node_genes.get(n, '')).split(',') if g})
        return ','.join(names)

    rows = []
    for bait in sorted(bait_nodes):
        in_a, in_b = bait in of_a, bait in of_b
        if not (in_a or in_b):
            continue
        members_a = communities_a[of_a[bait]] - {bait} if in_a else set()
        members_b = communities_b[of_b[bait]] - {bait} if in_b else set()
        union = members_a | members_b
        jaccard = len(members_a & members_b) / len(union) if union else 1.0
        if in_a and in_b:
            status = 'switch' if jaccard < switch_jaccard else 'kept'
        else:
            status = 'lost' if in_a else 'gained'
        row = {'node': bait, 'genes': node_genes.get(bait, ''), 'status': status,
               'community_a': of_a.get(bait, -1), 'community_b': of_b.get(bait, -1),
               'size_a': len(members_a) + in_a, 'size_b': len(members_b) + in_b, 'member_jaccard': jaccard,
               'community_a_transition': transition_of.get(of_a.get(bait, -1), ''),
               'joined_genes': genes_of(members_b - members_a), 'left_genes': genes_of(members_a - members_b)}
        if expression is not None:
            row['rna_log2FC'] = expression.get(bait, np.nan)
        if differential is not None:
            row['differential'] = bait in differential
        rows.append(row)
    table = pd.DataFrame(rows)
    test = {}
    if differential is not None and len(table):
        both = table[table['status'].isin(['switch', 'kept'])]
        if expression is not None:
            both = both[both['rna_log2FC'].notna()]
        switch = both['status'] == 'switch'
        de = both['differential'].astype(bool)
        counts = [[int((switch & de).sum()), int((switch & ~de).sum())],
                  [int((~switch & de).sum()), int((~switch & ~de).sum())]]
        odds, p = stats.fisher_exact(counts)
        test = {'baits_in_both': int(len(both)), 'switching': int(switch.sum()),
                'switching_differential': counts[0][0], 'kept_differential': counts[1][0],
                'fraction_de_switching': counts[0][0] / max(1, int(switch.sum())),
                'fraction_de_kept': counts[1][0] / max(1, int((~switch).sum())),
                'odds_ratio': float(odds), 'fisher_p': float(p)}
        if expression is not None and switch.any() and (~switch).any():
            test['mannwhitney_p_abs_log2FC'] = float(stats.mannwhitneyu(
                both.loc[switch, 'rna_log2FC'].abs(), both.loc[~switch, 'rna_log2FC'].abs(),
                alternative='greater').pvalue)
    return table, test


def _sample_without_replacement(rng, n, k, rows):
    """`rows` x `k` indices into range(n), each row without repeats: rows with a repeat are drawn again (vectorised
    rejection); for k > n / 4 the k smallest of n random keys per row."""
    if k > n / 4:
        return np.argpartition(rng.random((rows, n)), k - 1, axis=1)[:, :k]
    index = rng.integers(0, n, (rows, k))
    while True:
        ordered = np.sort(index, axis=1)
        repeat = (ordered[:, 1:] == ordered[:, :-1]).any(axis=1)
        if not repeat.any():
            return index
        index[repeat] = rng.integers(0, n, (int(repeat.sum()), k))


def _community_null(k):
    """Mean expression change of n_permutations random gene sets of size k, from generator (seed, k)."""
    state = _SHARED
    values, n_permutations = state['values'], state['n_permutations']
    rng = np.random.default_rng(np.random.SeedSequence([state['seed'], k]))
    rows = max(1, int(32e6 // (8 * (len(values) if k > len(values) / 4 else k))))
    means = []
    for start in range(0, n_permutations, rows):
        index = _sample_without_replacement(rng, len(values), k, min(rows, n_permutations - start))
        means.append(values[index].mean(axis=1))
    return np.concatenate(means) if means else np.zeros(0)


def community_expression_changes(communities_a, communities_b, main, transitions, edges, expression, min_genes=3,
                                 n_permutations=1000, seed=0, threads=1):
    """Per community of A and its main B community: contact rewiring and coordinated expression change.

    Rewiring counts the union-graph edges inside the pair (both ends in the A community or both in its main B
    community) by status. The expression change of a community is the mean RNA log2FC of its bait genes (A and B
    members); its coherence is tested against random sets of the same number of genes drawn from all genes with an
    expression change (two-sided, empirical p, BH q). Communities with fewer than `min_genes` genes are skipped; A
    communities that continue in the same B community are tested once (column communities_a lists them). The null of
    each gene-set size k is drawn from its own generator (seed, k), so sizes are computed in parallel over `threads`
    processes with a result independent of their number.
    """
    values = np.array(list(expression.values()), dtype=float)
    of_a = {n: i for i, members in enumerate(communities_a) for n in members}
    of_b = {n: i for i, members in enumerate(communities_b) for n in members}
    b_to_a = defaultdict(set)
    for index_a, index_b in main.items():
        b_to_a[index_b].add(index_a)
    edge_counts = defaultdict(lambda: defaultdict(int))
    for u, v, status in edges[['node_a', 'node_b', 'status']].itertuples(index=False):
        targets = set()
        if u in of_a and of_a.get(u) == of_a.get(v):
            targets.add(of_a[u])
        if u in of_b and of_b.get(u) == of_b.get(v):
            targets |= b_to_a.get(of_b[u], set())
        for index_a in targets:
            edge_counts[index_a][status] += 1
    transition_of = dict(zip(transitions.loc[transitions['side'] == 'A', 'community_a'],
                             transitions.loc[transitions['side'] == 'A', 'transition']))
    rows = []
    selected = []
    for index_a, members_a in enumerate(communities_a):
        index_b = main.get(index_a, -1)
        members_b = communities_b[index_b] if index_b >= 0 else set()
        genes = sorted(n for n in members_a | members_b if n in expression)
        if len(genes) >= min_genes:
            selected.append((index_a, index_b, genes))
    sizes = sorted({len(genes) for _, _, genes in selected})
    nulls = dict(zip(sizes, _pool_map(_community_null, sizes, threads,
                                      dict(values=values, n_permutations=n_permutations, seed=seed))))
    for index_a, index_b, genes in selected:
        counts = edge_counts[index_a]
        total = sum(counts.values())
        mean_change = float(np.mean([expression[n] for n in genes]))
        k = len(genes)
        null = nulls[k]
        p = (np.sum(np.abs(null) >= abs(mean_change)) + 1) / (n_permutations + 1)
        rows.append({'community_a': index_a, 'community_b': index_b, 'transition': transition_of.get(index_a, ''),
                     'n_genes': k, 'edges_shared': counts['shared'], 'edges_lost': counts['lost'],
                     'edges_gained': counts['gained'],
                     'rewired_fraction': (counts['lost'] + counts['gained']) / total if total else np.nan,
                     'mean_rna_log2FC': mean_change, 'coherence_p': float(p)})
    table = pd.DataFrame(rows)
    if len(table):
        # A communities with the same main B community form one unit: keep the pair with most genes, list all A ids;
        # A communities without a main B community (-1) stay separate units
        matched = table[table['community_b'] >= 0]
        sources = matched.groupby('community_b')['community_a'].apply(lambda v: ','.join(str(x) for x in sorted(v)))
        matched = matched.sort_values(['n_genes', 'community_a'], ascending=[False, True], kind='mergesort')
        matched = matched.drop_duplicates('community_b')
        matched = matched.assign(communities_a=matched['community_b'].map(sources))
        unmatched = table[table['community_b'] < 0]
        unmatched = unmatched.assign(communities_a=unmatched['community_a'].astype(str))
        table = pd.concat([matched, unmatched], ignore_index=True)
        table['coherence_q'] = _benjamini_hochberg(table['coherence_p'].to_numpy())
        table = table.sort_values('coherence_p').reset_index(drop=True)
    return table


def match_communities(communities_a, communities_b, bait_nodes):
    """Best matching community in B for every community in A by node and bait Jaccard index."""
    node_to_b = {}
    for index, members in enumerate(communities_b):
        for node in members:
            node_to_b[node] = index
    rows = []
    for index_a, members in enumerate(communities_a):
        overlap = defaultdict(int)
        for node in members:
            if node in node_to_b:
                overlap[node_to_b[node]] += 1
        best, best_jaccard = -1, 0.0
        for index_b, shared in overlap.items():
            jaccard = shared / len(members | communities_b[index_b])
            if jaccard > best_jaccard:
                best, best_jaccard = index_b, jaccard
        bait_jaccard = np.nan
        if best >= 0:
            baits_a = {n for n in members if n in bait_nodes}
            baits_b = {n for n in communities_b[best] if n in bait_nodes}
            union = baits_a | baits_b
            bait_jaccard = len(baits_a & baits_b) / len(union) if union else np.nan
        rows.append({'community_a': index_a, 'size_a': len(members), 'best_community_b': best,
                     'size_b': len(communities_b[best]) if best >= 0 else 0, 'jaccard': best_jaccard,
                     'bait_jaccard': bait_jaccard})
    return pd.DataFrame(rows)
