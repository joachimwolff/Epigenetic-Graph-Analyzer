"""Differential analysis of two capture Hi-C graphs with chromatin and expression changes.

Edges
    The union graph holds the interactions called (e.g. by CHiCAGO) in either
    condition. The strength of an edge in a condition is its score divided by
    the calling threshold (1 at the threshold, 0 where the pair was not
    called). A Chicdiff results table, when given, adds the tested change of
    each interaction (log2 fold change and adjusted p-value).

Genes
    For every bait gene the chromatin change around it is summarised on three
    levels of the union graph, each weighting a node's mark signal by the edge
    strength r in the respective condition:

    * promoter: the mark at the bait fragment itself,
    * direct:   sum over neighbours n of r(bait, n) * signal(n),
    * indirect: sum over nodes m at exactly two hops of
                sum over mediators n of r(bait, n) * r(n, m) * signal(m).

    The feature is log2(activity_B / activity_A). Expression changes are
    related to these features by nested linear models and by a hierarchical
    classification (promoter, direct, indirect, unexplained). The contribution
    of the indirect level is assessed by permuting node signals.
"""

import logging
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


class ActivityOperators:
    """Sparse operators mapping a node signal vector to direct and indirect activities of source nodes.

    direct_c = D_c @ signal_c and indirect_c = I_c @ signal_c with
    D_c[s, n] = r_c(s, n) and I_c[s, m] = sum_n r_c(s, n) r_c(n, m) restricted to
    nodes m at exactly two hops from s.
    """

    def __init__(self, graph, sources):
        import scipy.sparse as sp
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
        reach2 = (A_s @ A).tocsr()
        reach2.data[:] = 1.0
        # remove direct neighbours and the source itself from the two-hop set
        self_mask = sp.csr_matrix((np.ones(len(src)), (np.arange(len(src)), src)), shape=(len(src), n))
        excluded = A_s + self_mask
        excluded.data[:] = 1.0
        mask = reach2 - reach2.multiply(excluded)
        mask.eliminate_zeros()
        self.direct_a = W_a[src]
        self.direct_b = W_b[src]
        self.indirect_a = (W_a[src] @ W_a).multiply(mask).tocsr()
        self.indirect_b = (W_b[src] @ W_b).multiply(mask).tocsr()
        self.n_direct = np.asarray((self.direct_a > 0).sum(axis=1)).ravel()
        self.n_indirect = np.asarray((mask > 0).sum(axis=1)).ravel()
        self.source_index = src

    def vector(self, signal):
        return np.array([signal.get(n, 0.0) for n in self.nodes], dtype=float)


def gene_features(graph, marks, expression_change, pseudo_count=1.0, operators=None, mark_vectors=None):
    """Per bait node: expression change and promoter/direct/indirect log2 activity changes per mark.

    `marks` maps mark name -> (signal_a dict, signal_b dict);
    `expression_change` maps bait node -> RNA log2 fold change. Precomputed
    `operators` and `mark_vectors` (name -> (vector_a, vector_b)) are used when given.
    """
    if operators is None:
        operators = ActivityOperators(graph, list(expression_change))
    if mark_vectors is None:
        mark_vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    sources = operators.sources
    features = pd.DataFrame({
        'node': sources,
        'genes': [graph.nodes[n].get('genes', '') for n in sources],
        'rna_log2FC': [expression_change[n] for n in sources],
        'n_direct': operators.n_direct,
        'n_indirect': operators.n_indirect,
    })
    ones = np.ones(len(operators.nodes))
    features['contact_direct'] = np.log2((operators.direct_b @ ones + pseudo_count) /
                                         (operators.direct_a @ ones + pseudo_count))
    features['contact_indirect'] = np.log2((operators.indirect_b @ ones + pseudo_count) /
                                           (operators.indirect_a @ ones + pseudo_count))
    src = operators.source_index
    for mark, (vector_a, vector_b) in mark_vectors.items():
        features[mark + '_promoter'] = np.log2((vector_b[src] + pseudo_count) / (vector_a[src] + pseudo_count))
        features[mark + '_direct'] = np.log2((operators.direct_b @ vector_b + pseudo_count) /
                                             (operators.direct_a @ vector_a + pseudo_count))
        features[mark + '_indirect'] = np.log2((operators.indirect_b @ vector_b + pseudo_count) /
                                               (operators.indirect_a @ vector_a + pseudo_count))
    return features


def mirrored_features(graph, operators, tracks_a, tracks_b, pseudo_count=1.0):
    """Linear-proximity control of the direct and indirect features.

    Every neighbour m of a source s is replaced by the interval at the same
    distance on the opposite side of s (mirrored around the midpoint of s),
    keeping the contact weights. The mark signal of the mirrored intervals is
    read from the bigWig files, so the control has the genomic distance profile
    of the real neighbourhood without its 3D topology.
    """
    from .annotation import bigwig_interval_signal
    src_nodes = operators.sources
    columns = {}
    for level, (W_a, W_b) in (('direct', (operators.direct_a, operators.direct_b)),
                              ('indirect', (operators.indirect_a, operators.indirect_b))):
        coo_a = W_a.tocoo()
        coo_b = W_b.tocsr()
        rows, cols = coo_a.row, coo_a.col
        weight_a = coo_a.data
        weight_b = np.asarray(coo_b[rows, cols]).ravel()
        chroms, starts, ends = [], [], []
        for i, j in zip(rows, cols):
            source = graph.nodes[src_nodes[i]]
            node = graph.nodes[operators.nodes[j]]
            centre = source['start'] + source['end']
            chroms.append(source['chrom'])
            starts.append(centre - node['end'])
            ends.append(centre - node['start'])
        for mark in tracks_a:
            signal_a = np.mean([bigwig_interval_signal(p, chroms, starts, ends) for p in tracks_a[mark]], axis=0) \
                if len(rows) else np.zeros(0)
            signal_b = np.mean([bigwig_interval_signal(p, chroms, starts, ends) for p in tracks_b[mark]], axis=0) \
                if len(rows) else np.zeros(0)
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


def nested_models(features, marks):
    """R^2 of RNA change explained by promoter, + direct, + indirect features; partial F tests."""
    levels = ['promoter', 'direct', 'indirect']
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


def mirror_comparison(features, marks, n_bootstrap=1000, seed=0):
    """Variance in RNA change gained by graph neighbours versus their mirrored loci.

    Returns the R^2 of the nested models and the gains of the direct and the
    indirect level over their baselines for the graph and the mirror features,
    with paired bootstrap intervals (over genes) of graph minus mirror.
    """
    promoter = [m + '_promoter' for m in marks]
    direct = [m + '_direct' for m in marks]
    indirect = [m + '_indirect' for m in marks]
    direct_m = [m + '_direct_mirror' for m in marks]
    indirect_m = [m + '_indirect_mirror' for m in marks]

    def gains(df):
        r_p = model_r2(df, promoter)
        r_pd = model_r2(df, promoter + direct)
        r_pdi = model_r2(df, promoter + direct + indirect)
        r_pdm = model_r2(df, promoter + direct_m)
        r_pdii = model_r2(df, promoter + direct + indirect_m)
        return {'r2_promoter': r_p, 'r2_promoter_direct': r_pd, 'r2_promoter_direct_indirect': r_pdi,
                'gain_direct_graph': r_pd - r_p, 'gain_direct_mirror': r_pdm - r_p,
                'gain_indirect_graph': r_pdi - r_pd, 'gain_indirect_mirror': r_pdii - r_pd}

    observed = gains(features)
    rng = np.random.default_rng(seed)
    boot_direct, boot_indirect = [], []
    for _ in range(n_bootstrap):
        sample = features.iloc[rng.integers(0, len(features), len(features))]
        g = gains(sample)
        boot_direct.append(g['gain_direct_graph'] - g['gain_direct_mirror'])
        boot_indirect.append(g['gain_indirect_graph'] - g['gain_indirect_mirror'])
    boot_direct = np.array(boot_direct)
    boot_indirect = np.array(boot_indirect)
    observed.update({
        'n_genes': len(features),
        'direct_graph_minus_mirror': observed['gain_direct_graph'] - observed['gain_direct_mirror'],
        'direct_ci_low': float(np.quantile(boot_direct, 0.025)), 'direct_ci_high': float(np.quantile(boot_direct, 0.975)),
        'direct_fraction_bootstrap_above_0': float(np.mean(boot_direct > 0)),
        'indirect_graph_minus_mirror': observed['gain_indirect_graph'] - observed['gain_indirect_mirror'],
        'indirect_ci_low': float(np.quantile(boot_indirect, 0.025)),
        'indirect_ci_high': float(np.quantile(boot_indirect, 0.975)),
        'indirect_fraction_bootstrap_above_0': float(np.mean(boot_indirect > 0)),
        'n_bootstrap': n_bootstrap,
    })
    return observed


def feature_correlations(features, marks):
    rows = []
    levels = ['promoter', 'direct', 'indirect']
    if any(m + '_direct_mirror' in features for m in marks):
        levels += ['direct_mirror', 'indirect_mirror']
    for level in levels:
        for mark in marks:
            column = mark + '_' + level
            data = features[['rna_log2FC', column]].dropna()
            if len(data) < 3:
                continue
            rho, p = stats.spearmanr(data['rna_log2FC'], data[column])
            rows.append({'mark': mark, 'level': level, 'n': len(data), 'spearman_rho': rho, 'pvalue': p})
    for level in ('direct', 'indirect'):
        column = 'contact_' + level
        if column in features:
            data = features[['rna_log2FC', column]].dropna()
            if len(data) >= 3:
                rho, p = stats.spearmanr(data['rna_log2FC'], data[column])
                rows.append({'mark': 'contact', 'level': level, 'n': len(data), 'spearman_rho': rho, 'pvalue': p})
    return pd.DataFrame(rows)


def classify_genes(features, marks, mark_signs, rna_threshold=1.0, chromatin_threshold=1.0):
    """Hierarchical explanation of expression changes: promoter > direct > indirect > unexplained.

    A level explains a gene if at least one mark changes by at least
    `chromatin_threshold` (log2) on that level in the direction expected from
    the mark sign (+1 activating, -1 repressive) and the RNA change.
    """
    de = features[features['rna_log2FC'].abs() >= rna_threshold].copy()
    direction = np.sign(de['rna_log2FC'].to_numpy())
    explanation = np.array(['unexplained'] * len(de), dtype=object)
    explaining_marks = [''] * len(de)
    for level in ('indirect', 'direct', 'promoter'):
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
    return de


def _concordant(values, mark_sign, direction, threshold):
    values = np.nan_to_num(np.asarray(values, dtype=float))
    return (np.sign(values) * mark_sign == direction) & (np.abs(values) >= threshold)


def indirect_only_genes(graph, operators, mark_vectors, features, marks, mark_signs, rna_threshold=1.0,
                        chromatin_threshold=1.0, genes=None, top=3):
    """Differential genes explained at indirect neighbours only.

    A gene qualifies if its expression changes at least `rna_threshold`
    (log2), no mark changes concordantly at its promoter or its direct
    neighbours, and at least one mark changes concordantly on the indirect
    level. `beyond_mirror` marks genes whose explaining marks show no
    concordant change at the mirrored loci. For each explaining mark the
    two-hop nodes with the largest concordant contribution are listed with
    their annotated genes and the mediators that connect them to the gene.
    """
    data = features if genes is None else features[features['node'].isin(genes)]
    data = data[data['rna_log2FC'].abs() >= rna_threshold]
    row_of = {node: i for i, node in enumerate(operators.sources)}
    rows = []
    for _, gene in data.iterrows():
        direction = np.sign(gene['rna_log2FC'])
        promoter = any(_concordant(gene[m + '_promoter'], mark_signs.get(m, 1), direction, chromatin_threshold)
                       for m in marks)
        direct = any(_concordant(gene[m + '_direct'], mark_signs.get(m, 1), direction, chromatin_threshold)
                     for m in marks)
        explaining = [m for m in marks
                      if _concordant(gene[m + '_indirect'], mark_signs.get(m, 1), direction, chromatin_threshold)]
        if promoter or direct or not explaining:
            continue
        mirrored = [m for m in explaining if m + '_indirect_mirror' in gene and
                    _concordant(gene[m + '_indirect_mirror'], mark_signs.get(m, 1), direction, chromatin_threshold)]
        i = row_of[gene['node']]
        drivers = []
        for mark in explaining:
            vector_a, vector_b = mark_vectors[mark]
            row_a = operators.indirect_a.getrow(i)
            row_b = operators.indirect_b.getrow(i)
            columns = np.union1d(row_a.indices, row_b.indices)
            weight_a = np.asarray(row_a[:, columns].todense()).ravel()
            weight_b = np.asarray(row_b[:, columns].todense()).ravel()
            contribution = (weight_b * vector_b[columns] - weight_a * vector_a[columns]) * mark_signs.get(mark, 1) * direction
            for k in np.argsort(contribution)[::-1][:top]:
                if contribution[k] <= 0:
                    break
                node = operators.nodes[columns[k]]
                mediators = set(graph.neighbors(gene['node'])) & set(graph.neighbors(node))
                drivers.append('{}:{}[{}] via {}'.format(
                    mark, node, graph.nodes[node].get('genes', '') or ('promoter' if graph.nodes[node].get('bait') else 'OE'),
                    '+'.join(sorted(graph.nodes[m].get('genes', '') or m for m in mediators))))
        row = {'node': gene['node'], 'genes': gene['genes'], 'rna_log2FC': gene['rna_log2FC'],
               'n_direct': gene['n_direct'], 'n_indirect': gene['n_indirect'],
               'explaining_marks': ','.join(explaining), 'mirror_concordant_marks': ','.join(mirrored),
               'beyond_mirror': len(mirrored) < len(explaining), 'drivers': '; '.join(drivers)}
        for mark in marks:
            for level in ('promoter', 'direct', 'indirect', 'indirect_mirror'):
                column = '{}_{}'.format(mark, level)
                if column in gene:
                    row[column] = gene[column]
        rows.append(row)
    result = pd.DataFrame(rows)
    if len(result):
        result = result.sort_values(['beyond_mirror', 'rna_log2FC'], ascending=[False, False],
                                    key=lambda c: c.abs() if c.name == 'rna_log2FC' else c)
    return result


def permutation_test(graph, marks, expression_change, mark_signs, n_permutations=100, seed=0,
                     rna_threshold=1.0, chromatin_threshold=1.0):
    """Null distribution of the indirect level's contribution.

    Mark signals of all non-source nodes are permuted jointly (A and B values of
    a node move together, so the distribution of changes is preserved) while the
    graph and the signals at the source promoters stay fixed.
    """
    rng = np.random.default_rng(seed)
    operators = ActivityOperators(graph, list(expression_change))
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    observed = gene_features(graph, marks, expression_change, operators=operators, mark_vectors=vectors)
    observed_models = nested_models(observed, list(marks))
    observed_class = classify_genes(observed, list(marks), mark_signs, rna_threshold, chromatin_threshold)
    obs_delta = float(observed_models['delta_r2'].iloc[-1])
    obs_fraction = float((observed_class['explanation'] == 'indirect').mean()) if len(observed_class) else np.nan
    is_source = np.zeros(len(operators.nodes), dtype=bool)
    is_source[operators.source_index] = True
    others = np.nonzero(~is_source)[0]
    null_delta, null_fraction = [], []
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
        null_delta.append(float(nested_models(features, list(marks))['delta_r2'].iloc[-1]))
        classes = classify_genes(features, list(marks), mark_signs, rna_threshold, chromatin_threshold)
        null_fraction.append(float((classes['explanation'] == 'indirect').mean()) if len(classes) else np.nan)
    null_delta = np.array(null_delta)
    null_fraction = np.array(null_fraction)
    return {
        'n_genes': int(len(observed)),
        'n_differential_genes': int(len(observed_class)),
        'observed_indirect_delta_r2': obs_delta,
        'null_indirect_delta_r2_mean': float(np.mean(null_delta)),
        'null_indirect_delta_r2_sd': float(np.std(null_delta)),
        'empirical_p_delta_r2': float((np.sum(null_delta >= obs_delta) + 1) / (len(null_delta) + 1)),
        'observed_fraction_indirect': obs_fraction,
        'null_fraction_indirect_mean': float(np.nanmean(null_fraction)),
        'null_fraction_indirect_sd': float(np.nanstd(null_fraction)),
        'empirical_p_fraction': float((np.sum(null_fraction >= obs_fraction) + 1) / (len(null_fraction) + 1)),
        'n_permutations': int(n_permutations),
    }


def gene_level_indirect_test(operators, mark_vectors, features, marks, mark_signs, rna_threshold=1.0,
                             n_permutations=1000, seed=0, pseudo_count=1.0):
    """Per-gene test of the indirect level for differential genes.

    The statistic of a gene is the largest concordant indirect change over the marks,
    max_m sign_m * sign(RNA) * indirect_m (log2 B/A). Its null distribution comes from permuting the mark signals
    of all non-source nodes jointly (as in permutation_test), which keeps the graph and the distribution of
    changes but breaks which nodes carry them. Returns node, statistic, empirical p-value and BH q-value.
    """
    de = features[features['rna_log2FC'].abs() >= rna_threshold]
    row_of = {node: i for i, node in enumerate(operators.sources)}
    rows = np.array([row_of[n] for n in de['node']], dtype=np.int64)
    direction = np.sign(de['rna_log2FC'].to_numpy())
    ind_a = operators.indirect_a[rows]
    ind_b = operators.indirect_b[rows]

    def statistic(vectors):
        best = np.full(len(rows), -np.inf)
        for mark in marks:
            vector_a, vector_b = vectors[mark]
            change = np.log2((ind_b @ vector_b + pseudo_count) / (ind_a @ vector_a + pseudo_count))
            best = np.maximum(best, change * mark_signs.get(mark, 1) * direction)
        return best

    observed = statistic(mark_vectors)
    rng = np.random.default_rng(seed)
    is_source = np.zeros(len(operators.nodes), dtype=bool)
    is_source[operators.source_index] = True
    others = np.nonzero(~is_source)[0]
    exceed = np.zeros(len(rows))
    for _ in range(n_permutations):
        donors = others[rng.permutation(len(others))]
        permuted = {}
        for mark, (vector_a, vector_b) in mark_vectors.items():
            new_a, new_b = vector_a.copy(), vector_b.copy()
            new_a[others] = vector_a[donors]
            new_b[others] = vector_b[donors]
            permuted[mark] = (new_a, new_b)
        exceed += statistic(permuted) >= observed
    pvalues = (exceed + 1) / (n_permutations + 1)
    # genes without two-hop neighbours carry no indirect information
    pvalues[operators.n_indirect[rows] == 0] = 1.0
    return pd.DataFrame({'node': de['node'].to_numpy(), 'indirect_statistic': observed,
                         'indirect_p': pvalues, 'indirect_q': _benjamini_hochberg(pvalues)})


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
            best_a_for_b[index_b] = max(counts, key=lambda k: counts[k] / len(members | communities_a[k]))
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
            best = max(counts, key=lambda k: counts[k] / len(members | communities_b[k]))
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


def community_expression_changes(communities_a, communities_b, main, transitions, edges, expression, min_genes=3,
                                 n_permutations=1000, seed=0):
    """Per community of A and its main B community: contact rewiring and coordinated expression change.

    Rewiring counts the union-graph edges inside the pair (both ends in the A community or both in its main B
    community) by status. The expression change of a community is the mean RNA log2FC of its bait genes (A and B
    members); its coherence is tested against random sets of the same number of genes drawn from all genes with an
    expression change (two-sided, empirical p, BH q). Communities with fewer than `min_genes` genes are skipped.
    """
    rng = np.random.default_rng(seed)
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
    null_cache = {}
    for index_a, members_a in enumerate(communities_a):
        index_b = main.get(index_a, -1)
        members_b = communities_b[index_b] if index_b >= 0 else set()
        nodes = members_a | members_b
        genes = [n for n in nodes if n in expression]
        if len(genes) < min_genes:
            continue
        counts = edge_counts[index_a]
        total = sum(counts.values())
        mean_change = float(np.mean([expression[n] for n in genes]))
        k = len(genes)
        if k not in null_cache:
            null_cache[k] = np.array([rng.choice(values, k, replace=False).mean() for _ in range(n_permutations)])
        null = null_cache[k]
        p = (np.sum(np.abs(null) >= abs(mean_change)) + 1) / (n_permutations + 1)
        rows.append({'community_a': index_a, 'community_b': index_b, 'transition': transition_of.get(index_a, ''),
                     'n_genes': k, 'edges_shared': counts['shared'], 'edges_lost': counts['lost'],
                     'edges_gained': counts['gained'],
                     'rewired_fraction': (counts['lost'] + counts['gained']) / total if total else np.nan,
                     'mean_rna_log2FC': mean_change, 'coherence_p': float(p)})
    table = pd.DataFrame(rows)
    if len(table):
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
