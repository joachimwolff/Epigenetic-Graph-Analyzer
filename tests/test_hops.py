"""k-hop levels of the differential mode: operators, classification, candidates and per-level tests."""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer.differential import (ActivityOperators, classify_genes, gene_features,
                                                  gene_level_hop_test, gene_level_indirect_test, hop_levels,
                                                  indirect_only_genes, level_hops, level_name, mirror_comparison,
                                                  nested_models, permutation_test, proximity_columns,
                                                  read_proximity_support)


def weighted_graph(edges):
    """Union graph from (u, v, r_a, r_b) with the node attributes of union_graph."""
    graph = nx.Graph()
    for k, (u, v, r_a, r_b) in enumerate(edges):
        graph.add_edge(u, v, r_a=r_a, r_b=r_b, log2FC=0.0, qvalue=np.nan, status='shared')
    for i, n in enumerate(sorted(graph.nodes)):
        graph.nodes[n].update(chrom='chr1', start=100000 * (i + 1), end=100000 * (i + 1) + 4999, bait=True,
                              genes='G' + n)
    return graph


def row(operators, level, source, condition='a'):
    """{node: weight} of one source at one level."""
    H = operators.hop_a if condition == 'a' else operators.hop_b
    i = operators.sources.index(source)
    m = H[level_hops(level) if isinstance(level, str) else level].getrow(i)
    return {operators.nodes[j]: v for j, v in zip(m.indices, m.data)}


def test_level_names():
    assert [level_name(k) for k in (1, 2, 3, 5)] == ['direct', 'indirect', 'hop3', 'hop5']
    assert [level_hops(n) for n in ('promoter', 'direct', 'indirect', 'hop3', 'hop12')] == [0, 1, 2, 3, 12]
    with pytest.raises(ValueError):
        ActivityOperators(weighted_graph([('a', 'b', 1, 1)]), ['a'], max_hops=1)


def test_path_graph_exact_distance_and_path_products():
    # a - b - c - d - e with distinct strengths
    r = {('a', 'b'): 2.0, ('b', 'c'): 3.0, ('c', 'd'): 5.0, ('d', 'e'): 7.0}
    graph = weighted_graph([(u, v, w, 2 * w) for (u, v), w in r.items()])
    operators = ActivityOperators(graph, ['a', 'c'], max_hops=4)
    assert operators.levels == ['direct', 'indirect', 'hop3', 'hop4']
    assert row(operators, 1, 'a') == {'b': 2.0}
    assert row(operators, 2, 'a') == pytest.approx({'c': 2.0 * 3.0})
    assert row(operators, 3, 'a') == pytest.approx({'d': 2.0 * 3.0 * 5.0})
    assert row(operators, 4, 'a') == pytest.approx({'e': 2.0 * 3.0 * 5.0 * 7.0})
    # condition B doubles every strength: the k-hop weight scales by 2^k
    assert row(operators, 4, 'a', 'b') == pytest.approx({'e': 16 * 2.0 * 3.0 * 5.0 * 7.0})
    # from the middle: both sides at the same distance, nothing beyond the ends; walks back to the source or to
    # closer nodes never count
    assert row(operators, 1, 'c') == {'b': 3.0, 'd': 5.0}
    assert row(operators, 2, 'c') == pytest.approx({'a': 3.0 * 2.0, 'e': 5.0 * 7.0})
    assert row(operators, 3, 'c') == {} and row(operators, 4, 'c') == {}
    i = operators.sources.index('c')
    assert [operators.n_hop[k][i] for k in (1, 2, 3, 4)] == [2, 2, 0, 0]
    assert operators.n_indirect[i] == 2 and operators.indirect_a is operators.hop_a[2]


def test_triangle_and_square_masking():
    # triangle a, b, c with a tail c - d: c is a direct neighbour of a, so the walk a-b-c (two steps) does not
    # put c on the 2-hop level; d is at exactly two steps through c only
    graph = weighted_graph([('a', 'b', 2.0, 2.0), ('b', 'c', 3.0, 3.0), ('a', 'c', 5.0, 5.0), ('c', 'd', 7.0, 7.0)])
    operators = ActivityOperators(graph, ['a'], max_hops=3)
    assert set(row(operators, 1, 'a')) == {'b', 'c'}
    assert row(operators, 2, 'a') == pytest.approx({'d': 5.0 * 7.0})
    assert row(operators, 3, 'a') == {}
    # square a - b - c - d - a: c at two steps over two shortest paths, the weights add up
    graph = weighted_graph([('a', 'b', 2.0, 1.0), ('b', 'c', 3.0, 1.0), ('c', 'd', 5.0, 1.0), ('d', 'a', 7.0, 1.0),
                            ('c', 'e', 11.0, 1.0)])
    operators = ActivityOperators(graph, ['a'], max_hops=3)
    assert row(operators, 2, 'a') == pytest.approx({'c': 2.0 * 3.0 + 7.0 * 5.0})
    assert row(operators, 3, 'a') == pytest.approx({'e': (2.0 * 3.0 + 7.0 * 5.0) * 11.0})
    # the first two levels do not depend on max_hops
    two = ActivityOperators(graph, ['a'], max_hops=2)
    for k in (1, 2):
        assert (two.hop_a[k] != operators.hop_a[k]).nnz == 0


def features_table():
    """Five differential genes explained at promoter, 1, 2, 3 hops and not at all (H3K27me3 repressive)."""
    return pd.DataFrame({
        'node': ['p', 'd', 'i', 'h', 'u'], 'genes': ['P', 'D', 'I', 'H', 'U'],
        'rna_log2FC': [2.0, -2.0, 2.0, -2.0, 2.0],
        'n_direct': [1, 2, 2, 3, 1], 'n_indirect': [3, 3, 4, 5, 0], 'n_hop3': [6, 7, 8, 9, 0],
        'H3K27ac_promoter': [1.5, 0, 0, 0, 0], 'H3K27ac_direct': [0, -1.2, 0, 0, 0],
        'H3K27ac_indirect': [0, 0, 1.0, 0.0, 0], 'H3K27ac_hop3': [3, 0, 0, 0, -2.0],
        'H3K27me3_promoter': [0, 0, 0, 0, 0], 'H3K27me3_direct': [0, 0, 0, 0, 0],
        'H3K27me3_indirect': [0, 0, 0, -1.5, 0], 'H3K27me3_hop3': [0, 0, 0, 1.1, 0]})


def test_classification_with_three_hops(tmp_path):
    features = features_table()
    marks = ['H3K27ac', 'H3K27me3']
    assert hop_levels(features, marks) == ['direct', 'indirect', 'hop3']
    curve = tmp_path / 'hops.tsv'
    pd.DataFrame({'steps': [1, 2, 3], 'concordance': [0.66, 0.58, 0.54], 'ci_low': [0.64, 0.56, 0.52],
                  'ci_high': [0.68, 0.60, 0.56]}).to_csv(curve, sep='\t', index=False)
    support = read_proximity_support(str(curve))
    classes = classify_genes(features, marks, {'H3K27me3': -1}, proximity=support).set_index('node')
    # h: H3K27me3 -1.5 at two hops is the wrong direction for a down-regulated gene, +1.1 at three hops is right
    assert list(classes['explanation']) == ['promoter', 'direct', 'indirect', 'hop3', 'unexplained']
    assert list(classes['explaining_marks']) == ['H3K27ac', 'H3K27ac', 'H3K27ac', 'H3K27me3', '']
    assert classes['explaining_hops'].tolist()[:4] == [0, 1, 2, 3] and pd.isna(classes.loc['u', 'explaining_hops'])
    assert classes['n_hop_nodes_tested'].tolist()[:4] == [1, 2, 4, 9]
    assert classes.loc['p', 'proximity_support'] == 1.0 and classes.loc['p', 'proximity_support_ci'] == 'NA'
    assert classes.loc['h', 'proximity_support'] == 0.54 and classes.loc['h', 'proximity_support_ci'] == '0.52-0.56'
    assert np.isnan(classes.loc['u', 'proximity_support'])
    # levels beyond the step curve have no support value
    values, intervals = proximity_columns([4], support)
    assert np.isnan(values[0]) and intervals == ['']
    # without the 3-hop columns the classification is the former promoter / direct / indirect one
    two = classify_genes(features.drop(columns=['H3K27ac_hop3', 'H3K27me3_hop3']), marks, {'H3K27me3': -1})
    assert list(two['explanation']) == ['promoter', 'direct', 'indirect', 'unexplained', 'unexplained']
    models = nested_models(features.assign(rna_log2FC=[2.0, -2.0, 2.0, -2.0, 1.0]), marks)
    assert list(models['model']) == ['promoter', 'promoter+direct', 'promoter+direct+indirect',
                                     'promoter+direct+indirect+hop3']


def test_explanation_cascade():
    from epigeneticgraphanalyzer.differential import cascade_table, explanation_cascade
    marks = ['H3K27ac', 'H3K27me3']
    features = features_table()
    # own interactions: d (down) loses strength -> explained at the gene by interactions only (before its direct
    # mark change); p (up) gains -> marks and interactions; for i (up) and h (down) the change is discordant
    features['contact_direct'] = [1.2, -1.0, -1.5, 0.5, 0.0]
    # a sixth gene explained at three hops only by H3K27ac, with a concordant mirror (not beyond the mirror)
    extra = features.iloc[[4]].assign(node='m', genes='M', rna_log2FC=-2.0, H3K27ac_hop3=-2.0, n_hop3=5)
    features = pd.concat([features, extra], ignore_index=True)
    for mark in marks:
        for level in ('direct', 'indirect', 'hop3'):
            features['{}_{}_mirror'.format(mark, level)] = 0.0
    features.loc[features['node'] == 'm', 'H3K27ac_hop3_mirror'] = -1.5
    classes = classify_genes(features, marks, {'H3K27me3': -1})
    cascade = explanation_cascade(classes, marks, {'H3K27me3': -1})
    assert list(cascade['explained_at']) == ['gene', 'gene', 'step2', 'step3', 'unexplained', 'step3']
    assert list(cascade['explained_by'][:3]) == ['both', 'interactions', 'H3K27ac']
    assert cascade.loc['h', 'explained_by'] == 'H3K27me3' and cascade.loc['m', 'explained_by'] == 'H3K27ac'
    assert bool(cascade.loc['h', 'beyond_mirror']) and not bool(cascade.loc['m', 'beyond_mirror'])
    assert cascade['nodes_tested'].tolist()[:4] == [1, 1, 4, 9] and cascade.loc['m', 'nodes_tested'] == 5
    assert pd.isna(cascade.loc['u', 'nodes_tested']) and pd.isna(cascade.loc['u', 'beyond_mirror'])
    # without the interaction rule the cascade is the classification
    plain = explanation_cascade(classes, marks, {'H3K27me3': -1}, contact_threshold=np.inf)
    assert list(plain['explained_at'][:4]) == ['gene', 'step1', 'step2', 'step3']
    support = {1: (0.66, 0.64, 0.68), 2: (0.58, 0.56, 0.60), 3: (0.54, 0.52, 0.56)}
    table = cascade_table(cascade, dict(zip(classes['node'], classes['genes'])),
                          {'i': 0.01, 'h': 0.001, 'm': 0.002}, support, max_hops=3).set_index('level')
    assert list(table.index) == ['gene', 'gene_marks_only', 'gene_interactions_only', 'gene_both', 'step1', 'step2',
                                 'step3', 'unexplained']
    assert table['genes'].tolist() == [2, 0, 1, 1, 0, 1, 2, 1]
    assert table['genes_reaching'].tolist() == [6, 6, 6, 6, 4, 4, 3, 1]
    assert table.loc['step3', 'fraction_of_remaining'] == 2 / 3
    assert table.loc['step3', 'beyond_mirror'] == 1 and table.loc['step3', 'p_lt_0.05'] == 2
    assert table.loc['step3', 'beyond_mirror_p_lt_0.05'] == 1
    # BH across the two genes at step 3: q = 0.002 for both
    assert table.loc['step3', 'q_lt_0.05'] == 2 and table.loc['step3', 'genes_q_lt_0.05'] == 'H;M'
    assert table.loc['step3', 'proximity_support'] == 0.54 and table.loc['step3', 'proximity_ci'] == '0.52-0.56'
    assert table.loc['step3', 'median_nodes_tested'] == 7.0


def chain_design(n_units=40, n_planted=24, seed=0):
    """Units S - O1 - O2 - T (T at three hops of S); the RNA change of S follows the H3K27ac change at T."""
    rng = np.random.default_rng(seed)
    edges, units = [], []
    for k in range(n_units):
        s, o1, o2, t = ('{}{}'.format(x, k) for x in ('S', 'O', 'Q', 'T'))
        edges += [(s, o1, 2.0, 2.0), (o1, o2, 3.0, 3.0), (o2, t, 1.5, 1.5)]
        delta = rng.normal(0, 2)
        rna = delta + rng.normal(0, 0.3) if k < n_planted else rng.normal(0, 2)
        units.append(dict(S=s, O=o1, Q=o2, T=t, delta=delta, rna=rna, planted=k < n_planted))
    graph = weighted_graph(edges)
    signal_a = {n: 4.0 for n in graph.nodes}
    signal_b = dict(signal_a)
    for u in units:
        signal_b[u['T']] = 4.0 * 2 ** u['delta']
    return graph, {'H3K27ac': (signal_a, signal_b)}, {u['S']: u['rna'] for u in units}, units


def test_three_hop_candidates_and_per_level_test():
    graph, marks, change, units = chain_design()
    operators = ActivityOperators(graph, list(change), max_hops=3)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=vectors)
    assert np.allclose(features[['H3K27ac_promoter', 'H3K27ac_direct', 'H3K27ac_indirect']], 0)
    assert (features['n_hop3'] == 1).all()
    signs = {'H3K27ac': 1}
    classes = classify_genes(features, ['H3K27ac'], signs)
    assert set(classes['explanation']) <= {'hop3', 'unexplained'}
    report = indirect_only_genes(graph, operators, vectors, features, ['H3K27ac'], signs)
    assert set(report['explaining_hops']) == {3}
    assert set(report['node']) == set(classes.loc[classes['explanation'] == 'hop3', 'node'])
    first = report.iloc[0]
    unit = next(u for u in units if u['S'] == first['node'])
    # the driver is T, reached through O and Q
    assert '{}[G{}] via G{}>G{}'.format(unit['T'], unit['T'], unit['O'], unit['Q']) in first['drivers']
    # with only two hops nothing is explained beyond the promoter and direct levels
    two = gene_features(graph, marks, change, max_hops=2)
    assert len(indirect_only_genes(graph, ActivityOperators(graph, list(change)), vectors, two, ['H3K27ac'],
                                   signs)) == 0
    tests = gene_level_hop_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                                n_permutations=200, genes={'indirect': None, 'hop3': list(report['node'])})
    hop3 = tests['hop3'].set_index('node')
    assert set(hop3.index) == set(report['node'])
    strong = [u['S'] for u in units if u['planted'] and abs(u['delta']) >= 2 and u['S'] in hop3.index]
    assert hop3.loc[strong, 'hop3_p'].median() < 0.05
    # nobody has a 2-hop signal: the 2-hop P-values are 1 or large, and the same permutations give the same
    # 2-hop result as the single-level test
    single = gene_level_indirect_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                                      n_permutations=200)
    pd.testing.assert_frame_equal(single, tests['indirect'])
    result = permutation_test(graph, marks, change, signs, n_permutations=20, max_hops=3)
    assert result['observed_hop3_delta_r2'] > result['null_hop3_delta_r2_mean']
    assert result['empirical_p_delta_r2_hop3'] < 0.1
    assert 'observed_indirect_delta_r2' in result and result['n_permutations'] == 20


def test_mirror_comparison_levels():
    rng = np.random.default_rng(0)
    n = 80
    features = pd.DataFrame({'rna_log2FC': rng.normal(0, 1, n)})
    for level in ('promoter', 'direct', 'indirect', 'hop3'):
        features['M_' + level] = rng.normal(0, 1, n)
        if level != 'promoter':
            features['M_{}_mirror'.format(level)] = rng.normal(0, 1, n)
    features['M_hop3'] = features['rna_log2FC'] + rng.normal(0, 0.2, n)
    result = mirror_comparison(features, ['M'], n_bootstrap=50)
    assert result['gain_hop3_graph'] > 0.5 and result['hop3_ci_low'] > 0
    assert set(k for k in result if k.startswith('r2_')) == {
        'r2_promoter', 'r2_promoter_direct', 'r2_promoter_direct_indirect', 'r2_promoter_direct_indirect_hop3'}


def test_parallel_blocks_and_vectorised_statistic(monkeypatch):
    from epigeneticgraphanalyzer import differential
    graph, marks, change, units = chain_design(n_units=30, n_planted=20)
    operators = ActivityOperators(graph, list(change), max_hops=3)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=vectors)
    signs = {'H3K27ac': 1}
    levels = {'indirect': None, 'hop3': list(features['node'][:12])}
    # 2,500 permutations = three blocks; the result does not depend on the number of worker processes
    one = gene_level_hop_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                              n_permutations=2500, seed=3, genes=levels, threads=1)
    four = gene_level_hop_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                               n_permutations=2500, seed=3, genes=levels, threads=4)
    for level in levels:
        pd.testing.assert_frame_equal(one[level], four[level])
    # nor on the chunk size within a block
    monkeypatch.setattr(differential, 'BLOCK_MEMORY', 24 * 200 * 7)
    small = gene_level_hop_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                                n_permutations=2500, seed=3, genes=levels, threads=2)
    for level in levels:
        pd.testing.assert_frame_equal(one[level], small[level])
    # the vectorised statistic equals the per-permutation loop over the same permutations (block 0)
    n = 300
    result = gene_level_hop_test(operators, vectors, features, ['H3K27ac'], signs, rna_threshold=0.5,
                                 n_permutations=n, seed=5, genes={'hop3': None})['hop3']
    de = features[features['rna_log2FC'].abs() >= 0.5]
    rows = np.array([operators.sources.index(x) for x in de['node']])
    direction = np.sign(de['rna_log2FC'].to_numpy())
    H_a, H_b = operators.hop_a[3][rows], operators.hop_b[3][rows]
    vector_a, vector_b = vectors['H3K27ac']

    def statistic(a, b):
        return np.log2((H_b @ b + 1.0) / (H_a @ a + 1.0)) * direction

    observed = statistic(vector_a, vector_b)
    assert np.allclose(result['hop3_statistic'], observed, rtol=0, atol=1e-12)
    is_source = np.zeros(len(operators.nodes), dtype=bool)
    is_source[operators.source_index] = True
    others = np.nonzero(~is_source)[0]
    rng = np.random.default_rng(np.random.SeedSequence([5, 0]))
    exceed = np.zeros(len(rows))
    for _ in range(n):
        donors = others[rng.permutation(len(others))]
        a, b = vector_a.copy(), vector_b.copy()
        a[others], b[others] = vector_a[donors], vector_b[donors]
        exceed += statistic(a, b) >= observed
    assert np.array_equal(result['hop3_p'].to_numpy(), (exceed + 1) / (n + 1))


def test_community_null_threads_and_sampling():
    from epigeneticgraphanalyzer.differential import _sample_without_replacement, community_expression_changes
    rng = np.random.default_rng(0)
    for k in (3, 30, 90):
        index = _sample_without_replacement(rng, 100, k, 500)
        assert index.shape == (500, k) and index.min() >= 0 and index.max() < 100
        assert all(len(set(r)) == k for r in index)
    nodes = ['n{}'.format(i) for i in range(400)]
    expression = dict(zip(nodes, np.random.default_rng(1).normal(0, 1, 400)))
    communities_a = [set(nodes[i:i + s]) for i, s in ((0, 5), (10, 12), (40, 30), (100, 7))]
    communities_b = [set(c) for c in communities_a]
    transitions = pd.DataFrame({'side': ['A'] * 4, 'community_a': [0, 1, 2, 3], 'transition': ['stable'] * 4})
    edges = pd.DataFrame({'node_a': [], 'node_b': [], 'status': []})
    main = {i: i for i in range(4)}
    one = community_expression_changes(communities_a, communities_b, main, transitions, edges, expression,
                                       n_permutations=3000, threads=1)
    four = community_expression_changes(communities_a, communities_b, main, transitions, edges, expression,
                                        n_permutations=3000, threads=4)
    pd.testing.assert_frame_equal(one, four)
    assert len(one) == 4 and one['coherence_p'].between(1 / 3001, 1).all()


def test_mirror_includes_neighbours_gained_in_b(tmp_path):
    """A node reached in condition B only (through a gained edge) must enter the mirrored B activity."""
    from epigeneticgraphanalyzer.differential import mirrored_features
    from tests.synthetic import write_bigwig
    # a at 100 kb, b at 200 kb, c at 300 kb; source b; b-c shared, c-a gained in B, so a is 2 hops away in B only
    graph = weighted_graph([('b', 'c', 1.0, 1.0), ('c', 'a', 0.0, 1.0)])
    operators = ActivityOperators(graph, ['b'], max_hops=2)
    # the mirror of a around b lies at 300-305 kb (fragment 60)
    track_a = write_bigwig(str(tmp_path / 'a.bw'), {}, default=1.0)
    track_b = write_bigwig(str(tmp_path / 'b.bw'), {('chr1', 60): 10.0}, default=1.0)
    mirror = mirrored_features(graph, operators, {'m': [track_a]}, {'m': [track_b]})
    assert mirror.loc['b', 'm_indirect_mirror'] == pytest.approx(np.log2(11.0))
