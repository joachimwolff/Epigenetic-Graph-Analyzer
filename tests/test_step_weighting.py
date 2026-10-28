"""Weighting of the step feature (--stepWeighting product | unweighted | normalised) on a tiny graph."""

import networkx as nx
import numpy as np
import pytest

from epigeneticgraphanalyzer.differential import ActivityOperators, gene_features, gene_level_hop_test, weight_paths

# p - a - m and p - b - m: two shortest paths to m; p - c - n: one path to n
EDGES = [('p', 'a'), ('p', 'b'), ('a', 'm'), ('b', 'm'), ('p', 'c'), ('c', 'n')]
R_A = dict(zip(EDGES, (1.0, 2.0, 3.0, 1.0, 1.0, 2.0)))
R_B = dict(zip(EDGES, (2.0, 2.0, 1.0, 1.0, 4.0, 0.5)))
X_A = {'p': 1.0, 'a': 2.0, 'b': 4.0, 'c': 6.0, 'm': 10.0, 'n': 2.0}
X_B = {'p': 1.0, 'a': 3.0, 'b': 4.0, 'c': 2.0, 'm': 30.0, 'n': 8.0}


def tiny_graph(r_a=R_A, r_b=R_B):
    graph = nx.Graph()
    for (u, v), a in r_a.items():
        graph.add_edge(u, v, r_a=a, r_b=r_b[(u, v)], log2FC=0.0, qvalue=np.nan, status='shared')
    for i, n in enumerate(sorted(graph.nodes)):
        graph.nodes[n].update(chrom='chr1', start=100000 * (i + 1), end=100000 * (i + 1) + 4999, bait=n == 'p',
                              genes='G' + n if n == 'p' else '')
    return graph


def row(graph, step_feature, step_weighting, x_a=X_A, x_b=X_B):
    operators = ActivityOperators(graph, ['p'], max_hops=2, step_feature=step_feature, step_weighting=step_weighting)
    table = gene_features(graph, {'K27ac': (x_a, x_b)}, {'p': 3.0}, operators=operators)
    return table.set_index('node').loc['p']


def log2(b, a):
    return np.log2((b + 1.0) / (a + 1.0))


def mean_r(edge):
    return (R_A[edge] + R_B[edge]) / 2


def path_weights(r):
    """Step 1 and step 2 path products of p (by hand)."""
    step1 = {'a': r[('p', 'a')], 'b': r[('p', 'b')], 'c': r[('p', 'c')]}
    step2 = {'m': r[('p', 'a')] * r[('a', 'm')] + r[('p', 'b')] * r[('b', 'm')], 'n': r[('p', 'c')] * r[('c', 'n')]}
    return step1, step2


R_BAR = {e: mean_r(e) for e in EDGES}


def weighted(weights, x, mode):
    if mode == 'product':
        return sum(w * x[m] for m, w in weights.items())
    if mode == 'unweighted':
        return sum(x[m] for m, w in weights.items() if w > 0)
    total = sum(weights.values())
    return sum(w * x[m] for m, w in weights.items()) / total


@pytest.mark.parametrize('mode', ['product', 'unweighted', 'normalised'])
def test_combined_feature_by_hand(mode):
    r = row(tiny_graph(), 'combined', mode)
    a1, a2 = path_weights(R_A)
    b1, b2 = path_weights(R_B)
    assert r['K27ac_direct'] == pytest.approx(log2(weighted(b1, X_B, mode), weighted(a1, X_A, mode)))
    assert r['K27ac_indirect'] == pytest.approx(log2(weighted(b2, X_B, mode), weighted(a2, X_A, mode)))


@pytest.mark.parametrize('mode', ['product', 'unweighted', 'normalised'])
def test_decomposed_terms_by_hand(mode):
    r = row(tiny_graph(), 'decomposed', mode)
    _, bar2 = path_weights(R_BAR)
    _, a2 = path_weights(R_A)
    _, b2 = path_weights(R_B)
    x_mean = {n: (X_A[n] + X_B[n]) / 2 for n in X_A}
    mark = log2(weighted(bar2, X_B, mode), weighted(bar2, X_A, mode))
    assert r['K27ac_indirect'] == pytest.approx(mark)              # main column: the mark term
    assert r['K27ac_indirect_marks'] == pytest.approx(mark)
    assert r['K27ac_indirect_contacts'] == pytest.approx(log2(weighted(b2, x_mean, mode), weighted(a2, x_mean, mode)))
    assert r['K27ac_indirect_combined'] == pytest.approx(log2(weighted(b2, X_B, mode), weighted(a2, X_A, mode)))


def test_explicit_values():
    # step 2, mark term: r_bar paths to m = 1.5 * 2 + 2 * 1 = 5, to n = 2.5 * 1.25 = 3.125
    unweighted = row(tiny_graph(), 'decomposed', 'unweighted')
    assert unweighted['K27ac_indirect'] == pytest.approx(np.log2((30 + 8 + 1) / (10 + 2 + 1)))
    normalised = row(tiny_graph(), 'decomposed', 'normalised')
    mean_b = (5 * 30 + 3.125 * 8) / 8.125
    mean_a = (5 * 10 + 3.125 * 2) / 8.125
    assert normalised['K27ac_indirect'] == pytest.approx(np.log2((mean_b + 1) / (mean_a + 1)))
    product = row(tiny_graph(), 'decomposed', 'product')
    assert product['K27ac_indirect'] == pytest.approx(np.log2((5 * 30 + 3.125 * 8 + 1) / (5 * 10 + 3.125 * 2 + 1)))


def test_default_is_product():
    graph = tiny_graph()
    default = ActivityOperators(graph, ['p'], max_hops=2, step_feature='decomposed')
    product = ActivityOperators(graph, ['p'], max_hops=2, step_feature='decomposed', step_weighting='product')
    assert default.step_weighting == 'product'
    a = gene_features(graph, {'K27ac': (X_A, X_B)}, {'p': 3.0}, operators=default)
    b = gene_features(graph, {'K27ac': (X_A, X_B)}, {'p': 3.0}, operators=product)
    assert a.equals(b)


def test_unweighted_ignores_strengths_and_multiplicity():
    # any positive strengths give the same unweighted features; the two paths to m count once
    other_a = {e: 7.0 for e in EDGES}
    other_b = {e: 0.1 for e in EDGES}
    first = row(tiny_graph(), 'combined', 'unweighted')
    second = row(tiny_graph(other_a, other_b), 'combined', 'unweighted')
    for level in ('direct', 'indirect'):
        assert first['K27ac_' + level] == pytest.approx(second['K27ac_' + level])
    # the contact term vanishes when the same nodes are reached in both conditions
    assert row(tiny_graph(), 'decomposed', 'unweighted')['K27ac_indirect_contacts'] == pytest.approx(0.0)


def test_normalised_is_scale_free_and_a_mean():
    # scaling every strength by a constant leaves the normalised features unchanged
    scaled_a = {e: 10 * v for e, v in R_A.items()}
    scaled_b = {e: 10 * v for e, v in R_B.items()}
    base = row(tiny_graph(), 'combined', 'normalised')
    scaled = row(tiny_graph(scaled_a, scaled_b), 'combined', 'normalised')
    for level in ('direct', 'indirect'):
        assert base['K27ac_' + level] == pytest.approx(scaled['K27ac_' + level])
    # a constant signal gives the constant itself (a weighted mean) whatever the number of nodes
    operators = ActivityOperators(tiny_graph(), ['p'], max_hops=2, step_weighting='normalised')
    constant = np.full(len(operators.nodes), 5.0)
    for level in ('direct', 'indirect'):
        H_a, H_b, _ = operators.level(level)
        assert (H_a @ constant)[0] == pytest.approx(5.0)
        assert (H_b @ constant)[0] == pytest.approx(5.0)


def test_contact_columns_and_node_counts_unchanged():
    rows = {mode: row(tiny_graph(), 'decomposed', mode) for mode in ('product', 'unweighted', 'normalised')}
    for mode in ('unweighted', 'normalised'):
        for column in ('contact_direct', 'contact_indirect', 'n_direct', 'n_indirect', 'K27ac_promoter'):
            assert rows[mode][column] == pytest.approx(rows['product'][column])
    assert rows['product']['n_direct'] == 3 and rows['product']['n_indirect'] == 2


def test_weight_paths():
    import scipy.sparse as sp
    H = sp.csr_matrix(np.array([[2.0, 0.0, 6.0], [0.0, 0.0, 0.0]]))
    assert weight_paths(H, 'product') is H
    assert np.allclose(weight_paths(H, 'unweighted').toarray(), [[1, 0, 1], [0, 0, 0]])
    assert np.allclose(weight_paths(H, 'normalised').toarray(), [[0.25, 0, 0.75], [0, 0, 0]])
    with pytest.raises(ValueError):
        weight_paths(H, 'mean')
    with pytest.raises(ValueError):
        ActivityOperators(tiny_graph(), ['p'], max_hops=2, step_weighting='mean')


@pytest.mark.parametrize('mode', ['unweighted', 'normalised'])
def test_hop_test_uses_the_weighting(mode):
    graph = tiny_graph()
    operators = ActivityOperators(graph, ['p'], max_hops=2, step_feature='decomposed', step_weighting=mode)
    vectors = {'K27ac': (operators.vector(X_A), operators.vector(X_B))}
    table = gene_features(graph, {'K27ac': (X_A, X_B)}, {'p': 3.0}, operators=operators, mark_vectors=vectors)
    test = gene_level_hop_test(operators, vectors, table, ['K27ac'], {'K27ac': 1}, n_permutations=20, seed=0,
                               genes={'indirect': None})['indirect'].set_index('node')
    assert test.loc['p', 'indirect_statistic'] == pytest.approx(table.set_index('node').loc['p', 'K27ac_indirect'])


def test_cli_option():
    from epigeneticgraphanalyzer import egaDifferential
    base = ['--interactionsA', 'a', '--interactionsB', 'b', '--tss', 't', '--tracksA', 'K=a.bw', '--tracksB',
            'K=b.bw', '--expression', 'e.tsv', '-o', 'out']
    assert egaDifferential.parse_arguments(base).stepWeighting == 'product'
    assert egaDifferential.parse_arguments(base + ['--stepWeighting', 'normalised']).stepWeighting == 'normalised'
    with pytest.raises(SystemExit):
        egaDifferential.parse_arguments(base + ['--stepWeighting', 'mean'])
