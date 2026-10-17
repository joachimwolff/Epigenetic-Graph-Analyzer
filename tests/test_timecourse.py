"""Time course of the explanation cascade (timecourse, egaTimecourse) on toy trajectories."""

import itertools

import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer.timecourse import (classify_order, day_operators, event_patterns, expression_pattern,
                                                find_viewpoints, first_interval, knockdown_test, node_gene_choice,
                                                order_test, profile_at, shuffle_within_groups, timecourse_graph,
                                                viewpoint_profile)


def test_first_interval_threshold_and_significance():
    changes = np.array([[0.2, 1.5], [1.2, 2.0], [-1.1, 0.0], [0.5, -0.5], [np.nan, 1.0]])
    assert first_interval(changes, 1.0).tolist() == [1, 0, 0, -1, 1]
    significant = np.array([[True, True], [False, True], [True, True], [True, True], [True, False]])
    assert first_interval(changes, 1.0, significant).tolist() == [1, 1, 0, -1, -1]


def test_classify_order_two_intervals():
    # all 8 patterns of three layers over two intervals
    expected = {(0, 0, 0): 'simultaneous', (1, 1, 1): 'simultaneous', (0, 0, 1): 'forward', (0, 1, 1): 'forward',
                (1, 1, 0): 'reversed', (1, 0, 0): 'reversed', (0, 1, 0): 'other', (1, 0, 1): 'other'}
    names = ('forward', 'simultaneous', 'reversed', 'other')
    times = np.array(list(itertools.product([0, 1], repeat=3)))
    codes = classify_order(times)
    for t, code in zip(times, codes):
        assert names[code] == expected[tuple(t)]
    assert classify_order(np.array([[0, 1], [1, 1], [1, 0]])).tolist() == [0, 1, 2]
    # three intervals: strict cascade is forward
    assert classify_order(np.array([[0, 1, 2], [2, 1, 0], [0, 2, 1]])).tolist() == [0, 2, 3]


def test_shuffle_within_groups_keeps_labels_of_every_group():
    rng = np.random.default_rng(1)
    groups = np.repeat(np.arange(50), 4)
    labels = rng.integers(0, 2, len(groups))
    shuffled = shuffle_within_groups(labels, groups, rng)
    for g in range(50):
        assert sorted(shuffled[groups == g]) == sorted(labels[groups == g])
    assert not np.array_equal(shuffled, labels)


def toy_units(n_genes, pattern_of_gene):
    rows = []
    for g in range(n_genes):
        for e in range(2):
            c, m, x = pattern_of_gene(g)
            rows.append({'gene_node': 'g{}'.format(g), 'node': 'e{}_{}'.format(g, e), 't_contact': c, 't_mark': m,
                         't_expression': x})
    return pd.DataFrame(rows)


def test_order_test_detects_chromatin_before_expression():
    # every gene: contact and mark change in the first interval, expression in the second
    units = toy_units(60, lambda g: (0, 0, 1))
    res = order_test(units, ('contact', 'mark', 'expression'), n_permutations=500, n_bootstrap=200, seed=0,
                     gene_column='gene_node').set_index('category')
    assert res.loc['forward', 'observed'] == 120
    assert res.loc['forward', 'null_mean'] < 120
    assert res.loc['forward', 'p'] < 0.01
    assert res.loc['reversed', 'observed'] == 0
    assert res.loc['forward', 'boot_low'] == res.loc['forward', 'boot_high'] == 120


def test_order_test_exchangeable_layers_not_significant():
    rng = np.random.default_rng(3)
    patterns = {g: tuple(rng.integers(0, 2, 3)) for g in range(80)}
    rows = []
    for g in range(80):
        expression = int(rng.integers(0, 2))
        for e in range(3):
            c, m = rng.integers(0, 2, 2)
            rows.append({'gene_node': g, 't_contact': c, 't_mark': m, 't_expression': expression})
    units = pd.DataFrame(rows)
    res = order_test(units, ('contact', 'mark', 'expression'), n_permutations=500, n_bootstrap=100, seed=1,
                     gene_column='gene_node')
    assert res['observed'].sum() == len(units)
    assert (res['p'] > 0.01).all()
    assert patterns


def test_order_test_ignores_units_without_events_and_two_layers():
    units = pd.DataFrame({'gene_node': ['a', 'a', 'b', 'b'], 't_contact': [0, -1, 1, 0],
                          't_expression': [1, 1, 0, 0]})
    res = order_test(units, ('contact', 'expression'), n_permutations=50, n_bootstrap=20,
                     gene_column='gene_node').set_index('category')
    assert res['units'].iloc[0] == 3
    assert res.loc['before', 'observed'] == 1 and res.loc['after', 'observed'] == 1
    assert res.loc['same', 'observed'] == 1
    patterns = event_patterns(units, ('contact', 'expression'))
    assert patterns['units'].sum() == 3


def test_timecourse_graph_and_day_operators():
    def table(rows):
        return pd.DataFrame(rows, columns=['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end', 'score'])
    tables = {
        'D0': table([('chr1', 0, 10, 100, 110, 10.0), ('chr1', 100, 110, 200, 210, 5.0)]),
        'D3': table([('chr1', 0, 10, 100, 110, 20.0)]),
        'D6': table([('chr1', 0, 10, 100, 110, 5.0), ('chr1', 100, 110, 200, 210, 15.0)]),
    }
    graph = timecourse_graph(tables, reference_score=5.0)
    assert graph.number_of_edges() == 2
    assert graph.edges['chr1:100-110', 'chr1:200-210']['r_D3'] == 0.0
    ops, hops = day_operators(graph, ['chr1:0-10'], ['D0', 'D3', 'D6'], max_hops=2)
    i = ops.sources.index('chr1:0-10')
    mid, far = ops.index['chr1:100-110'], ops.index['chr1:200-210']
    assert hops['D0'][1][i, mid] == pytest.approx(2.0)
    assert hops['D3'][1][i, mid] == pytest.approx(4.0)
    assert hops['D6'][1][i, mid] == pytest.approx(1.0)
    # two steps: product of the edge strengths of the day
    assert hops['D0'][2][i, far] == pytest.approx(2.0)
    assert hops['D3'][2][i, far] == pytest.approx(0.0)
    assert hops['D6'][2][i, far] == pytest.approx(3.0)
    assert 'r_a' not in graph.edges['chr1:0-10', 'chr1:100-110']


def test_expression_pattern_and_gene_choice():
    lfc = np.array([[1.5, 0.2], [-2.0, 1.2], [0.3, -0.1]])
    padj = np.array([[0.01, 0.5], [0.001, 0.01], [0.01, 0.01]])
    assert expression_pattern(lfc, padj) == ['+0', '-+', '00']
    tables = {('D0', 'D3'): pd.DataFrame({'log2FoldChange': [1.5, 3.0, 4.0], 'padj': [0.01, 0.001, 0.5]},
                                         index=['A', 'B', 'C'])}
    choice = node_gene_choice({'n1': 'A,B', 'n2': 'C,A', 'n3': 'X'}, tables)
    assert choice == {'n1': 'B', 'n2': 'A'}


def test_knockdown_test_shift():
    rng = np.random.default_rng(0)
    cand = pd.DataFrame({'gene_node': np.repeat(['c{}'.format(i) for i in range(30)], 3),
                         'node': ['cn{}'.format(i) for i in range(90)]})
    ctrl = pd.DataFrame({'gene_node': np.repeat(['k{}'.format(i) for i in range(100)], 3),
                         'node': ['kn{}'.format(i) for i in range(300)]})
    values = {n: rng.normal(-1.0, 0.3) for n in cand['node']}
    values.update({n: rng.normal(0.0, 0.3) for n in ctrl['node']})
    res = knockdown_test(values, cand, ctrl, n_bootstrap=300, seed=0)
    assert res['nodes_candidate'] == 90 and res['nodes_control'] == 300
    assert res['median_difference'] == pytest.approx(-1.0, abs=0.2)
    assert res['md_high'] < 0
    assert res['cliffs_delta'] < -0.8 and res['delta_high'] < 0
    assert res['p_less'] < 1e-10
    # a node shared with a candidate gene leaves the control group
    shared = pd.concat([ctrl, pd.DataFrame({'gene_node': ['k0'], 'node': ['cn0']})])
    assert knockdown_test(values, cand, shared, n_bootstrap=10)['nodes_control'] == 300


def test_umi4c_viewpoints_and_profile():
    rows = []
    for i in range(1, 40):
        rows.append(('chr2', 100000, 105000, 'chr2', 100000 + i * 5000, 105000 + i * 5000, 50 if i == 3 else 5))
    rows.append(('chr2', 300000, 305000, 'chr3', 10000, 15000, 1))
    rows.append(('chr5', 0, 5000, 'chr5', 5000, 10000, 1))
    table = pd.DataFrame(rows, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count'])
    from epigeneticgraphanalyzer.timecourse import bin_totals
    viewpoints = find_viewpoints(bin_totals(table), min_fraction=0.05)
    assert viewpoints == [('chr2', 100000, 105000)]
    profile = viewpoint_profile(table, viewpoints[0], window=1000000)
    assert profile.sum() == 50 + 38 * 5
    assert profile_at(profile, 114000, 116000) == 50 + 5
    from epigeneticgraphanalyzer.timecourse import profile_over
    # overlapping intervals count each bin once
    assert profile_over(profile, [(114000, 116000), (115000, 115500)]) == 50 + 5


def test_layer_marginal_null_absorbs_a_global_layer_shift():
    # the mark changes in the first interval everywhere, expression in either interval at random: the within-gene
    # null reports mark before expression, the layer-marginal null (same per-layer interval frequencies) does not
    rng = np.random.default_rng(5)
    units = pd.DataFrame([{'gene_node': g, 't_mark': 0, 't_expression': int(rng.integers(0, 2))} for g in range(200)])
    res = order_test(units, ('mark', 'expression'), n_permutations=500, n_bootstrap=50, seed=2,
                     gene_column='gene_node').set_index('category')
    assert res.loc['before', 'p'] < 0.01
    assert res.loc['before', 'marginal_p'] > 0.05
    assert res.loc['before', 'marginal_ratio'] == pytest.approx(1.0)
