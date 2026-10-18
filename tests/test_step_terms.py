"""Decomposition of the step feature into a mark term (fixed strengths r_bar) and a contact term (fixed marks x_bar)."""

import json
import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import egaDifferential, egaExplain
from epigeneticgraphanalyzer.differential import (ActivityOperators, explanation_levels, gene_features,
                                                  gene_level_hop_test, step_terms, term_driver)

from synthetic import simulate_capture, write_bigwig, write_tss_bed


def tiny_graph(r_a, r_b):
    """Promoter p - enhancer e1 - enhancer e2 (two steps), and p - e3 (one step); strengths per condition."""
    graph = nx.Graph()
    for (u, v), a in r_a.items():
        graph.add_edge(u, v, r_a=a, r_b=r_b[(u, v)], log2FC=0.0, qvalue=np.nan, status='shared')
    for i, n in enumerate(sorted(graph.nodes)):
        graph.nodes[n].update(chrom='chr1', start=100000 * (i + 1), end=100000 * (i + 1) + 4999, bait=n == 'p',
                              genes='G' + n if n == 'p' else '')
    return graph


EDGES = [('p', 'e1'), ('e1', 'e2'), ('p', 'e3')]


def features(graph, signal_a, signal_b, step_feature, rna=3.0):
    operators = ActivityOperators(graph, ['p'], max_hops=2, step_feature=step_feature)
    marks = {'K27ac': (signal_a, signal_b)}
    table = gene_features(graph, marks, {'p': rna}, operators=operators)
    return operators, table.set_index('node').loc['p']


def test_only_contacts_change():
    # every contact strengthens 4-fold in B, the marks are the same: the change is a contact change only
    r_a = dict(zip(EDGES, (1.0, 1.0, 1.0)))
    r_b = dict(zip(EDGES, (4.0, 4.0, 4.0)))
    x = {'p': 5.0, 'e1': 10.0, 'e2': 20.0, 'e3': 8.0}
    graph = tiny_graph(r_a, r_b)
    _, combined = features(graph, x, x, 'combined')
    _, row = features(graph, x, x, 'decomposed')
    for level in ('direct', 'indirect'):
        assert row['K27ac_' + level] == pytest.approx(0.0)          # main column: the mark term
        assert row['K27ac_{}_marks'.format(level)] == pytest.approx(0.0)
        assert row['K27ac_{}_contacts'.format(level)] > 1.0
        # with fixed marks the contact term equals the combined feature
        assert row['K27ac_{}_contacts'.format(level)] == pytest.approx(combined['K27ac_' + level])
        assert row['K27ac_{}_combined'.format(level)] == pytest.approx(combined['K27ac_' + level])
    # contact_<level> keeps the strengths of the two conditions
    assert row['contact_direct'] == pytest.approx(combined['contact_direct'])
    assert row['contact_direct'] == pytest.approx(np.log2((8 + 1) / (2 + 1)))
    # the explicit formula: direct = 4 * (10 + 8) vs 1 * (10 + 8)
    assert row['K27ac_direct_contacts'] == pytest.approx(np.log2((4 * 18 + 1) / (18 + 1)))
    assert row['K27ac_indirect_contacts'] == pytest.approx(np.log2((16 * 20 + 1) / (20 + 1)))
    _, contacts = features(graph, x, x, 'contacts')
    assert contacts['K27ac_indirect'] == pytest.approx(row['K27ac_indirect_contacts'])
    _, marks = features(graph, x, x, 'marks')
    assert marks['K27ac_indirect'] == pytest.approx(0.0)


def test_only_marks_change():
    # same contacts in A and B, the mark rises at the enhancers
    r = dict(zip(EDGES, (2.0, 3.0, 1.0)))
    x_a = {'p': 5.0, 'e1': 4.0, 'e2': 2.0, 'e3': 1.0}
    x_b = {'p': 5.0, 'e1': 16.0, 'e2': 32.0, 'e3': 1.0}
    graph = tiny_graph(r, dict(r))
    _, combined = features(graph, x_a, x_b, 'combined')
    _, row = features(graph, x_a, x_b, 'decomposed')
    for level in ('direct', 'indirect'):
        assert row['K27ac_{}_contacts'.format(level)] == pytest.approx(0.0)
        assert row['K27ac_{}_marks'.format(level)] == pytest.approx(combined['K27ac_' + level])
        assert row['K27ac_' + level] == pytest.approx(combined['K27ac_' + level])
        assert row['K27ac_' + level] > 1.0
    assert row['K27ac_indirect'] == pytest.approx(np.log2((6 * 32 + 1) / (6 * 2 + 1)))
    # the promoter level is the same under every step feature
    assert row['K27ac_promoter'] == pytest.approx(combined['K27ac_promoter'])


def test_both_change_and_half_strength_of_uncalled_edges():
    # p - e1 gained in B (uncalled in A), e1 - e2 shared, p - e3 lost in B; the marks rise at e2
    r_a = dict(zip(EDGES, (0.0, 2.0, 3.0)))
    r_b = dict(zip(EDGES, (4.0, 2.0, 0.0)))
    x_a = {'p': 1.0, 'e1': 2.0, 'e2': 3.0, 'e3': 6.0}
    x_b = {'p': 1.0, 'e1': 2.0, 'e2': 12.0, 'e3': 6.0}
    graph = tiny_graph(r_a, r_b)
    operators, row = features(graph, x_a, x_b, 'decomposed')
    # r_bar: a gained / lost edge enters with half its strength, in both conditions
    i = operators.sources.index('p')
    direct = operators.hop_a[1].getrow(i)
    weights = {operators.nodes[j]: v for j, v in zip(direct.indices, direct.data)}
    assert weights == pytest.approx({'e1': 2.0, 'e3': 1.5})
    assert (operators.hop_a[2] != operators.hop_b[2]).nnz == 0
    # step 2: path weight r_bar(p,e1) * r_bar(e1,e2) = 2 * 2 = 4 in both conditions
    assert row['K27ac_indirect_marks'] == pytest.approx(np.log2((4 * 12 + 1) / (4 * 3 + 1)))
    # contact term at step 2: x_bar(e2) = 7.5, path weight 0 in A, 8 in B
    assert row['K27ac_indirect_contacts'] == pytest.approx(np.log2((8 * 7.5 + 1) / (0 + 1)))
    assert row['K27ac_indirect_combined'] == pytest.approx(np.log2((8 * 12 + 1) / (0 + 1)))
    assert row['K27ac_indirect_marks'] > 1.0 and row['K27ac_indirect_contacts'] > 1.0
    # step 1: the marks of e1 and e3 do not change (mark term 0), the contacts move from e3 to e1
    assert row['K27ac_direct_marks'] == pytest.approx(0.0)
    assert row['K27ac_direct_contacts'] == pytest.approx(np.log2((4 * 2 + 1) / (3 * 6 + 1)))
    # the two terms against the formula of step_terms
    C_a, C_b, _ = operators.condition_level('indirect')
    H_bar = operators.level('indirect')[0]
    va, vb = operators.vector(x_a), operators.vector(x_b)
    terms = step_terms(C_a, C_b, H_bar, va, vb)
    assert terms['marks'][i] == pytest.approx(row['K27ac_indirect_marks'])
    assert terms['contacts'][i] == pytest.approx(row['K27ac_indirect_contacts'])


def test_term_driver():
    assert term_driver(True, True, True) == 'both'
    assert term_driver(True, True, False) == 'marks'
    assert term_driver(True, False, True) == 'contacts'
    assert term_driver(True, False, False) == 'neither'
    assert term_driver(False, False, False) == ''
    assert term_driver(False, True, False) == 'marks'


def chain_graph(n_units, gained, seed=0):
    """n_units promoters p<i> - e<i> - f<i>; in the `gained` units the contacts double in B and the marks stay, in
    the others the marks at f<i> rise 8-fold and the contacts stay; a pool of unconnected nodes holds the other
    signals."""
    rng = np.random.default_rng(seed)
    graph = nx.Graph()
    x_a, x_b = {}, {}
    for i in range(n_units):
        p, e, f = 'p{}'.format(i), 'e{}'.format(i), 'f{}'.format(i)
        r_b = 4.0 if i in gained else 1.0
        graph.add_edge(p, e, r_a=1.0, r_b=r_b, log2FC=0.0, qvalue=np.nan, status='shared')
        graph.add_edge(e, f, r_a=1.0, r_b=r_b, log2FC=0.0, qvalue=np.nan, status='shared')
        x_a[p] = x_b[p] = 1.0
        x_a[e] = x_b[e] = float(rng.uniform(1, 2))
        x_a[f] = 4.0
        x_b[f] = 4.0 if i in gained else 32.0
    for j in range(400):
        a, b = 'z{}'.format(j), 'w{}'.format(j)
        graph.add_edge(a, b, r_a=1.0, r_b=1.0, log2FC=0.0, qvalue=np.nan, status='shared')
        x_a[a] = x_b[a] = x_a[b] = x_b[b] = float(rng.uniform(1, 5))
    for i, n in enumerate(sorted(graph.nodes)):
        graph.nodes[n].update(chrom='chr1', start=10000 * (i + 1), end=10000 * (i + 1) + 999,
                              bait=n.startswith('p'), genes=n if n.startswith('p') else '')
    return graph, x_a, x_b


def test_permutation_test_on_the_mark_term():
    gained = set(range(0, 10))
    graph, x_a, x_b = chain_graph(20, gained)
    change = {'p{}'.format(i): 3.0 for i in range(20)}
    results = {}
    for feature in ('combined', 'decomposed'):
        operators = ActivityOperators(graph, list(change), max_hops=2, step_feature=feature)
        vectors = {'K27ac': (operators.vector(x_a), operators.vector(x_b))}
        table = gene_features(graph, {'K27ac': (x_a, x_b)}, change, operators=operators, mark_vectors=vectors)
        results[feature] = (table.set_index('node'), gene_level_hop_test(
            operators, vectors, table, ['K27ac'], {'K27ac': 1}, rna_threshold=1.0, n_permutations=500, seed=1,
            genes={'indirect': None})['indirect'].set_index('node'))
    combined_table, combined_test = results['combined']
    decomposed_table, decomposed_test = results['decomposed']
    contact_genes = ['p{}'.format(i) for i in gained]
    mark_genes = ['p{}'.format(i) for i in range(10, 20)]
    # combined: contact gains alone give a concordant step-2 change; the mark term is 0 for these genes
    assert (combined_table.loc[contact_genes, 'K27ac_indirect'] > 1).all()
    assert decomposed_table.loc[contact_genes, 'K27ac_indirect'].abs().max() == pytest.approx(0.0)
    assert (decomposed_table.loc[contact_genes, 'K27ac_indirect_contacts'] > 1).all()
    assert (decomposed_table.loc[mark_genes, 'K27ac_indirect'] > 1).all()
    assert decomposed_table.loc[mark_genes, 'K27ac_indirect_contacts'].abs().max() == pytest.approx(0.0)
    # the permutation test on the mark term: genes with mark gains are significant, contact-only genes are not
    assert (decomposed_test.loc[mark_genes, 'indirect_p'] < 0.05).all()
    assert (decomposed_test.loc[contact_genes, 'indirect_p'] > 0.5).all()
    assert (combined_test.loc[contact_genes, 'indirect_statistic'] > 1).all()


def test_all_levels_report_the_terms():
    r_a = dict(zip(EDGES, (1.0, 1.0, 1.0)))
    r_b = dict(zip(EDGES, (4.0, 4.0, 1.0)))
    x = {'p': 5.0, 'e1': 10.0, 'e2': 20.0, 'e3': 8.0}
    graph = tiny_graph(r_a, r_b)
    operators = ActivityOperators(graph, ['p'], max_hops=2, step_feature='decomposed')
    vectors = {'K27ac': (operators.vector(x), operators.vector(x))}
    table = gene_features(graph, {'K27ac': (x, x)}, {'p': 3.0}, operators=operators, mark_vectors=vectors)
    levels = explanation_levels(table, graph, operators, vectors, ['K27ac'], {'declared': {'K27ac': 1}},
                                support=False)
    step2 = levels[levels['level'] == 'step2'].iloc[0]
    assert step2['driver_declared'] == 'contacts'
    assert step2['concordant_contact_term_declared'] == 'K27ac'
    assert step2['concordant_mark_term_declared'] == ''
    assert step2['concordant_combined_declared'] == 'K27ac'
    # the cascade uses the mark term: no concordant change at step 2
    assert step2['concordant_marks_declared'] == ''
    assert step2['mark_term_changes'] == 'K27ac:+0.00'
    assert 'driver_declared' in levels.columns and 'contact_term_changes' in levels.columns
    # without the decomposition the columns are absent
    plain = ActivityOperators(graph, ['p'], max_hops=2)
    plain_table = gene_features(graph, {'K27ac': (x, x)}, {'p': 3.0}, operators=plain)
    plain_levels = explanation_levels(plain_table, graph, plain, vectors, ['K27ac'], {'declared': {'K27ac': 1}},
                                      support=False)
    assert 'driver_declared' not in plain_levels.columns
    assert plain_levels.loc[plain_levels['level'] == 'step2', 'concordant_marks_declared'].iloc[0] == 'K27ac'
    assert isinstance(plain_levels, pd.DataFrame)


def test_unknown_step_feature():
    graph = tiny_graph(dict(zip(EDGES, (1.0, 1.0, 1.0))), dict(zip(EDGES, (1.0, 1.0, 1.0))))
    with pytest.raises(ValueError):
        ActivityOperators(graph, ['p'], max_hops=2, step_feature='both')


def test_differential_cli_decomposed(tmp_path):
    d = str(tmp_path)
    da, db = os.path.join(d, 'a'), os.path.join(d, 'b')
    os.makedirs(da)
    os.makedirs(db)
    planted = sorted(simulate_capture(da, seed=7)[1])
    gained = {(c, i, j + 3) for c, i, j in planted[:20]}
    ibed_a = simulate_capture(da, seed=7, planted=planted, n_noise=0, name='a7')[0]
    ibed_b, _, baits = simulate_capture(db, seed=8, planted=set(planted[20:]) | gained, n_noise=0)
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    rng = np.random.default_rng(0)
    bw_a = write_bigwig(os.path.join(d, 'a.bw'), {(c, j): 2.0 for c, _, j in planted})
    bw_b = write_bigwig(os.path.join(d, 'b.bw'), {(c, j): 4.0 for c, _, j in gained})
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)), 'B': rng.uniform(5, 50, len(baits))})
    expression_path = os.path.join(d, 'expression.tsv')
    expression.to_csv(expression_path, sep='\t', index=False)
    prefix = os.path.join(d, 'out', 'diff')
    egaDifferential.main(['--interactionsA', ibed_a, '--interactionsB', ibed_b, '-f', 'ibed', '--names', 'A', 'B',
                          '--tss', tss, '--tracksA', 'H3K27ac=' + bw_a, '--tracksB', 'H3K27ac=' + bw_b,
                          '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                          '--permutations', '3', '--bootstrap', '10', '--genePermutations', '50', '--roleBootstrap',
                          '10', '--threads', '1', '--roles', 'both', '--stepFeature', 'decomposed', '-o', prefix])
    features = pd.read_csv(prefix + '_gene_features.tsv', sep='\t')
    for term in ('combined', 'marks', 'contacts'):
        assert 'H3K27ac_indirect_' + term in features
    assert np.allclose(features['H3K27ac_indirect'], features['H3K27ac_indirect_marks'])
    assert 'H3K27ac_direct_combined_mirror' in features and 'H3K27ac_direct_contacts_mirror' in features
    with open(prefix + '_permutation.json') as handle:
        assert json.load(handle)['step_feature'] == 'decomposed'
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', keep_default_na=False)
    for column in ('mark_term_changes', 'contact_term_changes', 'combined_changes', 'driver_declared',
                   'driver_learned'):
        assert column in levels
    steps = levels[levels['level'].str.startswith('step') & ~levels['level'].str.endswith('state')]
    assert set(steps['driver_declared']) <= {'', 'marks', 'contacts', 'both', 'neither'}
    explain = os.path.join(d, 'out', 'explain')
    egaExplain.main(['-p', prefix, '-o', explain])
    genes = pd.read_csv(explain + '_genes.tsv', sep='\t')
    assert {'step_terms', 'first_step_driver'} <= set(genes.columns)
