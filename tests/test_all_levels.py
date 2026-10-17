"""Every level of the explanation cascade reported for every differential gene (explanation_levels) on toy data."""

import networkx as nx
import numpy as np

from epigeneticgraphanalyzer.differential import (ActivityOperators, LevelSigns, cascade_first_level, classify_genes,
                                                  explanation_cascade, explanation_combinations,
                                                  explanation_level_counts, explanation_level_names,
                                                  explanation_level_summary, explanation_levels, gene_features,
                                                  gene_level_hop_test)
from epigeneticgraphanalyzer.egaExplain import all_level_text
from epigeneticgraphanalyzer.states import gene_states


def toy():
    """Gene g (up) with a promoter bivalent in both conditions (fragment mean) whose H3K27me3 sits in a narrow peak
    only; its partner m1 carries no change, the node m2 two steps away gains the active mark P. Gene h (down) has no
    change anywhere."""
    graph = nx.Graph()
    g, m1, m2, h, k1 = 'chr1:0-1000', 'chr1:2000-3000', 'chr1:4000-5000', 'chr2:0-1000', 'chr2:2000-3000'
    for u, v in ((g, m1), (m1, m2), (h, k1)):
        graph.add_edge(u, v, r_a=1.0, r_b=1.0, status='shared', log2FC=0.0, qvalue=np.nan)
    for n in graph.nodes:
        chrom, rest = n.split(':')
        start, end = rest.split('-')
        graph.nodes[n].update(chrom=chrom, start=int(start), end=int(end), bait=n in (g, h),
                              genes={g: 'G', h: 'H'}.get(n, ''))
    # mark signals (changes): P gained at m2 only
    p_a = {g: 4.0, m1: 2.0, m2: 1.0, h: 3.0, k1: 2.0}
    p_b = dict(p_a, **{m2: 15.0})
    r_a = {n: 2.0 for n in graph.nodes}
    marks = {'P': (p_a, p_b), 'R': (r_a, dict(r_a))}
    expression = {g: 3.0, h: -2.0}
    # presence statistics of the states: mean (R present at g in both) and max1kb (R absent)
    mean = {'P': ({g: 4.0, m1: 0.0, m2: 0.0, h: 0.0, k1: 0.0}, {g: 4.0, m1: 0.0, m2: 3.0, h: 0.0, k1: 0.0}),
            'R': ({g: 3.0, m1: 0.0, m2: 0.0, h: 0.0, k1: 0.0}, {g: 3.0, m1: 0.0, m2: 0.0, h: 0.0, k1: 0.0})}
    max1kb = {'P': mean['P'], 'R': ({n: 0.0 for n in graph.nodes}, {n: 0.0 for n in graph.nodes})}
    return graph, marks, expression, {'mean': mean, 'max1kb': max1kb}


def test_promoter_state_and_step2_both_reported_in_rank_order():
    graph, marks, expression, presence = toy()
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    signs = {'P': 1, 'R': -1}
    classes = classify_genes(features, ['P', 'R'], signs)
    roles = {'P': 'active', 'R': 'repressive'}
    states = {('declared', s): gene_states(classes, graph, presence[s], roles, 2.0, names=('A', 'B'))
              for s in presence}
    tests = {'declared': gene_level_hop_test(operators, vectors, features, ['P', 'R'], signs, n_permutations=99,
                                             genes={'direct': None, 'indirect': None})}
    levels = explanation_levels(classes, graph, operators, vectors, ['P', 'R'], {'declared': signs}, states=states,
                                tests=tests, names=('A', 'B'))
    assert list(levels.loc[levels['node'] == 'chr1:0-1000', 'level']) == explanation_level_names(2)
    assert list(levels['rank'].unique()) == [1, 2, 3, 4, 5, 6]
    gene = levels[levels['node'] == 'chr1:0-1000'].set_index('level')
    # promoter state (poised in both under the fragment mean) and step 2 (P gained two steps away): both reported
    assert gene.loc['promoter_state', 'explains_declared_mean']
    assert gene.loc['promoter_state', 'result_declared_mean'] == 'poised in both'
    assert 'bivalent in A and B (P 4.0/4.0, R 3.0/3.0; mean)' in gene.loc['promoter_state', 'detail_declared_mean']
    assert gene.loc['step2', 'explains_declared_mean'] and gene.loc['step2', 'explains_declared_max1kb']
    assert gene.loc['step2', 'concordant_marks_declared'] == 'P'
    assert gene.loc['step2', 'strongest_node_declared'] == 'chr1:4000-5000'
    assert gene.loc['step2', 'n_concordant_nodes_declared'] == 1
    assert 0 < gene.loc['step2', 'p_declared'] <= 1
    assert not gene.loc['promoter_marks', 'explains_declared_mean']
    assert not gene.loc['step1_marks', 'explains_declared_mean']
    # under the max1kb statistic the promoter is active in both: the state level does not explain the gene
    assert gene.loc['promoter_state', 'result_declared_max1kb'] == 'unchanged'
    summary = explanation_level_summary(levels)
    assert summary.loc['chr1:0-1000', 'explaining_levels_declared_mean'] == 'promoter_state;step2'
    assert summary.loc['chr1:0-1000', 'first_level_declared_mean'] == 'promoter_state'
    assert summary.loc['chr1:0-1000', 'n_levels_explaining_declared_mean'] == 2
    assert summary.loc['chr1:0-1000', 'explaining_levels_declared_max1kb'] == 'step2'
    assert summary.loc['chr2:0-1000', 'first_level_declared_mean'] == 'none'
    # the first level is the level of the cascade
    cascade = explanation_cascade(classes, ['P', 'R'], signs, states=states[('declared', 'mean')])
    assert dict(cascade_first_level(cascade)) == dict(summary['first_level_declared_mean'])
    counts = explanation_level_counts(levels).set_index(['presence', 'level'])
    assert counts.loc[('mean', 'step2'), 'genes_any'] == 1 and counts.loc[('mean', 'step2'), 'genes_first'] == 0
    assert counts.loc[('max1kb', 'step2'), 'genes_first'] == 1
    assert counts.loc[('mean', 'promoter_state'), 'any_poised_in_both'] == 1
    assert counts.loc[('mean', 'none'), 'genes_any'] == 1
    combos = explanation_combinations(levels).set_index(['presence', 'combination'])
    assert combos.loc[('mean', 'promoter_state;step2'), 'genes'] == 1
    assert combos.loc[('mean', 'none'), 'genes'] == 1
    # per-gene text: every level in rank order
    text = all_level_text('G', 3.0, levels[levels['node'] == 'chr1:0-1000'], 'declared')
    order = [text.index(label) for label in ('promoter marks:', 'contacts:', 'promoter state:', 'step 1 marks:',
                                             'step 1 states:', 'step 2:')]
    assert order == sorted(order)
    assert text.startswith('G (up, log2 3.0): promoter marks: no concordant change')
    assert 'explaining levels: promoter_state;step2 (mean), step2 (max1kb)' in text


def test_learned_roles_reported_beside_declared():
    graph, marks, expression, presence = toy()
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    classes = classify_genes(features, ['P', 'R'], {'P': 1, 'R': -1})
    # learned: P without a consistent effect at step 2, so step 2 does not explain under the learned roles
    learned = LevelSigns({'promoter': {'P': 1, 'R': -1}, 'direct': {'P': 1}, 'indirect': {'P': 0}})
    levels = explanation_levels(classes, graph, operators, vectors, ['P', 'R'],
                                {'declared': {'P': 1, 'R': -1}, 'learned': learned})
    gene = levels[levels['node'] == 'chr1:0-1000'].set_index('level')
    assert gene.loc['step2', 'explains_declared_mean'] and not gene.loc['step2', 'explains_learned_mean']
    assert gene.loc['step2', 'result_learned_mean'] == 'no consistent role'
    assert gene.loc['promoter_state', 'result_declared_mean'] == 'not callable'
