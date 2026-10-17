"""In-silico perturbation of explaining nodes (design, egaDesign) on toy graphs."""

import json
import os

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer import egaDesign
from epigeneticgraphanalyzer.design import (GraphArrays, gene_design, level_nulls, level_signs, neighbourhood,
                                            null_p, rank_candidates, single_condition_contributions)
from epigeneticgraphanalyzer.differential import ActivityOperators

SIGNS = 'H3K27ac:+1;H3K27me3:-1'
LEVELS = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state', 'step2', 'step3']


def add(graph, name, chrom='chr1', start=0, a=None, b=None, bait=False, genes=''):
    a = a or {}
    b = b or {}
    graph.add_node(name, chrom=chrom, start=start, end=start + 4000, bait=bait, genes=genes,
                   H3K27ac_A=a.get('H3K27ac', 2.0), H3K27ac_B=b.get('H3K27ac', a.get('H3K27ac', 2.0)),
                   H3K27me3_A=a.get('H3K27me3', 1.0), H3K27me3_B=b.get('H3K27me3', a.get('H3K27me3', 1.0)))


def node_id(i):
    return 'chr1:{}-{}'.format(i * 10000, i * 10000 + 4000)


def star_graph():
    """Gene G with four step-1 partners; only E1 gains H3K27ac (2 -> 40). A second step behind E2."""
    g = nx.Graph()
    add(g, node_id(0), start=0, bait=True, genes='G')
    for i in range(1, 5):
        add(g, node_id(i), start=i * 10000, a={'H3K27ac': 2.0}, b={'H3K27ac': 40.0 if i == 1 else 2.0})
        g.add_edge(node_id(0), node_id(i), r_a=1.5, r_b=1.5, status='shared')
    add(g, node_id(5), start=50000)
    g.add_edge(node_id(2), node_id(5), r_a=1.0, r_b=1.0, status='shared')
    return g


def chain_graph():
    """G - M - E; only E (two steps away) gains H3K27ac."""
    g = nx.Graph()
    add(g, node_id(0), start=0, bait=True, genes='G')
    add(g, node_id(1), start=10000)
    add(g, node_id(2), start=20000, a={'H3K27ac': 1.0}, b={'H3K27ac': 60.0})
    g.add_edge(node_id(0), node_id(1), r_a=2.0, r_b=2.0, status='shared')
    g.add_edge(node_id(1), node_id(2), r_a=2.0, r_b=2.0, status='shared')
    return g


def levels_table(gene, explaining, concordant='H3K27ac', rna=2.0, max_hops=3, genes='G'):
    rows = []
    for rank, level in enumerate(LEVELS[:5] + ['step{}'.format(k) for k in range(2, max_hops + 1)], start=1):
        rows.append(dict(node=gene, genes=genes, rna_log2FC=rna, rank=rank, level=level,
                         signs_declared='contact:+1' if level == 'promoter_contacts' else
                         (SIGNS if 'state' not in level else np.nan),
                         concordant_marks_declared=concordant if level in explaining else '',
                         explains_declared_mean=level in explaining, detail_declared_mean='',
                         state_fragment_declared_mean='', p_declared=0.01, q_declared=0.05,
                         beyond_mirror_declared=True, proximity_support=0.6, proximity_support_ci='0.55-0.65'))
    return pd.DataFrame(rows)


def arrays_of(graph):
    sources = [n for n, d in graph.nodes(data=True) if d.get('bait')]
    return GraphArrays(graph, {'H3K27ac': ('H3K27ac_A', 'H3K27ac_B'), 'H3K27me3': ('H3K27me3_A', 'H3K27me3_B')},
                       sources=sources)


def test_neighbourhood_matches_activity_operators():
    rng = np.random.default_rng(3)
    g = nx.gnm_random_graph(60, 140, seed=4)
    g = nx.relabel_nodes(g, {i: node_id(i) for i in g.nodes})
    for i, n in enumerate(g.nodes):
        add(g, n, start=i * 10000, bait=i % 7 == 0)
    for u, v in g.edges:
        g.edges[u, v].update(r_a=float(rng.uniform(0, 3)), r_b=float(rng.uniform(0, 3)))
    sources = [n for n, d in g.nodes(data=True) if d['bait']]
    ops = ActivityOperators(g, sources, max_hops=3)
    arrays = arrays_of(g)
    for row, s in enumerate(ops.sources):
        hops = neighbourhood(arrays, arrays.index[s], 3)
        for k in (1, 2, 3):
            expected_a = ops.hop_a[k].getrow(row).toarray().ravel()
            expected_b = ops.hop_b[k].getrow(row).toarray().ravel()
            got_a = np.zeros(len(ops.nodes))
            got_b = np.zeros(len(ops.nodes))
            cols, w_a, w_b = hops[k]
            for c, a, b in zip(cols, w_a, w_b):
                got_a[ops.index[arrays.nodes[c]]] = a
                got_b[ops.index[arrays.nodes[c]]] = b
            assert np.allclose(got_a, expected_a)
            assert np.allclose(got_b, expected_b)


def test_removing_single_explaining_node_removes_explanation():
    g = star_graph()
    arrays = arrays_of(g)
    levels = levels_table(node_id(0), {'step1_marks'})
    signs = level_signs(levels)
    candidates, rows, hops = gene_design(arrays, levels.set_index('level'), 3, signs, threshold=1.0)
    before = {r['level']: r for r in rows}
    assert before['step1_marks']['statistic'] >= 1.0 and before['step1_marks']['explains_recomputed']
    assert [c['node'] for c in candidates] == [node_id(1)]
    c = candidates[0]
    assert c['level_lost_flatten'] and c['level_lost_silence'] and c['level_lost_remove']
    assert c['statistic_flatten'] < 1e-9 < c['statistic_before']
    assert 1 in hops
    ranked = rank_candidates(pd.DataFrame(candidates).assign(**{'ev_own_change': True}))
    assert ranked.iloc[0]['node'] == node_id(1) and ranked.iloc[0]['effect_fraction'] == 1.0


def test_two_step_explanation_removed_by_target_and_by_mediator():
    g = chain_graph()
    arrays = arrays_of(g)
    levels = levels_table(node_id(0), {'step2'})
    signs = level_signs(levels)
    candidates, rows, _ = gene_design(arrays, levels.set_index('level'), 3, signs, threshold=1.0)
    assert [c['node'] for c in candidates] == [node_id(2)]
    assert candidates[0]['level_lost_flatten'] and candidates[0]['level_lost_remove']
    # deleting the mediator removes the only path of two steps
    without_mediator = neighbourhood(arrays, arrays.index[node_id(0)], 3, removed=arrays.index[node_id(1)])
    assert len(without_mediator[1][0]) == 0 and len(without_mediator[2][0]) == 0


def test_own_contacts_and_own_promoter_candidates():
    g = nx.Graph()
    add(g, node_id(0), start=0, bait=True, genes='G', a={'H3K27ac': 1.0}, b={'H3K27ac': 30.0})
    add(g, node_id(1), start=10000)
    add(g, node_id(2), start=20000)
    g.add_edge(node_id(0), node_id(1), r_a=0.0, r_b=5.0, status='gained')
    g.add_edge(node_id(0), node_id(2), r_a=1.0, r_b=1.0, status='shared')
    arrays = arrays_of(g)
    levels = levels_table(node_id(0), {'promoter_marks', 'promoter_contacts'})
    candidates, _, _ = gene_design(arrays, levels.set_index('level'), 3, level_signs(levels), threshold=1.0)
    by_level = {c['level']: c for c in candidates}
    assert by_level['promoter_contacts']['node'] == node_id(1)
    assert by_level['promoter_contacts']['level_lost_silence']
    assert by_level['promoter_marks']['own_promoter']
    ranked = rank_candidates(pd.DataFrame(candidates))
    # the gene's own promoter is listed after every other node
    assert ranked.iloc[-1]['own_promoter'] and ranked.iloc[0]['node'] == node_id(1)


def test_null_of_unpermuted_neighbourhood_gives_p_one():
    g = star_graph()
    # every step-1 partner is a source (bait with expression): nothing is permuted, P = 1
    sources = list(g.nodes)
    arrays = GraphArrays(g, {'H3K27ac': ('H3K27ac_A', 'H3K27ac_B')}, sources=sources)
    s = arrays.index[node_id(0)]
    hops = neighbourhood(arrays, s, 1)[1]
    nulls = level_nulls(arrays, {node_id(0): hops}, {node_id(0): 1.0}, {'H3K27ac': 1}, n_permutations=50)
    observed = np.log2((hops[2] @ arrays.vectors['H3K27ac'][1][hops[0]] + 1) /
                       (hops[1] @ arrays.vectors['H3K27ac'][0][hops[0]] + 1))
    assert null_p(nulls[node_id(0)], observed) == 1.0
    assert null_p(np.sort(np.arange(99.0)), 1000.0) == 0.01


def test_single_condition_ranking():
    g = star_graph()
    for e in g.edges:
        g.edges[e]['score'] = 10.0
    arrays = GraphArrays(g, {'H3K27ac': ('H3K27ac_B', 'H3K27ac_B')}, weights=('score', 'score'), reference_score=5)
    table = single_condition_contributions(arrays, [node_id(0)], 2)
    top = table[table['rank'] == 1].iloc[0]
    assert top['node'] == node_id(1) and top['hops'] == 1
    assert set(table['node']) == {node_id(i) for i in range(1, 6)}


def write_prefix(directory, graph, levels):
    prefix = os.path.join(directory, 'A_vs_B')
    nx.write_graphml(graph, prefix + '_union.graphml')
    levels.to_csv(prefix + '_explanations_all_levels.tsv', sep='\t', index=False)
    baits = [n for n, d in graph.nodes(data=True) if d.get('bait')]
    pd.DataFrame({'node': baits}).to_csv(prefix + '_gene_features.tsv', sep='\t', index=False)
    rows = []
    for n, d in graph.nodes(data=True):
        for condition in ('A', 'B'):
            rows.append(dict(node=n, condition=condition, element_type='promoter' if d.get('bait') else
                             'non-promoter element', tss_genes=d.get('genes', ''),
                             class_declared_mean='active enhancer' if n == node_id(1) else 'unmarked element'))
    pd.DataFrame(rows).to_csv(prefix + '_elements.tsv', sep='\t', index=False)
    return prefix


def test_cli_outputs(tmp_path):
    g = star_graph()
    prefix = write_prefix(str(tmp_path), g, levels_table(node_id(0), {'step1_marks'}))
    out = str(tmp_path / 'design' / 'A_vs_B')
    egaDesign.main(['-p', prefix, '--permutations', '20', '--threads', '1', '-o', out])
    design = pd.read_csv(out + '_design.tsv', sep='\t')
    first = design[design['rank'] == 1].iloc[0]
    assert first['node'] == node_id(1) and bool(first['level_lost_flatten'])
    assert 'active enhancer' in first['element_class'] and first['expected_direction'] == 'down in B'
    assert 'remove the 1-step explanation of G' in first['sentence']
    bed = pd.read_csv(out + '_design_top.bed', sep='\t', header=None)
    assert bed.iloc[0, 2] - bed.iloc[0, 1] == 1000
    summary = json.load(open(out + '_design_summary.json'))
    assert summary['level_agreement'] == 1.0 or np.isnan(summary['level_agreement'])
    with open(out + '_design.txt') as handle:
        assert 'G (' in handle.read()
