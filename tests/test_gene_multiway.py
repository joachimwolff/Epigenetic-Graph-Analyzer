"""Constellation counting of analysis/figure_gene_multiway.py on a tiny synthetic incidence matrix."""

import os
import sys

import networkx as nx
import numpy as np
import scipy.sparse as sp

from epigeneticgraphanalyzer.genome import BinIndex, IntervalIndex
from epigeneticgraphanalyzer.multiway import MultiwayData

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from figure_gene_multiway import (ChromGroups, analyse, count_constellations, empirical_p,  # noqa: E402
                                  shifted_controls, split_count, unit_groups)


def incidence(rows, n_bins, n_groups):
    """MultiwayData-style matrix from {bin: [groups]}."""
    r, c = [], []
    for b, groups in rows.items():
        r += [b] * len(groups)
        c += list(groups)
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(n_bins, n_groups))


def test_count_constellations_by_hand():
    # groups 0..5; promoter in 0-4, mediator 1 in {0, 1}, mediator 2 in {2}, enhancers e1 {0, 2, 5}, e2 {0, 3}
    promoter = np.array([0, 1, 2, 3, 4])
    mediators = [np.array([0, 1]), np.array([2])]
    enhancers = [np.array([0, 2, 5]), np.array([0, 3])]
    counts = count_constellations(promoter, mediators, enhancers)
    assert counts['promoter'] == 5
    assert counts['promoter_enhancer'] == 3            # groups 0, 2, 3
    assert counts['promoter_mediator_enhancer'] == 2   # groups 0 (m1, e1, e2) and 2 (m2, e1)
    assert counts['promoter_step1_enhancer'] == 1      # group 0
    assert counts['promoter_step2_enhancer'] == 1      # group 2
    assert counts['all_units'] == 0                    # no group holds m1 and m2
    counts = count_constellations(promoter, mediators[:1], enhancers)
    assert counts['all_units'] == 1                    # group 0 holds p, m1, e1, e2


def test_empirical_p():
    assert empirical_p(5, [1, 2, 5, 6]) == 3 / 5
    assert empirical_p(10, [1, 2]) == 1 / 3


def test_bins_units_and_shifts():
    # one chromosome of 20 bins of 10 bp, 4 groups; graph nodes on bins 2 (promoter), 3 (mediator), 5 and 5
    # (two enhancers in one bin)
    data = MultiwayData(BinIndex({'chr1': 200}, 10),
                        incidence({2: [0, 1, 2], 3: [0, 1], 5: [0, 3], 7: [1, 2], 8: [1, 2], 0: [3],
                                   12: [2, 3], 13: [2], 15: [2]}, 20, 4),
                        np.ones(4), 'gam')
    graph = nx.Graph()
    for node, (s, e) in dict(p=(20, 29), m=(30, 39), e1=(50, 52), e2=(55, 58)).items():
        graph.add_node(node, chrom='chr1', start=s, end=e)
    path_nodes = [('promoter', ['p']), ('mediator', ['m']), ('enhancer', ['e1', 'e2'])]
    counts, table, controls, detail = analyse('toy', 'GAM', data, graph, path_nodes, 'chr1', n_shifts=20, seed=0)
    assert table.loc[table.node == 'e2', 'note'].iloc[0] == '2 nodes in one bin'
    assert counts.units_enhancer.iloc[0] == 1
    observed = counts.set_index('constellation')['observed']
    assert observed['promoter_enhancer'] == 1 and observed['all_units'] == 1   # group 0
    # covered bins 0, 2, 3, 5, 7, 8, 12, 13, 15: the only shift beyond the span (3 bins) with bins 2+k, 3+k and 5+k
    # covered is k = 10 (bins 12, 13, 15), whose constellation is held by group 2
    assert controls['offset'].tolist() == [10]
    assert controls['promoter_enhancer'].tolist() == [1] and controls['all_units'].tolist() == [1]
    assert counts.set_index('constellation').loc['all_units', 'p_empirical'] == 1.0
    chrom_groups = ChromGroups(data, 'chr1')
    shifted = shifted_controls(chrom_groups, [[2], [3], [5]], ['promoter', 'mediator', 'enhancer'], 5, 0, span=0)
    assert set(shifted['offset']) == {10}       # k = 0 lies within the span and is excluded
    for _, row in shifted.iterrows():
        sets = unit_groups(chrom_groups, [[2 + row.offset], [3 + row.offset], [5 + row.offset]])
        assert split_count(sets, ['promoter', 'mediator', 'enhancer'])['promoter_enhancer'] == \
            row.promoter_enhancer


def test_fragment_units_union_of_fragments():
    # fragments of 10 bp; node p spans fragments 1 and 2, groups hold either of them
    starts = np.arange(0, 200, 10)
    bins = IntervalIndex(['chr1'] * len(starts), starts, starts + 10)
    data = MultiwayData(bins, incidence({1: [0], 2: [1], 4: [0, 1], 6: [1], 11: [2], 12: [2], 14: [2], 16: [2]},
                                        len(starts), 3), np.ones(3), 'sprite')
    graph = nx.Graph()
    for node, (s, e) in dict(p=(10, 30), m=(40, 50), e=(60, 70)).items():
        graph.add_node(node, chrom='chr1', start=s, end=e)
    path_nodes = [('promoter', ['p']), ('mediator', ['m']), ('enhancer', ['e'])]
    counts, table, controls, _ = analyse('toy', 'fragments', data, graph, path_nodes, 'chr1', n_shifts=10, seed=1)
    assert table.set_index('node')['fragments'].to_dict() == {'p': 2, 'm': 1, 'e': 1}
    observed = counts.set_index('constellation')['observed']
    assert observed['promoter_enhancer'] == 1          # group 1 holds fragment 2 (p) and 6 (e)
    assert observed['all_units'] == 1
    # covered shifts beyond the span (58 bp) put m onto fragment 14 and e onto 16 (offsets 92..108), p onto fragments
    # among 10-13: every such constellation is held by group 2 alone
    assert len(controls) == 10
    assert controls['offset'].between(92, 108).all()
    assert (controls['all_units'] == 1).all() and (controls['promoter'] == 1).all()
