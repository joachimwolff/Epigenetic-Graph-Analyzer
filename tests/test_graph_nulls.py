"""Pair sets, degree classes and the distance-preserving rewiring of graph_nulls."""

import networkx as nx
import numpy as np

from epigeneticgraphanalyzer.graph_nulls import (_length_bin, degree_bin_classes, graph_distance_pairs,
                                                 node_midpoints, rewire_by_distance)


def _node(graph, chrom, position, bait):
    name = '{}:{}-{}'.format(chrom, position, position + 999)
    graph.add_node(name, chrom=chrom, start=position, end=position + 999, bait=bait)
    return name


def random_capture_graph(seed=0, n_baits=200, per_bait=6):
    rng = np.random.default_rng(seed)
    graph = nx.Graph()
    for k in range(n_baits):
        chrom = 'chr1' if k % 2 else 'chr2'
        position = int(rng.integers(1, 500)) * 10000
        bait = _node(graph, chrom, position, True)
        for _ in range(int(rng.integers(1, per_bait + 1))):
            other = _node(graph, chrom, max(0, position + int(rng.choice([-1, 1]) *
                                                                 10 ** rng.uniform(4, 5))) // 1000 * 1000, False)
            if other != bait:
                graph.add_edge(bait, other)
    return graph


def test_graph_distance_pairs_path():
    graph = nx.Graph()
    names = [_node(graph, 'chr1', 100000 * (i + 1), i % 2 == 0) for i in range(5)]
    graph.add_edges_from(zip(names[:-1], names[1:]))
    other = _node(graph, 'chr2', 100000, True)
    graph.add_edge(names[0], other)                     # trans edge: its pairs are not cis
    pairs = graph_distance_pairs(graph, max_distance=250000, max_steps=5)
    assert sorted(pairs['steps'].tolist()) == [1, 1, 1, 1, 2, 2, 2]
    assert (pairs['chrom'] == 'chr1').all() and pairs['distance'].max() <= 250000
    everything = graph_distance_pairs(graph, max_distance=None, max_steps=5)
    assert len(everything) == 10                         # all pairs of the five chr1 nodes
    assert everything.set_index(['node_a', 'node_b']).loc[(names[0], names[4]), 'steps'] == 4


def test_degree_bin_classes_by_kind():
    graph = nx.Graph()
    hub = _node(graph, 'chr1', 0, True)
    small = _node(graph, 'chr1', 1000000, True)
    others = [_node(graph, 'chr1', 10000 * (i + 1), False) for i in range(8)]
    graph.add_edges_from((hub, o) for o in others)
    graph.add_edge(small, others[0])
    node_bins = {n: i for i, n in enumerate(graph.nodes)}
    classes = degree_bin_classes(graph, node_bins, len(node_bins) + 3, n_classes=2)
    assert (classes[len(node_bins):] == -1).all()        # bins without nodes
    assert classes[node_bins[hub]] > classes[node_bins[small]] >= 2      # bait classes 2..3
    assert classes[node_bins[others[0]]] > classes[node_bins[others[1]]]  # degree 2 against degree 1
    assert all(0 <= classes[node_bins[o]] < 2 for o in others)


def test_rewiring_keeps_degrees_kinds_chromosomes_and_length_bins():
    graph = random_capture_graph()
    rewired, statistics = rewire_by_distance(graph, swaps_per_edge=5, log_bin_width=0.2, seed=3)
    assert statistics['accepted'] > 0 and statistics['fraction_edges_changed'] > 0
    assert rewired.number_of_edges() == graph.number_of_edges()
    assert dict(rewired.degree()) == dict(graph.degree())
    middle = node_midpoints(graph)

    def profile(g):
        rows = []
        for u, v in g.edges():
            assert g.nodes[u]['chrom'] == g.nodes[v]['chrom'] and u != v
            kinds = tuple(sorted((g.nodes[u]['bait'], g.nodes[v]['bait'])))
            rows.append((g.nodes[u]['chrom'], kinds, int(_length_bin(abs(middle[u] - middle[v]), 0.2))))
        return sorted(rows)
    assert profile(rewired) == profile(graph)
    # per node, the number of bait and of other-end neighbours is kept
    for node in graph.nodes:
        def kinds(g):
            return sorted(g.nodes[n]['bait'] for n in g.neighbors(node))
        assert kinds(rewired) == kinds(graph)
    again, _ = rewire_by_distance(graph, swaps_per_edge=5, log_bin_width=0.2, seed=3)
    assert set(map(frozenset, again.edges())) == set(map(frozenset, rewired.edges()))


def test_rewiring_wide_bins_and_length_limit():
    graph = random_capture_graph(seed=1)
    rewired, statistics = rewire_by_distance(graph, swaps_per_edge=5, log_bin_width=7.0, seed=0, max_length=50000)
    assert statistics['fraction_edges_changed'] > 0.3
    assert dict(rewired.degree()) == dict(graph.degree())
    middle = node_midpoints(graph)
    new = set(map(frozenset, rewired.edges())) - set(map(frozenset, graph.edges()))
    assert new and all(abs(middle[u] - middle[v]) <= 50000 for u, v in map(tuple, new))


def test_class_matched_controls_record_fallbacks():
    import pandas as pd

    from epigeneticgraphanalyzer.validation import sample_matched_controls

    class Pools:
        bait_array = np.arange(0, 100, 2)
        other_array = np.arange(1, 100, 2)
        chrom_index = np.zeros(200, dtype=np.int64)

        def pool(self, kind):
            return self.bait_array if kind == 'bait' else self.other_array

    bin_class = np.zeros(200, dtype=np.int64)
    bin_class[[10, 15]] = 7                       # a class no candidate at offset 5 shares
    pairs = pd.DataFrame({'bin_a': [20, 10], 'bin_b': [25, 15], 'offset': [5, 5],
                          'pair_class': ['bait-oe', 'bait-oe']})
    controls = sample_matched_controls(pairs, Pools(), {(20, 25), (10, 15)}, n_controls=4, seed=0,
                                       coverage_class=bin_class)
    assert controls.attrs['n_fallback'] == 1 and controls.attrs['fallback_index'].tolist() == [1]
    first = controls[controls['test_index'] == 0]
    assert (bin_class[first['bin_a']] == 0).all() and (bin_class[first['bin_b']] == 0).all()
    assert ((controls['bin_b'] - controls['bin_a']) == 5).all()


def test_graph_statistics():
    from epigeneticgraphanalyzer.graph_nulls import graph_statistics
    graph = nx.Graph()
    b = _node(graph, 'chr1', 1000000, True)
    o1 = _node(graph, 'chr1', 1100000, False)
    o2 = _node(graph, 'chr1', 4000000, False)
    b2 = _node(graph, 'chr2', 500000, True)
    t = _node(graph, 'chr3', 500000, False)
    graph.add_edges_from([(b, o1), (b, o2), (b2, t)])
    stats = graph_statistics(graph, max_distance=2000000)
    assert (stats['nodes'], stats['baits'], stats['edges'], stats['cis_edges']) == (5, 2, 3, 2)
    assert stats['cis_edges_within_limit'] == 1 and stats['median_edge_length'] == 0.5 * (100000 + 3000000)
    assert stats['mean_bait_degree'] == 1.5 and stats['largest_component_fraction'] == 3 / 5
