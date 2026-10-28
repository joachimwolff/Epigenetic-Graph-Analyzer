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
    eligible = controls.attrs['eligible']
    assert eligible[0] > 0 and eligible[1] > eligible[0]     # the fallback draws from all offset-matched pairs


def test_balance_helpers():
    from epigeneticgraphanalyzer.graph_nulls import local_edge_density, per_bin
    from epigeneticgraphanalyzer.validation import standardized_mean_difference
    assert standardized_mean_difference([1, 2, 3], [1, 2, 3]) == 0
    assert np.isclose(standardized_mean_difference([2, 4], [0, 2]), 2 / np.sqrt(2))
    graph = nx.Graph()
    names = [_node(graph, 'chr1', p, i == 0) for i, p in enumerate((0, 50000, 120000, 500000, 520000))]
    graph.add_edges_from([(names[0], names[1]), (names[0], names[2]), (names[3], names[4])])
    density = local_edge_density(graph, window=100000)
    assert density[names[0]] == 2 and density[names[1]] == 2 and density[names[3]] == 1
    values = per_bin({names[0]: 0, names[1]: 0, names[3]: 2}, {n: graph.degree(n) for n in names}, 3)
    assert values[0] == 2 and np.isnan(values[1]) and values[2] == 1


def test_walk_counts_and_personalized_pagerank():
    from epigeneticgraphanalyzer.graph_nulls import adjacency, pair_ppr, personalized_pagerank, walk_counts
    graph = nx.Graph([('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'a'), ('c', 'e'), ('e', 'f')])
    A, nodes = adjacency(graph)
    i = {n: k for k, n in enumerate(nodes)}
    # a - c: two shortest paths (via b and via d), b - e one (via c); a - e: a-b-c-e and a-d-c-e, d - f: d-c-e-f
    assert walk_counts(A, [i['a'], i['b']], [i['c'], i['e']], 2).tolist() == [2.0, 1.0]
    assert walk_counts(A, [i['a'], i['d']], [i['e'], i['f']], 3).tolist() == [2.0, 1.0]
    x = personalized_pagerank(A, [i['a']], damping=0.85, iterations=200)[:, 0]
    reference = nx.pagerank(graph, alpha=0.85, personalization={'a': 1.0}, weight=None, tol=1e-12, max_iter=5000)
    assert np.allclose(x, [reference[n] for n in nodes], atol=1e-8)
    score = pair_ppr(A, [i['a']], [i['f']], iterations=200, batch=1)
    reverse = nx.pagerank(graph, alpha=0.85, personalization={'f': 1.0}, weight=None, tol=1e-12, max_iter=5000)
    assert np.isclose(score[0], 0.5 * (reference['f'] + reverse['a']), atol=1e-8)


def test_thinning_and_top_degree_nodes():
    from epigeneticgraphanalyzer.graph_nulls import thin_pairs, top_degree_nodes
    graph = nx.Graph()
    for name, bait in (('p1', True), ('p2', True), ('m', False), ('o1', False), ('o2', False), ('o3', False)):
        graph.add_node(name, bait=bait)
    graph.add_edges_from([('p1', 'm'), ('p2', 'm'), ('m', 'o1'), ('m', 'o2'), ('p1', 'o3'), ('o3', 'o2')])
    assert top_degree_nodes(graph, 0.1) == {'m'}
    a = np.array(['p1', 'p1', 'p2', 'o1'])
    b = np.array(['p2', 'o1', 'o2', 'o2'])
    for seed in range(5):
        keep = thin_pairs(graph, a, b, by_bait=True, seed=seed)
        baits = [n for k in np.nonzero(keep)[0] for n in (a[k], b[k]) if graph.nodes[n]['bait']]
        assert len(baits) == len(set(baits)) and keep[3]          # o1 - o2 has no bait: always kept
        keep = thin_pairs(graph, a, b, by_bait=False, by_mediator=True, seed=seed)
        # p1-p2, p1-o1, p2-o2 and o1-o2 all have m as a mediator: one pair is kept
        assert keep.sum() == 1


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


def test_rewiring_keeps_edge_attributes_with_the_anchor():
    graph = random_capture_graph(seed=2)
    for k, (u, v) in enumerate(graph.edges()):
        graph.edges[u, v]['score'] = float(k)
    plain, _ = rewire_by_distance(graph, swaps_per_edge=5, log_bin_width=0.2, seed=4)
    rewired, _ = rewire_by_distance(graph, swaps_per_edge=5, log_bin_width=0.2, seed=4, keep_attributes=True)
    assert set(map(frozenset, plain.edges())) == set(map(frozenset, rewired.edges()))
    assert sorted(d['score'] for _, _, d in rewired.edges(data=True)) == \
        sorted(d['score'] for _, _, d in graph.edges(data=True))
    # per bait, the scores of its edges to other ends move with it
    for node in graph.nodes:
        if not graph.nodes[node]['bait']:
            continue

        def scores(g):
            return sorted(g.edges[node, n]['score'] for n in g.neighbors(node) if not g.nodes[n]['bait'])
        assert scores(rewired) == scores(graph)
