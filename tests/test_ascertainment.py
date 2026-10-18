"""Node mapping status of the pair validation (validation.node_mapping_status)."""

import networkx as nx
import numpy as np

from epigeneticgraphanalyzer.validation import map_nodes_to_bins, node_mapping_status


class Bins:
    def bin_of(self, chrom, position):
        return position // 1000 if chrom == 'chr1' else -1


class Data:
    bins = Bins()
    coverage_raw = np.array([3, 0, 5, 1, 2], dtype=float)


class Lift:
    def position(self, chrom, position):
        return None if position == 3499 else (chrom, position)


def test_node_mapping_status_explains_map_nodes_to_bins():
    graph = nx.Graph()
    for name, chrom, start, bait in (('a', 'chr1', 0, True), ('b', 'chr1', 1000, False), ('c', 'chr2', 0, False),
                                     ('d', 'chr1', 3000, True), ('e', 'chr1', 4000, False)):
        graph.add_node(name, chrom=chrom, start=start, end=start + 999, bait=bait)
    graph.add_edges_from([('a', 'b'), ('a', 'c'), ('a', 'd'), ('d', 'e')])
    table = node_mapping_status(graph, Data(), Lift()).set_index('node')
    assert table['status'].to_dict() == {'a': 'mapped', 'b': 'zero_coverage', 'c': 'no_bin',
                                         'd': 'liftover_failed', 'e': 'mapped'}
    assert table.loc['a', 'degree'] == 3 and table.loc['e', 'bin'] == 4
    mapped = map_nodes_to_bins(graph, Data(), Lift())
    assert set(mapped) == set(table.index[table['status'] == 'mapped'])
