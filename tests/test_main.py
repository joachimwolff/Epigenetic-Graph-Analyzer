import unittest
import networkx as nx

from epigeneticgraphanalyzer.egaValidateWithSPRITE import build_bin_offsets, build_graph_triplets

class TestMain(unittest.TestCase):
    def test_build_bin_offsets(self):
        offsets = build_bin_offsets({'chr1': 1000, 'chr2': 1200}, 200)
        self.assertEqual(offsets['chr1'], 0)
        self.assertEqual(offsets['chr2'], 5)

    def test_build_graph_triplets_empty(self):
        graph = nx.Graph()
        graph.add_node('A', position='chr1:10-20')
        graph.add_node('B', position='chr1:30-40')
        graph.add_edge('A', 'B')
        triplets = build_graph_triplets(graph, {'A': 0, 'B': 1})
        self.assertEqual(triplets, [])

if __name__ == '__main__':
    unittest.main()