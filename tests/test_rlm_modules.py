import importlib
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

import networkx as nx
import pandas as pd


class DummyBW:
    def __init__(self, path):
        self.path = path

    def stats(self, chrom, start, end, type='mean'):
        del chrom, start, end, type
        if 'rep1' in self.path:
            return [2.0]
        return [4.0]

    def close(self):
        return None


def import_with_fake_pybigwig(module_name):
    fake = types.ModuleType('pyBigWig')
    fake.open = lambda path: DummyBW(path)
    with mock.patch.dict(sys.modules, {'pyBigWig': fake}):
        if module_name in sys.modules:
            del sys.modules[module_name]
        return importlib.import_module(module_name)


class TestRlmModules(unittest.TestCase):
    def test_compute_rlm_main(self):
        compute = import_with_fake_pybigwig('epigeneticgraphanalyzer.egaComputeRLM')

        with tempfile.TemporaryDirectory() as tmpdir:
            graph_path = os.path.join(tmpdir, 'g.graphml')
            out_dir = os.path.join(tmpdir, 'out')
            deseq = os.path.join(tmpdir, 'deseq.tsv')

            g = nx.Graph()
            g.add_node('A', position='chr1:0-10', name='GeneA')
            g.add_node('B', position='chr1:10-20', name='GeneB')
            g.add_edge('A', 'B')
            nx.write_graphml(g, graph_path)

            pd.DataFrame([
                {'Gene': 'GeneA', 'log2FoldChange': 1.5},
                {'Gene': 'GeneB', 'log2FoldChange': -1.0},
            ]).to_csv(deseq, sep='\t', index=False)

            argv = [
                'egaComputeRLM',
                '--graph', graph_path,
                '--conditions', 'condA',
                '--bigwigFiles', 'condA_rep1.bw,condA_rep2.bw',
                '--deseqCounts', deseq,
                '--outputFolder', out_dir,
            ]
            with mock.patch('sys.argv', argv):
                compute.main()

            self.assertTrue(os.path.exists(os.path.join(out_dir, 'rlm_table.csv')))
            self.assertTrue(os.path.exists(os.path.join(out_dir, 'graph_annotated.graphml')))

    def test_merge_rlm_main(self):
        merge = import_with_fake_pybigwig('epigeneticgraphanalyzer.egaMergeRLM')

        with tempfile.TemporaryDirectory() as tmpdir:
            graph_path = os.path.join(tmpdir, 'g.graphml')
            out_dir = os.path.join(tmpdir, 'out')

            g = nx.Graph()
            g.add_node('A', position='chr1:0-10', name='GeneA')
            g.add_node('B', position='chr1:10-20', name='GeneB')
            g.add_edge('A', 'B')
            nx.write_graphml(g, graph_path)

            argv = [
                'egaMergeRLM',
                '--graph', graph_path,
                '--conditions', 'condA', 'condB',
                '--bigwigFiles', 'condA_rep1.bw,condA_rep2.bw', 'condB_rep1.bw,condB_rep2.bw',
                '--outputFolder', out_dir,
            ]
            with mock.patch('sys.argv', argv):
                merge.main()

            out_table = os.path.join(out_dir, 'merged_rlm_table.csv')
            self.assertTrue(os.path.exists(out_table))
            df = pd.read_csv(out_table)
            self.assertIn('condB_vs_condA_log2FC', df.columns)


if __name__ == '__main__':
    unittest.main()
