import importlib
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

import networkx as nx


class FakeScatter:
    def __init__(self, *args, **kwargs):
        del args
        self.kwargs = kwargs
        self.marker = types.SimpleNamespace(color=None)
        self.text = None


class FakeLayout:
    def __init__(self, *args, **kwargs):
        del args, kwargs


class FakeFigure:
    def __init__(self, data=None, layout=None):
        self.data = data
        self.layout = layout

    def write_html(self, path):
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write('<html></html>')


def import_with_fake_plotly():
    fake_graph_objs = types.ModuleType('plotly.graph_objs')
    fake_graph_objs.Scatter = FakeScatter
    fake_graph_objs.Layout = FakeLayout
    fake_graph_objs.Figure = FakeFigure

    fake_plotly = types.ModuleType('plotly')
    fake_plotly.graph_objs = fake_graph_objs

    with mock.patch.dict(sys.modules, {'plotly': fake_plotly, 'plotly.graph_objs': fake_graph_objs}):
        module_name = 'epigeneticgraphanalyzer.egaPlotlyGraph'
        if module_name in sys.modules:
            del sys.modules[module_name]
        return importlib.import_module(module_name)


class TestPlotlyCli(unittest.TestCase):
    def test_plotly_cli_main(self):
        module = import_with_fake_plotly()

        with tempfile.TemporaryDirectory() as tmpdir:
            graph_path = os.path.join(tmpdir, 'g.graphml')
            out_dir = os.path.join(tmpdir, 'out')
            graph = nx.Graph()
            graph.add_edge('A', 'B')
            nx.write_graphml(graph, graph_path)

            argv = ['egaPlotlyGraph', graph_path, out_dir, 'prefix_', '--partitionId', '7']
            with mock.patch('sys.argv', argv):
                module.main()

            expected = os.path.join(out_dir, 'prefix_7.html')
            self.assertTrue(os.path.exists(expected))


if __name__ == '__main__':
    unittest.main()
