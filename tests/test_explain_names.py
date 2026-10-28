import json

import networkx as nx

from epigeneticgraphanalyzer.egaExplain import condition_names


def condition_graph():
    """A union graph whose first node is a differentially expressed promoter (rna_log2FC) and whose condition B sorts
    before condition A, as naive B cells (A) against GM12878 (B)."""
    graph = nx.Graph()
    graph.add_node('chr1:1-10', chrom='chr1', start=1, end=10, rna_log2FC=2.0, H3K27ac_naiveB=1.0,
                   H3K27ac_GM12878=3.0)
    graph.add_node('chr1:11-20', chrom='chr1', start=11, end=20, H3K27ac_naiveB=2.0, H3K27ac_GM12878=1.0)
    graph.add_edge('chr1:1-10', 'chr1:11-20')
    return graph


def test_condition_names_from_the_support_file(tmp_path):
    prefix = str(tmp_path / 'naiveB_vs_GM12878')
    with open(prefix + '_support.json', 'w') as handle:
        json.dump({'names': ['naiveB', 'GM12878']}, handle)
    assert condition_names(condition_graph(), prefix) == ['naiveB', 'GM12878']


def test_condition_names_from_the_graph_keep_the_order_and_skip_expression():
    assert condition_names(condition_graph()) == ['naiveB', 'GM12878']
