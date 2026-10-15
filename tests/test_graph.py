import networkx as nx
import pandas as pd

from epigeneticgraphanalyzer.graph import (annotate_genes, build_graph, community_table, detect_communities,
                                           neighbourhood_layers, read_graph, write_graph)
from epigeneticgraphanalyzer.indirect import direct_contacts, indirect_contacts, pair_class, two_hop_triplets


def _table(rows):
    """rows: (bait_start, other_start, other_is_bait, score)."""
    return pd.DataFrame([{
        'chrom': 'chr1', 'bait_start': b, 'bait_end': b + 999, 'other_start': o, 'other_end': o + 999,
        'other_is_bait': ob, 'distance': abs(o - b), 'count': 10.0, 'expected': 2.0, 'score': s,
        'significant': True} for b, o, ob, s in rows])


def toy_graph():
    # B1 = 100000, B2 = 300000 (baits); O1 = 50000, O2 = 200000, O3 = 400000
    rows = [(100000, 50000, False, 5.0), (100000, 200000, False, 4.0),
            (300000, 200000, False, 6.0), (300000, 400000, False, 3.0)]
    tss = pd.DataFrame({'chrom': ['chr1', 'chr1'], 'tss': [100500, 300500], 'strand': ['+', '-'],
                        'gene_id': ['ENSG1', 'ENSG2'], 'gene_name': ['B1', 'B2'], 'gene_type': ['pc', 'pc']})
    return build_graph(_table(rows), tss=tss, weight='score')


def test_nodes_are_shared_between_bait_and_other_end():
    rows = [(100000, 300000, True, 5.0), (300000, 500000, False, 5.0)]
    graph = build_graph(_table(rows))
    assert graph.number_of_nodes() == 3
    assert graph.nodes['chr1:300000-300999']['bait']
    assert nx.shortest_path_length(graph, 'chr1:100000-100999', 'chr1:500000-500999') == 2


def test_gene_annotation_and_attributes():
    graph = toy_graph()
    assert graph.nodes['chr1:100000-100999']['genes'] == 'B1'
    assert graph.nodes['chr1:300000-300999']['gene_ids'] == 'ENSG2'
    assert graph.nodes['chr1:50000-50999']['bait'] is False
    assert graph['chr1:100000-100999']['chr1:200000-200999']['weight'] == 4.0


def test_indirect_contacts_classes_and_scores():
    graph = toy_graph()
    pairs = indirect_contacts(graph)
    found = {(a, b): c for a, b, c in zip(pairs['node_a'], pairs['node_b'], pairs['pair_class'])}
    o1, o2, o3 = 'chr1:50000-50999', 'chr1:200000-200999', 'chr1:400000-400999'
    b1, b2 = 'chr1:100000-100999', 'chr1:300000-300999'
    assert found[(o1, o2)] == 'oe-oe'
    assert found[(o2, o3)] == 'oe-oe'
    assert found[(b1, b2)] == 'bait-bait'
    assert len(found) == 3
    row = pairs[(pairs['node_a'] == b1) & (pairs['node_b'] == b2)].iloc[0]
    assert row['n_mediators'] == 1
    assert row['path_score'] == 4.0  # min(4, 6)
    assert row['distance'] == 200000


def test_indirect_contacts_respect_distance_and_existing_edges():
    graph = toy_graph()
    graph.add_edge('chr1:50000-50999', 'chr1:200000-200999', weight=1.0)
    pairs = indirect_contacts(graph, max_distance=180000)
    assert len(pairs) == 0 or not ((pairs['node_a'] == 'chr1:50000-50999') &
                                   (pairs['node_b'] == 'chr1:200000-200999')).any()
    assert (pairs['distance'] <= 180000).all()


def test_triplets_and_direct_contacts():
    graph = toy_graph()
    triplets = two_hop_triplets(graph)
    assert len(triplets) == 3
    assert set(triplets['mediator']) == {'chr1:100000-100999', 'chr1:200000-200999', 'chr1:300000-300999'}
    direct = direct_contacts(graph)
    assert len(direct) == 4
    assert set(direct['pair_class']) == {'bait-oe'}


def test_pair_class():
    assert pair_class(True, True) == 'bait-bait'
    assert pair_class(True, False) == 'bait-oe'
    assert pair_class(False, False) == 'oe-oe'


def test_communities_split_two_cliques():
    graph = nx.Graph()
    for offset in (0, 10):
        for i in range(5):
            for j in range(i + 1, 5):
                graph.add_edge('chr1:{}-{}'.format((offset + i) * 1000, (offset + i) * 1000 + 999),
                               'chr1:{}-{}'.format((offset + j) * 1000, (offset + j) * 1000 + 999), weight=1.0)
    graph.add_edge('chr1:4000-4999', 'chr1:10000-10999', weight=0.01)
    for node in graph.nodes:
        start = int(node.split(':')[1].split('-')[0])
        graph.nodes[node].update(chrom='chr1', start=start, end=start + 999, bait=False, genes='')
    for method in ('leiden', 'louvain'):
        communities = detect_communities(graph, method=method)
        assert len(communities) == 2
        assert sorted(len(c) for c in communities) == [5, 5]
    table = community_table(graph, communities)
    assert list(table['n_nodes']) == [5, 5]


def test_graphml_roundtrip(tmp_path):
    graph = toy_graph()
    detect_communities(graph)
    path = str(tmp_path / 'g.graphml')
    write_graph(graph, path)
    loaded = read_graph(path)
    assert loaded.number_of_edges() == graph.number_of_edges()
    assert loaded.nodes['chr1:100000-100999']['bait'] is True
    assert loaded.nodes['chr1:100000-100999']['community'] == graph.nodes['chr1:100000-100999']['community']


def test_neighbourhood_layers():
    graph = toy_graph()
    layers = neighbourhood_layers(graph, 'chr1:100000-100999', 2)
    assert layers[0] == {'chr1:50000-50999', 'chr1:200000-200999'}
    assert layers[1] == {'chr1:300000-300999'}


def test_annotate_genes_handles_multiple_tss():
    fragments = pd.DataFrame({'node': ['chr1:0-999'], 'chrom': ['chr1'], 'start': [0], 'end': [999]})
    tss = pd.DataFrame({'chrom': ['chr1'] * 3, 'tss': [10, 20, 5000], 'strand': ['+'] * 3,
                        'gene_id': ['E1', 'E2', 'E3'], 'gene_name': ['A', 'B', 'C'], 'gene_type': [''] * 3})
    genes = annotate_genes(fragments, tss)
    assert genes['chr1:0-999'] == (['A', 'B'], ['E1', 'E2'])


def test_read_expression_with_trailing_tabs(tmp_path):
    from epigeneticgraphanalyzer.io import read_expression
    path = tmp_path / 'rpkm.tsv'
    path.write_text('gene_id\tE003\tE004\nENSG1.5\t1.5\t2.0\t\nENSG2\t0.0\t3.0\t\n')
    table = read_expression(str(path))
    assert list(table.index) == ['ENSG1', 'ENSG2']
    assert table.loc['ENSG1', 'E004'] == 2.0


def test_read_homer_with_bait_design(tmp_path):
    import numpy as np
    from epigeneticgraphanalyzer.io import read_interactions
    header = ('InteractionID\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)\tPeakID(2)\tchr(2)\t'
              'start(2)\tend(2)\tstrand(2)\tTotal Reads(2)\tDistance\tInteraction Reads\tExpected Reads\tZ-score\t'
              'LogP\tFDR(Benjamini)\tCircos Thickness')
    rows = ['i1\tchr1-50000\tchr1\t50312\t55312\t+\t100\tchr1-0\tchr1\t130\t5130\t+\t90\t50000\t40\t5\t3\t-46.05\t0.0\t10',
            'i2\tchr1-300000\tchr1\t300000\t305000\t+\t100\tchr1-400000\tchr1\t400000\t405000\t+\t90\t100000\t3\t2\t'
            '0.5\t-1.0\t0.9\t1']
    path = tmp_path / 'homer.txt'
    path.write_text(header + '\n' + '\n'.join(rows) + '\n')
    baits = pd.DataFrame({'chrom': ['chr1'], 'start': [1000], 'end': [2000]})
    table = read_interactions([str(path)], 'homer', baits=baits)
    assert len(table) == 1  # the second pair has no baited window
    row = table.iloc[0]
    assert row['bait_start'] == 0 and row['bait_end'] == 4999 and row['other_start'] == 50000  # snapped to the grid
    assert np.isclose(row['score'], 20.0, atol=0.01) and row['count'] == 40 and row['fdr'] == 0.0


def test_bait_design_orients_symmetric_calls():
    import numpy as np
    from epigeneticgraphanalyzer.io import apply_baits
    rows = [('chr1', 0, 999, 50000, 50999, 8.0), ('chr1', 50000, 50999, 0, 999, 8.0),  # one pair, both orientations
            ('chr1', 90000, 90999, 0, 999, 6.0),                                     # other end is the bait
            ('chr1', 0, 999, 200000, 200999, 7.0), ('chr1', 200000, 200999, 0, 999, 7.5),  # bait-bait
            ('chr1', 300000, 300999, 400000, 400999, 9.0)]                           # no bait
    table = pd.DataFrame(rows, columns=['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end', 'score'])
    table['other_is_bait'] = False
    baits = pd.DataFrame({'chrom': ['chr1', 'chr1'], 'start': [0, 200000], 'end': [999, 200999]})
    out = apply_baits(table, baits)
    assert len(out) == 3
    assert (out['bait_start'] == 0).all()
    assert list(out['other_start']) == [50000, 90000, 200000]
    assert list(out['other_is_bait']) == [False, False, True]
    assert np.isclose(out.loc[out['other_start'] == 200000, 'score'].iloc[0], 7.5)
