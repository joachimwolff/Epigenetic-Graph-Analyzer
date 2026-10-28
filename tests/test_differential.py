import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.differential import (ActivityOperators, annotate_chicdiff, classify_genes,
                                                  differential_edges_from_scores, gene_features, match_communities, nested_models,
                                                  permutation_test, union_graph)

UNIT = 150000


def _row(chrom, bait, other, other_is_bait, score):
    return {'chrom': chrom, 'bait_start': bait, 'bait_end': bait + 4999, 'other_start': other,
            'other_end': other + 4999, 'other_is_bait': other_is_bait, 'distance': abs(other - bait),
            'count': np.nan, 'expected': np.nan, 'score': float(score), 'significant': True}


def node(chrom, start):
    return '{}:{}-{}'.format(chrom, start, start + 4999)


def design(n_units=46, n_planted=30, seed=0):
    """Units S - O - B2 - O2: the RNA change of S follows the H3K27ac change at the two-hop promoter B2."""
    rng = np.random.default_rng(seed)
    rows_a, rows_b, units = [], [], []
    for k in range(n_units):
        chrom = 'chr1' if k < 26 else 'chr2'
        x = (k if k < 26 else k - 26) * UNIT + 10000
        s, o, b2, o2 = x, x + 20000, x + 40000, x + 60000
        for bait, other, ob in ((s, o, False), (b2, o, False), (b2, o2, False)):
            rows_a.append(_row(chrom, bait, other, ob, 17.5))
            rows_b.append(_row(chrom, bait, other, ob, 17.5))
        delta = rng.normal(0, 2)
        planted = k < n_planted
        rna = delta + rng.normal(0, 0.3) if planted else rng.normal(0, 2)
        units.append({'chrom': chrom, 'S': node(chrom, s), 'O': node(chrom, o), 'B2': node(chrom, b2),
                      'O2': node(chrom, o2), 'delta': delta, 'rna': rna, 'planted': planted})
    return pd.DataFrame(rows_a), pd.DataFrame(rows_b), units


def build(units, table_a, table_b):
    edges = differential_edges_from_scores(table_a, table_b, 'A', 'B')
    graph = union_graph(edges, 'A', 'B')
    signal_a = {n: 4.0 for n in graph.nodes}
    signal_b = {n: 4.0 for n in graph.nodes}
    for u in units:
        signal_b[u['B2']] = 4.0 * 2 ** u['delta']
    marks = {'H3K27ac': (signal_a, signal_b)}
    change = {u['S']: u['rna'] for u in units}
    for u in units:
        graph.nodes[u['S']]['genes'] = 'G_' + u['S']
    return edges, graph, marks, change


def test_score_edges_and_chicdiff_annotation():
    a = pd.DataFrame([_row('chr1', 0, 50000, False, 8), _row('chr1', 0, 80000, False, 6)])
    b = pd.DataFrame([_row('chr1', 0, 50000, False, 8), _row('chr1', 0, 120000, False, 12)])
    edges = differential_edges_from_scores(a, b)
    status = dict(zip(edges['node_b'], edges['status']))
    assert status == {'chr1:50000-54999': 'shared', 'chr1:80000-84999': 'lost', 'chr1:120000-124999': 'gained'}
    gained = edges[edges['node_b'] == 'chr1:120000-124999'].iloc[0]
    assert gained['enrichment_A'] == 0 and np.isclose(gained['enrichment_B'], 12 / 5)
    # Chicdiff regions may span several fragments; the gained edge lies inside a tested region
    chicdiff = pd.DataFrame({'chrom': ['chr1', 'chr1'], 'bait_start': [0, 0], 'bait_end': [4999, 4999],
                             'other_start': [110000, 40000], 'other_end': [130000, 60000],
                             'log2FC': [2.5, 0.1], 'pvalue': [1e-5, 0.8], 'padj': [1e-3, 0.9]})
    annotated = annotate_chicdiff(edges, chicdiff)
    row = annotated.set_index('node_b')
    assert row.loc['chr1:120000-124999', 'chicdiff_significant']
    assert row.loc['chr1:120000-124999', 'log2FC'] == 2.5
    assert row.loc['chr1:50000-54999', 'chicdiff_tested'] and not row.loc['chr1:50000-54999', 'chicdiff_significant']
    assert not row.loc['chr1:80000-84999', 'chicdiff_tested']


def test_activity_operators_two_hop_sets():
    table_a, table_b, units = design(n_units=3)
    _, graph, _, change = build(units, table_a, table_b)
    operators = ActivityOperators(graph, list(change))
    u = units[0]
    i = operators.sources.index(u['S'])
    direct = {operators.nodes[j] for j in operators.direct_a[i].indices}
    indirect = {operators.nodes[j] for j in operators.indirect_a[i].indices}
    assert direct == {u['O']}
    assert indirect == {u['B2']}
    # weights are products of enrichments (21/6)^2
    assert np.isclose(operators.indirect_a[i].data[0], (21 / 6) ** 2)


def test_indirect_level_explains_planted_genes():
    table_a, table_b, units = design()
    _, graph, marks, change = build(units, table_a, table_b)
    features = gene_features(graph, marks, change)
    assert np.allclose(features['H3K27ac_promoter'], 0)
    assert np.allclose(features['H3K27ac_direct'], 0)
    models = nested_models(features, ['H3K27ac'])
    assert models['delta_r2'].iloc[-1] > 0.15
    assert models['f_pvalue'].iloc[-1] < 1e-2
    classes = classify_genes(features, ['H3K27ac'], {'H3K27ac': 1}, rna_threshold=1.0, chromatin_threshold=1.0)
    planted = {u['S'] for u in units if u['planted']}
    explained = classes[classes['explanation'] == 'indirect']
    assert len(explained) > 0
    assert explained['node'].isin(planted).mean() > 0.8
    result = permutation_test(graph, marks, change, {'H3K27ac': 1}, n_permutations=30)
    assert result['empirical_p_delta_r2'] < 0.05
    assert result['observed_indirect_delta_r2'] > result['null_indirect_delta_r2_mean']


def test_repressive_mark_sign():
    features = pd.DataFrame({'node': ['a', 'b'], 'rna_log2FC': [2.0, 2.0], 'H3K27me3_promoter': [-2.0, 2.0],
                             'H3K27me3_direct': [0.0, 0.0], 'H3K27me3_indirect': [0.0, 0.0]})
    classes = classify_genes(features, ['H3K27me3'], {'H3K27me3': -1})
    assert list(classes['explanation']) == ['promoter', 'unexplained']


def test_match_communities():
    a = [{'x', 'y', 'z'}, {'p', 'q'}]
    b = [{'x', 'y'}, {'p', 'q', 'r'}]
    table = match_communities(a, b, bait_nodes={'x', 'p'})
    assert list(table['best_community_b']) == [0, 1]
    assert np.isclose(table['jaccard'].iloc[0], 2 / 3)
    assert table['bait_jaccard'].iloc[0] == 1.0


def test_score_based_differential_and_indirect_only_report(tmp_path):
    from epigeneticgraphanalyzer.differential import (ActivityOperators, consensus_scores,
                                                      differential_edges_from_scores, indirect_only_genes)
    from epigeneticgraphanalyzer.io import read_washu
    # gene S (bait) - shared enhancer O - promoter B2 with an H3K27ac gain; S itself has no change
    lines_a = ['chr1,10000,14999\tchr1,30000,34999\t8.0', 'chr1,50000,54999\tchr1,30000,34999\t9.0',
               'chr1,10000,14999\tchr1,90000,94999\t6.0']
    lines_b = ['chr1,10000,14999\tchr1,30000,34999\t8.5', 'chr1,50000,54999\tchr1,30000,34999\t9.5']
    paths = []
    for name, lines in (('a1', lines_a), ('a2', lines_a[:2]), ('b1', lines_b), ('b2', lines_b)):
        path = tmp_path / (name + '.txt')
        path.write_text('\n'.join(lines) + '\n')
        paths.append(str(path))
    table = read_washu(paths[0])
    assert set(table['other_is_bait']) == {False}
    # the WashU browser layout of the same calls gives the same table
    browser = tmp_path / 'browser.txt'
    browser.write_text('\n'.join('{}\t{}\t{}\t{}:{}-{},{}\t1\t.'.format(*b.split(','), *o.split(','), s)
                                 for b, o, s in (line.split('\t') for line in lines_a)) + '\n')
    assert read_washu(str(browser)).equals(table)
    a = consensus_scores([read_washu(paths[0]), read_washu(paths[1])], 'intersect')
    assert len(a) == 2
    assert len(consensus_scores([read_washu(paths[0]), read_washu(paths[1])], 'union')) == 3
    b = consensus_scores([read_washu(paths[2]), read_washu(paths[3])])
    edges = differential_edges_from_scores(a, b, 'A', 'B')
    assert set(edges['status']) == {'shared'}
    graph = union_graph(edges, 'A', 'B')
    s, o, b2 = 'chr1:10000-14999', 'chr1:30000-34999', 'chr1:50000-54999'
    for n in graph.nodes:
        graph.nodes[n]['genes'] = {s: 'GeneS', b2: 'GeneB2'}.get(n, '')
    signal_a = {n: 4.0 for n in graph.nodes}
    signal_b = dict(signal_a)
    signal_b[b2] = 32.0
    marks = {'H3K27ac': (signal_a, signal_b)}
    change = {s: 3.0}
    operators = ActivityOperators(graph, list(change))
    vectors = {m: (operators.vector(x), operators.vector(y)) for m, (x, y) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=vectors)
    report = indirect_only_genes(graph, operators, vectors, features, ['H3K27ac'], {'H3K27ac': 1})
    assert list(report['genes']) == ['GeneS']
    assert 'GeneB2' in report['drivers'].iloc[0]
    assert o in report['drivers'].iloc[0]


def test_gene_level_indirect_test_finds_planted_genes():
    from epigeneticgraphanalyzer.differential import gene_level_indirect_test
    table_a, table_b, units = design(n_units=46, n_planted=30)
    _, graph, marks, change = build(units, table_a, table_b)
    operators = ActivityOperators(graph, list(change))
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=vectors)
    result = gene_level_indirect_test(operators, vectors, features, ['H3K27ac'], {'H3K27ac': 1}, rna_threshold=0.5,
                                      n_permutations=200)
    planted = {u['S'] for u in units if u['planted']}
    strong = {u['S'] for u in units if u['planted'] and abs(u['delta']) >= 2}
    p = result.set_index('node')['indirect_p']
    assert p[[n for n in p.index if n in strong]].median() < 0.05
    assert p[[n for n in p.index if n in planted]].median() < 0.5 * p[[n for n in p.index if n not in planted]].median()
    assert ((result['indirect_q'] >= result['indirect_p']) | np.isclose(result['indirect_q'], result['indirect_p'])).all()


def test_community_transitions_and_gene_switches():
    from epigeneticgraphanalyzer.differential import (community_expression_changes, community_transitions,
                                                      gene_community_changes)
    p_nodes = {'p{}'.format(i) for i in range(1, 11)}
    n_nodes = {'n{}'.format(i) for i in range(1, 6)}
    a = [{'s1', 's2', 's3', 's4'}, p_nodes, {'m1', 'm2', 'm3'}, n_nodes, {'d1', 'd2', 'd3', 'd4'}]
    b = [{'s1', 's2', 's3', 's4'}, {'p1', 'p2'}, p_nodes - {'p1', 'p2'}, {'m1', 'm2', 'm3'} | n_nodes,
         {'f1', 'f2', 'f3'}, {'d1'}]
    nodes_a = set().union(*a)
    nodes_b = set().union(*b)
    table, main = community_transitions(a, b, nodes_a, nodes_b)
    kind = dict(zip(table.loc[table['side'] == 'A', 'community_a'], table.loc[table['side'] == 'A', 'transition']))
    assert kind == {0: 'stable', 1: 'split', 2: 'merged', 3: 'merged', 4: 'dissolved'}
    assert list(table.loc[table['side'] == 'B', 'community_b']) == [4]
    genes = {n: n.upper() for n in nodes_a | nodes_b}
    baits = {'s1', 'p1', 'm1', 'd2'}
    expression = {'s1': 0.1, 'p1': 2.0, 'm1': 3.0, 'd2': -1.0}
    changes, test = gene_community_changes(a, b, table, genes, baits, expression=expression,
                                           differential={'p1', 'm1'})
    status = dict(zip(changes['node'], changes['status']))
    assert status == {'s1': 'kept', 'p1': 'switch', 'm1': 'switch', 'd2': 'lost'}
    assert 'N1' in changes.set_index('node').loc['m1', 'joined_genes']
    assert test['switching'] == 2 and test['switching_differential'] == 2
    edges = pd.DataFrame({'node_a': ['s1', 's1', 'm1'], 'node_b': ['s2', 's3', 'n1'],
                          'status': ['shared', 'lost', 'gained']})
    expression = {**{n: 3.0 for n in ['m1', 'm2', 'm3', 'n1', 'n2', 'n3', 'n4', 'n5']},
                  **{n: float(v) for n, v in zip(['s1', 's2', 's3', 's4'], [0.1, -0.2, 0.0, 0.1])},
                  **{'x{}'.format(i): v for i, v in enumerate(np.random.default_rng(1).normal(0, 1, 200))}}
    result = community_expression_changes(a, b, main, table, edges, expression, n_permutations=500)
    by_a = result.set_index('community_a')
    assert by_a.loc[0, 'edges_shared'] == 1 and by_a.loc[0, 'edges_lost'] == 1
    assert by_a.loc[2, 'coherence_p'] < 0.01 and by_a.loc[0, 'coherence_p'] > 0.2
    # A communities without a main B community are tested separately, A communities sharing one B are one unit
    a = [{'u1', 'u2', 'u3'}, {'v1', 'v2', 'v3'}, {'w1', 'w2', 'w3'}, {'y1', 'y2', 'y3'}]
    b = [{'w1', 'w2', 'w3', 'y1', 'y2', 'y3'}]
    expression = {n: 1.0 for n in set().union(*a)}
    expression.update({'x{}'.format(i): v for i, v in enumerate(np.random.default_rng(2).normal(0, 1, 100))})
    empty = pd.DataFrame({'node_a': [], 'node_b': [], 'status': []})
    transitions = pd.DataFrame({'side': ['A'] * 4, 'community_a': [0, 1, 2, 3],
                                'transition': ['dissolved', 'dissolved', 'merged', 'merged']})
    result = community_expression_changes(a, b, {2: 0, 3: 0}, transitions, empty, expression, n_permutations=100)
    assert sorted(result.loc[result['community_b'] < 0, 'community_a']) == [0, 1]
    assert list(result.loc[result['community_b'] == 0, 'communities_a']) == ['2,3']


def test_replicate_rule_present_in_all_called_in_one():
    from epigeneticgraphanalyzer.differential import consensus_scores
    # pair x: 7 and 3 (present in both, called in one) -> kept; y: 8 and missing -> dropped; z: 3 and 4 -> dropped
    r1 = pd.DataFrame([_row('chr1', 0, 50000, False, 7), _row('chr1', 0, 80000, False, 8),
                       _row('chr1', 0, 90000, False, 3)])
    r2 = pd.DataFrame([_row('chr1', 0, 50000, False, 3), _row('chr1', 0, 90000, False, 4)])
    kept = consensus_scores([r1, r2], 'all', threshold=5)
    assert list(kept['other_start']) == [50000]
    assert kept['score'].iloc[0] == 5.0 and kept['n_replicates_called'].iloc[0] == 1
    assert len(consensus_scores([r1, r2], 'intersect', threshold=5)) == 0
    assert set(consensus_scores([r1, r2], 'union', threshold=5)['other_start']) == {50000, 80000}
    # without a threshold (published, already thresholded files) 'all' requires the pair in every replicate
    assert set(consensus_scores([r1, r2], 'all')['other_start']) == {50000, 90000}


def test_edge_chromatin_association():
    from epigeneticgraphanalyzer.egaExplain import edge_chromatin, matched_concordance
    rng = np.random.default_rng(0)
    n = 300
    edges = pd.DataFrame({'status': ['shared'] * n + ['gained'] * n + ['lost'] * n,
                          'bin': rng.integers(1, 4, 3 * n)})
    # active mark: up at gained, down at lost anchors; repressive mark the opposite
    edges['H3K27ac'] = np.concatenate([rng.normal(0, 1, n), rng.normal(1, 1, n), rng.normal(-1, 1, n)])
    edges['H3K27me3'] = np.concatenate([rng.normal(0, 1, n), rng.normal(-1, 1, n), rng.normal(1, 1, n)])
    table = edge_chromatin(edges, ['H3K27ac', 'H3K27me3'], {'H3K27me3'}).set_index(['mark', 'interactions'])
    assert table.loc[('H3K27ac', 'gained'), 'concordance_vs_shared'] > 0.7
    assert table.loc[('H3K27ac', 'lost'), 'concordance_vs_shared'] < 0.3
    assert table['concordant'].all()
    assert table['mannwhitney_p_expected_direction'].max() < 1e-6
    assert matched_concordance(np.array([1.0]), np.array([1.0]), np.array([1]), np.array([1])) == 0.5
