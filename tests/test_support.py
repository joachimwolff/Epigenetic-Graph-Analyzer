"""Support classes of gained and lost edges (stable, threshold-sensitive, not testable), their propagation to the
all-level explanations, provenance tags, and the earlier columns unchanged (toy data)."""

import os
import re

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import egaDifferential, egaExplain, egaReport, explorer
from epigeneticgraphanalyzer.differential import (ActivityOperators, classify_genes, consensus_scores,
                                                  differential_edges_from_scores, explanation_level_counts,
                                                  explanation_level_summary, explanation_levels, gene_features,
                                                  union_graph)
from epigeneticgraphanalyzer.support import (PathSupport, all_levels_provenance, edge_support, input_detection,
                                             level_provenance, support_summary)
from synthetic import fragments, simulate_capture, write_bigwig, write_tss_bed


def pairs(rows):
    """Interaction table of (bait start, other start, score) on chr1, fragments of 1 kb."""
    return pd.DataFrame([dict(chrom='chr1', bait_start=b, bait_end=b + 1000, other_start=o, other_end=o + 1000,
                              other_is_bait=False, distance=abs(o - b), count=np.nan, expected=np.nan, score=s)
                         for b, o, s in rows])


def toy_edges(b_floor_files=True):
    """Condition A: two replicates exported down to score 1; condition B: two replicates holding called pairs only
    (score >= 5). Pairs (bait 0, other ...):
      10000 shared, 20000 gained with A 4 / 4.5 (sensitive), 30000 gained with A 1.5 / 2 (stable),
      40000 gained, absent from A (stable: A lists scores down to 1), 50000 lost, absent from B (not testable),
      60000 lost, called in one B replicate only (sensitive: replicate)."""
    a1 = pairs([(0, 10000, 9), (0, 20000, 4), (0, 30000, 1.5), (0, 50000, 8), (0, 60000, 7), (0, 90000, 1.0)])
    a2 = pairs([(0, 10000, 8), (0, 20000, 4.5), (0, 30000, 2), (0, 50000, 6), (0, 60000, 9), (0, 90000, 1.2)])
    b1 = pairs([(0, 10000, 7), (0, 20000, 8), (0, 30000, 9), (0, 40000, 6), (0, 60000, 6)])
    b2 = pairs([(0, 10000, 9), (0, 20000, 7), (0, 30000, 8), (0, 40000, 7)])
    tables = {'A': [a1, a2], 'B': [b1, b2]}
    table_a = consensus_scores(tables['A'], 'all', 5)
    table_b = consensus_scores(tables['B'], 'intersect')
    edges = differential_edges_from_scores(table_a, table_b, 'A', 'B', 5.0)
    detection = input_detection(tables, {'A': ['a1', 'a2'], 'B': ['b1', 'b2']}, ('A', 'B'), 5.0, 3.0)
    return edges, tables, detection


def test_input_detection_from_lowest_score():
    _, _, detection = toy_edges()
    assert list(detection['min_score']) == [1.0, 1.2, 6.0, 7.0]
    assert list(detection['resolves_floor']) == [True, True, False, False]
    assert list(detection['absence']) == ['below floor', 'below floor', 'not testable', 'not testable']


def test_support_classes():
    edges, tables, detection = toy_edges()
    support = edge_support(edges, tables, detection, ('A', 'B'), 5.0, 3.0)
    got = dict(zip(edges['node_b'].str.split(':').str[1], zip(edges['status'], support['support'])))
    assert got['10000-11000'] == ('shared', 'stable')
    assert got['20000-21000'] == ('gained', 'threshold-sensitive')
    assert got['30000-31000'] == ('gained', 'stable')
    assert got['40000-41000'] == ('gained', 'stable')
    assert got['50000-51000'] == ('lost', 'not testable')
    assert got['60000-61000'] == ('lost', 'threshold-sensitive')
    reasons = dict(zip(edges['node_b'].str.split(':').str[1], support['support_reason']))
    assert 'sub-threshold score 4.5 in A' in reasons['20000-21000']
    assert 'called in 1 of 2 replicates of B' in reasons['60000-61000']
    assert 'holding scores >= 6' in reasons['50000-51000']
    other = dict(zip(edges['node_b'].str.split(':').str[1], support['other_max_score']))
    assert other['20000-21000'] == 4.5 and np.isnan(other['40000-41000']) and np.isnan(other['10000-11000'])
    summary = support_summary(pd.concat([edges, support], axis=1)).set_index(['status', 'support'])
    assert summary.loc[('gained', 'stable'), 'edges'] == 2
    assert summary.loc[('gained+lost', 'threshold-sensitive'), 'edges'] == 2


def test_appended_edge_columns_keep_the_earlier_ones():
    edges, tables, detection = toy_edges()
    annotated = pd.concat([edges, edge_support(edges, tables, detection, ('A', 'B'), 5.0, 3.0)], axis=1)
    assert list(annotated.columns[:len(edges.columns)]) == list(edges.columns)
    pd.testing.assert_frame_equal(annotated[edges.columns], edges)
    graph = union_graph(annotated, 'A', 'B')
    plain = union_graph(edges, 'A', 'B')
    for u, v, d in plain.edges(data=True):
        extra = dict(graph.edges[u, v])
        assert extra.pop('support') in ('stable', 'threshold-sensitive', 'not testable')
        assert np.isnan(extra.pop('qvalue')) and np.isnan(d['qvalue'])
        assert extra == {k: v for k, v in d.items() if k != 'qvalue'}


def chain_graph():
    """Gene g (up): g -shared- m1 -gained, threshold-sensitive- m2, m2 gains the active mark P (step 2). Gene h (up):
    h -shared- k1 -gained, stable- k2, k2 gains P. Gene f (up): f -lost, threshold-sensitive- f1 (own contact lost,
    expression up: no contact explanation) and f -gained, threshold-sensitive- f2 with P gained at f2 (step 1). H also
    reaches k4 (P gained) over a threshold-sensitive contact h - k3: its step 2 explains without it."""
    graph = nx.Graph()
    g, m1, m2 = 'chr1:0-1000', 'chr1:2000-3000', 'chr1:4000-5000'
    h, k1, k2, k3, k4 = 'chr2:0-1000', 'chr2:2000-3000', 'chr2:4000-5000', 'chr2:6000-7000', 'chr2:8000-9000'
    f, f1, f2 = 'chr3:0-1000', 'chr3:2000-3000', 'chr3:4000-5000'
    for u, v, status, support in ((g, m1, 'shared', 'stable'), (m1, m2, 'gained', 'threshold-sensitive'),
                                  (h, k1, 'shared', 'stable'), (k1, k2, 'gained', 'stable'),
                                  (h, k3, 'gained', 'threshold-sensitive'), (k3, k4, 'shared', 'stable'),
                                  (f, f1, 'lost', 'threshold-sensitive'), (f, f2, 'gained', 'threshold-sensitive')):
        graph.add_edge(u, v, r_a=0.0 if status == 'gained' else 1.0, r_b=0.0 if status == 'lost' else 1.0,
                       status=status, log2FC=0.0, qvalue=np.nan, support=support)
    for n in graph.nodes:
        chrom, rest = n.split(':')
        start, end = rest.split('-')
        graph.nodes[n].update(chrom=chrom, start=int(start), end=int(end), bait=n in (g, h, f),
                              genes={g: 'G', h: 'H', f: 'F'}.get(n, ''))
    p_a = {n: 1.0 for n in graph.nodes}
    p_b = dict(p_a, **{m2: 20.0, k2: 20.0, k4: 20.0, f2: 20.0})
    marks = {'P': (p_a, p_b)}
    return graph, marks, {g: 3.0, h: 3.0, f: 3.0}


def chain_levels(support=None):
    graph, marks, expression = chain_graph()
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    classes = classify_genes(features, ['P'], {'P': 1})
    return explanation_levels(classes, graph, operators, vectors, ['P'], {'declared': {'P': 1}},
                              names=('A', 'B'), support=support), graph


def test_path_support_weakest_edge_best_path():
    graph, _, _ = chain_graph()
    paths = PathSupport(graph, 'chr1:0-1000', 2)
    assert paths.node_class('chr1:2000-3000') == 'stable'
    assert paths.node_class('chr1:4000-5000') == 'threshold-sensitive'
    summary = paths.summary(['chr1:4000-5000'])
    assert summary['support'] == 'threshold-sensitive'
    assert summary['edges'] == {'stable': 1, 'threshold-sensitive': 1, 'not testable': 0}


def test_support_propagates_to_levels():
    levels, _ = chain_levels()
    by = levels.set_index(['genes', 'level'])
    g = by.loc[('G', 'step2')]
    assert g['explains_declared_mean'] and g['support_declared_mean'] == 'threshold-sensitive'
    assert g['threshold_sensitive_only_declared_mean']
    assert g['edges_stable_declared_mean'] == 1 and g['edges_threshold_sensitive_declared_mean'] == 1
    assert g['support_detail_declared_mean'] == 'rests on 1 threshold-sensitive contact'
    h = by.loc[('H', 'step2')]
    assert h['explains_declared_mean'] and h['support_declared_mean'] == 'stable'
    assert not h['threshold_sensitive_only_declared_mean']
    assert h['edges_threshold_sensitive_declared_mean'] == 1 and h['edges_stable_declared_mean'] == 3
    assert h['support_detail_declared_mean'].endswith('explains without the threshold-sensitive and untested ones')
    f = by.loc[('F', 'step1_marks')]
    assert f['explains_declared_mean'] and f['threshold_sensitive_only_declared_mean']
    # the promoter marks carry no contact; the own contacts of F: one gained, one lost (summed strength unchanged)
    assert by.loc[('G', 'promoter_marks'), 'support_declared_mean'] == ''
    assert by.loc[('F', 'promoter_contacts'), 'edges_threshold_sensitive_declared_mean'] == 2
    summary = explanation_level_summary(levels).set_index(levels.drop_duplicates('node').set_index('node')['genes'])
    assert summary.loc['G', 'only_threshold_sensitive_declared_mean']
    assert summary.loc['G', 'first_level_support_declared_mean'] == 'threshold-sensitive'
    assert not summary.loc['H', 'only_threshold_sensitive_declared_mean']
    assert summary.loc['F', 'threshold_sensitive_levels_declared_mean'] == 'step1_marks'
    counts = explanation_level_counts(levels).set_index('level')
    assert counts.loc['step2', 'any_threshold_sensitive_only'] == 1
    assert counts.loc['step2', 'any_support_stable'] == 1


def test_provenance_tags():
    levels, _ = chain_levels()
    by = levels.set_index(['genes', 'level'])
    assert by.loc[('G', 'promoter_marks'), 'provenance'] == 'measured'
    assert by.loc[('G', 'promoter_contacts'), 'provenance'] == 'called'
    assert by.loc[('G', 'step2'), 'provenance_chain'] == 'MC>D>I'
    assert level_provenance('promoter_state') == ('called', 'M>C>I')
    table = all_levels_provenance(levels.columns).set_index('column')
    assert table.loc['mark_changes', 'provenance'] == 'measured'
    assert table.loc['explains_declared_mean', 'provenance'] == 'interpretation'
    assert table.loc['contacts_gained', 'provenance'] == 'called'
    assert table.loc['support_declared_mean', 'provenance'] == 'derived'
    assert set(table['tag']) <= {'M', 'C', 'D', 'I', ''}


def test_levels_backward_compatible():
    """The columns written before the support and provenance columns are identical, in the same order."""
    with_support, _ = chain_levels()
    without, _ = chain_levels(support=False)
    graph, _, _ = chain_graph()
    for _, _, d in graph.edges(data=True):
        d.pop('support')
    assert 'support_declared_mean' not in without
    earlier = [c for c in without.columns if not c.startswith('provenance')]
    assert list(with_support.columns[:len(earlier)]) == earlier
    pd.testing.assert_frame_equal(with_support[earlier], without[earlier])
    new = list(with_support.columns[len(earlier):])
    assert new[:2] == ['provenance', 'provenance_chain'] and all(
        c.startswith(('support_', 'edges_', 'threshold_sensitive_only_')) for c in new[2:])
    summary_a = explanation_level_summary(with_support)
    summary_b = explanation_level_summary(without)
    pd.testing.assert_frame_equal(summary_a[summary_b.columns], summary_b)
    assert list(summary_a.columns[:len(summary_b.columns)]) == list(summary_b.columns)
    counts_a = explanation_level_counts(with_support)
    counts_b = explanation_level_counts(without)
    assert list(counts_a.columns[:len(counts_b.columns)]) == list(counts_b.columns)
    pd.testing.assert_frame_equal(counts_a[counts_b.columns], counts_b)


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    """egaDifferential on synthetic data: condition A exported down to score 1 (two replicates, gained pairs listed
    at sub-threshold scores 4 or 1.5), condition B called pairs only."""
    d = str(tmp_path_factory.mktemp('support'))
    da, db = os.path.join(d, 'a'), os.path.join(d, 'b')
    os.makedirs(da)
    os.makedirs(db)
    planted = sorted(simulate_capture(da, seed=7)[1])
    gained = sorted({(c, i, j + 3) for c, i, j in planted[:20]})
    ibed_a = [simulate_capture(da, seed=s, planted=planted, n_noise=0, name='a{}'.format(s))[0] for s in (7, 9)]
    # sub-threshold scores of the gained pairs in both A replicates: half at 4 (threshold-sensitive), half at 1.5
    for path in ibed_a:
        table = pd.read_csv(path, sep='\t')
        extra = []
        for k, (chrom, i, j) in enumerate(gained):
            bait, other = fragments(chrom)[i], fragments(chrom)[j]
            extra.append((chrom, bait[1], bait[2], 'GENE_{}'.format(i), chrom, other[1], other[2], '.', 5,
                          4.0 if k % 2 == 0 else 1.5))
        table = pd.concat([table, pd.DataFrame(extra, columns=table.columns)], ignore_index=True)
        table.to_csv(path, sep='\t', index=False)
    ibed_b, _, baits = simulate_capture(db, seed=8, planted=set(planted[20:]) | set(gained), n_noise=0)
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    rng = np.random.default_rng(0)
    bait_set = set(baits)
    up = {(c, i) for c, i in baits[::3]}
    ac_a = write_bigwig(os.path.join(d, 'ac_a.bw'), {})
    ac_b_values = {(c, j): 4.0 for c, _, j in gained}
    ac_b_values.update({b: 6.0 for b in up})
    ac_b = write_bigwig(os.path.join(d, 'ac_b.bw'), ac_b_values)
    me_a = write_bigwig(os.path.join(d, 'me_a.bw'), {b: 4.0 for b in bait_set})
    me_b = write_bigwig(os.path.join(d, 'me_b.bw'), {b: (1.0 if b in up else 4.0) for b in bait_set})
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)),
                               'B': [rng.uniform(150, 300) if b in up else rng.uniform(5, 50) for b in baits]})
    expression_path = os.path.join(d, 'expression.tsv')
    expression.to_csv(expression_path, sep='\t', index=False)
    prefix = os.path.join(d, 'out', 'cond1_vs_cond2')
    egaDifferential.main(['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--names', 'cond1',
                          'cond2', '--scoreThreshold', '5', '--replicates', 'all', '--supportFloor', '3',
                          '--tss', tss, '--tracksA', 'H3K27ac=' + ac_a, 'H3K27me3=' + me_a,
                          '--tracksB', 'H3K27ac=' + ac_b, 'H3K27me3=' + me_b,
                          '--markRoles', 'H3K27ac=active', 'H3K27me3=repressive',
                          '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                          '--permutations', '5', '--bootstrap', '20', '--roleBootstrap', '20',
                          '--genePermutations', '200', '--communityPermutations', '200', '--threads', '1',
                          '-o', prefix])
    explain = os.path.join(d, 'out', 'explain', 'cond1_vs_cond2')
    egaExplain.main(['-p', prefix, '-o', explain])
    return dict(dir=d, prefix=prefix, explain=explain)


def test_run_writes_support(run):
    prefix = run['prefix']
    edges = pd.read_csv(prefix + '_edges.tsv.gz', sep='\t')
    assert {'support', 'support_reason', 'other_max_score', 'replicates_floor_cond1', 'replicates_floor_cond2'} <= \
        set(edges.columns)
    assert list(edges.columns[-5:]) == ['support', 'support_reason', 'other_max_score', 'replicates_floor_cond1',
                                        'replicates_floor_cond2']
    gained = edges[edges['status'] == 'gained']
    assert (gained['support'] == 'threshold-sensitive').any() and (gained['support'] == 'stable').any()
    lost = edges[edges['status'] == 'lost']
    assert lost.empty or (lost['support'] == 'not testable').all()
    assert (edges.loc[edges['status'] == 'shared', 'support'] == 'stable').all()
    inputs = pd.read_csv(prefix + '_support_inputs.tsv', sep='\t')
    assert list(inputs['absence']) == ['below floor', 'below floor', 'not testable']
    summary = pd.read_csv(prefix + '_support.tsv', sep='\t')
    assert {'status', 'support', 'edges', 'fraction', 'threshold', 'floor'} <= set(summary.columns)
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t')
    assert {'provenance', 'provenance_chain', 'support_declared_mean', 'threshold_sensitive_only_declared_mean'} <= \
        set(levels.columns)
    provenance = pd.read_csv(prefix + '_explanations_provenance.tsv', sep='\t')
    assert list(provenance['column']) == list(levels.columns)
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    assert 'only_threshold_sensitive_declared_mean' in classes
    genes = pd.read_csv(run['explain'] + '_genes.tsv', sep='\t')
    assert {'first_level_support', 'only_threshold_sensitive', 'support_note', 'provenance'} <= set(genes.columns)


def test_report_flags_provenance_and_overlay(run):
    prefix = run['prefix']
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t')
    explained = levels[levels['explains_declared_mean'].astype(bool) & (levels['support_declared_mean'].fillna('') != '')]
    gene = str(explained['genes'].iloc[0]).split(',')[0] if len(explained) else \
        str(levels['genes'].iloc[0]).split(',')[0]
    out = os.path.join(run['dir'], 'reports')
    outputs = egaReport.main(['-p', prefix, '--gene', gene, '--overlay', '--pdf', 'none', '-o', out])
    text = open(outputs[0][0]).read()
    assert 'class="prov prov-M"' in text or 'class="prov prov-D"' in text
    assert 'Provenance of the statements' in text and 'measured' in text and 'interpretation' in text
    assert 'Plain overlay and EGA' in text
    assert 'Support of the contacts' in text
    assert '<th class="">Contacts</th>' in text
    stem = os.path.join(out, 'cond1_vs_cond2_gene_{}_overlay'.format(gene))
    assert os.path.getsize(stem + '_inputs.pdf') > 0 and os.path.getsize(stem + '_ega.pdf') > 0
    assert not re.search(r'(?:src|href)\s*=\s*["\']\s*(?:https?:)?//', text)


def test_explorer_carries_support(run):
    path = explorer.main(['-p', run['prefix'], '-o', os.path.join(run['dir'], 'explorer')])
    with open(path) as handle:
        text = handle.read()
    data = explorer.decode(text)
    edges = data['graph']['edges']
    assert len(edges['t']) == len(edges['u']) and set(edges['t']) <= {0, 1, 2} and 1 in edges['t']
    assert data['meta']['supportDefinition'] and 'prov prov-M' in data['meta']['provenance']
    gene = next(g for g, e in data['genes'].items() if e.get('evidence'))
    assert 'Source' in data['genes'][gene]['evidence'][0]['cols']
    assert 'sup-ts' in text and 'lad-ts' in text
