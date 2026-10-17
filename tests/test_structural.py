"""Structural variants and copy number as a level of the explanation cascade (structural.py) on toy data: a planted
breakpoint explains a gene, without breakpoint input the outputs are unchanged."""

import os

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer import egaDifferential, egaExplain
from epigeneticgraphanalyzer.differential import (ActivityOperators, cascade_first_level, cascade_table,
                                                  classify_genes, explanation_cascade, explanation_details,
                                                  explanation_level_names, explanation_level_summary,
                                                  explanation_levels, gene_features)
from epigeneticgraphanalyzer.structural import (add_sv_edges, flag_sv_edges, read_copy_number,
                                                read_structural_variants, structural_genes, sv_edge_table, sv_label,
                                                sv_status)

from synthetic import fragments, simulate_capture, write_bigwig, write_tss_bed


def toy():
    """Gene G (up) on chr8 near a breakend whose partner on chr21 holds an active element in B only; gene H (up) on
    chr8 far from it; gene K (down) whose relative copy number drops. No mark changes anywhere."""
    graph = nx.Graph()
    g, g1 = 'chr8:1000000-1005000', 'chr8:1050000-1055000'
    h, h1 = 'chr8:9000000-9005000', 'chr8:9050000-9055000'
    k, k1 = 'chr8:20000000-20005000', 'chr8:20050000-20055000'
    p, p1 = 'chr21:5000000-5005000', 'chr21:5100000-5105000'
    for u, v in ((g, g1), (h, h1), (k, k1), (p, p1)):
        graph.add_edge(u, v, r_a=1.0, r_b=1.0, status='shared', log2FC=0.0, qvalue=np.nan)
    for n in graph.nodes:
        chrom, rest = n.split(':')
        start, end = rest.split('-')
        graph.nodes[n].update(chrom=chrom, start=int(start), end=int(end), bait=n in (g, h, k, p),
                              genes={g: 'G', h: 'H', k: 'K', p: 'RUNX1'}.get(n, ''))
    values = {n: 2.0 for n in graph.nodes}
    marks = {'P': (values, dict(values))}
    expression = {g: 3.0, h: 2.0, k: -2.0}
    elements = pd.DataFrame([(n, c, 'active enhancer' if (n == p1 and c == 'B') else 'unmarked element')
                             for n in graph.nodes for c in ('A', 'B')],
                            columns=['node', 'condition', 'class_declared_mean'])
    svs = pd.DataFrame([dict(chrom1='chr8', start1=1300000, end1=1400000, chrom2='chr21', start2=5050000,
                             end2=5150000, name='t(8;21)', present_A=False, present_B=True, pairs_A=0.0,
                             pairs_B=400.0, obs_exp_A=np.nan, obs_exp_B=150.0)])
    return graph, marks, expression, elements, svs, (g, h, k, p)


def run_levels(graph, marks, expression, structural=None, structural_levels=()):
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    classes = classify_genes(features, list(marks), {'P': 1})
    levels = explanation_levels(classes, graph, operators, vectors, list(marks), {'declared': {'P': 1}},
                                names=('A', 'B'), structural=structural, structural_levels=structural_levels)
    cascade = explanation_cascade(classes, list(marks), {'P': 1}, structural=structural)
    return classes, levels, cascade


def test_sv_label_and_status():
    assert sv_label('chr21', 'chr8') == 't(8;21)'
    assert sv_label('chrX', 'chr3') == 't(3;X)'
    assert sv_label('chr8', 'chr8', 12300000, 45600000) == 'chr8 join 12.3-45.6 Mb'
    svs = pd.DataFrame([dict(present_A=True, present_B=True, obs_exp_A=200.0, obs_exp_B=210.0, pairs_A=1, pairs_B=1),
                        dict(present_A=True, present_B=True, obs_exp_A=20.0, obs_exp_B=200.0, pairs_A=1, pairs_B=1),
                        dict(present_A=True, present_B=False, obs_exp_A=50.0, obs_exp_B=1.0, pairs_A=1, pairs_B=1),
                        dict(present_A=False, present_B=False, obs_exp_A=1.0, obs_exp_B=1.0, pairs_A=1, pairs_B=1),
                        dict(present_A=False, present_B=True, obs_exp_A=18.0, obs_exp_B=22.0, pairs_A=45, pairs_B=55)])
    status = sv_status(svs, ('A', 'B'))
    assert list(status['status']) == ['present in both', 'stronger in B', 'lost', 'absent', 'similar support in both']
    assert list(status['direction']) == [0, 1, -1, 0, 0]


def test_planted_breakpoint_explains_gene():
    graph, marks, expression, elements, svs, (g, h, k, p) = toy()
    svs = svs.assign(sv_id='sv1', type='inter', pos1=1350000, pos2=5100000)
    copy_number = {'A': pd.DataFrame({'chrom': ['chr8'], 'start': [19000000], 'end': [21000000], 'value': [1.0]}),
                   'B': pd.DataFrame({'chrom': ['chr8'], 'start': [19000000], 'end': [21000000], 'value': [0.5]})}
    genes, structural = structural_genes(graph, expression, svs, ('A', 'B'), window=1000000, elements=elements,
                                         marks=marks, copy_number=copy_number, differential=set(expression))
    # G: promoter 0.35 Mb from the breakend, join present in B only, partner with an active element in B
    assert list(genes['node']) == [g]
    row = genes.iloc[0]
    assert row['sv_status'] == 'gained' and row['hijacking_candidate']
    assert row['partner_active_B'] == 1 and row['partner_active_A'] == 0
    assert 'RUNX1' in row['partner_genes']
    assert structural.loc[g, 'sv_explains'] and structural.loc[g, 'sv_result'] == 'concordant'
    text = structural.loc[g, 'sv_detail']
    assert text.startswith('promoter 0.3 Mb from the t(8;21) breakpoint (chr8:1.35 Mb / chr21:5.10 Mb)')
    assert 'present in B only' in text and 'chromosome 21 near RUNX1' in text
    assert not structural.loc[h, 'sv_explains'] and structural.loc[h, 'sv_result'] == 'no breakpoint nearby'
    # K: copy number halved, expression down
    assert structural.loc[k, 'cn_explains'] and abs(structural.loc[k, 'cn_log2'] + 1.0) < 1e-9
    classes, levels, cascade = run_levels(graph, marks, expression, structural,
                                          ('structural_variant', 'copy_number'))
    assert list(levels.loc[levels['node'] == g, 'level']) == \
        explanation_level_names(2, ('structural_variant', 'copy_number'))
    assert explanation_level_names(2, ('structural_variant', 'copy_number'))[3:5] == ['structural_variant',
                                                                                        'copy_number']
    at = levels.set_index(['node', 'level'])
    assert at.loc[(g, 'structural_variant'), 'explains_declared_mean']
    assert at.loc[(g, 'structural_variant'), 'sv_name'] == 't(8;21)'
    assert not at.loc[(h, 'structural_variant'), 'explains_declared_mean']
    assert at.loc[(k, 'copy_number'), 'explains_declared_mean']
    assert cascade.loc[g, 'explained_at'] == 'structural_variant'
    assert cascade.loc[g, 'explained_type'] == 'structural variant'
    assert cascade.loc[k, 'explained_type'] == 'copy number'
    assert cascade.loc[h, 'explained_at'] == 'unexplained'
    summary = explanation_level_summary(levels)
    assert dict(cascade_first_level(cascade)) == dict(summary['first_level_declared_mean'])
    table = cascade_table(cascade, max_hops=2, structural=True).set_index('level')
    assert table.loc['structural_variant', 'genes'] == 2
    assert table.loc['structural_variant_breakpoint', 'genes'] == 1
    assert table.loc['structural_variant_copy_number', 'genes'] == 1
    assert table.loc['step1', 'genes_reaching'] == 1
    details = explanation_details(classes, cascade, graph, marks, {'P': 1}, structural=structural)
    assert details.loc[g, 'explanation_detail'].startswith('explained by a structural variant: promoter 0.3 Mb')
    assert details.loc[k, 'explanation_detail'].startswith('explained by a copy number change')
    # the same join in both conditions: reported, explains nothing
    both = svs.assign(present_A=True, pairs_A=380.0, obs_exp_A=140.0)
    _, structural_both = structural_genes(graph, expression, both, ('A', 'B'), window=1000000, elements=elements)
    assert structural_both.loc[g, 'sv_result'] == 'present in both'
    assert 'present in A and B' in structural_both.loc[g, 'sv_detail']
    assert not structural_both['sv_explains'].any()
    # discordant: the join is lost in B while G goes up
    lost = svs.assign(present_A=True, present_B=False)
    _, structural_lost = structural_genes(graph, expression, lost, ('A', 'B'), window=1000000, elements=elements)
    assert structural_lost.loc[g, 'sv_result'] == 'discordant'


def test_without_structural_input_unchanged():
    graph, marks, expression, _, _, _ = toy()
    classes, levels, cascade = run_levels(graph, marks, expression)
    assert list(levels['level'].unique()) == explanation_level_names(2)
    assert not ({'structural_variant', 'copy_number', 'sv_name'} & (set(levels.columns) | set(levels['level'])))
    assert set(cascade['explained_at']) <= {'gene', 'step1', 'step2', 'unexplained'}
    table = cascade_table(cascade, max_hops=2)
    assert 'structural_variant' not in set(table['level'])


def test_edges_flagged_and_junction_edges():
    graph, _, expression, _, svs, (g, h, k, p) = toy()
    svs = svs.assign(sv_id='sv1', type='inter', pos1=1350000, pos2=5100000)
    graph.add_edge(g, p, r_a=1.0, r_b=1.0, status='gained', log2FC=1.0, qvalue=np.nan)
    assert flag_sv_edges(graph, svs, 1000000) == 1
    assert graph.edges[g, p]['sv'] and graph.edges[g, p]['sv_name'] == 't(8;21)'
    assert 'sv' not in graph.edges[h, 'chr8:9050000-9055000']
    contacts = pd.DataFrame([dict(join='t(8;21)', bait_chrom='chr8', bait_start=1000000, bait_end=1005000,
                                  partner_chrom='chr21', partner_start=5000000, partner_end=5100000, pairs_A=0,
                                  pairs_B=12, inter_total_A=1e6, inter_total_B=1e6)])
    edges = sv_edge_table(graph, contacts, svs, ('A', 'B'))
    assert list(edges['status']) == ['gained'] and edges['sv'].all()
    combined = add_sv_edges(graph, edges)
    assert combined.edges[g, 'chr21:5000000-5100000']['sv']
    assert combined.nodes['chr21:5000000-5100000']['sv_partner']
    assert 'chr21:5000000-5100000' not in graph


def test_read_structural_variant_formats(tmp_path):
    d = str(tmp_path)
    hot = {'A': [('chr8', 92000000, 'chr21', 34000000, 1000, 5.0, 200.0), ('chr1', 0, 'chr2', 0, 10, 5.0, 2.0)],
           'B': [('chr8', 93000000, 'chr21', 34000000, 20, 5.0, 4.0)]}
    items = []
    for c, rows in hot.items():
        path = os.path.join(d, c + '_trans_hotspots.tsv')
        pd.DataFrame(rows, columns=['chrom1', 'start1', 'chrom2', 'start2', 'pairs', 'expected', 'obs_exp']) \
            .to_csv(path, sep='\t', index=False)
        items.append('{}={}'.format(c, path))
    svs = read_structural_variants(items, ('A', 'B'), min_pairs=50, min_obs_exp=20)
    assert len(svs) == 1
    r = svs.iloc[0]
    assert r['name'] == 't(8;21)' and r['present_A'] and not r['present_B']
    assert r['pairs_A'] == 1000 and r['start1'] == 92000000
    # BEDPE with a conditions column, plain BEDPE with an 11th column, condition=file
    bedpe = os.path.join(d, 'calls.bedpe')
    with open(bedpe, 'w') as handle:
        handle.write('#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tconditions\t'
                     'pairs_A\tobs_exp_A\tpairs_B\tobs_exp_B\n')
        handle.write('chr8\t92000000\t92100000\tchr21\t34800000\t34900000\tt(8;21)\t200\t.\t.\tA,B\t1000\t200\t900'
                     '\t180\n')
    plain = os.path.join(d, 'plain.bedpe')
    with open(plain, 'w') as handle:
        handle.write('chr3\t1000\t2000\tchr3\t9000000\t9001000\t.\t0\t.\t.\tB\n')
    single = os.path.join(d, 'single.bedpe')
    with open(single, 'w') as handle:
        handle.write('chr5\t1000\t2000\tchr7\t5000\t6000\tmy_join\n')
    svs = read_structural_variants([bedpe, plain, 'A=' + single], ('A', 'B'))
    assert list(svs['name']) == ['t(8;21)', 'chr3 join 0.0-9.0 Mb', 'my_join']
    assert list(svs['present_A']) == [True, False, True]
    assert list(svs['present_B']) == [True, True, False]
    assert list(svs['type']) == ['inter', 'intra', 'inter']
    status = sv_status(svs, ('A', 'B'))
    assert list(status['status']) == ['present in both', 'gained', 'lost']
    cn = os.path.join(d, 'a.bedGraph')
    with open(cn, 'w') as handle:
        handle.write('track type=bedGraph\nchr1\t0\t1000000\t1.0\n')
    assert read_copy_number(['A=' + cn], ('A', 'B'))['A']['value'].tolist() == [1.0]


def _differential_inputs(d):
    da, db = os.path.join(d, 'a'), os.path.join(d, 'b')
    os.makedirs(da)
    os.makedirs(db)
    planted = sorted(simulate_capture(da, seed=7)[1])
    gained = {(c, i, j + 3) for c, i, j in planted[:20]}
    ibed_a = [simulate_capture(da, seed=s, planted=planted, n_noise=0, name='a{}'.format(s))[0] for s in (7, 9)]
    ibed_b, _, baits = simulate_capture(db, seed=8, planted=set(planted[20:]) | gained, n_noise=0)
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    rng = np.random.default_rng(0)
    bw_a = write_bigwig(os.path.join(d, 'a.bw'), {})
    bw_b = write_bigwig(os.path.join(d, 'b.bw'), {(c, j): 4.0 for c, _, j in gained})
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)), 'B': rng.uniform(5, 50, len(baits))})
    path = os.path.join(d, 'expression.tsv')
    expression.to_csv(path, sep='\t', index=False)
    common = ['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--names', 'A', 'B', '--tss',
                                            tss, '--tracksA', 'H3K27ac=' + bw_a, '--tracksB', 'H3K27ac=' + bw_b,
                                            '--expression', path, '--expressionColumns', 'A', 'B', '--geneKey',
                                            'genes', '--permutations', '5', '--bootstrap', '20', '--roleBootstrap',
                                            '20', '--genePermutations', '50', '--threads', '1']
    return common, gained


def test_differential_cli_with_planted_breakpoint(tmp_path):
    d = str(tmp_path)
    common, gained = _differential_inputs(d)
    plain = os.path.join(d, 'out', 'plain')
    egaDifferential.main(common + ['-o', plain])
    classes = pd.read_csv(plain + '_gene_classes.tsv', sep='\t')
    # an up-regulated gene on chr2 that the promoter levels leave unexplained, and a partner element that
    # gains the mark in B (active in B only) on chr1
    candidates = classes[(classes['rna_log2FC'] > 0) & (classes['explained_at'] != 'gene') &
                         classes['node'].str.startswith('chr2:')]
    gene = candidates.iloc[0]
    start = int(gene['node'].split(':')[1].split('-')[0])
    partner = sorted((c, j) for c, _, j in gained if c == 'chr1')[0]
    other = fragments('chr1')[partner[1]]
    bedpe = os.path.join(d, 'sv.bedpe')
    with open(bedpe, 'w') as handle:
        handle.write('#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tconditions\t'
                     'pairs_A\tobs_exp_A\tpairs_B\tobs_exp_B\n')
        handle.write('chr2\t{}\t{}\tchr1\t{}\t{}\t.\t100\t.\t.\tB\t0\t0\t300\t100\n'.format(
            start + 20000, start + 30000, other[1], other[2] + 1))
    cn_a, cn_b = os.path.join(d, 'A.bedGraph'), os.path.join(d, 'B.bedGraph')
    with open(cn_a, 'w') as a, open(cn_b, 'w') as b:
        for chrom, size in (('chr1', 4000000), ('chr2', 3000000)):
            a.write('{}\t0\t{}\t1.0\n'.format(chrom, size))
            b.write('{}\t0\t{}\t1.0\n'.format(chrom, size))
    sv = os.path.join(d, 'out', 'sv')
    egaDifferential.main(common + ['--structuralVariants', bedpe, '--copyNumber', 'A=' + cn_a, 'B=' + cn_b,
                                   '--svWindow', '0.2', '-o', sv])
    sv_classes = pd.read_csv(sv + '_gene_classes.tsv', sep='\t').set_index('node')
    assert sv_classes.loc[gene['node'], 'sv_explains']
    assert sv_classes.loc[gene['node'], 'explained_at'] == 'structural_variant'
    assert 't(1;2) breakpoint' in sv_classes.loc[gene['node'], 'explanation_detail']
    genes = pd.read_csv(sv + '_sv_genes.tsv', sep='\t')
    assert gene['node'] in set(genes['node'])
    assert genes.set_index('node').loc[gene['node'], 'hijacking_candidate']
    levels = pd.read_csv(sv + '_explanations_all_levels.tsv', sep='\t')
    names = list(levels.loc[levels['node'] == gene['node'], 'level'])
    assert names[3:5] == ['structural_variant', 'copy_number']
    assert len(levels) == 9 * len(classes)
    cascade = pd.read_csv(sv + '_explanation_cascade.tsv', sep='\t').set_index('level')
    assert cascade.loc['structural_variant', 'genes'] >= 1
    assert cascade.loc['structural_variant_copy_number', 'genes'] == 0
    assert os.path.exists(sv + '_sv_breakpoints.tsv')
    # the outputs of the other levels do not change with the structural variant input
    plain_levels = pd.read_csv(plain + '_explanations_all_levels.tsv', sep='\t')
    shared = levels[~levels['level'].isin(['structural_variant', 'copy_number'])].drop(columns='rank')
    shared = shared[[c for c in plain_levels.columns if c != 'rank']].reset_index(drop=True)
    pd.testing.assert_frame_equal(shared, plain_levels.drop(columns='rank'), check_dtype=False)
    explained = sv_classes['explained_at'] != 'structural_variant'
    pd.testing.assert_series_equal(sv_classes.loc[explained, 'explained_at'],
                                   classes.set_index('node').loc[explained[explained].index, 'explained_at'])
    egaExplain.main(['-p', sv, '-o', os.path.join(d, 'out', 'explain')])
    report = pd.read_csv(os.path.join(d, 'out', 'explain_genes.tsv'), sep='\t').set_index('node')
    assert 'structural variant: promoter' in report.loc[gene['node'], 'all_levels']
    assert 'explained by a structural variant' in report.loc[gene['node'], 'summary']
    # the report names the structural variant in plain words, in the gene report and in the main report
    from epigeneticgraphanalyzer import egaReport
    name = str(gene['genes']).split(',')[0]
    out = os.path.join(d, 'out', 'report')
    outputs = egaReport.main(['-p', sv, '--explainPrefix', os.path.join(d, 'out', 'explain'), '--gene', name,
                              '--pdf', 'none', '-o', out])
    page = next(o[0] for o in outputs if o[0].endswith('_gene_{}.html'.format(name)))
    text = open(page).read()
    assert 'from the t(1;2) breakpoint; the rearrangement is present in B only' in text
    assert 'Promoter: structural variant nearby' in text
    outputs = egaReport.main(['-p', sv, '--explainPrefix', os.path.join(d, 'out', 'explain'), '--top', '1',
                              '--pdf', 'none', '-o', os.path.join(d, 'out', 'main')])
    main = open(outputs[0][0]).read()
    assert 'a structural variant near the promoter' in main.lower()
