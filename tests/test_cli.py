"""End-to-end runs of every command line tool on synthetic data."""

import json
import os

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer import (egaAnnotateGraph, egaBuildMultiway, egaCreateGraph, egaDifferential, egaExplain, egaIndirectContacts, egaPlotLocus, egaPlotNeighbourhood,
                                     egaValidateMultiway)
from epigeneticgraphanalyzer.graph import read_graph

from synthetic import (FRAGMENT, fragments, simulate_capture, simulate_gam, simulate_sprite, write_bigwig,
                       write_chrom_sizes, write_cool, write_tss_bed)


def test_full_pipeline(tmp_path):
    d = str(tmp_path)
    ibed, planted, baits = simulate_capture(d, seed=5)
    sizes = write_chrom_sizes(os.path.join(d, 'sizes.txt'))
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    table = pd.read_csv(ibed, sep='\t')

    prefix = os.path.join(d, 'graph', 'sample')
    egaCreateGraph.main(['-i', ibed, '-f', 'ibed', '--scoreThreshold', '5', '--tss', tss, '-o', prefix,
                         '--resolution', '1.0'])
    graph = read_graph(prefix + '.graphml')
    assert graph.number_of_edges() == len(table)
    assert sum(1 for _, _, d in graph.edges(data=True) if d['weight'] >= 7) >= len(planted)
    communities = pd.read_csv(prefix + '_communities.tsv', sep='\t')
    assert communities['n_nodes'].sum() == graph.number_of_nodes()
    assert any(graph.nodes[n]['genes'].startswith('GENE') for n in graph.nodes if graph.nodes[n]['bait'])

    egaIndirectContacts.main(['-g', prefix + '.graphml', '-o', prefix, '--triplets'])
    indirect = pd.read_csv(prefix + '_indirect.tsv', sep='\t')
    assert len(indirect) > 0
    assert set(indirect['pair_class']) <= {'oe-oe', 'bait-oe', 'bait-bait'}

    bigwig = write_bigwig(os.path.join(d, 'k27ac.bw'), {(c, i): 5.0 for c, i in baits})
    egaAnnotateGraph.main(['-g', prefix + '.graphml', '-t', 'H3K27ac=' + bigwig, '-o', prefix + '_annotated'])
    annotated = read_graph(prefix + '_annotated.graphml')
    bait_values = [annotated.nodes[n]['H3K27ac'] for n in annotated.nodes if annotated.nodes[n]['bait']]
    assert np.allclose(bait_values, 5.0)

    # multi-way data with the planted triplets of the first baits as hubs
    hubs = []
    for chrom, i in baits[:10]:
        partners = sorted(j for c, b, j in planted if c == chrom and b == i)
        if len(partners) >= 2:
            hubs.append((chrom, [i * FRAGMENT, partners[0] * FRAGMENT, partners[1] * FRAGMENT]))
    clusters = simulate_sprite(os.path.join(d, 'clusters.txt'), hubs, n_background=5000)
    sprite_npz = os.path.join(d, 'sprite.npz')
    egaBuildMultiway.main(['--type', 'sprite', '--input', clusters, '--chromSizes', sizes, '--binSize', '25000',
                           '-o', sprite_npz])
    gam = simulate_gam(os.path.join(d, 'gam.tsv'), hubs, n_tubes=100)
    gam_npz = os.path.join(d, 'gam.npz')
    egaBuildMultiway.main(['--type', 'gam', '--input', gam, '--chromSizes', sizes, '-o', gam_npz])
    for name, npz in (('sprite', sprite_npz), ('gam', gam_npz)):
        out = os.path.join(d, 'validation', name)
        egaValidateMultiway.main(['-g', prefix + '.graphml', '-m', npz, '-o', out, '--controls', '3'])
        summary = pd.read_csv(out + '_summary.tsv', sep='\t')
        assert {'direct', 'indirect'} <= set(summary['set'])
        assert os.path.exists(out + '_indirect_pairs.tsv.gz')

    chrom, positions = hubs[0]
    locus = os.path.join(d, 'locus.png')
    egaPlotLocus.main(['-g', prefix + '.graphml', '-m', sprite_npz, '--style', 'map', '--region',
                       '{}:{}-{}'.format(chrom, max(0, min(positions) - 200000), max(positions) + 200000), '-o', locus])
    assert os.path.getsize(locus) > 0
    community = max(range(graph.number_of_nodes()),
                    key=lambda c: sum(1 for _, x in graph.nodes(data='community') if x == c))
    egaPlotLocus.main(['-g', prefix + '.graphml', '-m', gam_npz, '--highlightCommunity', str(community),
                       '--shifts', '5', '--minLoci', '2', '-o', locus.replace('.png', '_gam.png')])

    egaPlotLocus.main(['-g', prefix + '.graphml', '-m', sprite_npz, '--style', 'communities', '--region',
                       '{}:0-3000000'.format(chrom), '--minLoci', '2', '-o', locus.replace('.png', '_communities.png')])
    assert os.path.getsize(locus.replace('.png', '_communities.png')) > 0

    gene = annotated.nodes[next(n for n in annotated.nodes if annotated.nodes[n]['genes'])]['genes'].split(',')[0]
    plot = os.path.join(d, 'plot.png')
    egaPlotNeighbourhood.main(['-g', prefix + '_annotated.graphml', '--gene', gene, '--colorBy', 'H3K27ac',
                               '-o', plot])
    assert os.path.getsize(plot) > 0


def test_differential_cli(tmp_path):
    d = str(tmp_path)
    da, db = os.path.join(d, 'a'), os.path.join(d, 'b')
    os.makedirs(da)
    os.makedirs(db)
    planted = sorted(simulate_capture(da, seed=7)[1])
    gained = {(c, i, j + 3) for c, i, j in planted[:20]}
    ibed_a = [simulate_capture(da, seed=s, planted=planted, n_noise=0, name='a{}'.format(s))[0] for s in (7, 9)]
    ibed_b, _, baits = simulate_capture(db, seed=8, planted=set(planted[20:]) | gained, n_noise=0)
    # Chicdiff results for the gained interactions (other-end regions of three fragments)
    rows = []
    for c, i, j in sorted(gained):
        bait, other = fragments(c)[i], fragments(c)[j]
        rows.append({'baitchr': c, 'baitstart': bait[1], 'baitend': bait[2], 'OEchr': c, 'OEstart': other[1] - FRAGMENT,
                     'OEend': other[2] + FRAGMENT, 'log2FoldChange': 3.0, 'pvalue': 1e-6, 'weighted_pvalue': 1e-6,
                     'weighted_padj': 1e-4})
    chicdiff = os.path.join(d, 'chicdiff_results.tsv')
    pd.DataFrame(rows).to_csv(chicdiff, sep='\t', index=False)
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    rng = np.random.default_rng(0)
    bw_a = write_bigwig(os.path.join(d, 'a.bw'), {})
    bw_b = write_bigwig(os.path.join(d, 'b.bw'), {(c, j): 4.0 for c, _, j in gained})
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)), 'B': rng.uniform(5, 50, len(baits))})
    expression_path = os.path.join(d, 'expression.tsv')
    expression.to_csv(expression_path, sep='\t', index=False)
    prefix = os.path.join(d, 'out', 'diff')
    egaDifferential.main(['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--chicdiff', chicdiff,
                          '--names', 'A', 'B',
                          '--tss', tss, '--tracksA', 'H3K27ac=' + bw_a, '--tracksB', 'H3K27ac=' + bw_b,
                          '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                          '--permutations', '5', '--bootstrap', '20', '-o', prefix])
    edges = pd.read_csv(prefix + '_edges.tsv.gz', sep='\t')
    assert {'gained', 'lost', 'shared'} & set(edges['status'])
    assert (edges.loc[edges['status'] == 'gained', 'log2FC'] > 0).mean() > 0.8
    planted_gained = {'{}:{}-{}|{}:{}-{}'.format(*fragments(c)[i], *fragments(c)[j]) for c, i, j in gained}
    is_gained = (edges['node_a'] + '|' + edges['node_b']).isin(planted_gained)
    assert is_gained.sum() == len(gained)
    assert edges.loc[is_gained, 'chicdiff_significant'].all()
    assert not edges.loc[edges['status'] == 'shared', 'chicdiff_significant'].any()
    features = pd.read_csv(prefix + '_gene_features.tsv', sep='\t')
    assert len(features) > 0
    with open(prefix + '_permutation.json') as handle:
        result = json.load(handle)
    assert result['n_permutations'] == 5
    assert os.path.exists(prefix + '_community_matching.tsv')
    cascade = pd.read_csv(prefix + '_explanation_cascade.tsv', sep='\t')
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    assert cascade.loc[cascade['level'].isin(['gene', 'step1', 'step2', 'step3', 'unexplained']), 'genes'].sum() == \
        len(classes)
    assert set(classes['explained_at']) <= {'gene', 'step1', 'step2', 'step3', 'unexplained'}
    for name in ('_mark_roles.tsv', '_state_transitions.tsv', '_bivalent_transitions.tsv', '_explanation_types.tsv',
                 '_state_thresholds.tsv', '_state_marks.tsv'):
        assert os.path.exists(prefix + name)
    # one active mark without a repressive one: no bivalent states
    assert set(classes['promoter_transition']) == {'states not callable'}
    # --markRoles gives the result of --markSign
    other = os.path.join(d, 'out', 'roles')
    egaDifferential.main(['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--chicdiff', chicdiff,
                          '--names', 'A', 'B', '--tss', tss, '--tracksA', 'H3K27ac=' + bw_a, '--tracksB',
                          'H3K27ac=' + bw_b, '--markSign', '--markRoles', 'H3K27ac=active',
                          '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                          '--permutations', '5', '--bootstrap', '20', '--roleBootstrap', '20', '-o', other])
    pd.testing.assert_frame_equal(classes, pd.read_csv(other + '_gene_classes.tsv', sep='\t'))
    # roles learned per level: long table, learned cascade and per-gene text
    for mode in ('learned', 'both'):
        learned = os.path.join(d, 'out', mode)
        egaDifferential.main(['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--names', 'A',
                              'B', '--tss', tss, '--tracksA', 'H3K27ac=' + bw_a, '--tracksB', 'H3K27ac=' + bw_b,
                              '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                              '--permutations', '5', '--bootstrap', '20', '--roleBootstrap', '20', '--roles', mode,
                              '-o', learned])
        roles = pd.read_csv(learned + '_mark_roles.tsv', sep='\t')
        assert {'promoter', 'step1', 'step2', 'step3', 'promoter_vs_contact'} <= set(roles['level'])
        assert 'contact' in set(roles['mark'])
        learned_classes = pd.read_csv(learned + '_gene_classes.tsv', sep='\t')
        assert learned_classes['explanation_detail'].astype(str).str.len().gt(0).all()
    assert os.path.exists(learned + '_explanation_cascade_learned.tsv')
    assert 'explanation_detail_learned' in learned_classes
    egaExplain.main(['-p', learned, '-o', os.path.join(d, 'out', 'explain')])
    report = pd.read_csv(os.path.join(d, 'out', 'explain_genes.tsv'), sep='\t')
    assert len(report) == len(learned_classes)
    # every level for every differential gene, both role sets and both presence statistics
    levels = pd.read_csv(learned + '_explanations_all_levels.tsv', sep='\t')
    assert len(levels) == 7 * len(learned_classes)   # --maxHops 3
    for column in ('explains_declared_mean', 'explains_declared_max1kb', 'explains_learned_mean',
                   'explains_learned_max1kb', 'p_declared', 'q_learned', 'detail_learned_max1kb'):
        assert column in levels
    assert {'first_level_declared_mean', 'explaining_levels_learned_max1kb'} <= set(learned_classes.columns)
    for name in ('_explanation_levels.tsv', '_explanation_combinations.tsv', '_explanation_cascade_max1kb.tsv',
                 '_explanation_cascade_learned_max1kb.tsv'):
        assert os.path.exists(learned + name)
    assert report['all_levels'].str.contains('promoter marks: .*step 3: ', regex=True).all()
    assert 'all_levels_learned' in report
    from epigeneticgraphanalyzer import egaPlotDifferential
    figure = os.path.join(d, 'out', 'differential.png')
    gene = 'GENE{}'.format(baits.index(next(b for b in baits if b[0] == sorted(gained)[0][0] and b[1] == sorted(gained)[0][1])))
    egaPlotDifferential.main(['-g', prefix + '_union.graphml', '--gene', gene, '--steps', '2', '--sideBySide',
                              '--candidates', prefix + '_indirect_only_genes.tsv', '-o', figure])
    assert os.path.getsize(figure) > 0


def test_graph_from_cool(tmp_path):
    d = str(tmp_path)
    ibed, planted, baits = simulate_capture(d, seed=3)
    cool = write_cool(os.path.join(d, 'sample.cool'), ibed)
    bait_bed = os.path.join(d, 'baits.bed')
    pd.DataFrame([fragments(c)[i] for c, i in baits]).to_csv(bait_bed, sep='\t', header=False, index=False)
    prefix = os.path.join(d, 'cool_graph')
    egaCreateGraph.main(['-i', cool, '-f', 'cool', '--baits', bait_bed, '--scoreThreshold', '10', '-o', prefix])
    graph = read_graph(prefix + '.graphml')
    table = pd.read_csv(ibed, sep='\t')
    assert graph.number_of_edges() == int((table['N_reads'] >= 10).sum())
    assert sum(1 for _, b in graph.nodes(data='bait') if b) > 0
