"""eQTL mapping, control matching and graph path features on toy data."""

import gzip
import os

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer import egaEqtl
from epigeneticgraphanalyzer.eqtl import (FragmentGraph, FragmentMap, PairSet, add_features, distance_bin,
                                          gene_promoters, lead_per_gene, map_pairs, odds_table, overlay_nodes,
                                          path_features, sample_controls, significant_pairs, tile_bins,
                                          variant_fragments)

SIZE = 10000


def toy_map(n=100, chroms=('chr1', 'chr2')):
    """Fragments of 10 kb, 1-based closed: fragment k of chr1 is [k*10000+1, (k+1)*10000]."""
    return FragmentMap(tile_bins({c: n * SIZE for c in chroms}, SIZE))


def toy_tss(genes):
    """genes: {gene_id: (chrom, [0-based TSS])} as io.read_tss_from_gtf returns them."""
    rows = [(c, t, '+', g, g, 'protein_coding') for g, (c, ts) in genes.items() for t in ts]
    return pd.DataFrame(rows, columns=['chrom', 'tss', 'strand', 'gene_id', 'gene_name', 'gene_type'])


def frag_pos(k, offset=5000):
    return k * SIZE + offset


def test_fragment_map_locates_one_based_closed():
    fmap = toy_map()
    assert fmap.locate(['chr1'], [1])[0] == 0
    assert fmap.locate(['chr1'], [10000])[0] == 0
    assert fmap.locate(['chr1'], [10001])[0] == 1
    assert fmap.locate(['chr2'], [5])[0] == 100
    assert fmap.locate(['chr3'], [5])[0] == -1
    assert fmap.locate(['chr1'], [100 * SIZE + 1])[0] == -1


def test_overlay_nodes_replaces_overlapped_fragments():
    fragments = tile_bins({'chr1': 100000}, SIZE)
    nodes = pd.DataFrame({'chrom': ['chr1'], 'start': [15001], 'end': [32000]})
    fmap = FragmentMap(overlay_nodes(fragments, nodes))
    # fragments 2 (10001-20000), 3, 4 (30001-40000) overlap the node and are replaced by it
    assert fmap.n == 10 - 3 + 1
    k = fmap.locate(['chr1'], [20000])[0]
    assert (fmap.start[k], fmap.end[k]) == (15001, 32000)
    assert fmap.locate(['chr1'], [12000])[0] == -1  # the uncovered rest of a replaced fragment


def test_gene_promoters_keep_baited_tss_only():
    fmap = toy_map()
    bait = np.zeros(fmap.n, dtype=bool)
    bait[[10, 12]] = True
    tss = toy_tss({'A': ('chr1', [frag_pos(10), frag_pos(12), frag_pos(20)]), 'B': ('chr1', [frag_pos(30)])})
    promoters = gene_promoters(tss, fmap, bait)
    assert set(promoters) == {'A'}
    chrom, positions, fragments = promoters['A']
    assert chrom == 'chr1' and list(fragments) == [10, 12]
    assert list(positions) == [frag_pos(10) + 1, frag_pos(12) + 1]


def test_map_pairs_exclusions_and_distance():
    fmap = toy_map()
    bait = np.zeros(fmap.n, dtype=bool)
    bait[[10, 50]] = True
    tss = toy_tss({'A': ('chr1', [frag_pos(10)]), 'B': ('chr1', [frag_pos(50)]), 'C': ('chr1', [frag_pos(70)])})
    promoters = gene_promoters(tss, fmap, bait)
    eqtls = pd.DataFrame([
        ('chr1', frag_pos(13), 'A'),          # kept, 30 kb
        ('chr1', frag_pos(13, 6000), 'A'),     # same fragment and gene: duplicate
        ('chr1', frag_pos(10, 100), 'A'),      # in the promoter fragment
        ('chr1', frag_pos(70), 'C'),           # gene without baited promoter
        ('chr1', frag_pos(20), 'Z'),           # gene not annotated
        ('chr2', frag_pos(5), 'B'),            # other chromosome
        ('chrX', 100, 'A'),                    # not autosomal
        ('chr1', 100 * SIZE + 50, 'A'),        # beyond the map
        ('chr1', frag_pos(99), 'A'),           # beyond max distance
    ], columns=['chrom', 'pos', 'gene_id'])
    pairs, counts = map_pairs(eqtls, fmap, promoters, all_genes={'A', 'B', 'C'}, max_distance=500000)
    counts = dict(zip(counts['step'], counts['n']))
    assert len(pairs) == 1
    row = pairs.iloc[0]
    assert (row['var_frag'], row['gene_id'], row['distance']) == (13, 'A', 29999)
    assert counts['variant_in_promoter_fragment'] == 1
    assert counts['gene_promoter_not_baited'] == 1
    assert counts['gene_not_annotated'] == 1
    assert counts['other_chromosome'] == 1
    assert counts['not_autosomal'] == 1
    assert counts['variant_not_on_fragment'] == 1
    assert counts['beyond_max_distance'] == 1
    assert counts['duplicate_fragment_gene'] == 1


def test_path_features_and_mirror():
    # 0 - 1 - 2 - 3 - 4 ; 5 isolated ; 6 - 7 in another component
    graph = FragmentGraph(8, [0, 1, 2, 3, 6], [1, 2, 3, 4, 7])
    assert path_features(graph, 0, [1]) == (True, False, False, '1')
    assert path_features(graph, 0, [2]) == (False, True, False, '2')
    assert path_features(graph, 0, [3]) == (False, False, True, '3')
    assert path_features(graph, 0, [4]) == (False, False, False, '>3')
    assert path_features(graph, 0, [6]) == (False, False, False, 'none')
    assert path_features(graph, 5, [0]) == (False, False, False, 'none')
    # several promoter fragments: the closest counts
    assert path_features(graph, 0, [4, 2])[3] == '2'
    # a triangle has a direct edge and a two-step path; a square adds three-step paths
    tri = FragmentGraph(3, [0, 1, 0], [1, 2, 2])
    assert path_features(tri, 0, [2]) == (True, True, False, '1')
    square = FragmentGraph(4, [0, 1, 2, 3], [1, 2, 3, 0])
    assert path_features(square, 0, [1]) == (True, False, True, '1')

    fmap = toy_map(n=20, chroms=('chr1',))
    bait = np.zeros(fmap.n, dtype=bool)
    bait[10] = True
    graph = FragmentGraph(fmap.n, [10, 10], [13, 8])
    mirrored = graph.mirrored(fmap, bait)
    assert sorted(mirrored.neighbours(10).tolist()) == [7, 12]
    assert mirrored.degree[10] == 2 and mirrored.n_edges == 2


def test_pairset():
    s = PairSet([1, 1, 2], ['A', 'B', 'A'])
    assert list(s.contains([1, 2, 2, 3], ['B', 'A', 'B', 'Q'])) == [True, True, False, False]


def control_setup():
    fmap = toy_map()
    baited = [10, 20, 30, 40, 60]
    bait = np.zeros(fmap.n, dtype=bool)
    bait[baited] = True
    genes = {'G{}'.format(k): ('chr1', [frag_pos(k)]) for k in baited}
    promoters = gene_promoters(toy_tss(genes), fmap, bait)
    # eQTL: variant on fragment 25 for G20 (50 kb); significant also: fragment 25 - G30 (eGene of the fragment),
    # fragment 15 - G20 (LD proxy); other genes' variants on fragments 15, 35, 45 and 75
    leads = pd.DataFrame([('chr1', frag_pos(25), 'G20')], columns=['chrom', 'pos', 'gene_id'])
    significant = pd.DataFrame([('chr1', frag_pos(25), 'G20'), ('chr1', frag_pos(25), 'G30'),
                                ('chr1', frag_pos(15), 'G20'), ('chr1', frag_pos(15), 'G10'),
                                ('chr1', frag_pos(35), 'G40'), ('chr1', frag_pos(45), 'G40'),
                                ('chr1', frag_pos(75), 'G60')], columns=['chrom', 'pos', 'gene_id'])
    return fmap, bait, promoters, leads, significant


def test_sample_controls_match_distance_and_exclude_eqtls():
    fmap, bait, promoters, leads, significant = control_setup()
    pairs, _ = map_pairs(leads, fmap, promoters)
    assert pairs['distance'].tolist() == [49999]
    sig = significant_pairs(significant, fmap)
    universe = sorted(set(significant['gene_id']) & set(promoters))
    candidates = variant_fragments(significant, fmap)
    table = sample_controls(pairs, promoters, universe, candidates, sig, n_controls=5, window=0.1)
    gene_controls = table[table['group'] == 'control_gene']
    # genes at 50 kb +- 26% from fragment 25: G30 (50 kb, but an eGene of fragment 25 -> excluded),
    # G10 (150 kb), G40 (150 kb), G60 (350 kb) are out of the window -> no gene control except none
    assert gene_controls.empty
    frag_controls = table[table['group'] == 'control_fragment']
    # fragments with significant variants at ~50 kb from G20's TSS: 15 (LD proxy of G20 -> excluded), 25 (the
    # eQTL fragment itself -> excluded) -> none; 35 is 150 kb away
    assert frag_controls.empty

    # widen the window: G30 stays excluded, G10 and G40 (150 kb) are accepted; fragment 35 and 45 too
    table = sample_controls(pairs, promoters, universe, candidates, sig, n_controls=5, window=0.5)
    gene_controls = table[table['group'] == 'control_gene']
    assert sorted(gene_controls['gene_id']) == ['G10', 'G40']
    assert (gene_controls['var_frag'] == 25).all()
    frag_controls = table[table['group'] == 'control_fragment']
    assert sorted(frag_controls['var_frag']) == [35]
    assert (frag_controls['gene_id'] == 'G20').all()
    assert table.loc[table['group'] == 'eqtl', 'set_id'].tolist() == [0]
    # at most n_controls per type
    table = sample_controls(pairs, promoters, universe, candidates, sig, n_controls=1, window=0.5, seed=3)
    assert (table['group'] == 'control_gene').sum() == 1


def test_features_and_odds_on_planted_contacts():
    """eQTL pairs with a two-step path, controls without: the two-step odds ratio is large, direct is null."""
    fmap = toy_map(n=400, chroms=('chr1', 'chr2', 'chr3'))
    bait = np.zeros(fmap.n, dtype=bool)
    rows, a, b = [], [], []
    promoters = {}
    for c, chrom in enumerate(('chr1', 'chr2', 'chr3')):
        for k in range(20, 380, 40):
            gene = '{}_{}'.format(chrom, k)
            g = c * 400 + k
            bait[g] = True
            promoters[gene] = (chrom, np.array([frag_pos(k)]), np.array([g]))
            # eQTL fragment k+5 connects over mediator k+2; control fragment k-5 has no path
            a += [g, g + 2]
            b += [g + 2, g + 5]
            rows.append((k // 40, 'eqtl', chrom, g + 5, frag_pos(k + 5), gene, 50000))
            rows.append((k // 40, 'control_fragment', chrom, g - 5, frag_pos(k - 5), gene, 50000))
    graph = FragmentGraph(fmap.n, a, b)
    table = pd.DataFrame(rows, columns=['set_id', 'group', 'chrom', 'var_frag', 'var_pos', 'gene_id', 'distance'])
    features = add_features(table, graph, promoters, fmap, bait)
    assert (features.loc[features['group'] == 'eqtl', 'shortest'] == '2').all()
    assert (features.loc[features['group'] != 'eqtl', 'shortest'] == 'none').all()
    assert (features['distance_bin'] == distance_bin([50000])[0]).all()
    odds = odds_table(features, n_boot=50)
    two = odds[(odds['control'] == 'control') & (odds['contrast'] == 'two_step_any') &
               (odds['distance_bin'] == 'all')].iloc[0]
    assert two['eqtl_fraction'] == 1.0 and two['control_fraction'] == 0.0 and two['odds_ratio'] > 100
    direct = odds[(odds['control'] == 'control') & (odds['contrast'] == 'direct') &
                  (odds['distance_bin'] == 'all')].iloc[0]
    assert direct['odds_ratio'] == 1.0


def test_lead_per_gene():
    df = pd.DataFrame({'chrom': 'chr1', 'pos': [1, 2, 3], 'gene_id': ['A', 'A', 'B'], 'pvalue': [1e-5, 1e-9, 0.01],
                       'abs_z': [4, 6, 2]})
    lead = lead_per_gene(df)
    assert dict(zip(lead['gene_id'], lead['pos'])) == {'A': 2, 'B': 3}


def test_cli_end_to_end(tmp_path):
    d = str(tmp_path)
    rng = np.random.default_rng(1)
    chroms = ['chr{}'.format(i) for i in range(1, 5)]
    n = 300
    with open(os.path.join(d, 'frag.rmap'), 'w') as handle:
        for c in chroms:
            for k in range(n):
                handle.write('{}\t{}\t{}\t{}\n'.format(c[3:], k * SIZE + 1, (k + 1) * SIZE, k))
    gtf_lines, baits, nodes, edges, leads, sig = [], [], {}, [], [], []
    for c in chroms:
        for k in range(10, n - 10, 12):
            gene = 'ENSG{}{:04d}'.format(c[3:], k)
            tss = k * SIZE + 5000
            gtf_lines.append('{}\tx\ttranscript\t{}\t{}\t.\t+\t.\tgene_id "{}.1"; gene_type "protein_coding"; '
                             'gene_name "{}";'.format(c, tss, tss + 1000, gene, gene))
            baits.append((c, k))
            # eQTL fragment k+3..k+5 (30-50 kb) with a direct edge for half of the genes
            step = int(rng.integers(3, 6))
            leads.append((c, (k + step) * SIZE + 100, gene))
            sig.append((c, (k + step) * SIZE + 200, gene))
            if rng.random() < 0.5:
                edges.append(((c, k), (c, k + step)))
            if rng.random() < 0.3:
                edges.append(((c, k), (c, k - 3)))
    with gzip.open(os.path.join(d, 'genes.gtf.gz'), 'wt') as handle:
        handle.write('\n'.join(gtf_lines) + '\n')
    with open(os.path.join(d, 'baits.baitmap'), 'w') as handle:
        for c, k in baits:
            handle.write('{}\t{}\t{}\t{}\tg\n'.format(c[3:], k * SIZE + 1, (k + 1) * SIZE, k))

    def node(x):
        return '{}:{}-{}'.format(x[0], x[1] * SIZE + 1, (x[1] + 1) * SIZE)
    for x, y in edges:
        nodes[node(x)] = (x[0], x[1] * SIZE + 1, (x[1] + 1) * SIZE, (x[0], x[1]) in set(baits))
        nodes[node(y)] = (y[0], y[1] * SIZE + 1, (y[1] + 1) * SIZE, (y[0], y[1]) in set(baits))
    pd.DataFrame([(k,) + v for k, v in nodes.items()], columns=['node', 'chrom', 'start', 'end', 'bait']).to_csv(
        os.path.join(d, 'g_nodes.tsv'), sep='\t', index=False)
    pd.DataFrame([(node(x), node(y), 10.0) for x, y in edges], columns=['node_a', 'node_b', 'weight']).to_csv(
        os.path.join(d, 'g_edges.tsv'), sep='\t', index=False)
    pd.DataFrame(leads, columns=['chrom', 'pos', 'gene_id']).to_csv(os.path.join(d, 'leads.tsv'), sep='\t',
                                                                     index=False)
    pd.DataFrame(sig, columns=['chrom', 'pos', 'gene_id']).to_csv(os.path.join(d, 'sig.tsv'), sep='\t', index=False)
    prefix = os.path.join(d, 'out', 'toy')
    egaEqtl.main(['--leads', os.path.join(d, 'leads.tsv'), '--significant', os.path.join(d, 'sig.tsv'),
                  '--nodes', os.path.join(d, 'g_nodes.tsv'), '--edges', os.path.join(d, 'g_edges.tsv'),
                  '--rmap', os.path.join(d, 'frag.rmap'), '--baits', os.path.join(d, 'baits.baitmap'),
                  '--gtf', os.path.join(d, 'genes.gtf.gz'), '--controls', '3', '--window', '0.5',
                  '--bootstrap', '50', '--modelBootstrap', '10', '--mirror', '-o', prefix])
    pairs = pd.read_csv(prefix + '_pairs.tsv.gz', sep='\t')
    observed = pairs[pairs['graph'] == 'observed']
    eqtl = observed[observed['group'] == 'eqtl']
    assert len(eqtl) == len(leads)
    assert 0.3 < eqtl['direct'].mean() < 0.7
    # controls are other variant-bearing fragments of the same gene or other genes of the fragment; never eQTLs
    assert observed.loc[observed['group'] != 'eqtl', 'direct'].mean() < eqtl['direct'].mean()
    odds = pd.read_csv(prefix + '_odds.tsv', sep='\t')
    assert set(odds['graph']) == {'observed', 'mirrored'}
    models = pd.read_csv(prefix + '_models.tsv', sep='\t')
    assert {'all_pairs', 'without_direct'} <= set(models['model'])
    mapping = pd.read_csv(prefix + '_mapping.tsv', sep='\t')
    assert dict(zip(mapping['step'], mapping['n']))['pairs'] == len(leads)


def test_sample_controls_match_fragment_length_and_bait_status():
    fmap, bait, promoters, leads, significant = control_setup()
    pairs, _ = map_pairs(leads, fmap, promoters)
    sig = significant_pairs(significant, fmap)
    universe = sorted(set(significant['gene_id']) & set(promoters))
    candidates = variant_fragments(significant, fmap)
    length = fmap.length(np.arange(fmap.n)).astype(float)
    table = sample_controls(pairs, promoters, universe, candidates, sig, window=0.5, fragment_length=length,
                            bait=bait)
    assert sorted(table.loc[table['group'] == 'control_fragment', 'var_frag']) == [35]
    # a control fragment twice as long as the eQTL fragment is not a match
    longer = length.copy()
    longer[35] = 2 * length[25]
    table = sample_controls(pairs, promoters, universe, candidates, sig, window=0.5, fragment_length=longer,
                            length_window=0.15)
    assert (table['group'] == 'control_fragment').sum() == 0
    # nor is a baited fragment for an eQTL fragment without bait
    other = bait.copy()
    other[35] = True
    table = sample_controls(pairs, promoters, universe, candidates, sig, window=0.5, bait=other)
    assert (table['group'] == 'control_fragment').sum() == 0
    # gene controls share the eQTL fragment and are unaffected
    assert sorted(table.loc[table['group'] == 'control_gene', 'gene_id']) == ['G10', 'G40']
