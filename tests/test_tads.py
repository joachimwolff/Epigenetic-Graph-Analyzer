"""TAD / boundary readers, crossing logic, rankings and the --tads displays (egaTadCrossing) on toy data."""

import os

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from epigeneticgraphanalyzer import egaExplain, egaReport, egaTadCrossing
from epigeneticgraphanalyzer import tads as td
from test_report import run  # noqa: F401  (module fixture: a synthetic egaDifferential run)


def write(path, text):
    with open(path, 'w') as handle:
        handle.write(text)
    return str(path)


# ----------------------------------------------------------------------------------------------------------------
# readers

def hicexplorer_files(d, chrom='6'):
    prefix = os.path.join(str(d), 'cond_25kb')
    write(prefix + '_boundaries.bed', '{c}\t975000\t1025000\tB00001\t-1.5\t.\n'
                                      '{c}\t1975000\t2025000\tB00002\t-0.2\t.\n'.format(c=chrom))
    write(prefix + '_domains.bed', '{c}\t0\t1000000\tID_0.01_1\t-1.5\t.\t0\t1000000\t31,120,180\n'
                                   '{c}\t1000000\t2000000\tID_0.01_2\t-0.2\t.\t1000000\t2000000\t51,160,44\n'
                                   '{c}\t2000000\t3000000\tID_0.01_3\t-0.2\t.\t2000000\t3000000\t31,120,180\n'
          .format(c=chrom))
    write(prefix + '_boundaries.gff', '{c}\tHiCExplorer\tboundary\t975000\t1025000\t-1.5\t.\t.\tID=B00001;delta=0.9;'
                                      'pvalue=0.001;tad_sep=-1.5\n'.format(c=chrom))
    write(prefix + '_score.bedgraph', ''.join('{}\t{}\t{}\t{}\n'.format(chrom, s, s + 25000, v) for s, v in
                                              zip(range(0, 3000000, 25000), np.sin(np.arange(120) / 6.0))))
    write(prefix + '_tad_score.bm', '#{"step":25000}\n' + ''.join(
        '{}\t{}\t{}\t{}\t{}\n'.format(chrom, s, s + 25000, 0.1, 0.3) for s in range(0, 3000000, 25000)))
    return prefix


def test_read_hicexplorer(tmp_path):
    prefix = hicexplorer_files(tmp_path)
    for path in (prefix, prefix + '_domains.bed', prefix + '_boundaries.bed', prefix + '_tad_score.bm'):
        assert td.detect_format(path) == 'hicexplorer'
        tads = td.read_tads(path)
        assert len(tads) == 2 and tads.fmt == 'hicexplorer'
        b = tads.boundaries
        assert list(b['pos']) == [1000000, 2000000]
        # strength = -TAD-separation score: the boundary at 1 Mb is the stronger one
        assert list(b['strength']) == [1.5, 0.2] and list(b['strength_rank']) == [1.0, 2.0]
        assert b['delta'].iloc[0] == 0.9 and np.isnan(b['delta'].iloc[1])
        assert len(tads.domains) == 3 and tads.track is not None


def test_read_hicexplorer_score_only(tmp_path):
    prefix = os.path.join(str(tmp_path), 'x')
    write(prefix + '_score.bedgraph', ''.join('chr1\t{}\t{}\t{}\n'.format(s, s + 10000, v) for s, v in
                                              zip(range(0, 1000000, 10000), np.cos(np.arange(100) / 5.0))))
    tads = td.read_tads(prefix + '_score.bedgraph')
    assert len(tads) >= 2  # local minima of the score
    assert (tads.boundaries['strength'] > 0).all()


def test_read_arrowhead(tmp_path):
    path = write(tmp_path / 'domains.txt',
                 'chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcorner_score\tUvar\tLvar\tUsign\tLsign\n'
                 '1\t100000\t500000\t1\t100000\t500000\t255,255,0\t0.9\t0.1\t0.1\t0.5\t0.5\n'
                 '1\t500000\t900000\t1\t500000\t900000\t255,255,0\t0.6\t0.1\t0.1\t0.5\t0.5\n'
                 '1\t100000\t900000\t1\t100000\t900000\t255,255,0\t0.3\t0.1\t0.1\t0.5\t0.5\n')
    assert td.detect_format(path) == 'arrowhead'
    tads = td.read_tads(path, targets={'chr1'})
    b = tads.boundaries
    assert set(b['chrom']) == {'chr1'}
    assert list(b['pos']) == [100000, 500000, 900000]
    # strength: the largest corner score of the domains ending at the boundary
    assert list(b['strength']) == [0.9, 0.9, 0.6]
    assert len(tads.domains) == 3


def test_read_cooltools(tmp_path):
    rows = []
    for k, s in enumerate(range(0, 1000000, 50000)):
        rows.append(dict(chrom='chr2', start=s, end=s + 50000, is_bad_bin=False,
                         log2_insulation_score_100000=-1.0 if s == 500000 else 0.2,
                         boundary_strength_100000=0.8 if s == 500000 else np.nan,
                         is_boundary_100000=s == 500000,
                         log2_insulation_score_200000=0.0, boundary_strength_200000=np.nan,
                         is_boundary_200000=s in (300000, 700000)))
    path = str(tmp_path / 'ins.tsv')
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    assert td.detect_format(path) == 'cooltools'
    small = td.read_tads(path, window=100000)
    assert list(small.boundaries['pos']) == [525000] and small.boundaries['strength'].iloc[0] == 0.8
    large = td.read_tads(path)  # default: the largest window
    assert list(large.boundaries['pos']) == [325000, 725000]
    with pytest.raises(ValueError):
        td.read_tads(path, window=12345)


def test_read_bed(tmp_path):
    boundaries = write(tmp_path / 'b.bed', 'chr1\t990000\t1010000\tb1\t3.0\nchr1\t1990000\t2010000\tb2\t1.0\n')
    domains = write(tmp_path / 'd.bed', 'chr1\t0\t1000000\nchr1\t1000000\t2000000\n')
    assert td.detect_format(boundaries) == 'bed'
    b = td.read_tads(boundaries)
    assert b.fmt == 'bed-boundaries' and list(b.boundaries['strength']) == [3.0, 1.0] and b.domains is None
    d = td.read_tads(domains)
    assert d.fmt == 'bed-domains' and list(d.boundaries['pos']) == [0, 1000000, 2000000]
    forced = td.read_tads(domains, fmt='bed-boundaries')
    assert list(forced.boundaries['pos']) == [500000, 1500000]


def test_read_ega(tmp_path):
    prefix = str(tmp_path / 'ega')
    bins = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(0, 500000, 25000)})
    bins['end'] = bins['start'] + 25000
    bins['capture_insulation'] = 0.0
    bins['score'] = np.where(bins['start'] == 250000, 0.7, np.nan)
    bins['boundary'] = bins['start'] == 250000
    bins.to_csv(prefix + '_bins.tsv', sep='\t', index=False)
    write(prefix + '_boundaries.bed', 'chr1\t250000\t275000\tboundary\t0.7\n')
    write(prefix + '_domains.bed', 'chr1\t0\t262500\nchr1\t262500\t500000\n')
    for path in (prefix, prefix + '_bins.tsv', prefix + '_boundaries.bed'):
        assert td.detect_format(path) == 'ega'
        tads = td.read_tads(path)
        assert list(tads.boundaries['pos']) == [262500] and tads.boundaries['strength'].iloc[0] == 0.7


# ----------------------------------------------------------------------------------------------------------------
# chromosome names, liftover

def test_chrom_harmonisation(tmp_path):
    prefix = hicexplorer_files(tmp_path, chrom='6')
    tads = td.read_tads(prefix, targets={'chr6', 'chr7'})
    assert set(tads.boundaries['chrom']) == {'chr6'} and set(tads.domains['chrom']) == {'chr6'}
    mapper = td.chrom_mapper({'6', 'X', 'MT'})
    assert mapper('chr6') == '6' and mapper('chrX') == 'X' and mapper('chrM') == 'MT' and mapper('chr9') is None
    assert td.chrom_mapper({'chrM'})('MT') == 'chrM'
    assert td.chrom_mapper()('6') == 'chr6'
    # chromosomes absent from the graph are dropped
    assert len(td.read_tads(prefix, targets={'chr1'})) == 0


def test_liftover(tmp_path):
    chain = write(tmp_path / 'toy.chain', 'chain 1000 chr6 4000000 + 0 4000000 chr6 4000000 + 100 4000100 1\n'
                                          '4000000\n\n')
    prefix = hicexplorer_files(tmp_path, chrom='6')
    tads = td.read_tads(prefix, chain=chain, targets={'chr6'})
    assert list(tads.boundaries['pos']) == [1000100, 2000100]
    assert list(tads.domains['start']) == [100, 1000100, 2000100]


# ----------------------------------------------------------------------------------------------------------------
# crossing

def toy_tads(domains=True):
    b = pd.DataFrame({'chrom': 'chr1', 'start': [990000, 1990000], 'end': [1010000, 2010000],
                      'strength': [2.0, 0.5]})
    d = pd.DataFrame({'chrom': 'chr1', 'start': [0, 1000000, 2000000], 'end': [1000000, 2000000, 2800000]}) \
        if domains else None
    return td.Tads(b, d)


def edges_frame(pairs, status=None):
    return pd.DataFrame({'node_a': [a for a, _ in pairs], 'node_b': [b for _, b in pairs],
                         'status': status or ['shared'] * len(pairs)})


@pytest.mark.parametrize('domains', [True, False])
def test_edge_crossings(domains):
    tads = toy_tads(domains)
    edges = edges_frame([('chr1:500000-505000', 'chr1:1500000-1505000'),     # spans the boundary at 1 Mb
                         ('chr1:1100000-1105000', 'chr1:1800000-1805000'),   # same domain
                         ('chr1:2600000-2605000', 'chr1:500000-505000'),     # spans both, reversed order
                         ('chr1:2900000-2905000', 'chr1:2100000-2105000'),   # one end outside the domains
                         ('chr1:100000-105000', 'chr2:100000-105000')])      # trans
    out = td.edge_crossings(edges, tads)
    assert list(out['crosses']) == [True, False, True, False, False]
    assert list(out['n_boundaries']) == [1, 0, 2, 0, 0]
    assert out['strongest_boundary'].iloc[2] == 'chr1:1000000'
    assert out['boundaries'].iloc[2] == 'chr1:1000000,chr1:2000000'
    relation = list(out['domain_relation'])
    assert relation[:3] == ['cross domain', 'same domain', 'cross domain'] and relation[4] == 'trans'
    assert relation[3] == ('outside domains' if domains else 'same domain')
    # strength filter: only the boundary at 1 Mb (above the median) counts
    strong = td.edge_crossings(edges, tads, min_strength_pct=0.75)
    assert list(strong['n_boundaries']) == [1, 0, 1, 0, 0]


def test_edge_on_boundary_midpoint_not_crossing():
    tads = toy_tads()
    out = td.edge_crossings(edges_frame([('chr1:999000-1001000', 'chr1:1500000-1505000')]), tads)
    assert not out['crosses'].iloc[0]


def test_nested_domains():
    b = pd.DataFrame({'chrom': 'chr1', 'start': [100000, 500000, 900000], 'end': [100001, 500001, 900001]})
    d = pd.DataFrame({'chrom': 'chr1', 'start': [100000, 500000, 100000], 'end': [500000, 900000, 900000]})
    out = td.edge_crossings(edges_frame([('chr1:200000-200001', 'chr1:700000-700001')]), td.Tads(b, d))
    # crosses the nested boundary at 500 kb but both ends lie in the outer domain
    assert out['crosses'].iloc[0] and out['domain_relation'].iloc[0] == 'same domain'


def test_sum_binomial_sf():
    p = td.sum_binomial_sf(5, [10], [0.2])
    assert np.isclose(p, stats.binom.sf(4, 10, 0.2))
    p2 = td.sum_binomial_sf(3, [4, 6], [0.1, 0.3])
    sim = stats.binom.pmf(np.arange(5)[:, None], 4, 0.1) * stats.binom.pmf(np.arange(7)[None, :], 6, 0.3)
    total = np.add.outer(np.arange(5), np.arange(7))
    assert np.isclose(p2, sim[total >= 3].sum())
    assert td.sum_binomial_sf(0, [3], [0.5]) == 1.0


def test_rankings():
    tads = toy_tads()
    bait = 'chr1:900000-905000'
    pairs, status = [], []
    for k in range(8):  # gained contacts across the boundary at 1 Mb
        pairs.append((bait, 'chr1:{}-{}'.format(1200000 + k * 20000, 1205000 + k * 20000)))
        status.append('gained')
    for k in range(8):  # lost contacts inside its domain, of the same lengths
        pairs.append((bait, 'chr1:{}-{}'.format(600000 - k * 20000, 605000 - k * 20000)))
        status.append('lost')
    other = 'chr1:2500000-2505000'
    for k in range(6):
        pairs.append((other, 'chr1:{}-{}'.format(2100000 + k * 50000, 2105000 + k * 50000)))
        status.append('shared')
    ann = td.annotate_edges(edges_frame(pairs, status), tads, names=['A', 'B'])
    promoters = td.promoter_ranking(ann, {bait, other}, {bait: 'G1', other: 'G2'})
    top = promoters.iloc[0]
    assert top['node'] == bait and top['genes'] == 'G1'
    assert top['gained_crossing'] == 8 and top['lost_within'] == 8
    assert np.isclose(top['gained_p'], 1 / stats.binom(16, 0.5).pmf(8) / 2 ** 16)  # 1 / C(16, 8)
    assert np.isnan(promoters.set_index('node').loc[other, 'gained_p'])
    bounds = td.boundary_ranking(ann, tads).set_index('id')
    assert bounds.loc['chr1:1000000', 'gained_crossing'] == 8
    assert bounds.loc['chr1:1000000', 'gained_p'] < 0.05
    summary = td.condition_summary(ann, ['A', 'B']).set_index('group')
    assert summary.loc['gained', 'crossing'] == 8 and summary.loc['lost', 'crossing'] == 0
    # presence per condition from the status when no significant_ columns exist
    assert ann['in_A'].sum() == 14 and ann['in_B'].sum() == 14


def test_gene_sentence():
    tads = toy_tads()
    bait = 'chr1:900000-905000'
    pairs = [(bait, 'chr1:{}-{}'.format(1200000 + k * 20000, 1205000 + k * 20000)) for k in range(4)]
    pairs += [(bait, 'chr1:{}-{}'.format(100000 + k * 50000, 105000 + k * 50000)) for k in range(3)]
    ann = td.annotate_edges(edges_frame(pairs, ['gained'] * 4 + ['lost'] * 3), tads, names=['A', 'B'])
    regions = pd.DataFrame({'chrom': ['chr1'], 'start': [1190000], 'end': [1250000], 'name': ['ENH']})
    crossing = td.GeneCrossing(ann, tads, ['A', 'B'], regions=regions)
    row = crossing.gene_row(bait, 'G1')
    text = row['sentence']
    assert text.startswith('G1: in B, 4 of its 4 gained contacts cross a boundary')
    assert 'chr1:1.000 Mb' in text and 'strength rank 1 of 2' in text and 'ENH: 3 partners' in text
    assert 'its 3 lost contacts (A only) stay within its own domain' in text
    assert row['gained_crossing'] == 4 and row['lost_within'] == 3


# ----------------------------------------------------------------------------------------------------------------
# command line and displays on a synthetic run

@pytest.fixture(scope='module')
def tads_file(run, tmp_path_factory):  # noqa: F811
    d = tmp_path_factory.mktemp('tads')
    # one boundary per Mb on chr1 / chr2 of the synthetic genome, Ensembl names ('1') to exercise harmonisation
    text = ''.join('{}\t{}\t{}\tB{}\t{}\t.\n'.format(c, p - 12500, p + 12500, k, -1.0 - k / 10.0)
                   for k, (c, p) in enumerate([('1', 1000000), ('1', 2000000), ('1', 3000000), ('2', 1500000)]))
    return write(d / 'toy_boundaries.bed', text)


def test_tad_crossing_cli(run, tads_file, tmp_path):  # noqa: F811
    out = str(tmp_path / 'tc')
    analysis = egaTadCrossing.main(['-p', run['prefix'], '--tads', tads_file, '--tadFormat', 'bed-boundaries',
                                    '-o', out])
    for suffix in ('_edges.tsv.gz', '_summary.tsv', '_promoters.tsv', '_boundaries.tsv', '_genes.tsv'):
        assert os.path.exists(out + suffix), suffix
    edges = pd.read_csv(out + '_edges.tsv.gz', sep='\t')
    union = pd.read_csv(run['prefix'] + '_edges.tsv.gz', sep='\t')
    assert len(edges) == len(union) and edges['crosses'].any() and (~edges['crosses']).any()
    assert set(analysis.tads.boundaries['chrom']) == {'chr1', 'chr2'}
    genes = pd.read_csv(out + '_genes.tsv', sep='\t')
    assert genes['sentence'].fillna('').str.len().gt(0).any()


def test_explain_backward_compatible(run, tads_file, tmp_path):  # noqa: F811
    plain = str(tmp_path / 'plain')
    egaExplain.main(['-p', run['prefix'], '-o', plain])
    base = pd.read_csv(plain + '_genes.tsv', sep='\t')
    reference = pd.read_csv(os.path.join(os.path.dirname(run['prefix']), 'explain', 'cond1_vs_cond2_genes.tsv'),
                            sep='\t')
    pd.testing.assert_frame_equal(base, reference)
    assert 'tad_crossing' not in base
    with_tads = str(tmp_path / 'tads')
    egaExplain.main(['-p', run['prefix'], '--tads', tads_file, '-o', with_tads])
    extended = pd.read_csv(with_tads + '_genes.tsv', sep='\t')
    assert {'tad_crossing', 'gained_crossing', 'lost_crossing', 'shared_crossing'} <= set(extended)
    pd.testing.assert_frame_equal(extended[[c for c in base if c != 'summary']],
                                  base[[c for c in base if c != 'summary']])


def test_report_with_tads(run, tads_file, tmp_path):  # noqa: F811
    classes = run['classes']
    gene = str(classes['genes'].iloc[0]).split(',')[0]
    plain = egaReport.main(['-p', run['prefix'], '--gene', gene, '--pdf', 'none', '-o', str(tmp_path / 'a')])
    assert 'id="boundaries"' not in open(plain[0][0]).read()
    out = egaReport.main(['-p', run['prefix'], '--gene', gene, '--pdf', 'none', '--tads', tads_file,
                          '-o', str(tmp_path / 'b')])
    text = open(out[0][0]).read()
    assert 'id="boundaries"' in text and 'Domain boundaries' in text
    assert 'Boundaries: toy (hicexplorer, 4 boundaries' in text
