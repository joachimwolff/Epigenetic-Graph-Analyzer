"""Compartment scores, eigenvectors and domain statistics on toy data."""

import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import compartments as cp


def _toy_compartments(n=120, block=10, seed=0):
    """Alternating A/B blocks; returns the +1/-1 label per bin."""
    rng = np.random.default_rng(seed)
    sizes = rng.integers(block // 2, block * 2, size=n)
    label = np.concatenate([np.full(s, 1.0 if i % 2 == 0 else -1.0) for i, s in enumerate(sizes)])[:n]
    return label


def _toy_hic(label, seed=0, decay=1.0):
    """Contact matrix with power-law decay and a compartment term: same compartment contacts x3."""
    rng = np.random.default_rng(seed)
    n = len(label)
    i, j = np.indices((n, n))
    dist = np.abs(i - j) + 1.0
    lam = 200.0 / dist ** decay * np.where(np.outer(label, label) > 0, 3.0, 1.0)
    m = rng.poisson(lam).astype(float)
    m = np.triu(m) + np.triu(m, 1).T
    return m


def test_log_enrichment_and_density():
    x = np.array([0.0, 1.0, 2.0, 1.0])
    e = cp.log_enrichment(x)
    assert e[1] == pytest.approx(0.0) and e[2] > 0 > e[0]
    bins = cp.make_bins({'chr1': 250, 'chr2': 100}, 100)
    assert list(bins.end) == [100, 200, 250, 100]
    d = cp.interval_density(bins, ['chr1', 'chr1', 'chr1', 'chr2'], [5, 150, 160, 99])
    assert list(d) == [1, 2, 0, 1]


def test_chromatin_pc1_orientation_and_role_score():
    rng = np.random.default_rng(1)
    label = _toy_compartments(300)
    active = label[:, None] * np.array([1.0, 0.8]) + rng.normal(0, 0.5, (300, 2))
    repressive = -0.6 * label[:, None] + rng.normal(0, 0.5, (300, 1))
    feats = np.hstack([active, repressive])
    score, loadings, explained = cp.chromatin_pc1(feats, active=np.array([True, True, False]))
    assert np.corrcoef(score, label)[0, 1] > 0.8          # oriented: positive = A
    assert loadings[0] > 0 and loadings[2] < 0 and 0 < explained <= 1
    # the sign must not depend on the column order
    flipped, _, _ = cp.chromatin_pc1(feats[:, ::-1], active=np.array([False, True, True]))
    assert np.corrcoef(flipped, score)[0, 1] > 0.999
    role = cp.role_score(feats, ['active', 'enhancer', 'repressive'])
    assert np.corrcoef(role, label)[0, 1] > 0.8
    with pytest.raises(ValueError):
        cp.role_score(feats, ['repressive'] * 3)
    # a row with a missing track gets NaN, the others are still scored
    feats[5, 1] = np.nan
    score, _, _ = cp.chromatin_pc1(feats, active=np.array([True, True, False]))
    assert np.isnan(score[5]) and np.isfinite(score[6])


def test_contact_eigenvector_recovers_compartments():
    label = _toy_compartments(150)
    m = _toy_hic(label)
    vec, index = cp.contact_eigenvector(m, orient_by=label)
    assert index == 0
    assert np.mean(np.sign(vec) == label) > 0.95
    oe = cp.observed_expected(m)
    # expected removed: diagonals average to 1
    assert np.nanmean(np.diagonal(oe, 3)) == pytest.approx(1.0)


def test_contact_eigenvector_selects_compartments_over_arms():
    label = np.tile(np.repeat([1.0, -1.0], 10), 8)          # balanced in both arms, uncorrelated with the arms
    m = _toy_hic(label, seed=7)
    arm = np.where(np.arange(160) < 80, 1, -1)
    m = m * np.where(np.outer(arm, arm) > 0, 10.0, 0.1)     # strong arm pattern dominates the first eigenvector
    first, _ = cp.contact_eigenvector(m, orient_by=label)
    chosen, index = cp.contact_eigenvector(m, orient_by=label, n_candidates=3)
    assert abs(np.corrcoef(first, arm)[0, 1]) > 0.7 and np.mean(np.sign(first) == label) < 0.8
    assert index > 0 and np.mean(np.sign(chosen) == label) > 0.9


def test_contact_eigenvector_sparse_coverage():
    label = _toy_compartments(150, seed=3)
    m = _toy_hic(label, seed=3)
    covered = np.zeros(150, dtype=bool)
    covered[::3] = True                                   # only every third bin baited
    m[~covered] = 0
    m[:, ~covered] = 0
    vec, _ = cp.contact_eigenvector(m, orient_by=label)
    assert np.all(np.isnan(vec[~covered]))
    assert np.mean(np.sign(vec[covered]) == label[covered]) > 0.9


def test_profile_eigenvector_and_per_chrom():
    rng = np.random.default_rng(2)
    label = _toy_compartments(200, seed=2)
    profiles = label[:, None] * rng.uniform(0.5, 1.5, 6) + rng.normal(0, 0.4, (200, 6))
    vec = cp.profile_eigenvector(profiles, orient_by=label)
    assert np.mean(np.sign(vec) == label) > 0.9
    with pytest.raises(ValueError):
        cp.profile_eigenvector(profiles[:, :2])
    bins = cp.make_bins({'chr1': 100 * 100, 'chr2': 100 * 100}, 100)
    both = cp.per_chrom(cp.profile_eigenvector, profiles, bins, orient_by=label)
    assert np.mean(np.sign(both) == label) >= 0.85


def test_compare_to_reference_metrics():
    ref = np.array([1, 2, -1, -2, 0, 1, -1, 2, -2, 1] * 4, dtype=float)
    chroms = np.array(['chr1'] * 20 + ['chr2'] * 20)
    table = cp.compare_to_reference(ref * 2 + 5, ref, chroms, threshold='median', min_bins=5)
    genome = table[table.chrom == 'genome'].iloc[0]
    assert genome.pearson == pytest.approx(1.0) and genome.auc == pytest.approx(1.0)
    assert genome.n_bins == 36                            # reference zeros skipped
    # with threshold zero the shifted estimate calls every bin A
    zero = cp.compare_to_reference(ref * 2 + 5, ref, chroms, threshold='zero', min_bins=5)
    assert zero[zero.chrom == 'genome'].sign_agreement.iloc[0] < 0.6
    inverted = cp.compare_to_reference(-ref, ref, chroms, min_bins=5)
    assert inverted[inverted.chrom == 'genome'].auc.iloc[0] == pytest.approx(0.0)


def test_loco_regression_is_out_of_chromosome():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(400, 3))
    y = x @ np.array([1.0, -0.5, 0.0]) + rng.normal(0, 0.1, 400)
    chroms = np.repeat(['chr1', 'chr2', 'chr3', 'chr4'], 100)
    pred = cp.loco_regression(x, y, chroms)
    assert np.corrcoef(pred, y)[0, 1] > 0.95


def test_community_boundaries_and_enrichment():
    nodes = pd.DataFrame({'chrom': 'chr1', 'start': [0, 100, 200, 300, 400], 'end': [50, 150, 250, 350, 450],
                          'community': [1, 1, 2, 2, 1]})
    b = cp.community_boundaries(nodes)
    assert list(b.position) == [175, 375]
    reference = {'chr1': np.array([170, 380])}
    out = cp.boundary_enrichment({'chr1': b.position.values}, reference, {'chr1': (0, 10000)}, window=10,
                                 n_shifts=200, null_positions={'chr1': np.arange(0, 10000, 50)})
    assert out['observed_fraction'] == 1.0
    assert out['shift_enrichment'] > 10 and out['shift_p'] < 0.05
    assert out['matched_p'] < 0.05


def test_domain_statistics():
    domains = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(0, 100000, 10000),
                            'end': np.arange(10000, 100001, 10000)})
    assert list(cp.domain_id(domains, 'chr1', [5, 15000, 200000])) == [0, 1, -1]
    # irregular domains, spans and edges inside them
    rng = np.random.default_rng(5)
    edges_pos = np.concatenate([[0], np.cumsum(rng.integers(5000, 40000, 40))])
    irregular = pd.DataFrame({'chrom': 'chr1', 'start': edges_pos[:-1], 'end': edges_pos[1:]})
    hi = int(edges_pos[-1])
    spans = pd.DataFrame({'chrom': 'chr1', 'start': irregular.start + 1000, 'end': irregular.end - 1000})
    out = cp.spans_within_domain(spans, irregular, {'chr1': (0, hi)}, n_shifts=200)
    assert out['observed_fraction'] == 1.0 and out['shift_fraction'] < 0.5 and out['shift_p'] < 0.05
    edges = pd.DataFrame({'chrom': 'chr1', 'mid_a': irregular.start + 1000, 'mid_b': irregular.end - 1000})
    e = cp.edges_crossing(edges, irregular, {'chr1': (0, hi)}, n_shifts=200)
    assert e['observed_fraction'] == 0.0 and e['shift_fraction'] > 0.3 and e['p_deficit'] < 0.05
    nested = pd.DataFrame({'chrom': ['chr1', 'chr1'], 'start': [0, 100], 'end': [1000, 200]})
    assert list(cp.domain_id(nested, 'chr1', [50, 150])) == [0, 1]


def test_compartment_switches():
    a = np.array([1.0, -1.0, 0.2, -0.1, np.nan, 1.0])
    b = np.array([-1.0, 1.0, -0.1, 0.2, 1.0, 1.2])
    calls = cp.compartment_switches(a, b, threshold=0.5)
    assert list(calls) == ['A->B', 'B->A', 'A', 'B', 'NA', 'A']
    z = cp.standardise_per_chrom(np.array([1.0, 2.0, 3.0, 10.0, 20.0, 30.0]), np.array(['c1'] * 3 + ['c2'] * 3))
    assert z[1] == 0 and z[4] == 0 and z[0] == pytest.approx(z[3])


def test_flatten_domains():
    nested = pd.DataFrame({'chrom': 'chr1', 'start': [0, 100, 2000], 'end': [1000, 200, 3000]})
    flat = cp.flatten_domains(nested)
    assert list(zip(flat.start, flat.end)) == [(0, 100), (100, 200), (200, 1000), (2000, 3000)]
    ids = cp.domain_id(flat, 'chr1', [50, 150, 500, 1500])
    assert ids[0] != ids[1] != ids[2] and ids[3] == -1


def test_ega_compartments_cli(tmp_path):
    import pyBigWig
    from epigeneticgraphanalyzer import egaCompartments
    size, res = 4000000, 100000
    (tmp_path / 'sizes').write_text('chr1\t{}\n'.format(size))
    label = np.repeat([1.0, -1.0, 1.0, -1.0], 10)          # 40 bins of 100 kb
    rng = np.random.default_rng(0)
    paths = {}
    for name, sign in (('K27ac', 1.0), ('K4me1', 1.0), ('K27me3', -0.5)):
        path = str(tmp_path / (name + '.bw'))
        bw = pyBigWig.open(path, 'w')
        bw.addHeader([('chr1', size)])
        values = np.exp(sign * label + rng.normal(0, 0.3, 40))
        bw.addEntries(['chr1'] * 40, list(range(0, size, res)), ends=list(range(res, size + 1, res)),
                      values=[float(v) for v in values])
        bw.close()
        paths[name] = path
    ref = tmp_path / 'ref.bedGraph'
    ref.write_text(''.join('chr1\t{}\t{}\t{}\n'.format(i * res, (i + 1) * res, v) for i, v in enumerate(label)))
    genes = tmp_path / 'genes.bed'
    genes.write_text(''.join('chr1\t{}\t{}\n'.format(i * res + 10, i * res + 20)
                             for i in range(40) if label[i] > 0))
    nodes = pd.DataFrame({'node': ['n{}'.format(i) for i in range(8)], 'chrom': 'chr1',
                          'start': np.arange(8) * 400000 + 1000, 'end': np.arange(8) * 400000 + 5000,
                          'bait': True, 'community': [0, 0, 1, 1, 2, 2, 3, 3]})
    nodes.to_csv(tmp_path / 'nodes.tsv', sep='\t', index=False)
    pd.DataFrame({'node_a': ['n0', 'n2'], 'node_b': ['n1', 'n3'], 'weight': [5, 5]}).to_csv(
        tmp_path / 'edges.tsv', sep='\t', index=False)
    (tmp_path / 'tads.bed').write_text(''.join('chr1\t{}\t{}\n'.format(i * 800000, (i + 1) * 800000)
                                               for i in range(5)))
    prefix = str(tmp_path / 'out' / 'toy')
    egaCompartments.main(['--chromSizes', str(tmp_path / 'sizes'), '--tracks'] +
                         ['{}={}'.format(k, v) for k, v in paths.items()] +
                         ['--markRoles', 'K27ac=active', 'K4me1=enhancer', 'K27me3=repressive',
                          '--genes', str(genes), '--reference', str(ref), '--nodes', str(tmp_path / 'nodes.tsv'),
                          '--edges', str(tmp_path / 'edges.tsv'), '--domains', str(tmp_path / 'tads.bed'),
                          '--shifts', '50', '-o', prefix])
    evaluation = pd.read_csv(prefix + '_evaluation.tsv', sep='\t')
    genome = evaluation[(evaluation.chrom == 'genome') & (evaluation.estimate == 'chromatin_pc1')].iloc[0]
    assert genome.pearson > 0.9 and genome.auc > 0.95
    domains = pd.read_csv(prefix + '_domains.tsv', sep='\t')
    assert set(domains.test) == {'community_boundaries_at_tad_boundaries', 'communities_within_one_tad',
                                 'edges_crossing_tad_boundaries'}


def test_eigenvector_selection_margin():
    rng = np.random.default_rng(8)
    phasing = rng.normal(size=500)
    noise = rng.normal(size=(500, 2))

    def mix(r):
        return r * phasing + np.sqrt(1 - r ** 2) * noise[:, 0]

    vecs = np.column_stack([mix(0.5), 0.6 * phasing + 0.8 * noise[:, 1]])
    assert cp._select(vecs, phasing)[1] == 0                # slightly better second eigenvector: keep the first
    vecs = np.column_stack([mix(0.2), 0.9 * phasing + 0.44 * noise[:, 1]])
    assert cp._select(vecs, phasing)[1] == 1
