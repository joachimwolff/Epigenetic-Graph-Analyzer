"""Boundary calling from capture insulation, CTCF / cohesin and motif orientation on toy data."""

import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import boundaries as bd

# JASPAR MA0139-like counts (forward orientation CCASYAGRKGGCRSY)
CONSENSUS = 'GCCACCAGGGGGCGC'


def _pfm():
    counts = np.full((len(CONSENSUS), 4), 2.0)
    for k, base in enumerate(CONSENSUS):
        counts[k, bd.BASES[base]] = 94.0
    return counts


def _revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def test_motif_orientation():
    pwm = bd.log_odds(_pfm())
    thr = bd.score_threshold(pwm, 1e-4)
    rng = np.random.default_rng(1)
    background = ''.join(rng.choice(list('ACGT'), size=301))
    fwd = background[:100] + CONSENSUS + background[100 + len(CONSENSUS):]
    rev = background[:100] + _revcomp(CONSENSUS) + background[100 + len(CONSENSUS):]
    both = background[:50] + CONSENSUS + background[50 + len(CONSENSUS):200] + _revcomp(CONSENSUS) + \
        background[200 + len(CONSENSUS):]
    codes = np.vstack([bd.encode(s) for s in (fwd, rev, both)])
    f, r = bd.scan_motif(codes, pwm)
    assert list(bd.orientation(f, r, thr)) == ['+', '-', 'both']


def test_score_threshold_tail():
    pwm = bd.log_odds(_pfm())
    rng = np.random.default_rng(2)
    codes = rng.integers(0, 4, size=(200000, len(CONSENSUS))).astype(np.uint8)
    f, r = bd.scan_motif(codes, pwm)
    for p in (1e-2, 1e-3):
        thr = bd.score_threshold(pwm, p)
        # single windows of random sequence: the false positive rate is at most p (discrete scores) and the next
        # lower attainable score would exceed p
        rate = np.mean(f >= thr)
        assert 0 < rate <= 1.2 * p
        assert np.mean(f >= thr - 5.5) > p


def test_classify_loops():
    sites = pd.DataFrame({'chrom': 'chr1', 'summit': [1000, 50000, 100000, 150000],
                          'ctcf': [10.0, 12.0, 8.0, 9.0], 'cohesin': 5.0, 'cohesin_bound': True,
                          'orientation': ['+', '-', '-', '+']})
    ends = pd.DataFrame({'chrom': 'chr1',
                         'left_start': [0, 49000, 0, 99000], 'left_end': [2000, 51000, 2000, 101000],
                         'right_start': [49000, 99000, 149000, 149000], 'right_end': [51000, 101000, 151000, 151000]})
    cls, strength = bd.classify_loops(ends, sites)
    assert list(cls) == ['convergent', 'tandem_reverse', 'tandem_forward', 'divergent']
    assert strength[0] == pytest.approx(22.0)
    # cohesin-free sites are ignored
    cls2, _ = bd.classify_loops(ends, sites.assign(cohesin_bound=False))
    assert set(cls2) == {'none'}


def _toy_graph(boundary=1_000_000, size=2_000_000, n_baits=80, seed=3):
    """Baits along one chromosome; edges of random length that never cross the planted boundary."""
    rng = np.random.default_rng(seed)
    baits = np.sort(rng.integers(20_000, size - 20_000, n_baits))
    rows = []
    for b in baits:
        for _ in range(15):
            length = int(rng.exponential(150_000)) + 5_000
            side = rng.choice([-1, 1])
            other = b + side * length
            if not (0 < other < size):
                continue
            if (b < boundary) != (other < boundary):
                other = b - side * length         # reflect away from the boundary
                if not (0 < other < size) or (b < boundary) != (other < boundary):
                    continue
            rows.append(('chr1', b, other, True, False, rng.gamma(2.0, 3.0)))
    return pd.DataFrame(rows, columns=['chrom', 'a', 'b', 'anchor_a', 'anchor_b', 'weight'])


def test_capture_insulation_dip_at_planted_boundary():
    edges = _toy_graph(boundary=1_012_500)     # the centre of bin 40
    bins = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(0, 2_000_000, 25_000)})
    bins['end'] = bins.start + 25_000
    ins = bd.capture_insulation(bins, edges, min_expected=3)
    k = 1_000_000 // 25_000
    assert ins.observed.iat[k] == 0
    assert ins.insulation.iat[k] == pytest.approx(np.nanmin(ins.insulation.values))
    calls = bd.local_minima(ins.insulation.values, bd.chrom_codes(bins), prominence=0.5)
    assert np.min(np.abs(calls - k)) <= 1
    # the expectation keeps the edge count: observed and expected sum to the same total span
    assert ins.observed.sum() == pytest.approx(ins.expected.sum(), rel=0.1)


def _toy_features(n_chrom=6, n_bins=200, res=25_000, seed=4):
    """Several chromosomes with planted boundaries: insulation dip + reverse CTCF left / forward CTCF right of the
    boundary, decoy CTCF sites elsewhere."""
    rng = np.random.default_rng(seed)
    bins, sites, truth = [], [], []
    for c in range(n_chrom):
        chrom = 'chr{}'.format(c + 1)
        b = pd.DataFrame({'chrom': chrom, 'start': np.arange(n_bins) * res})
        b['end'] = b.start + res
        bins.append(b)
        bounds = np.sort(rng.choice(np.arange(10, n_bins - 10, 20), size=6, replace=False) + rng.integers(-3, 4, 6))
        truth.append((chrom, bounds))
        for k in bounds:     # boundary: reverse site in the boundary bin, forward site in the next bin
            sites.append((chrom, k * res + 12_000, 20.0, 6.0, True, '-'))
            sites.append((chrom, (k + 1) * res + 12_000, 20.0, 6.0, True, '+'))
        for k in rng.integers(2, n_bins - 2, 6):      # decoys: same sites, inward-facing (forward left, reverse right)
            sites.append((chrom, k * res + 12_000, 20.0, 6.0, True, '+'))
            sites.append((chrom, (k + 1) * res + 12_000, 20.0, 6.0, True, '-'))
    bins = pd.concat(bins, ignore_index=True)
    sites = pd.DataFrame(sites, columns=['chrom', 'summit', 'ctcf', 'cohesin', 'cohesin_bound', 'orientation'])
    ins = pd.DataFrame({'observed': 20.0, 'expected': 20.0, 'expected_global': 20.0, 'observed_w': 20.0,
                        'expected_w': 20.0, 'covered': True}, index=bins.index)
    noise = rng.normal(0, 0.15, len(bins))
    ins['insulation'] = noise
    labels = np.zeros(len(bins), dtype=bool)
    idx = bd.bin_index(bins)
    for chrom, bounds in truth:
        first = idx[chrom][0]
        labels[first + bounds] = True
        ins.loc[first + bounds, 'insulation'] -= 1.0
    ins['insulation_w'] = ins.insulation
    ins['insulation_global'] = ins.insulation
    return bins, sites, ins, labels


def test_planted_boundary_found_on_held_out_chromosome():
    bins, sites, ins, labels = _toy_features()
    feats = bd.bin_features(bins, sites, ins, None, context_bp=25_000)
    # the orientation pattern is highest at the planted boundaries
    assert feats.flip.values[labels].mean() > 3 * feats.flip.values[~labels].mean()
    codes = bd.chrom_codes(bins)
    train = codes < 4
    test = ~train
    model = bd.LogisticModel(bd.METHODS['combined']).fit(feats[train], labels[train])
    score = model.predict(feats)
    thr, f1_train = bd.best_threshold(score[train], np.flatnonzero(labels[train]), codes[train])
    calls = bd.peak_calls(np.where(test, score, -np.inf), codes, thr)
    m = bd.match_boundaries(calls, np.flatnonzero(labels & test), codes, tolerance=1)
    assert m['recall'] >= 0.9 and m['precision'] >= 0.8
    # CTCF without orientation cannot tell the planted boundaries from the inward-facing decoys
    aps = {}
    for method in ('ctcf_cohesin', 'ctcf_orientation'):
        s = bd.LogisticModel(bd.METHODS[method]).fit(feats[train], labels[train]).predict(feats)
        aps[method] = bd.average_precision(s[test], labels[test])
    assert aps['ctcf_orientation'] > aps['ctcf_cohesin'] + 0.2


def test_orientation_pattern_sides():
    """Reverse site in the bin, forward one bin to the right: flip pattern at the bin; the mirrored arrangement
    (forward left, reverse right) is the inside-loop pattern."""
    res = 10_000
    bins = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(20) * res})
    bins['end'] = bins.start + res
    sites = pd.DataFrame({'chrom': 'chr1', 'summit': [50_500, 61_000], 'ctcf': 10.0, 'cohesin': 5.0,
                          'cohesin_bound': True, 'orientation': ['-', '+']})
    f = bd.bin_features(bins, sites, None, None, context_bp=res)
    assert f.flip.iat[5] > 0 and f.antiflip.iat[5] == 0
    swapped = sites.assign(orientation=['+', '-'])
    g = bd.bin_features(bins, swapped, None, None, context_bp=res)
    assert g.flip.iat[5] == 0 and g.antiflip.iat[5] > 0


def test_diamond_insulation_block_boundary():
    n = 60
    i, j = np.indices((n, n))
    m = 100.0 / (np.abs(i - j) + 1)
    m[(i < 30) != (j < 30)] *= 0.2
    ins = bd.diamond_insulation(bd.band_from_dense(m, 10), 5)
    assert np.nanargmin(ins) in (29, 30)


def test_metrics_and_partitions():
    codes = np.zeros(100, dtype=int)
    codes[50:] = 1
    m = bd.match_boundaries([10, 30, 49, 70], [11, 50, 90], codes, tolerance=1)
    # 49 (chromosome 0) must not match 50 (chromosome 1)
    assert m['precision'] == pytest.approx(0.25) and m['recall'] == pytest.approx(1 / 3)
    assert bd.average_precision([3, 2, 1, 0], [1, 1, 0, 0]) == pytest.approx(1.0)
    assert bd.average_precision([0, 1, 2, 3], [1, 1, 0, 0]) == pytest.approx(0.5, abs=0.1)
    a = bd.segment_labels(codes, [20, 70])
    assert bd.variation_of_information(a, a) == pytest.approx(0.0)
    agree, _, _ = bd.local_rand(a, a, codes, 10)
    assert agree == 1.0
    b = bd.segment_labels(codes, [25, 80])
    assert bd.variation_of_information(a, b) > 0


def test_cli_capture_only(tmp_path):
    from epigeneticgraphanalyzer import egaBoundaries
    edges = _toy_graph(boundary=1_012_500)

    def nid(p):
        return 'chr1:{}-{}'.format(p - 1000, p + 1000)

    baits = set(np.unique(edges.a))
    nodes = {nid(p): ('chr1', p - 1000, p + 1000, p in baits) for p in np.concatenate([edges.a, edges.b])}
    pd.DataFrame([(k, *v, 1, 0, '') for k, v in nodes.items()],
                 columns=['node', 'chrom', 'start', 'end', 'bait', 'degree', 'community', 'genes']).to_csv(
        tmp_path / 'n.tsv', sep='\t', index=False)
    pd.DataFrame({'node_a': [nid(a) for a in edges.a], 'node_b': [nid(b) for b in edges.b],
                  'weight': edges.weight}).to_csv(tmp_path / 'e.tsv', sep='\t', index=False)
    (tmp_path / 'sizes').write_text('chr1\t2000000\n')
    (tmp_path / 'ref.bed').write_text('chr1\t1000000\t1025000\n')
    prefix = str(tmp_path / 'out')
    egaBoundaries.main(['--nodes', str(tmp_path / 'n.tsv'), '--edges', str(tmp_path / 'e.tsv'), '--chromSizes',
                        str(tmp_path / 'sizes'), '--binSize', '25000', '--minExpected', '3', '--prominence', '0.5',
                        '--reference', str(tmp_path / 'ref.bed'), '-o', prefix])
    calls = pd.read_csv(prefix + '_boundaries.bed', sep='\t', header=None)
    assert np.min(np.abs(calls[1].values - 1_000_000)) <= 25_000
    ev = pd.read_csv(prefix + '_evaluation.tsv', sep='\t')
    assert ev[(ev.region == 'all') & (ev.tolerance_bins == 1)].recall.iat[0] == 1.0
