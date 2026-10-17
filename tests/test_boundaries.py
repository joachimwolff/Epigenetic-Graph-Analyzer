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


def _toy_evidence(n_chrom=4, n_bins=200, res=25_000, seed=4):
    """Planted boundaries: insulation dip + reverse CTCF in the boundary bin and forward CTCF in the next bin.
    Decoys: inward-facing CTCF pairs (forward left, reverse right) without a dip, and dips without CTCF."""
    rng = np.random.default_rng(seed)
    bins, sites, truth, dips = [], [], [], []
    for c in range(n_chrom):
        chrom = 'chr{}'.format(c + 1)
        b = pd.DataFrame({'chrom': chrom, 'start': np.arange(n_bins) * res})
        b['end'] = b.start + res
        bins.append(b)
        slots = rng.permutation(np.arange(8, n_bins - 8, 12))
        bounds, ctcf_decoys, dip_decoys = np.sort(slots[:5]), slots[5:9], slots[9:13]
        truth.append((chrom, bounds))
        dips.append((chrom, np.concatenate([bounds, dip_decoys])))
        for k in bounds:
            sites.append((chrom, k * res + 12_000, 20.0, 6.0, True, '-'))
            sites.append((chrom, (k + 1) * res + 12_000, 20.0, 6.0, True, '+'))
        for k in ctcf_decoys:
            sites.append((chrom, k * res + 12_000, 20.0, 6.0, True, '+'))
            sites.append((chrom, (k + 1) * res + 12_000, 20.0, 6.0, True, '-'))
    bins = pd.concat(bins, ignore_index=True)
    sites = pd.DataFrame(sites, columns=['chrom', 'summit', 'ctcf', 'cohesin', 'cohesin_bound', 'orientation'])
    ins = pd.DataFrame({'insulation': rng.normal(0, 0.05, len(bins)), 'covered': True}, index=bins.index)
    labels = np.zeros(len(bins), dtype=bool)
    idx = bd.bin_index(bins)
    for (chrom, bounds), (_, d) in zip(truth, dips):
        first = idx[chrom][0]
        labels[first + bounds] = True
        ins.loc[first + d, 'insulation'] -= 1.0
    return bins, sites, ins, labels


def test_planted_boundary_rules():
    bins, sites, ins, labels = _toy_evidence()
    codes = bd.chrom_codes(bins)
    ev = bd.bin_evidence(bins, sites, ins, None, dip_prominence=0.3, ctcf_min=10)
    truth = np.flatnonzero(labels)
    res = {}
    for rule in ('dip', 'ctcf', 'ctcf_pattern', 'dip_or_ctcf', 'dip_and_ctcf', 'dip_and_ctcf_pattern'):
        flags, score = bd.rule_flags(ev, rule)
        res[rule] = bd.match_boundaries(bd.calls_from_flags(flags, score, codes), truth, codes, tolerance=1)
    # the planted boundaries (dip + convergent-boundary CTCF pattern) are all found and nothing else
    assert res['dip_and_ctcf_pattern']['recall'] == 1.0 and res['dip_and_ctcf_pattern']['precision'] == 1.0
    # each evidence alone also finds them but picks up its decoys
    assert res['dip']['recall'] == 1.0 and res['dip']['precision'] < 0.7
    assert res['ctcf']['recall'] == 1.0 and res['ctcf']['precision'] < 0.7
    # the orientation pattern rejects the inward-facing decoys
    assert res['ctcf_pattern']['precision'] > res['ctcf']['precision']
    assert res['dip_or_ctcf']['n_called'] > res['dip_and_ctcf']['n_called']


def test_orientation_pattern_sides():
    """Reverse site in the bin, forward one bin to the right: boundary pattern at the bin; the mirrored
    arrangement (forward left, reverse right, the inside of a loop) is not."""
    res = 10_000
    bins = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(20) * res})
    bins['end'] = bins.start + res
    sites = pd.DataFrame({'chrom': 'chr1', 'summit': [50_500, 61_000], 'ctcf': 10.0, 'cohesin': 5.0,
                          'cohesin_bound': True, 'orientation': ['-', '+']})
    ev = bd.bin_evidence(bins, sites)
    assert ev.ctcf_pattern.iat[5] and ev.pattern_score.iat[5] == 10.0
    assert not ev.ctcf_pattern.iat[3] and not ev.ctcf_pattern.iat[8]
    swapped = bd.bin_evidence(bins, sites.assign(orientation=['+', '-']))
    assert not swapped.ctcf_pattern.any() and swapped.ctcf.iat[5]
    weak = bd.bin_evidence(bins, sites, ctcf_min=20)
    assert not weak.ctcf.any()
    unbound = bd.bin_evidence(bins, sites.assign(cohesin_bound=False))
    assert not unbound.ctcf.any()


def test_calls_from_flags_one_per_run():
    codes = np.array([0] * 6 + [1] * 4)
    flags = np.array([1, 1, 0, 1, 1, 1, 1, 1, 0, 0], dtype=bool)
    score = np.array([1, 2, 0, 5, 3, 9, 4, 1, 0, 0], dtype=float)
    # runs: [0, 1], [3, 4, 5] on chromosome 0 and [6, 7] on chromosome 1 (a run does not cross chromosomes)
    assert list(bd.calls_from_flags(flags, score, codes)) == [1, 5, 6]


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
                        str(tmp_path / 'sizes'), '--binSize', '25000', '--minExpected', '3', '--dipProminence', '0.5',
                        '--reference', str(tmp_path / 'ref.bed'), '-o', prefix])
    calls = pd.read_csv(prefix + '_boundaries.bed', sep='\t', header=None)
    assert np.min(np.abs(calls[1].values - 1_000_000)) <= 25_000
    ev = pd.read_csv(prefix + '_evaluation.tsv', sep='\t')
    assert ev[(ev.region == 'all') & (ev.tolerance_bins == 1)].recall.iat[0] == 1.0
