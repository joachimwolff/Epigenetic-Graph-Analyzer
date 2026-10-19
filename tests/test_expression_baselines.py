"""Window baseline of the expression models (analysis/expression_baselines.py) on a tiny synthetic example."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from expression_baselines import log2_change, nested_test, window_signals  # noqa: E402

from epigeneticgraphanalyzer.differential import nested_models  # noqa: E402

# promoter fragment chr1:10000-12000 (centre 11000)
SOURCES = pd.DataFrame({'chrom': ['chr1', 'chr3'], 'start': [10000, 10000], 'end': [12000, 12000]})
# elements: overlapping the promoter, at 10 kb (right), 20 kb, 10 kb (left), 50 kb, and one on chr2
ELEMENTS = pd.DataFrame({'chrom': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr2'],
                         'start': [10500, 20500, 30500, 500, 60500, 20500],
                         'end': [11500, 21500, 31500, 1500, 61500, 21500]})
A = np.array([100.0, 2.0, 4.0, 0.0, 50.0, 70.0])
B = np.array([100.0, 4.0, 4.0, 6.0, 90.0, 80.0])


def test_window_weighted_by_inverse_distance():
    S, n = window_signals(SOURCES, ELEMENTS, np.column_stack([A, B]), 0, 40000)
    # weights 1/10 kb, 1/20 kb, 1/10 kb -> 0.4, 0.2, 0.4; the overlapping element, the one at 50 kb and chr2 are out
    assert n.tolist() == [3, 0]
    assert S[0] == pytest.approx([0.4 * 2 + 0.2 * 4 + 0.4 * 0, 0.4 * 4 + 0.2 * 4 + 0.4 * 6])
    assert S[1].tolist() == [0.0, 0.0]
    assert log2_change(S[0, 0], S[0, 1])[()] == pytest.approx(np.log2(5.8 / 2.6))


def test_window_range_bounds_and_exponent():
    values = np.column_stack([A, B])
    # lo < d <= hi: d = 10 kb excluded at lo = 10 kb, d = 20 kb included at hi = 20 kb
    S, n = window_signals(SOURCES, ELEMENTS, values, 10000, 20000)
    assert n[0] == 1 and S[0].tolist() == [4.0, 4.0]
    # per-source ranges; hi <= lo gives no elements
    S, n = window_signals(SOURCES, ELEMENTS, values, np.array([0, 0]), np.array([0, 40000]))
    assert n.tolist() == [0, 0] and S[0].tolist() == [0.0, 0.0]
    # exponent 0: unweighted mean
    S, _ = window_signals(SOURCES, ELEMENTS, values, 0, 60000, exponent=0.0)
    assert S[0] == pytest.approx([np.mean([2, 4, 0, 50]), np.mean([4, 4, 6, 90])])
    # the 1/d weights are the same in both conditions: equal signals give a change of 0
    S, _ = window_signals(SOURCES, ELEMENTS, np.column_stack([A, A]), 0, 60000)
    assert log2_change(S[0, 0], S[0, 1])[()] == pytest.approx(0.0)


def test_nested_test_matches_nested_models():
    rng = np.random.default_rng(1)
    n = 200
    table = pd.DataFrame({'rna_log2FC': rng.normal(size=n)})
    for mark in ('K27ac', 'K4me3'):
        for level in ('promoter', 'direct', 'indirect'):
            table['{}_{}'.format(mark, level)] = rng.normal(size=n) + 0.3 * table['rna_log2FC']
    reference = nested_models(table, ['K27ac', 'K4me3'], levels=['direct', 'indirect']).set_index('model')
    result = nested_test(table, ['K27ac_promoter', 'K4me3_promoter'], ['K27ac_direct', 'K4me3_direct'])
    assert result['r2_base'] == pytest.approx(reference.loc['promoter', 'r2'])
    assert result['delta_r2'] == pytest.approx(reference.loc['promoter+direct', 'delta_r2'])
    assert result['f_pvalue'] == pytest.approx(reference.loc['promoter+direct', 'f_pvalue'])
