import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from crispr_head_to_head import METHODS, coverage, evaluate, lowest_for_missing, score_table, strata  # noqa: E402


def pairs():
    return pd.DataFrame(dict(
        name=['p1', 'p2', 'p3', 'p4', 'p5', 'p6'], gene=['A', 'A', 'B', 'B', 'C', 'C'],
        element=['e1', 'e2', 'e3', 'e4', 'e5', 'e6'],
        distanceToTSS=[5000, -30000, 100000, 15000, -250000, 40000], same_fragment=[False, False, True, False,
                                                                                    False, False],
        scorable=[True, True, False, True, False, True], prediction_found=['True', 'False', 'True', 'True', 'True',
                                                                           'True'],
        abc=[0.1, 0.0, 0.02, 0.3, 0.001, 0.05], re2g=[0.5, 0.0, 0.1, 0.7, 0.01, 0.2],
        powerlaw=[0.05, 0.0, 0.01, 0.03, 0.001, 0.02], direct=[8.0, 0.0, 0.0, 6.0, 0.0, 0.0],
        two_hop=[0.0, 5.0, 0.0, 0.0, 0.0, 3.0], ppr=[0.01, 0.002, 0.0, 0.02, 0.0, 0.001], y=[1, 0, 0, 1, 0, 0]))


def logistic():
    return pd.DataFrame({'M1_direct': [0.6, 0.1, 0.5, 0.2], 'M3_all_graph': [0.7, 0.2, 0.4, 0.1]},
                        index=['p1', 'p2', 'p4', 'p6'])


def test_score_table_marks_uncovered_pairs_missing():
    s = score_table(pairs(), logistic())
    assert list(s.columns) == METHODS
    assert s['distance'].notna().all() and s.loc[0, 'distance'] == -5000
    assert np.isnan(s.loc[1, 'encode_re2g']) and s.loc[3, 'encode_re2g'] == 0.7
    assert np.isnan(s.loc[2, 'ega_direct']) and np.isnan(s.loc[4, 'ega_logistic_M1'])
    assert s.loc[3, 'ega_logistic_M3'] == 0.4


def test_lowest_for_missing_ranks_below_all_observed():
    s = score_table(pairs(), logistic())
    f = lowest_for_missing(s)
    assert f.notna().all().all()
    assert f.loc[1, 'encode_re2g'] < s['encode_re2g'].min()
    assert f.loc[2, 'ega_direct'] < s['ega_direct'].min()


def test_coverage_and_strata():
    p = pairs()
    c = coverage(score_table(p, logistic()), p).set_index('method')
    assert c.loc['distance', 'covered'] == 6 and c.loc['encode_re2g', 'covered'] == 5
    assert c.loc['ega_direct', 'covered'] == 4 and c.loc['ega_direct', 'positives_covered'] == 2
    s = strata(p)
    assert list(s['element_off_tss_fragment']) == [True, True, False, True, True, True]
    assert list(s['distance_gt_20kb']) == [False, True, True, False, True, True]


def test_evaluate_intersection_uses_complete_pairs_only():
    p = pairs()
    m = evaluate(p, score_table(p, logistic()), 'toy', n_boot=10, seed=0, threads=1)
    inter = m[(m['stratum'] == 'all_pairs') & (m['evaluation'] == 'intersection')]
    assert (inter['pairs'] == 3).all()   # p1, p4, p6
    allp = m[(m['stratum'] == 'all_pairs') & (m['evaluation'] == 'missing_lowest')]
    assert (allp['pairs'] == 6).all()
    assert 'delta_vs_distance_auroc' in m.columns and 'delta_vs_encode_re2g_auprc' in m.columns
