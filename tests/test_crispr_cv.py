import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from crispr_benchmark import prediction_found  # noqa: E402
from crispr_cv import (MODELS, auroc, bootstrap_metrics, folds, gene_resamples, lr_test, out_of_fold,  # noqa: E402
                       split_overlap)


def synthetic(n_genes=60, per_gene=12, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        chrom = 'chr{}'.format(g % 6 + 1)
        for k in range(per_gene):
            distance = rng.integers(1000, 500000)
            direct = rng.choice([0.0, rng.uniform(5, 20)], p=[0.8, 0.2])
            logit = -1 - 1.5 * (np.log10(distance + 1000) - 4.5) + 0.2 * direct
            start = 1000000 * g + 1000 * k
            rows.append(dict(gene='G{}'.format(g), chrom=chrom, chromStart=start, chromEnd=start + 500,
                             element='{}:{}'.format(chrom, start), log_distance=np.log10(distance + 1000),
                             element_degree=rng.integers(0, 10), direct=direct, two_hop=rng.uniform(0, 5),
                             hops=rng.uniform(0, 1), community=rng.integers(0, 2), ppr=rng.uniform(0, 0.01),
                             y=int(rng.uniform() < 1 / (1 + np.exp(-logit)))))
    return pd.DataFrame(rows)


def test_auroc_matches_sklearn_with_ties():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 200)
    s = rng.integers(0, 5, 200).astype(float)
    assert np.isclose(auroc(y, s), roc_auc_score(y, s))
    assert np.isnan(auroc(np.zeros(5), np.arange(5)))


def test_gene_resamples_keep_genes_whole():
    genes = ['a', 'a', 'b', 'c', 'c', 'c']
    for index in gene_resamples(genes, 20, seed=3):
        drawn = pd.Series(np.array(genes)[index]).value_counts()
        size = {'a': 2, 'b': 1, 'c': 3}
        assert all(count % size[g] == 0 for g, count in drawn.items())
    a = gene_resamples(genes, 5, seed=3)
    b = gene_resamples(genes, 5, seed=3)
    assert all(np.array_equal(x, z) for x, z in zip(a, b))


def test_folds_never_share_genes_or_chromosomes():
    table = synthetic()
    for train, test in folds(table, 'gene'):
        assert not set(table['gene'].iloc[train]) & set(table['gene'].iloc[test])
    chrom_folds = folds(table, 'chromosome')
    assert len(chrom_folds) == table['chrom'].nunique()
    for train, test in chrom_folds:
        assert table['chrom'].iloc[test].nunique() == 1
        assert not set(table['chrom'].iloc[train]) & set(table['chrom'].iloc[test])


def test_out_of_fold_predicts_every_pair_and_direct_helps():
    table = synthetic(n_genes=120)
    splits = folds(table, 'gene')
    base = out_of_fold(table, MODELS['M0_distance_degree'], splits)
    direct = out_of_fold(table, MODELS['M1_direct'], splits)
    assert not np.isnan(base).any() and not np.isnan(direct).any()
    assert auroc(table['y'], direct) > auroc(table['y'], base)


def test_bootstrap_metrics_intervals_and_paired_delta():
    table = synthetic()
    scores = {'ref': table['log_distance'] * -1, 'same': table['log_distance'] * -1}
    m = bootstrap_metrics(table['y'], scores, table['gene'], n_boot=50, seed=0, reference='ref')
    for _, r in m.iterrows():
        assert r['auroc_lo'] <= r['auroc'] <= r['auroc_hi']
        assert r['auprc_lo'] <= r['auprc_hi']
    same = m.set_index('score').loc['same']
    assert same['delta_auroc'] == 0 and same['delta_auroc_lo'] == 0 and same['delta_auroc_hi'] == 0


def test_lr_test_detects_informative_feature():
    table = synthetic(n_genes=150)
    informative = lr_test(table, ['log_distance', 'element_degree'], ['direct'], 'x')
    assert informative['lr_p'] < 1e-3 and informative['coefficient'] > 0
    assert informative['df'] == 1 and 0 <= informative['wald_cluster_gene_p'] <= 1
    joint = lr_test(table, ['log_distance'], ['hops', 'ppr'], 'x')
    assert joint['df'] == 2 and np.isnan(joint['coefficient'])


def test_split_overlap_counts():
    a = pd.DataFrame(dict(gene=['g1', 'g2'], element=['chr1:100-200', 'chr1:500-600'], chrom=['chr1', 'chr1'],
                          chromStart=[100, 500], chromEnd=[200, 600]))
    b = pd.DataFrame(dict(gene=['g1', 'g3', 'g3'], element=['chr1:150-250', 'chr1:500-600', 'chr2:1-10'],
                          chrom=['chr1', 'chr1', 'chr2'], chromStart=[150, 500, 1], chromEnd=[250, 600, 10]))
    o = split_overlap(a, b)
    assert o['shared_genes'] == 1 and o['identical_elements'] == 1
    assert o['heldout_elements_overlapping_training'] == 2 and o['identical_pairs'] == 0
    assert o['heldout_pairs_of_shared_genes'] == 1


def test_prediction_found():
    pairs = pd.DataFrame(dict(measuredGeneSymbol=['A', 'A', 'B'], chrom=['chr1', 'chr1', 'chr1'],
                              chromStart=[100, 1000, 100], chromEnd=[200, 1100, 200]))
    predictions = pd.DataFrame(dict(TargetGene=['A'], chrom=['chr1'], start=[150], end=[300]))
    assert list(prediction_found(pairs, predictions)) == [True, False, False]
