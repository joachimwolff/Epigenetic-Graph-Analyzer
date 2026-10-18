"""Leakage-aware evaluation of graph features against CRISPRi-tested enhancer-gene pairs (revision analysis B2).

Input are the pair tables of crispr_benchmark.py for the two K562 pair sets of the ENCODE-rE2G CRISPRi benchmark
(EPCrisprBenchmark): the training pairs (EPCrisprBenchmark_combined_data.training_K562) and the K562 part of the
held-out pairs (EPCrisprBenchmark_combined_data.heldout_5_cell_types, CellType == K562). Only scorable pairs enter
the models: TSS fragment a bait with interactions in the graph, element not on the TSS fragment.

Nested logistic models, logit P(Regulated) = b0 + b1 log10(|d| + 1 kb) + b2 degree [+ b3 direct [+ b4 two_step
[+ b5 hops + b6 community + b7 ppr]]], see MODELS and FORMULAS. They are evaluated

  * out of fold on the training pairs with leave-one-chromosome-out CV (23 folds) and with GroupKFold by gene
    (5 folds; all pairs of a gene in one fold), so that no gene and no chromosome contributes to both fit and test;
  * externally: fitted on all training pairs, applied to the held-out pairs (other screens, no shared pairs).

AUROC and AUPRC with 95% percentile intervals from a bootstrap over genes (all pairs of a resampled gene enter
together); differences to the distance + degree model use the same resamples (paired). Coefficient tables with
model-based and gene-clustered standard errors; likelihood-ratio tests of each nested step, also restricted to
pairs without a direct edge.

Usage: python crispr_cv.py <training_pairs.tsv.gz> <heldout_pairs.tsv.gz> <output dir> [--threads 6]
       [--bootstrap 1000] [--seed 0]
"""

import argparse
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed
from scipy import stats
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

BASE = ['log_distance', 'element_degree']
MODELS = {
    'M0_distance_degree': BASE,
    'M1_direct': BASE + ['direct'],
    'M2_direct_two_step': BASE + ['direct', 'two_hop'],
    'M3_all_graph': BASE + ['direct', 'two_hop', 'hops', 'community', 'ppr'],
}
FORMULAS = {
    'M0_distance_degree': 'logit P(Regulated) = b0 + b1*log10(|d|+1000) + b2*degree',
    'M1_direct': 'logit P(Regulated) = b0 + b1*log10(|d|+1000) + b2*degree + b3*direct',
    'M2_direct_two_step': 'logit P(Regulated) = b0 + b1*log10(|d|+1000) + b2*degree + b3*direct + b4*two_step',
    'M3_all_graph': 'logit P(Regulated) = b0 + b1*log10(|d|+1000) + b2*degree + b3*direct + b4*two_step '
                    '+ b5*hops + b6*community + b7*ppr',
}
VARIABLES = {
    'log_distance': 'log10(|distanceToTSS| + 1000), distanceToTSS = element centre to TSS (bp)',
    'element_degree': 'largest number of graph edges of a fragment overlapping the element',
    'direct': 'largest CHiCAGO score of an edge element fragment - TSS fragment (0 without edge)',
    'two_hop': 'largest min(CHiCAGO score) over paths element - mediator - TSS fragment (0 without such path)',
    'hops': '1 / shortest path length element - TSS fragment (0 if > 6 edges or disconnected)',
    'community': '1 if element and TSS fragment are in the same Leiden community, else 0',
    'ppr': 'personalised PageRank of the element fragment, restart at the TSS fragment (score-weighted)',
}


def prepare(path):
    """Scorable pairs of a crispr_benchmark.py table with gene, element and chromosome identifiers."""
    table = pd.read_csv(path, sep='\t')
    for column in ('tss_in_graph', 'same_fragment', 'Regulated'):
        table[column] = table[column].astype(str).str.upper() == 'TRUE'
    table['gene'] = table['measuredGeneSymbol']
    table['element'] = table['chrom'] + ':' + table['chromStart'].astype(str) + '-' + table['chromEnd'].astype(str)
    table['y'] = table['Regulated'].astype(int)
    table['scorable'] = table['tss_in_graph'] & ~table['same_fragment']
    return table


def describe(table, label):
    """Counts of pairs, positives, negatives, independent genes and elements of a pair set."""
    return dict(set=label, pairs=len(table), positives=int(table['y'].sum()),
                negatives=int((1 - table['y']).sum()), genes=table['gene'].nunique(),
                genes_with_positive=table.loc[table['y'] == 1, 'gene'].nunique(),
                elements=table['element'].nunique(), chromosomes=table['chrom'].nunique(),
                datasets=';'.join('{}:{}'.format(k, v) for k, v in table['Dataset'].value_counts().items()))


def split_overlap(training, heldout):
    """Shared genes, identical elements, overlapping elements and identical pairs between two pair sets."""
    shared_genes = set(training['gene']) & set(heldout['gene'])
    identical_elements = set(training['element']) & set(heldout['element'])
    pairs_a = set(zip(training['element'], training['gene']))
    pairs_b = set(zip(heldout['element'], heldout['gene']))
    overlapping = 0
    by_chrom = {c: g for c, g in training.groupby('chrom')}
    for element in heldout.drop_duplicates('element').itertuples():
        t = by_chrom.get(element.chrom)
        if t is not None and np.any((t['chromEnd'] > element.chromStart) & (t['chromStart'] < element.chromEnd)):
            overlapping += 1
    return dict(heldout_genes=heldout['gene'].nunique(), shared_genes=len(shared_genes),
                heldout_pairs_of_shared_genes=int(heldout['gene'].isin(shared_genes).sum()),
                heldout_elements=heldout['element'].nunique(), identical_elements=len(identical_elements),
                heldout_elements_overlapping_training=overlapping, identical_pairs=len(pairs_a & pairs_b))


def fit_logit(table, features, groups=None):
    """Unpenalised logistic fit; with `groups` gene-clustered (sandwich) covariance."""
    x = sm.add_constant(table[features].astype(float).to_numpy(), has_constant='add')
    y = table['y'].to_numpy()
    if groups is None:
        return sm.Logit(y, x).fit(disp=0, maxiter=200)
    return sm.Logit(y, x).fit(disp=0, maxiter=200, cov_type='cluster',
                              cov_kwds={'groups': pd.factorize(table[groups])[0]})


def predict(fit, table, features):
    x = sm.add_constant(table[features].astype(float).to_numpy(), has_constant='add')
    return fit.predict(x)


def folds(table, scheme, seed=0):
    """Train/test index pairs: 'chromosome' = leave one chromosome out, 'gene' = 5-fold GroupKFold by gene."""
    if scheme == 'chromosome':
        return list(LeaveOneGroupOut().split(table, groups=table['chrom']))
    if scheme == 'gene':
        return list(GroupKFold(n_splits=5, shuffle=True, random_state=seed).split(table, groups=table['gene']))
    raise ValueError(scheme)


def out_of_fold(table, features, splits):
    """Out-of-fold predicted probabilities; each pair is predicted by a model fitted without its fold."""
    prediction = np.full(len(table), np.nan)
    for train, test in splits:
        fit = fit_logit(table.iloc[train], features)
        prediction[test] = predict(fit, table.iloc[test], features)
    return prediction


def auroc(y, score):
    """Area under the ROC curve (Mann-Whitney, ties counted half)."""
    y = np.asarray(y).astype(bool)
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = stats.rankdata(score)
    return (ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def auprc(y, score):
    y = np.asarray(y)
    if y.sum() == 0:
        return np.nan
    return average_precision_score(y, score)


def gene_resamples(genes, n, seed):
    """`n` bootstrap index vectors: genes drawn with replacement, every pair of a drawn gene included."""
    codes, unique = pd.factorize(pd.Series(genes))
    members = [np.flatnonzero(codes == i) for i in range(len(unique))]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        draw = rng.integers(0, len(unique), len(unique))
        out.append(np.concatenate([members[i] for i in draw]))
    return out


def _metrics_on(y, scores, index):
    yy = y[index]
    return {name: (auroc(yy, s[index]), auprc(yy, s[index])) for name, s in scores.items()}


def bootstrap_metrics(y, scores, genes, n_boot=1000, seed=0, reference=None, threads=1):
    """AUROC and AUPRC with gene-bootstrap 95% CIs for several scores of the same pairs.

    scores: dict name -> score vector (higher = more likely regulated). With `reference`, paired differences of
    every score to that score with CIs and the bootstrap fraction of resamples where the difference is <= 0.
    """
    y = np.asarray(y).astype(int)
    scores = {k: np.asarray(v, dtype=float) for k, v in scores.items()}
    point = _metrics_on(y, scores, np.arange(len(y)))
    resamples = gene_resamples(genes, n_boot, seed)
    chunks = np.array_split(np.arange(n_boot), max(1, threads))
    parts = Parallel(n_jobs=threads)(delayed(lambda ids: [_metrics_on(y, scores, resamples[i]) for i in ids])(c)
                                     for c in chunks)
    boot = [b for part in parts for b in part]
    rows = []
    for name in scores:
        row = dict(score=name)
        for k, metric in enumerate(('auroc', 'auprc')):
            values = np.array([b[name][k] for b in boot], dtype=float)
            row[metric] = point[name][k]
            row[metric + '_lo'], row[metric + '_hi'] = np.nanpercentile(values, [2.5, 97.5])
            if reference is not None and name != reference:
                ref = np.array([b[reference][k] for b in boot], dtype=float)
                diff = values - ref
                row['delta_' + metric] = point[name][k] - point[reference][k]
                row['delta_' + metric + '_lo'], row['delta_' + metric + '_hi'] = np.nanpercentile(diff, [2.5, 97.5])
                row['delta_' + metric + '_frac_le0'] = float(np.nanmean(diff <= 0))
        rows.append(row)
    return pd.DataFrame(rows)


def coefficient_table(table, model, label):
    """Coefficients, model-based and gene-clustered SEs, Wald P and 95% CI of one model."""
    features = MODELS[model]
    fit = fit_logit(table, features)
    robust = fit_logit(table, features, groups='gene')
    rows = []
    for i, name in enumerate(['intercept'] + features):
        rows.append(dict(set=label, model=model, term=name, coefficient=fit.params[i], se=fit.bse[i],
                         z=fit.tvalues[i], p=fit.pvalues[i], ci_lo=fit.conf_int()[i][0], ci_hi=fit.conf_int()[i][1],
                         se_cluster_gene=robust.bse[i], p_cluster_gene=robust.pvalues[i],
                         n=len(table), positives=int(table['y'].sum()), log_likelihood=fit.llf, aic=fit.aic))
    return rows


def lr_test(table, base, added, label):
    """Likelihood-ratio test of `added` features over `base`; with the gene-clustered Wald P of the added terms."""
    m0 = fit_logit(table, base)
    m1 = fit_logit(table, base + added)
    robust = fit_logit(table, base + added, groups='gene')
    lr = 2 * (m1.llf - m0.llf)
    k = len(added)
    restriction = np.zeros((k, len(base) + k + 1))
    for j in range(k):
        restriction[j, len(base) + 1 + j] = 1
    wald = robust.wald_test(restriction, scalar=True)
    return dict(set=label, base='+'.join(base), added='+'.join(added), n=len(table),
                positives=int(table['y'].sum()), genes=table['gene'].nunique(), lr_statistic=lr, df=k,
                lr_p=stats.chi2.sf(lr, k), wald_cluster_gene_statistic=float(wald.statistic),
                wald_cluster_gene_p=float(wald.pvalue),
                coefficient=m1.params[-1] if k == 1 else np.nan, se=m1.bse[-1] if k == 1 else np.nan)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('training')
    parser.add_argument('heldout')
    parser.add_argument('out')
    parser.add_argument('--threads', type=int, default=6)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    out = lambda name: os.path.join(args.out, name)  # noqa: E731

    training_all = prepare(args.training)
    heldout_all = prepare(args.heldout)
    training = training_all[training_all['scorable']].reset_index(drop=True)
    heldout = heldout_all[heldout_all['scorable']].reset_index(drop=True)

    # (a) what the two sets are
    description = pd.DataFrame([describe(training_all, 'training_all_pairs'), describe(training, 'training_scorable'),
                                describe(training[training['direct'] == 0], 'training_scorable_without_direct_edge'),
                                describe(heldout_all, 'heldout_all_pairs'), describe(heldout, 'heldout_scorable'),
                                describe(heldout[heldout['direct'] == 0], 'heldout_scorable_without_direct_edge')])
    description.to_csv(out('pair_sets.tsv'), sep='\t', index=False)
    print(description.to_string(), flush=True)
    overlap = pd.Series({**{'all_' + k: v for k, v in split_overlap(training_all, heldout_all).items()},
                         **{'scorable_' + k: v for k, v in split_overlap(training, heldout).items()}})
    overlap.to_csv(out('training_heldout_overlap.tsv'), sep='\t', header=False)
    print(overlap.to_string(), flush=True)

    with open(out('formulas.txt'), 'w') as handle:
        for model, formula in FORMULAS.items():
            handle.write('{}\t{}\n'.format(model, formula))
        for variable, meaning in VARIABLES.items():
            handle.write('{}\t{}\n'.format(variable, meaning))

    # (b) out-of-fold predictions on the training pairs and external predictions on the held-out pairs
    predictions = training[['name', 'gene', 'element', 'chrom', 'Dataset', 'y', 'direct']].copy()
    metrics = []
    for scheme in ('chromosome', 'gene'):
        splits = folds(training, scheme, args.seed)
        scores = {}
        for model, features in MODELS.items():
            scores[model] = out_of_fold(training, features, splits)
            predictions['{}_{}'.format(scheme, model)] = scores[model]
        for subset, mask in (('scorable', np.ones(len(training), bool)),
                             ('scorable_without_direct_edge', training['direct'].to_numpy() == 0)):
            sub = training[mask]
            m = bootstrap_metrics(sub['y'], {k: v[mask] for k, v in scores.items()}, sub['gene'], args.bootstrap,
                                  args.seed, reference='M0_distance_degree', threads=args.threads)
            m.insert(0, 'evaluation', 'training_oof_' + scheme)
            m.insert(1, 'subset', subset)
            for k, v in describe(sub, subset).items():
                if k not in ('set', 'datasets', 'chromosomes'):
                    m[k] = v
            m['folds'] = len(splits)
            metrics.append(m)
            print(m.to_string(), flush=True)
    predictions.to_csv(out('oof_predictions_training.tsv.gz'), sep='\t', index=False)

    external = heldout[['name', 'gene', 'element', 'chrom', 'Dataset', 'y', 'direct']].copy()
    scores = {}
    for model, features in MODELS.items():
        scores[model] = predict(fit_logit(training, features), heldout, features)
        external[model] = scores[model]
    external.to_csv(out('external_predictions_heldout.tsv.gz'), sep='\t', index=False)
    for subset, mask in (('scorable', np.ones(len(heldout), bool)),
                         ('scorable_without_direct_edge', heldout['direct'].to_numpy() == 0)):
        sub = heldout[mask]
        m = bootstrap_metrics(sub['y'], {k: v[mask] for k, v in scores.items()}, sub['gene'], args.bootstrap,
                              args.seed, reference='M0_distance_degree', threads=args.threads)
        m.insert(0, 'evaluation', 'fit_training_test_heldout')
        m.insert(1, 'subset', subset)
        for k, v in describe(sub, subset).items():
            if k not in ('set', 'datasets', 'chromosomes'):
                m[k] = v
        m['folds'] = 1
        metrics.append(m)
        print(m.to_string(), flush=True)
    pd.concat(metrics).to_csv(out('oof_metrics.tsv'), sep='\t', index=False)

    # (c) coefficients and likelihood-ratio tests
    coefficients = []
    for label, table in (('training_scorable', training), ('heldout_scorable', heldout)):
        for model in MODELS:
            coefficients += coefficient_table(table, model, label)
    coefficients = pd.DataFrame(coefficients)
    coefficients.to_csv(out('coefficients.tsv'), sep='\t', index=False)
    print(coefficients.to_string(), flush=True)

    tests = []
    for label, table in (('training_scorable', training), ('heldout_scorable', heldout)):
        tests.append(lr_test(table, BASE, ['direct'], label))
        tests.append(lr_test(table, BASE + ['direct'], ['two_hop'], label))
        tests.append(lr_test(table, BASE + ['direct', 'two_hop'], ['hops', 'community', 'ppr'], label))
        for feature in ('two_hop', 'hops', 'community', 'ppr'):
            tests.append(lr_test(table, BASE, [feature], label))
        no_direct = table[table['direct'] == 0]
        for feature in ('two_hop', 'hops', 'community', 'ppr'):
            tests.append(lr_test(no_direct, BASE, [feature], label + '_without_direct_edge'))
        tests.append(lr_test(no_direct, BASE, ['two_hop', 'hops', 'community', 'ppr'],
                             label + '_without_direct_edge'))
    tests = pd.DataFrame(tests)
    tests.to_csv(out('lr_tests.tsv'), sep='\t', index=False)
    print(tests.to_string(), flush=True)


if __name__ == '__main__':
    main()
