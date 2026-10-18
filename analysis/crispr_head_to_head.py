"""Head-to-head comparison of EGA graph scores with published enhancer-gene predictions on identical CRISPRi pairs
(revision analysis B3).

Scores (higher = predicted regulatory), all on the same pairs of crispr_benchmark.py tables:

  distance             -|distanceToTSS|
  abc                  ABC score (ABC.Score.Feature) of the published K562 ENCODE-rE2G predictions
  encode_re2g          ENCODE-rE2G score (Score) of the same predictions
  powerlaw             power-law 3D contact estimate (3DContact.Feature) of the same predictions
  ega_direct           CHiCAGO score of the direct edge element - TSS fragment (0 without edge)
  ega_two_step         strongest two-step path element - mediator - TSS fragment
  ega_ppr              personalised PageRank from the TSS fragment
  ega_logistic_M1      logistic model log distance + degree + direct (crispr_cv.py)
  ega_logistic_M3      logistic model log distance + degree + direct + two-step + hops + community + PPR

The logistic models are never evaluated on pairs they were fitted on: on the training pairs the out-of-fold
predictions of 5-fold GroupKFold by gene, on the held-out pairs the models fitted on all training pairs.
Published predictions are used as released (no tool is rerun). The ENCODE-rE2G model was trained on the K562
training pairs, so its scores there are in-sample; the held-out K562 pairs are new to every method.

Coverage: published scores exist for pairs whose element overlaps a predicted candidate element of the gene; EGA
scores exist for pairs with the TSS fragment a bait with interactions and the element off the TSS fragment.
Evaluation on the intersection of covered pairs and on all pairs with uncovered pairs given the lowest score, each
for all pairs, pairs with the element off the TSS fragment, and pairs more than 20 kb from the TSS. AUROC and
AUPRC with gene-bootstrap 95% CIs (crispr_cv.bootstrap_metrics), paired differences to distance and to ENCODE-rE2G.

Usage: python crispr_head_to_head.py <training_pairs.tsv.gz> <heldout_pairs.tsv.gz> <crispr_cv output dir>
       <output dir> [--threads 6] [--bootstrap 1000] [--seed 0]
"""

import argparse
import os

import numpy as np
import pandas as pd

from crispr_cv import bootstrap_metrics, prepare

PUBLISHED = {'abc': 'abc', 'encode_re2g': 're2g', 'powerlaw': 'powerlaw'}
GRAPH = {'ega_direct': 'direct', 'ega_two_step': 'two_hop', 'ega_ppr': 'ppr'}
LOGISTIC = {'ega_logistic_M1': 'M1_direct', 'ega_logistic_M3': 'M3_all_graph'}
METHODS = ['distance'] + list(PUBLISHED) + list(GRAPH) + list(LOGISTIC)
MIN_DISTANCE = 20000


def score_table(pairs, logistic):
    """Per-pair score of every method, NaN where the method cannot score the pair.

    pairs: prepared crispr_benchmark.py table; logistic: DataFrame indexed like the scorable pairs' 'name' with
    one column per logistic model (out-of-sample predictions).
    """
    scores = pd.DataFrame(index=pairs.index)
    scores['distance'] = -pairs['distanceToTSS'].abs().astype(float)
    found = pairs['prediction_found'].astype(str).str.upper() == 'TRUE'
    for method, column in PUBLISHED.items():
        scores[method] = pairs[column].astype(float).where(found)
    for method, column in GRAPH.items():
        scores[method] = pairs[column].astype(float).where(pairs['scorable'])
    for method, model in LOGISTIC.items():
        values = pairs['name'].map(logistic[model])
        scores[method] = values.where(pairs['scorable'])
    return scores


def lowest_for_missing(scores):
    """Missing scores replaced by a value below every observed score of the method."""
    out = scores.copy()
    for column in out:
        observed = out[column].dropna()
        floor = observed.min() - 1.0 if len(observed) else 0.0
        out[column] = out[column].fillna(floor)
    return out


def coverage(scores, pairs):
    rows = []
    y = pairs['y'].to_numpy()
    for method in scores:
        covered = scores[method].notna().to_numpy()
        rows.append(dict(method=method, pairs=len(pairs), covered=int(covered.sum()),
                         fraction=float(covered.mean()), positives=int(y.sum()),
                         positives_covered=int(y[covered].sum())))
    return pd.DataFrame(rows)


def strata(pairs):
    return {'all_pairs': np.ones(len(pairs), bool),
            'element_off_tss_fragment': ~pairs['same_fragment'].to_numpy(),
            'distance_gt_20kb': pairs['distanceToTSS'].abs().to_numpy() > MIN_DISTANCE}


def evaluate(pairs, scores, label, n_boot, seed, threads):
    rows = []
    complete = scores.notna().all(axis=1).to_numpy()
    filled = lowest_for_missing(scores)
    for stratum, mask in strata(pairs).items():
        for mode, use, table in (('intersection', mask & complete, scores), ('missing_lowest', mask, filled)):
            sub = pairs[use]
            values = {m: table.loc[use, m].to_numpy() for m in METHODS}
            parts = []
            for reference in ('distance', 'encode_re2g'):
                m = bootstrap_metrics(sub['y'], values, sub['gene'], n_boot, seed, reference=reference,
                                      threads=threads)
                m = m.rename(columns={c: c.replace('delta_', 'delta_vs_' + reference + '_')
                                      for c in m.columns if c.startswith('delta_')})
                parts.append(m.set_index('score'))
            m = parts[0].join(parts[1][[c for c in parts[1].columns if c.startswith('delta_')]]).reset_index()
            m.insert(0, 'pair_set', label)
            m.insert(1, 'stratum', stratum)
            m.insert(2, 'evaluation', mode)
            m['pairs'] = len(sub)
            m['positives'] = int(sub['y'].sum())
            m['negatives'] = int(len(sub) - sub['y'].sum())
            m['genes'] = sub['gene'].nunique()
            m['elements'] = sub['element'].nunique()
            rows.append(m)
            print(m[['pair_set', 'stratum', 'evaluation', 'score', 'auroc', 'auroc_lo', 'auroc_hi', 'auprc',
                     'auprc_lo', 'auprc_hi', 'pairs', 'positives']].to_string(), flush=True)
    return pd.concat(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('training')
    parser.add_argument('heldout')
    parser.add_argument('cv_dir')
    parser.add_argument('out')
    parser.add_argument('--threads', type=int, default=6)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    oof = pd.read_csv(os.path.join(args.cv_dir, 'oof_predictions_training.tsv.gz'), sep='\t')
    ext = pd.read_csv(os.path.join(args.cv_dir, 'external_predictions_heldout.tsv.gz'), sep='\t')
    logistic_training = oof.set_index('name')[['gene_' + m for m in LOGISTIC.values()]]
    logistic_training.columns = list(LOGISTIC.values())
    logistic_heldout = ext.set_index('name')[list(LOGISTIC.values())]
    assert logistic_training.index.is_unique and logistic_heldout.index.is_unique

    results, coverages = [], []
    for label, path, logistic in (('training_K562', args.training, logistic_training),
                                  ('heldout_K562', args.heldout, logistic_heldout)):
        pairs = prepare(path)
        scores = score_table(pairs, logistic)
        pd.concat([pairs[['name', 'gene', 'element', 'distanceToTSS', 'same_fragment', 'y']], scores],
                  axis=1).to_csv(os.path.join(args.out, label + '_scores.tsv.gz'), sep='\t', index=False)
        c = coverage(scores, pairs)
        c.insert(0, 'pair_set', label)
        coverages.append(c)
        print(c.to_string(), flush=True)
        results.append(evaluate(pairs, scores, label, args.bootstrap, args.seed, args.threads))
    pd.concat(coverages).to_csv(os.path.join(args.out, 'coverage.tsv'), sep='\t', index=False)
    pd.concat(results).to_csv(os.path.join(args.out, 'head_to_head.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
