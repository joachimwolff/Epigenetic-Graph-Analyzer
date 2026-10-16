"""Chromosome block bootstrap of three results whose items (hubs, genes, pairs) share nodes and chromosomes.

Whole chromosomes are drawn with replacement (1,000 resamples, percentile 95% interval), as for the step curve
(epigeneticgraphanalyzer.validation.block_bootstrap). Where a statistic has a null value (0 for differences and
assortativities, 1 for ratios), `fraction_wrong_side` is the fraction of resamples on the null side of it (at or
beyond the null value, opposite to the point estimate).

  hubs     promoter hubs with the promoter and >= 2 partners in one SPRITE cluster against shifted hubs
           (analysis/hub_multiway.py per-hub table), per hub size class and overall: fractions, their difference
           and ratio.
  nested   Delta R^2 of the direct and the indirect level of the nested expression models
           (differential.nested_models on the egaDifferential gene feature table, genes with direct and indirect
           neighbours), refitted in every resample; adjusted Delta R^2 alongside.
  chaser   ChAseR assortativity (Newman's numeric assortativity, every edge counted in both directions) of the
           pair networks at 1 to 5 steps (analysis/chaser_inputs.py) and of distance-matched randomisations of
           the step-1 network (analysis/chaser_random_networks.R), and the difference step 1 minus random.
           The reimplementation is checked against the stored ChAseR values to 3 decimals first.

Usage:
  python block_bootstrap_results.py hubs <hubs.tsv.gz> <output prefix>
  python block_bootstrap_results.py nested <gene_features.tsv> <nested_models.tsv> <output prefix>
  python block_bootstrap_results.py chaser <chaser input dir> <stored chaser tsv> <random networks dir> <output prefix>
Options: --bootstrap N (default 1000), --seed S (default 0).
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.differential import nested_models
from epigeneticgraphanalyzer.genome import parse_node_id
from epigeneticgraphanalyzer.validation import (assortativity_from_moments, assortativity_moments, block_resample_items,
                                                block_resamples, edge_assortativity)

SIZE_CLASSES = ((3, 3), (4, 5), (6, 10), (11, 1000))
STEPS = range(1, 6)


def interval(estimates):
    estimates = np.asarray(estimates, dtype=float)
    return float(np.nanpercentile(estimates, 2.5)), float(np.nanpercentile(estimates, 97.5))


def wrong_side(estimate, estimates, null):
    """Fraction of resamples at the null value or on the side of it opposite to the point estimate."""
    estimates = np.asarray(estimates, dtype=float)
    estimates = estimates[~np.isnan(estimates)]
    if not len(estimates) or np.isnan(estimate):
        return np.nan
    return float(np.mean(estimates <= null) if estimate > null else np.mean(estimates >= null))


def describe(name, estimate, estimates, null=None):
    low, high = interval(estimates)
    row = {name: estimate, name + '_ci_low': low, name + '_ci_high': high}
    if null is not None:
        row[name + '_fraction_wrong_side'] = wrong_side(estimate, estimates, null)
    return row


def chromosome_of(nodes):
    return np.array([parse_node_id(n)[0] for n in nodes])


# --------------------------------------------------------------------------------------------------- hubs

def size_class(size):
    for low, high in SIZE_CLASSES:
        if low <= size <= high:
            return '{}-{}'.format(low, high) if high < 1000 else '>={}'.format(low)
    return None


def run_hubs(args):
    table = pd.read_csv(args.inputs[0], sep='\t')
    table['chrom'] = chromosome_of(table['promoter'])
    table['class'] = [size_class(s) for s in table['hub_size']]
    table['observed'] = (table['groups_promoter_plus2'] > 0).astype(float)
    labels, resamples = block_resamples(table['chrom'].to_numpy(), args.bootstrap, args.seed)
    rows, replicates = [], {}
    classes = [c for c in pd.unique(table['class'].dropna())]
    order = [size_class(low) for low, _ in SIZE_CLASSES]
    for name in [c for c in order if c in classes] + ['all']:
        sub = table if name == 'all' else table[table['class'] == name]
        # additive per-chromosome sums: hubs, hubs supported, mean shifted support
        index = np.searchsorted(labels, sub['chrom'].to_numpy())
        n = np.bincount(index, minlength=len(labels)).astype(float)
        hit = np.bincount(index, weights=sub['observed'], minlength=len(labels))
        shifted = np.bincount(index, weights=sub['shifted_any_promoter_plus2'], minlength=len(labels))
        weights = np.array([np.bincount(chosen, minlength=len(labels)) for chosen in resamples], dtype=float)
        bn, bh, bs = weights @ n, weights @ hit, weights @ shifted
        with np.errstate(invalid='ignore', divide='ignore'):
            b_obs, b_shift = bh / bn, bs / bn
            b_diff, b_ratio = b_obs - b_shift, b_obs / b_shift
        obs, shift = hit.sum() / n.sum(), shifted.sum() / n.sum()
        row = {'hub_size': name, 'hubs': int(n.sum()), 'chromosomes': int((n > 0).sum())}
        row.update(describe('hubs_fraction', obs, b_obs))
        row.update(describe('shifted_fraction', shift, b_shift))
        row.update(describe('difference', obs - shift, b_diff, null=0.0))
        row.update(describe('ratio', obs / shift, b_ratio, null=1.0))
        rows.append(row)
        replicates[name + '_difference'] = b_diff
        replicates[name + '_ratio'] = b_ratio
    result = pd.DataFrame(rows)
    result.insert(len(result.columns), 'n_bootstrap', args.bootstrap)
    result.to_csv(args.prefix + '_hub_support_bootstrap.tsv', sep='\t', index=False)
    pd.DataFrame(replicates).to_csv(args.prefix + '_hub_support_bootstrap_replicates.tsv.gz', sep='\t', index=False)
    print(result.to_string(index=False))


# --------------------------------------------------------------------------------------------------- nested

def nested_statistics(features, marks):
    models = nested_models(features, marks).set_index('model')
    n = int(models['n_genes'].iloc[0])

    def adjusted(row):
        return 1 - (1 - row['r2']) * (n - 1) / (n - row['n_parameters'])

    p, pd_, pdi = (models.loc[m] for m in ('promoter', 'promoter+direct', 'promoter+direct+indirect'))
    return {'r2_promoter': p['r2'], 'r2_promoter_direct': pd_['r2'], 'r2_promoter_direct_indirect': pdi['r2'],
            'delta_r2_direct': pd_['delta_r2'], 'delta_r2_indirect': pdi['delta_r2'],
            'delta_adj_r2_direct': adjusted(pd_) - adjusted(p), 'delta_adj_r2_indirect': adjusted(pdi) - adjusted(pd_)}


def run_nested(args):
    features_path, stored_path = args.inputs[:2]
    features = pd.read_csv(features_path, sep='\t')
    marks = [c[:-len('_promoter')] for c in features.columns if c.endswith('_promoter')]
    informative = features[(features['n_direct'] > 0) & (features['n_indirect'] > 0)].dropna(subset=['rna_log2FC'])
    informative = informative.reset_index(drop=True)
    observed = nested_statistics(informative, marks)
    stored = pd.read_csv(stored_path, sep='\t').set_index('model')
    check = {'delta_r2_direct': stored.loc['promoter+direct', 'delta_r2'],
             'delta_r2_indirect': stored.loc['promoter+direct+indirect', 'delta_r2'],
             'r2_promoter': stored.loc['promoter', 'r2']}
    for key, value in check.items():
        if not np.isclose(observed[key], value, rtol=1e-6, atol=1e-9):
            raise SystemExit('{}: recomputed {} = {} differs from the stored {}'.format(features_path, key,
                                                                                       observed[key], value))
    chroms = chromosome_of(informative['node'])
    boot = [nested_statistics(informative.iloc[items], marks)
            for items in block_resample_items(chroms, args.bootstrap, args.seed)]
    boot = pd.DataFrame(boot)
    row = {'comparison': os.path.basename(features_path).replace('_gene_features.tsv', ''), 'marks': ','.join(marks),
           'n_genes': len(informative), 'chromosomes': len(np.unique(chroms))}
    for key in ('r2_promoter', 'r2_promoter_direct', 'r2_promoter_direct_indirect'):
        row.update(describe(key, observed[key], boot[key]))
    for key in ('delta_r2_direct', 'delta_r2_indirect', 'delta_adj_r2_direct', 'delta_adj_r2_indirect'):
        row.update(describe(key, observed[key], boot[key], null=0.0))
    row['n_bootstrap'] = args.bootstrap
    pd.DataFrame([row]).to_csv(args.prefix + '_nested_models_bootstrap.tsv', sep='\t', index=False)
    boot.to_csv(args.prefix + '_nested_models_bootstrap_replicates.tsv.gz', sep='\t', index=False)
    print(pd.Series(row).to_string())


# --------------------------------------------------------------------------------------------------- chaser

def read_network(path, features):
    """Unique undirected edges (as ChAseR's binary adjacency) with feature values at both ends, per chromosome.
    Nodes without features get 0 (missingv = 0 in analysis/chaser_assortativity.R)."""
    edges = pd.read_csv(path, sep='\t', header=None, dtype={0: str, 3: str})
    a = edges[0] + ':' + edges[1].astype(str) + '-' + edges[2].astype(str)
    b = edges[3] + ':' + edges[4].astype(str) + '-' + edges[5].astype(str)
    keep = (a != b).to_numpy()
    a, b, chrom = a[keep].to_numpy(), b[keep].to_numpy(), edges[0][keep].to_numpy()
    if np.any(chrom != edges[3][keep].to_numpy()):
        raise SystemExit('{}: trans edges are not expected'.format(path))
    low, high = np.where(a < b, a, b), np.where(a < b, b, a)
    unique = ~pd.DataFrame({'u': low, 'v': high}).duplicated().to_numpy()
    low, high, chrom = low[unique], high[unique], chrom[unique]
    values = features.reindex(np.concatenate([low, high])).fillna(0.0).to_numpy()
    return chrom, values[:len(low)], values[len(low):]


def network_moments(network, marks, labels):
    """moments[label index, mark, 4] of a network."""
    chrom, x_u, x_v = network
    moments = np.zeros((len(labels), len(marks), 4))
    index = np.searchsorted(labels, chrom)
    if np.any(labels[np.clip(index, 0, len(labels) - 1)] != chrom):
        raise SystemExit('chromosome missing from the label set')
    for k in np.unique(index):
        sel = index == k
        for j in range(len(marks)):
            moments[k, j] = assortativity_moments(x_u[sel, j], x_v[sel, j])
    return moments


def run_chaser(args):
    chaser_dir, stored_path, random_dir = args.inputs[:3]
    features = pd.read_csv(os.path.join(chaser_dir, 'features.txt'), sep='\t')
    marks = list(features.columns)
    stored = pd.read_csv(stored_path, sep='\t')
    networks = {'step {}'.format(s): read_network(os.path.join(chaser_dir, 'network_step{}.txt'.format(s)), features)
                for s in STEPS}
    random_paths = sorted(glob.glob(os.path.join(random_dir, 'random_step1_*.txt')))
    randoms = [read_network(p, features) for p in random_paths]
    # verification against ChAseR: stored values of the real networks and the per-randomisation R values
    checks = []
    for name, (chrom, x_u, x_v) in networks.items():
        for j, mark in enumerate(marks):
            ours = edge_assortativity(x_u[:, j], x_v[:, j])
            ref = stored[(stored['network'] == name) & (stored['mark'] == mark)]['assortativity'].iloc[0]
            checks.append(dict(network=name, mark=mark, python=ours, chaser=ref))
    chas_path = os.path.join(random_dir, 'random_step1_chas.tsv')
    random_chas = pd.read_csv(chas_path, sep='\t') if os.path.exists(chas_path) else None
    for r, (chrom, x_u, x_v) in enumerate(randoms):
        for j, mark in enumerate(marks):
            ours = edge_assortativity(x_u[:, j], x_v[:, j])
            ref = np.nan
            if random_chas is not None:
                ref = random_chas[(random_chas['randomisation'] == r + 1) &
                                  (random_chas['mark'] == mark)]['assortativity'].iloc[0]
            checks.append(dict(network='random {}'.format(r + 1), mark=mark, python=ours, chaser=ref))
    checks = pd.DataFrame(checks)
    checks['agree_3_decimals'] = np.round(checks['python'], 3) == np.round(checks['chaser'], 3)
    checks.to_csv(args.prefix + '_chaser_reimplementation_check.tsv', sep='\t', index=False)
    stored_random = stored[stored['network'].str.contains('random')].set_index('mark')
    print(checks.to_string(index=False))
    if not checks['agree_3_decimals'].all():
        raise SystemExit('reimplementation does not reproduce ChAseR to 3 decimals')

    labels = np.unique(np.concatenate([n[0] for n in networks.values()] + [n[0] for n in randoms]))
    real = {name: network_moments(net, marks, labels) for name, net in networks.items()}
    random_moments = np.stack([network_moments(net, marks, labels) for net in randoms])  # (R, labels, marks, 4)
    _, resamples = block_resamples(labels, args.bootstrap, args.seed)
    weights = np.array([np.bincount(chosen, minlength=len(labels)) for chosen in resamples], dtype=float)

    def assort(moments):            # (..., marks, 4) -> (..., marks)
        return np.apply_along_axis(lambda m: assortativity_from_moments(*m), -1, moments)

    boot_random = assort(np.einsum('bl,rlmk->brmk', weights, random_moments)).mean(axis=1)   # (B, marks)
    obs_random = assort(random_moments.sum(axis=1)).mean(axis=0)
    rows, replicates = [], {}
    for name, moments in real.items():
        observed = assort(moments.sum(axis=0))
        boot = assort(np.einsum('bl,lmk->bmk', weights, moments))
        for j, mark in enumerate(marks):
            row = {'network': name, 'mark': mark, 'edges': int(moments[:, j, 0].sum() / 2),
                   'chromosomes': int((moments[:, j, 0] > 0).sum())}
            row.update(describe('assortativity', observed[j], boot[:, j], null=0.0))
            replicates['{} {}'.format(name, mark)] = boot[:, j]
            if name == 'step 1':
                row.update(describe('random_assortativity', obs_random[j], boot_random[:, j]))
                row['random_assortativity_stored_chaser'] = stored_random.loc[mark, 'assortativity']
                row['random_sd_stored_chaser'] = stored_random.loc[mark, 'sd']
                row['random_sd_over_randomisations'] = float(np.std(assort(random_moments.sum(axis=1))[:, j], ddof=1))
                row.update(describe('difference_to_random', observed[j] - obs_random[j],
                                    boot[:, j] - boot_random[:, j], null=0.0))
                replicates['{} {} minus random'.format(name, mark)] = boot[:, j] - boot_random[:, j]
            rows.append(row)
    result = pd.DataFrame(rows)
    result['n_randomisations'] = len(randoms)
    result['n_bootstrap'] = args.bootstrap
    result.to_csv(args.prefix + '_chaser_bootstrap.tsv', sep='\t', index=False)
    pd.DataFrame(replicates).to_csv(args.prefix + '_chaser_bootstrap_replicates.tsv.gz', sep='\t', index=False)
    print(result.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('analysis', choices=['hubs', 'nested', 'chaser'])
    parser.add_argument('inputs', nargs='+')
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    args.prefix = args.inputs.pop()
    os.makedirs(os.path.dirname(os.path.abspath(args.prefix)), exist_ok=True)
    {'hubs': run_hubs, 'nested': run_nested, 'chaser': run_chaser}[args.analysis](args)


if __name__ == '__main__':
    main()
