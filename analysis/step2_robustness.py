"""Robustness of the step-2 (indirect) gain in R^2 of the mark term (egaDifferential --stepFeature decomposed).

Usage: python step2_robustness.py <prefix> <tracksA> <tracksB> <out.tsv> [--hox chr:start-end ...] [--threads N]
       [--bootstrap N] [--seed S]

<prefix> is an egaDifferential run with --stepFeature decomposed (<prefix>_union.graphml, _gene_features.tsv,
_elements.tsv); tracksA / tracksB are mark=file.bw lists (comma-separated, one argument each, as --tracksA / --tracksB)
for the mirrored loci. The gain of a level is the R^2 of the nested model (promoter, + direct, + indirect; marks of
the run) over the model without it, on the informative genes (nodes at one and at two steps), as nested_models; the
mirror gain replaces the indirect columns by those of the mirrored loci, as mirror_comparison. Rows:

  observed                the run's features (check: equals <prefix>_nested_models.tsv)
  non-promoter step 2     step-2 sum restricted to nodes without a TSS (element type 'non-promoter element' of
                          <prefix>_elements.tsv); mirror over the same restricted nodes
  loco                    leave one chromosome out (one row per chromosome), and the range
  no Hox genes            genes whose promoter fragment lies in a --hox window dropped
  no Hox genes or nodes   in addition the step-2 nodes in a --hox window carry no signal
  no Cook top 1%          genes of the 1% largest Cook's distances of the full model (promoter + direct + indirect)
                          dropped
  block bootstrap         chromosomes resampled with replacement (N reps): 2.5 / 97.5 percentiles and fraction > 0 of
                          the indirect gain and of graph minus mirror
"""

import argparse
import re

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.annotation import parse_track_arguments
from epigeneticgraphanalyzer.differential import ActivityOperators, _ols_r2, mirrored_features
from epigeneticgraphanalyzer.graph import read_graph


def parse():
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('tracksA')
    p.add_argument('tracksB')
    p.add_argument('out')
    p.add_argument('--hox', nargs='*', default=[])
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--bootstrap', type=int, default=1000)
    p.add_argument('--seed', type=int, default=0)
    return p.parse_args()


class Restricted:
    """ActivityOperators view whose step levels keep only the nodes of `mask` at the levels in `levels`."""

    def __init__(self, operators, mask, levels):
        import scipy.sparse as sp
        self.ops = operators
        self.diag = sp.diags(mask.astype(float))
        self.levels_masked = set(levels)

    def __getattr__(self, name):
        return getattr(self.ops, name)

    def level(self, level):
        H_a, H_b, n = self.ops.level(level)
        if level in self.levels_masked:
            return (H_a @ self.diag).tocsr(), (H_b @ self.diag).tocsr(), n
        return H_a, H_b, n

    def condition_level(self, level):
        H_a, H_b, n = self.ops.condition_level(level)
        if level in self.levels_masked:
            return (H_a @ self.diag).tocsr(), (H_b @ self.diag).tocsr(), n
        return H_a, H_b, n


def gains(df, marks, indirect_suffix='', mirror_suffix=None):
    """(gain of direct, gain of indirect, gain of the indirect mirror) on the rows of df."""
    y = df['rna_log2FC'].to_numpy(dtype=float)

    def r2(columns):
        return _ols_r2(y, np.nan_to_num(df[columns].to_numpy(dtype=float)))[0]
    promoter = ['{}_promoter'.format(m) for m in marks]
    direct = promoter + ['{}_direct'.format(m) for m in marks]
    r_p, r_d = r2(promoter), r2(direct)
    r_i = r2(direct + ['{}_indirect{}'.format(m, indirect_suffix) for m in marks])
    out = dict(n_genes=len(df), r2_promoter_direct=r_d, gain_direct=r_d - r_p, gain_indirect=r_i - r_d)
    if mirror_suffix is not None:
        r_m = r2(direct + ['{}_indirect{}'.format(m, mirror_suffix) for m in marks])
        out['gain_indirect_mirror'] = r_m - r_d
        out['graph_minus_mirror'] = out['gain_indirect'] - out['gain_indirect_mirror']
    return out


def indirect_columns(operators, vectors, marks, suffix, pseudo=1.0):
    H_a, H_b, _ = operators.level('indirect')
    return {'{}_indirect{}'.format(m, suffix):
            np.log2((H_b @ vectors[m][1] + pseudo) / (H_a @ vectors[m][0] + pseudo)) for m in marks}


def in_windows(nodes, windows):
    hit = np.zeros(len(nodes), dtype=bool)
    for i, n in enumerate(nodes):
        chrom, start, end = re.split('[:-]', n)
        for c, s, e in windows:
            if chrom == c and int(end) >= s and int(start) <= e:
                hit[i] = True
    return hit


def cooks_distance(df, columns):
    y = df['rna_log2FC'].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(y)), np.nan_to_num(df[columns].to_numpy(dtype=float))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    residual = y - X @ beta
    p = X.shape[1]
    q, _ = np.linalg.qr(X)
    leverage = (q ** 2).sum(axis=1)
    s2 = residual @ residual / (len(y) - p)
    return residual ** 2 / (p * s2) * leverage / (1 - leverage) ** 2


def main():
    args = parse()
    windows = []
    for w in args.hox:
        chrom, start, end = re.split('[:-]', w)
        windows.append((chrom, int(float(start)), int(float(end))))
    features = pd.read_csv(args.prefix + '_gene_features.tsv', sep='\t')
    marks = [c[:-len('_promoter')] for c in features if c.endswith('_promoter')]
    informative = features[(features['n_direct'] > 0) & (features['n_indirect'] > 0)].reset_index(drop=True)
    graph = read_graph(args.prefix + '_union.graphml')
    with open(args.prefix + '_support.json') as handle:
        import json
        name_a, name_b = json.load(handle)['names']
    for n in graph.nodes:
        graph.nodes[n]['start'] = int(graph.nodes[n]['start'])
        graph.nodes[n]['end'] = int(graph.nodes[n]['end'])
        graph.nodes[n]['chrom'] = n.split(':')[0]
    operators = ActivityOperators(graph, list(features['node']), max_hops=2, step_feature='marks')
    vectors = {m: (operators.vector({n: float(d.get(m + '_' + name_a, 0.0)) for n, d in graph.nodes(data=True)}),
                   operators.vector({n: float(d.get(m + '_' + name_b, 0.0)) for n, d in graph.nodes(data=True)}))
               for m in marks}
    row = {n: i for i, n in enumerate(operators.sources)}
    rows = np.array([row[n] for n in informative['node']])
    # check: the recomputed mark term equals the run's
    check = indirect_columns(operators, vectors, marks, '_check')
    for m in marks:
        informative[m + '_indirect_check'] = check[m + '_indirect_check'][rows]
        diff = np.nanmax(np.abs(informative[m + '_indirect_check'] - informative[m + '_indirect']))
        print('{} indirect recomputed vs run: max abs difference {:.2e}'.format(m, diff))
    results = []
    observed = gains(informative, marks, '', '_mirror')
    results.append(dict(analysis='observed', **observed))

    # (a) step-2 nodes without a TSS
    elements = pd.read_csv(args.prefix + '_elements.tsv', sep='\t', usecols=['node', 'condition', 'element_type'])
    elements = elements[elements['condition'] == name_a].set_index('node')['element_type']
    non_promoter = elements.reindex(operators.nodes).fillna('non-promoter element').to_numpy() == \
        'non-promoter element'
    restricted = Restricted(operators, non_promoter, ['indirect'])
    for m, values in indirect_columns(restricted, vectors, marks, '_nonpromoter').items():
        informative[m] = values[rows]
    tracks_a = parse_track_arguments(args.tracksA.split())
    tracks_b = parse_track_arguments(args.tracksB.split())
    mirror = mirrored_features(graph, restricted, tracks_a, tracks_b, threads=args.threads)
    for m in marks:
        informative[m + '_indirect_nonpromoter_mirror'] = mirror['{}_indirect_mirror'.format(m)].reindex(
            informative['node']).to_numpy()
    n_np = np.asarray((abs(restricted.level('indirect')[0]) > 0).sum(axis=1)).ravel()[rows]
    results.append(dict(analysis='non-promoter step 2', note='genes with non-promoter step-2 nodes: {}'.format(
        int((n_np > 0).sum())), **gains(informative, marks, '_nonpromoter', '_nonpromoter_mirror')))
    has = n_np > 0
    results.append(dict(analysis='non-promoter step 2, genes with such nodes',
                        **gains(informative[has].reset_index(drop=True), marks, '_nonpromoter',
                                '_nonpromoter_mirror')))
    # same gene subset with all step-2 nodes, for comparison
    results.append(dict(analysis='all step 2, genes with non-promoter step-2 nodes',
                        **gains(informative[has].reset_index(drop=True), marks, '', '_mirror')))

    # (b) leave one chromosome out
    chrom = informative['node'].str.split(':').str[0]
    loco = []
    for c in sorted(chrom.unique()):
        g = gains(informative[chrom != c].reset_index(drop=True), marks, '', '_mirror')
        loco.append(g['gain_indirect'])
        results.append(dict(analysis='loco', dropped=c, n_dropped=int((chrom == c).sum()), **g))
    i_min, i_max = int(np.argmin(loco)), int(np.argmax(loco))
    chroms = sorted(chrom.unique())
    results.append(dict(analysis='loco range', gain_indirect=np.nan,
                        note='min {:.5f} (without {}), max {:.5f} (without {})'.format(
                            loco[i_min], chroms[i_min], loco[i_max], chroms[i_max])))
    if windows:
        hox_gene = in_windows(list(informative['node']), windows)
        kept = informative[~hox_gene].reset_index(drop=True)
        results.append(dict(analysis='no Hox genes', n_dropped=int(hox_gene.sum()),
                            note=','.join(informative.loc[hox_gene, 'genes'].astype(str)),
                            **gains(kept, marks, '', '_mirror')))
        hox_node = in_windows(operators.nodes, windows)
        no_hox = Restricted(operators, ~hox_node, ['indirect'])
        for m, values in indirect_columns(no_hox, vectors, marks, '_nohox').items():
            informative[m] = values[rows]
        mirror = mirrored_features(graph, no_hox, tracks_a, tracks_b, threads=args.threads)
        for m in marks:
            informative[m + '_indirect_nohox_mirror'] = mirror['{}_indirect_mirror'.format(m)].reindex(
                informative['node']).to_numpy()
        kept = informative[~hox_gene].reset_index(drop=True)
        n_affected = int((np.asarray(operators.level('indirect')[0][:, np.nonzero(hox_node)[0]].sum(axis=1))
                          .ravel()[rows] > 0).sum())
        results.append(dict(analysis='no Hox genes or step-2 nodes', n_dropped=int(hox_gene.sum()),
                            note='genes with Hox step-2 nodes: {}'.format(n_affected),
                            **gains(kept, marks, '_nohox', '_nohox_mirror')))
    # Cook's distance of the full model
    full = ['{}_{}'.format(m, lv) for lv in ('promoter', 'direct', 'indirect') for m in marks]
    cook = cooks_distance(informative, full)
    top = cook >= np.quantile(cook, 0.99)
    results.append(dict(analysis="no Cook's top 1%", n_dropped=int(top.sum()),
                        note=','.join(informative.loc[top, 'genes'].astype(str).head(80)),
                        **gains(informative[~top].reset_index(drop=True), marks, '', '_mirror')))

    # (c) chromosome-block bootstrap
    rng = np.random.default_rng(args.seed)
    groups = {c: np.nonzero((chrom == c).to_numpy())[0] for c in chroms}
    boot = []
    for _ in range(args.bootstrap):
        pick = rng.choice(chroms, size=len(chroms), replace=True)
        index = np.concatenate([groups[c] for c in pick])
        boot.append(gains(informative.iloc[index].reset_index(drop=True), marks, '', '_mirror'))
    boot = pd.DataFrame(boot)
    for column in ('gain_indirect', 'graph_minus_mirror', 'gain_direct'):
        values = boot[column].to_numpy()
        results.append(dict(analysis='block bootstrap ' + column, n_genes=len(informative),
                            **{column: observed[column]},
                            note='2.5% {:.5f}, 97.5% {:.5f}, fraction > 0 {:.3f}, mean {:.5f}, {} reps, {} '
                                 'chromosomes'.format(np.quantile(values, 0.025), np.quantile(values, 0.975),
                                                      float(np.mean(values > 0)), float(values.mean()),
                                                      args.bootstrap, len(chroms))))
    table = pd.DataFrame(results)
    table.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 20, 'display.max_rows', 200,
                           'display.max_colwidth', 120):
        print(table.to_string(index=False))


if __name__ == '__main__':
    main()
