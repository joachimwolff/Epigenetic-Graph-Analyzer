"""Promoter hubs of the interaction graph against three-body contacts predicted by a polymer model from Hi-C.

Uses the GM12878 example of HLM-Nbody (Liu, Zhang and Hyeon 2021): for promoter viewpoints in two domains (5 kb
bins), the predicted three-body contact probabilities p(viewpoint, j, k) (<fit>.vpAAA-BBB.pm) and the observed SPRITE
three-body counts (sprite-data/TADxx.p3.vpAAA-BBB.txt). A bin pair (j, k) is a hub pair if both bins hold direct
graph partners of a viewpoint promoter fragment. Per viewpoint and pooled, the script reports the Spearman
correlation of prediction and observation and, in a linear model log(observed + 1) ~ log(predicted) + hub, whether
hub pairs co-occur more in SPRITE than the polymer model predicts.

Usage: python polymer_comparison.py <HLM-Nbody sprite dir> <graph.graphml> <output prefix>
"""

import glob
import os
import re
import sys

import networkx as nx
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from epigeneticgraphanalyzer.graph import read_graph

BIN = 5000


def domains(table):
    result = {}
    for line in open(table):
        match = re.match(r'#TAD(\d+): (chr\w+):(\d+)-(\d+) N: (\d+)', line)
        if match:
            result[int(match.group(1))] = (match.group(2), int(match.group(3)), int(match.group(5)))
    return result


def load_matrix(path):
    return np.loadtxt(path, comments='#')


def main():
    sprite_dir, graph_path, prefix = sys.argv[1:4]
    graph = read_graph(graph_path)
    info = domains(os.path.join(sprite_dir, 'sprite-data', 'CHROMATIX-suppTable3b.txt'))
    rows, pooled = [], []
    for observed_path in sorted(glob.glob(os.path.join(sprite_dir, 'sprite-data', 'TAD*.p3.vp*.txt'))):
        name = os.path.basename(observed_path)
        tad, first, last = map(int, re.match(r'TAD(\d+)\.p3\.vp(\d+)-(\d+)\.txt', name).groups())
        predicted_paths = glob.glob(os.path.join(sprite_dir, 'TAD{:02d}'.format(tad), '*.vp{:03d}-{:03d}.pm'.format(
            first, last)))
        if not predicted_paths:
            continue
        chrom, start, n = info[tad]
        observed = load_matrix(observed_path)
        predicted = load_matrix(predicted_paths[0])
        # graph nodes of every bin, and the direct partners of the viewpoint fragments
        bins = {}
        for node, data in graph.nodes(data=True):
            if data['chrom'] != chrom:
                continue
            middle = (int(data['start']) + int(data['end'])) // 2
            b = (middle - start) // BIN
            if 0 <= b < n:
                bins.setdefault(node, b)
        viewpoint = [v for v, b in bins.items() if first <= b <= last and graph.nodes[v].get('bait')]
        partner_bins = {bins[u] for v in viewpoint for u in graph.neighbors(v) if u in bins} - set(range(first, last + 1))
        hub = np.zeros((n, n), dtype=bool)
        for j in partner_bins:
            for k in partner_bins:
                if j != k:
                    hub[j, k] = True
        keep = np.triu(np.ones((n, n), dtype=bool), 1)
        for b in range(first, last + 1):
            keep[b, :] = keep[:, b] = False
        keep &= np.isfinite(observed) & np.isfinite(predicted) & (predicted > 0)
        o, p, h = observed[keep], predicted[keep], hub[keep]
        rho = stats.spearmanr(p, o).correlation
        frame = pd.DataFrame({'log_observed': np.log(o + 1), 'log_predicted': np.log(p), 'hub': h.astype(float),
                              'viewpoint': name})
        pooled.append(frame)
        row = dict(viewpoint=name, chrom=chrom, viewpoint_fragments=len(viewpoint), partner_bins=len(partner_bins),
                   pairs=int(keep.sum()), hub_pairs=int(h.sum()), spearman_polymer=rho)
        if h.sum() >= 3:
            fit = sm.OLS(frame['log_observed'], sm.add_constant(frame[['log_predicted', 'hub']])).fit()
            ratio = frame['log_observed'] - fit.params['log_predicted'] * frame['log_predicted']
            row.update(hub_coefficient=fit.params['hub'], hub_p=fit.pvalues['hub'],
                       spearman_hub=stats.spearmanr(h, o).correlation,
                       mannwhitney_p_residual=stats.mannwhitneyu(ratio[h], ratio[~h], alternative='greater').pvalue)
        rows.append(row)
        print(row, flush=True)
    table = pd.DataFrame(rows)
    if pooled:
        frame = pd.concat(pooled, ignore_index=True)
        design = pd.get_dummies(frame[['log_predicted', 'hub', 'viewpoint']], columns=['viewpoint'], drop_first=True,
                                dtype=float)
        fit = sm.OLS(frame['log_observed'], sm.add_constant(design)).fit()
        table = pd.concat([table, pd.DataFrame([dict(viewpoint='pooled (viewpoint as covariate)', pairs=len(frame),
                                                     hub_pairs=int(frame['hub'].sum()),
                                                     spearman_polymer=stats.spearmanr(frame['log_predicted'],
                                                                                      frame['log_observed']).correlation,
                                                     hub_coefficient=fit.params['hub'], hub_p=fit.pvalues['hub'])])])
    table.to_csv(prefix + '_polymer_vs_graph.tsv', sep='\t', index=False)
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
