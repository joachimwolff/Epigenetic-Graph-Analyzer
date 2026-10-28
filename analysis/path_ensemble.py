"""Multi-way support of indirect pairs by the number of shortest paths and by personalised PageRank.

For the tested pairs of a step curve (hop_distance.py --pairTable), the number of distinct shortest paths between
the two ends (for two steps the number of mediators) and a symmetric personalised PageRank score (unweighted graph,
damping 0.85, mean of PPR_a(b) and PPR_b(a)) are computed on the graph. The NPMI concordance against the matched
controls is reported per class of the number of shortest paths (1, 2-3, > 3) and per quartile of the PageRank score,
with a chromosome block bootstrap. Over all pairs at two steps within 2 Mb of the graph, the share of pairs that
depend on a single mediator is reported as well.

Usage: python path_ensemble.py <graph.graphml> <pairs.tsv.gz> <output prefix> [--steps 2 3] [--bootstrap 1000]
"""

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.graph_nulls import adjacency, graph_distance_pairs, pair_ppr, walk_counts
from epigeneticgraphanalyzer.validation import block_bootstrap


def path_class(n_paths):
    n_paths = np.asarray(n_paths)
    return np.where(n_paths <= 1, '1', np.where(n_paths <= 3, '2-3', '>3'))


def stratum_rows(table, column, steps, n_boot):
    rows = []
    for stratum, group in table.groupby(column):
        values = group['npmi_concordance'].to_numpy()
        low, high = block_bootstrap(values, group['chrom'].to_numpy(), n_boot=n_boot, seed=steps)
        rows.append(dict(steps=steps, stratifier=column, stratum=stratum, n=int(np.sum(~np.isnan(values))),
                         concordance=float(np.nanmean(values)), ci_low=low, ci_high=high,
                         mean_npmi=float(group['npmi'].mean()), mean_npmi_control=float(group['npmi_control'].mean()),
                         median_distance=float(group['distance'].median()),
                         median_paths=float(group['n_paths'].median()), median_ppr=float(group['ppr'].median())))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('pairs')
    parser.add_argument('prefix')
    parser.add_argument('--steps', type=int, nargs='+', default=[2, 3])
    parser.add_argument('--bootstrap', type=int, default=1000)
    args = parser.parse_args()
    graph = read_graph(args.graph)
    A, nodes = adjacency(graph)
    index = {n: k for k, n in enumerate(nodes)}
    pairs = pd.read_csv(args.pairs, sep='\t')
    rows, tables, overall = [], [], []
    for steps in args.steps:
        table = pairs[pairs['steps'] == steps].reset_index(drop=True)
        a = table['node_a'].map(index).to_numpy()
        b = table['node_b'].map(index).to_numpy()
        table['n_paths'] = walk_counts(A, a, b, steps)
        table['ppr'] = pair_ppr(A, a, b)
        table['path_class'] = path_class(table['n_paths'])
        table['ppr_quartile'] = pd.qcut(table['ppr'], 4, labels=['Q1', 'Q2', 'Q3', 'Q4']).astype(str)
        rows.extend(stratum_rows(table, 'path_class', steps, args.bootstrap))
        rows.extend(stratum_rows(table, 'ppr_quartile', steps, args.bootstrap))
        rho_paths = float(stats.spearmanr(table['n_paths'], table['distance'])[0])
        rho_ppr = float(stats.spearmanr(table['ppr'], table['distance'])[0])
        rho_conc = float(stats.spearmanr(table['ppr'], table['npmi_concordance'], nan_policy='omit')[0])
        overall.append(dict(steps=steps, set='tested pairs', n=len(table),
                            single_path_fraction=float((table['n_paths'] == 1).mean()),
                            spearman_paths_distance=rho_paths, spearman_ppr_distance=rho_ppr,
                            spearman_ppr_concordance=rho_conc))
        tables.append(table)
        print(overall[-1], flush=True)
    # all pairs at two steps within 2 Mb: how many depend on a single mediator
    every = graph_distance_pairs(graph, 2000000, 2)
    every = every[every['steps'] == 2]
    counts = walk_counts(A, every['node_a'].map(index).to_numpy(), every['node_b'].map(index).to_numpy(), 2)
    overall.append(dict(steps=2, set='all pairs at two steps within 2 Mb', n=len(every),
                        single_path_fraction=float((counts == 1).mean()),
                        mediators_2_3_fraction=float(((counts >= 2) & (counts <= 3)).mean()),
                        mediators_over_3_fraction=float((counts > 3).mean()),
                        median_mediators=float(np.median(counts))))
    print(overall[-1], flush=True)
    pd.DataFrame(rows).to_csv(args.prefix + '_strata.tsv', sep='\t', index=False)
    pd.DataFrame(overall).to_csv(args.prefix + '_overall.tsv', sep='\t', index=False)
    pd.concat(tables).to_csv(args.prefix + '_pairs.tsv.gz', sep='\t', index=False)


if __name__ == '__main__':
    main()
