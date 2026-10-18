"""Multi-way support of graph paths against degree-preserving rewired graphs.

The graph is rewired by double-edge swaps that keep every node degree, the bait / other-end kind of both ends of
every edge, the chromosome and the length bin of every edge (graph_nulls.rewire_by_distance). On the real graph and
on each rewired graph, pairs at graph distance 1, 2 and 3 within 2 Mb are evaluated against the same multi-way data
as in hop_distance.py (coverage-matched controls at the same offset, concordance of NPMI with a chromosome block
bootstrap). Controls of a rewired graph exclude the pairs connected within five steps in the rewired and in the real
graph. For a rewired graph, the pairs not connected within five steps in the real graph (novel paths) are also
evaluated on their own. Next to the concordance, absolute co-occurrence is reported for test and control pairs (unweighted groups):
`together` the share of pairs found together in at least one group, `share` the mean share of the groups holding
one end that also hold the other, `groups_both` the mean number of groups holding both ends. If paths of the real
graph only reflected the degree and the distance of their ends, the rewired graphs would reach the same values.

Usage: python rewired_null.py <graph.graphml> <multiway.npz> <output prefix> [--chain FILE] [--rewirings 10]
           [--pairs 50000] [--steps 1 2 3] [--swapsPerEdge 10] [--logBinWidth 0.1] [--maxLength 2000000]
           [--seed 0]
"""

import argparse
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hop_distance import MAX_DISTANCE, MAX_STEPS, evaluate_steps, step_pairs  # noqa: E402

from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402
from epigeneticgraphanalyzer.graph import read_graph  # noqa: E402
from epigeneticgraphanalyzer.graph_nulls import rewire_by_distance  # noqa: E402
from epigeneticgraphanalyzer.multiway import MultiwayData  # noqa: E402
from epigeneticgraphanalyzer.validation import LocusPools, map_nodes_to_bins  # noqa: E402


def cooccurrence_rates(data, bin_a, bin_b):
    """Unweighted groups holding both ends, whether there is one, and the mean share over both ends of the groups
    holding one end that also hold the other."""
    n_ab = data.cooccurrence(bin_a, bin_b, weighted=False)
    n_a = data.coverage(bin_a, weighted=False)
    n_b = data.coverage(bin_b, weighted=False)
    with np.errstate(divide='ignore', invalid='ignore'):
        share = 0.5 * (np.where(n_a > 0, n_ab / n_a, 0.0) + np.where(n_b > 0, n_ab / n_b, 0.0))
    return dict(groups_both=n_ab, together=(n_ab > 0).astype(float), share=share)


def _rates(rows, tables, data):
    for row in rows:
        table, controls = tables[row['steps']]
        test = cooccurrence_rates(data, table['bin_a'].to_numpy(), table['bin_b'].to_numpy())
        control = cooccurrence_rates(data, controls['bin_a'].to_numpy(), controls['bin_b'].to_numpy())
        for key in ('together', 'share', 'groups_both'):
            row[key] = float(np.mean(test[key]))
            row[key + '_control'] = float(np.mean(control[key]))
    return rows


def evaluate_graph(name, graph, node_bins, data, pools, real_keys, args):
    """Step rows of one graph; for a rewired graph also of its novel pairs (not connected within five steps in the
    real graph)."""
    pairs = step_pairs(graph, node_bins, data, MAX_DISTANCE, MAX_STEPS)
    keys = set(zip(pairs['bin_a'], pairs['bin_b']))
    excluded = keys | real_keys
    sets = [(name, pairs)]
    if real_keys:
        in_real = np.array([k in real_keys for k in zip(pairs['bin_a'], pairs['bin_b'])], dtype=bool)
        sets.append((name + '_novel', pairs[~in_real].reset_index(drop=True)))
    rows = []
    for label, subset_pairs in sets:
        results, tables = evaluate_steps(subset_pairs, data, pools, excluded, args.pairs, args.coverageClasses,
                                         n_boot=args.bootstrap, steps_list=args.steps, return_tables=True)
        results = _rates(results, tables, data)
        for row in results:
            row['graph'] = label
            subset = subset_pairs[subset_pairs['steps'] == row['steps']]
            row['shared_with_real'] = float(np.mean([k in real_keys for k in zip(subset['bin_a'], subset['bin_b'])])) \
                if real_keys else 1.0
        rows.extend(results)
    return rows, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('prefix')
    parser.add_argument('--chain')
    parser.add_argument('--rewirings', type=int, default=10)
    parser.add_argument('--pairs', type=int, default=50000)
    parser.add_argument('--steps', type=int, nargs='+', default=[1, 2, 3])
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--swapsPerEdge', type=float, default=10)
    parser.add_argument('--logBinWidth', type=float, default=0.1)
    parser.add_argument('--maxLength', type=int, default=MAX_DISTANCE, help='Longest new edge (bp)')
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
    liftover = LiftOver(args.chain) if args.chain else None
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    node_bins = map_nodes_to_bins(graph, data, liftover)
    pools = LocusPools(graph, node_bins, data)
    rows, rewiring = [], []
    real_rows, real_pairs = evaluate_graph('real', graph, node_bins, data, pools, set(), args)
    rows.extend(real_rows)
    real_keys = set(zip(real_pairs['bin_a'], real_pairs['bin_b']))
    for k in range(args.rewirings):
        start = time.time()
        rewired, statistics = rewire_by_distance(graph, swaps_per_edge=args.swapsPerEdge,
                                                 log_bin_width=args.logBinWidth, seed=args.seed + k,
                                                 max_length=args.maxLength)
        statistics.update(graph='rewired_{}'.format(k), seed=args.seed + k, seconds=time.time() - start)
        rewiring.append(statistics)
        print(statistics, flush=True)
        results, _ = evaluate_graph('rewired_{}'.format(k), rewired, node_bins, data, pools, real_keys, args)
        rows.extend(results)
        pd.DataFrame(rows).to_csv(args.prefix + '_steps.tsv', sep='\t', index=False)
        pd.DataFrame(rewiring).to_csv(args.prefix + '_rewiring.tsv', sep='\t', index=False)
    table = pd.DataFrame(rows)
    table.to_csv(args.prefix + '_steps.tsv', sep='\t', index=False)
    # real against the distribution over rewired graphs
    summary = []
    metrics = ['concordance', 'mean_npmi', 'mean_npmi_control', 'together', 'together_control', 'share',
               'share_control', 'groups_both', 'groups_both_control', 'median_distance', 'pairs_available']
    table['null'] = np.where(table['graph'] == 'real', 'real',
                             np.where(table['graph'].str.endswith('_novel'), 'rewired_novel', 'rewired'))
    for (steps, null_kind), null in table[table['null'] != 'real'].groupby(['steps', 'null']):
        real = table[(table['graph'] == 'real') & (table['steps'] == steps)].iloc[0]
        for metric in metrics:
            values = null[metric].to_numpy(dtype=float)
            summary.append(dict(steps=steps, null=null_kind, metric=metric, real=float(real[metric]), rewired_mean=values.mean()
                                if len(values) else np.nan, rewired_sd=values.std(ddof=1) if len(values) > 1
                                else np.nan, rewired_min=values.min() if len(values) else np.nan,
                                rewired_max=values.max() if len(values) else np.nan, n_rewired=len(values),
                                rewired_at_least_real=int((values >= real[metric]).sum())))
    pd.DataFrame(summary).to_csv(args.prefix + '_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
