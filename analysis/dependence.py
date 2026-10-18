"""Dependence between tested pairs of the step curves: thinned pair sets, hub exclusion, leave-one-chromosome-out
and a genomic block bootstrap.

Pairs of a step curve share nodes (a promoter with many partners) and, at two steps, mediators, so they are not
independent. On the same graph / multi-way pairing as hop_distance.py (coverage-matched controls, NPMI
concordance):
  (a) step 2 thinned to one pair per promoter (no bait end in two pairs), one pair per mediator (no mediator shared
      by two pairs) and both at once, greedy over the pairs in a fixed random order;
  (b) steps 1-5 without pairs that touch the 1% or 5% nodes of highest degree;
  (c) leave-one-chromosome-out range of the concordance of the full step curve;
  (d) a block bootstrap over 10 Mb genomic blocks next to the chromosome bootstrap.

Usage: python dependence.py <graph.graphml> <multiway.npz> <output.tsv> [--chain FILE] [--pairs 50000]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hop_distance import evaluate_steps, step_pairs  # noqa: E402

from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402
from epigeneticgraphanalyzer.graph import read_graph  # noqa: E402
from epigeneticgraphanalyzer.graph_nulls import node_midpoints, thin_pairs, top_degree_nodes  # noqa: E402
from epigeneticgraphanalyzer.multiway import MultiwayData  # noqa: E402
from epigeneticgraphanalyzer.validation import (LocusPools, block_bootstrap, map_nodes_to_bins,  # noqa: E402
                                                  validate_pairs)

BLOCK = 10000000


def genomic_blocks(table, middle, size=BLOCK):
    """Block label of every pair: chromosome and the size-bp window holding the midpoint between both ends."""
    centre = (table['node_a'].map(middle).to_numpy() + table['node_b'].map(middle).to_numpy()) // 2
    return np.array(['{}:{}'.format(c, m // size) for c, m in zip(table['chrom'], centre)])


def describe(table, label, steps, middle, n_boot, extra=None):
    values = table['npmi_concordance'].to_numpy()
    low, high = block_bootstrap(values, table['chrom'].to_numpy(), n_boot=n_boot, seed=steps)
    low10, high10 = block_bootstrap(values, genomic_blocks(table, middle), n_boot=n_boot, seed=steps)
    row = dict(analysis=label, steps=steps, n=int(np.sum(~np.isnan(values))), concordance=float(np.nanmean(values)),
               ci_low_chromosome=low, ci_high_chromosome=high, ci_low_10mb=low10, ci_high_10mb=high10,
               median_distance=float(table['distance'].median()), mean_npmi=float(table['npmi'].mean()),
               mean_npmi_control=float(table['npmi_control'].mean()))
    row.update(extra or {})
    print(row, flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('output')
    parser.add_argument('--chain')
    parser.add_argument('--pairs', type=int, default=50000)
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    node_bins = map_nodes_to_bins(graph, data, LiftOver(args.chain) if args.chain else None)
    pairs = step_pairs(graph, node_bins, data)
    pools = LocusPools(graph, node_bins, data)
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    max_offset = int(pairs['offset'].max())
    middle = node_midpoints(graph)
    rows = []
    # full step curve: chromosome and 10 Mb block bootstrap, leave-one-chromosome-out
    _, tables = evaluate_steps(pairs, data, pools, excluded, args.pairs, args.coverageClasses, n_boot=args.bootstrap,
                               return_tables=True)
    for steps, (table, _) in tables.items():
        values = table['npmi_concordance'].to_numpy()
        chroms = table['chrom'].to_numpy()
        loco = {c: float(np.nanmean(values[chroms != c])) for c in np.unique(chroms)}
        low_c, high_c = min(loco, key=loco.get), max(loco, key=loco.get)
        rows.append(describe(table, 'full', steps, middle, args.bootstrap,
                             dict(loco_min=loco[low_c], loco_min_chrom=low_c, loco_max=loco[high_c],
                                  loco_max_chrom=high_c, n_blocks_10mb=len(np.unique(genomic_blocks(table, middle))),
                                  n_chromosomes=len(loco))))
    # (b) without pairs touching the highest-degree nodes
    for share in (0.01, 0.05):
        hubs = top_degree_nodes(graph, share)
        touching = pairs['node_a'].isin(hubs).to_numpy() | pairs['node_b'].isin(hubs).to_numpy()
        _, tables = evaluate_steps(pairs[~touching].reset_index(drop=True), data, pools, excluded, args.pairs,
                                   args.coverageClasses, n_boot=args.bootstrap, return_tables=True)
        for steps, (table, _) in tables.items():
            rows.append(describe(table, 'without top {:g}% degree nodes'.format(100 * share), steps, middle,
                                 args.bootstrap, dict(excluded_nodes=len(hubs),
                                                      min_excluded_degree=min(graph.degree(n) for n in hubs),
                                                      pairs_available=int(((pairs['steps'] == steps) &
                                                                           ~touching).sum()))))
    # (a) step 2 thinned by promoter and by mediator
    two = pairs[pairs['steps'] == 2].reset_index(drop=True)
    for label, by_bait, by_mediator in (('one pair per promoter', True, False),
                                        ('one pair per mediator', False, True),
                                        ('one pair per promoter and mediator', True, True)):
        keep = thin_pairs(graph, two['node_a'].to_numpy(), two['node_b'].to_numpy(), by_bait, by_mediator,
                          seed=args.seed)
        subset = two[keep].reset_index(drop=True)
        n_available = len(subset)
        if len(subset) > args.pairs:
            subset = subset.sample(n=args.pairs, random_state=2).reset_index(drop=True)
        table, _ = validate_pairs(subset, data, pools, excluded, label, max_offset=max_offset, seed=2,
                                  coverage_matched=args.coverageClasses, with_expected=False)
        rows.append(describe(table, label, 2, middle, args.bootstrap, dict(pairs_available=n_available)))
    pd.DataFrame(rows).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
