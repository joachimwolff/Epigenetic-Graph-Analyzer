"""Multi-way support of fragment pairs as a function of their distance in the graph (number of steps).

Pairs at graph distance 1 (direct contacts), 2, 3, 4 and 5 or more steps, and pairs in the same connected component
range are sampled within 2 Mb; each set is compared with control pairs at the same fragment offset drawn from the
same pools of graph loci (the matched controls of the pair validation, concordance of NPMI).

Usage: python hop_distance.py <graph.graphml> <multiway.npz> <output.tsv> [chain file] [pairs per step=50000]
"""

import sys

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.indirect import pair_class
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import LocusPools, map_nodes_to_bins, pairs_to_bins, validate_pairs

MAX_DISTANCE = 2000000
MAX_STEPS = 5


def main():
    graph_path, multiway_path, output = sys.argv[1:4]
    liftover = LiftOver(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != '-' else None
    n_pairs = int(sys.argv[5]) if len(sys.argv) > 5 else 50000
    graph = read_graph(graph_path)
    data = MultiwayData.load(multiway_path)
    node_bins = map_nodes_to_bins(graph, data, liftover)
    middle = {n: (int(d['start']) + int(d['end'])) // 2 for n, d in graph.nodes(data=True)}
    rows = []
    for source in graph.nodes:
        lengths = nx.single_source_shortest_path_length(graph, source, cutoff=MAX_STEPS)
        chrom = graph.nodes[source]['chrom']
        for target, steps in lengths.items():
            if target <= source or graph.nodes[target]['chrom'] != chrom:
                continue
            distance = abs(middle[target] - middle[source])
            if distance > MAX_DISTANCE:
                continue
            rows.append((source, target, chrom, distance, steps))
    pairs = pd.DataFrame(rows, columns=['node_a', 'node_b', 'chrom', 'distance', 'steps'])
    pairs['pair_class'] = [pair_class(bool(graph.nodes[a].get('bait')), bool(graph.nodes[b].get('bait')))
                           for a, b in zip(pairs['node_a'], pairs['node_b'])]
    pairs = pairs_to_bins(pairs, node_bins)
    pools = LocusPools(graph, node_bins, data)
    # every graph-connected pair within MAX_STEPS is excluded from the controls
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    max_offset = int(pairs['offset'].max())
    results = []
    for steps in range(1, MAX_STEPS + 1):
        subset = pairs[pairs['steps'] == steps]
        if len(subset) > n_pairs:
            subset = subset.sample(n=n_pairs, random_state=steps)
        subset = subset.reset_index(drop=True)
        _, summary = validate_pairs(subset, data, pools, excluded, 'steps_{}'.format(steps), max_offset=max_offset,
                                    seed=steps)
        row = summary[(summary['metric'] == 'npmi') & (summary['pair_class'] == 'all') &
                      (summary['stratum'] == 'all')].iloc[0]
        results.append(dict(steps=steps, pairs_available=int((pairs['steps'] == steps).sum()), n=int(row['n']),
                            median_distance=float(subset['distance'].median()),
                            concordance=float(row['fraction_above_control']),
                            mean_npmi=float(row['test_mean']), mean_npmi_control=float(row['control_mean'])))
        print(results[-1], flush=True)
    pd.DataFrame(results).to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    main()
