"""Ascertainment of the step curves: graph nodes and node pairs lost on the way to the multi-way data.

A node is lost if its midpoint cannot be lifted to the assembly of the multi-way data, falls outside its bins, or
lies in a bin never seen in any group (zero coverage); this is the node mapping of the pair validation
(validation.map_nodes_to_bins). A pair at graph distance 1..5 within 2 Mb is lost if an end is lost, if both ends
fall in the same bin, or if the mapped ends lie on different chromosomes or more than 2 Mb apart. The lost share is
reported by graph distance, pair class (bait-bait, bait-OE, OE-OE), genomic distance class and degree class of the
higher-degree end (1, 2, 3-4, 5-9, >= 10; quartiles of the degree collapse to two classes because most nodes have
degree 1-3), and per node by bait status and degree class.

Usage: python ascertainment.py <graph.graphml> <multiway.npz> <output prefix> [--chain FILE]
"""

import argparse

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.graph_nulls import graph_distance_pairs
from epigeneticgraphanalyzer.indirect import pair_class
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import distance_stratum, node_mapping_status

MAX_DISTANCE = 2000000
REASONS = ['liftover_failed', 'no_bin', 'zero_coverage', 'same_bin', 'beyond_limit_after_mapping']


def degree_classes(degree):
    """Degree classes 1, 2, 3-4, 5-9, >=10 (quartiles of the node degree collapse: most nodes have degree 1-3)."""
    edges = np.array([1, 2, 3, 5, 10])
    labels = np.array(['1', '2', '3-4', '5-9', '>=10'])
    return labels[np.searchsorted(edges, np.asarray(degree), side='right') - 1]


def bin_midpoints(bins, ids):
    ids = np.asarray(ids, dtype=np.int64)
    return (np.asarray(bins.bin_start(ids)) + np.asarray(bins.bin_end(ids))) // 2


def summarise(table, dimension):
    rows = []
    for level, group in table.groupby(dimension):
        row = dict(dimension=dimension, level=level, n=len(group), kept=int((group['reason'] == 'kept').sum()))
        row['dropped_fraction'] = 1 - row['kept'] / len(group)
        for reason in REASONS:
            row[reason] = int((group['reason'] == reason).sum())
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('prefix')
    parser.add_argument('--chain')
    args = parser.parse_args()
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    liftover = LiftOver(args.chain) if args.chain else None
    nodes = node_mapping_status(graph, data, liftover)
    nodes['degree_class'] = degree_classes(nodes['degree'])
    nodes['kind'] = np.where(nodes['bait'], 'bait', 'other end')
    nodes['reason'] = np.where(nodes['status'] == 'mapped', 'kept', nodes['status'])
    node_rows = summarise(nodes.assign(all='all'), 'all') + summarise(nodes, 'kind') + \
        summarise(nodes, 'degree_class')
    pd.DataFrame(node_rows).to_csv(args.prefix + '_nodes.tsv', sep='\t', index=False)
    print(pd.DataFrame(node_rows).to_string(), flush=True)

    status = dict(zip(nodes['node'], nodes['status']))
    bin_of = dict(zip(nodes['node'], nodes['bin']))
    degree_class = dict(zip(nodes['node'], nodes['degree_class']))
    degree = dict(zip(nodes['node'], nodes['degree']))
    pairs = graph_distance_pairs(graph, MAX_DISTANCE, 5)
    pairs['pair_class'] = [pair_class(bool(graph.nodes[a].get('bait')), bool(graph.nodes[b].get('bait')))
                           for a, b in zip(pairs['node_a'], pairs['node_b'])]
    pairs['distance_class'] = distance_stratum(pairs['distance'].to_numpy())
    pairs['degree_class'] = [degree_class[a] if degree[a] >= degree[b] else degree_class[b]
                             for a, b in zip(pairs['node_a'], pairs['node_b'])]
    reason = np.full(len(pairs), 'kept', dtype=object)
    for k, (a, b) in enumerate(zip(pairs['node_a'], pairs['node_b'])):
        for node in (a, b):
            if status[node] != 'mapped':
                reason[k] = status[node]
                break
    mapped = reason == 'kept'
    bin_a = pairs['node_a'].map(bin_of).to_numpy()
    bin_b = pairs['node_b'].map(bin_of).to_numpy()
    same = mapped & (bin_a == bin_b)
    reason[same] = 'same_bin'
    rest = np.nonzero(reason == 'kept')[0]
    chrom_a = data.bins.chrom_index_of_bins(bin_a[rest])
    chrom_b = data.bins.chrom_index_of_bins(bin_b[rest])
    far = (chrom_a != chrom_b) | (np.abs(bin_midpoints(data.bins, bin_a[rest]) -
                                         bin_midpoints(data.bins, bin_b[rest])) > MAX_DISTANCE)
    reason[rest[far]] = 'beyond_limit_after_mapping'
    pairs['reason'] = reason
    pair_rows = []
    for dimension in ('steps', 'pair_class', 'distance_class', 'degree_class'):
        pair_rows.extend(summarise(pairs, dimension))
    for steps, group in pairs.groupby('steps'):
        for dimension in ('pair_class', 'distance_class', 'degree_class'):
            for row in summarise(group, dimension):
                row['dimension'] = 'steps={} {}'.format(steps, dimension)
                pair_rows.append(row)
    table = pd.DataFrame(pair_rows)
    table.to_csv(args.prefix + '_pairs.tsv', sep='\t', index=False)
    print(table[table['dimension'] == 'steps'].to_string(), flush=True)


if __name__ == '__main__':
    main()
