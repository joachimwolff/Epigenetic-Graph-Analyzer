"""Graph node pairs whose ends fall into the same or adjacent windows of a binned multi-way data set (GAM).

For node pairs at graph distance 1..5 within 2 Mb, the ends are mapped as in the pair validation (midpoint, liftover,
covered bins) and the pairs are counted by the window offset of their ends: same window (offset 0, never tested),
adjacent windows (offset 1) and further apart, per graph distance and pair class.

Usage: python gam_windows.py <graph.graphml> <multiway.npz> <output.tsv> [--chain FILE]
"""

import argparse

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.graph_nulls import graph_distance_pairs
from epigeneticgraphanalyzer.indirect import pair_class
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import map_nodes_to_bins


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('output')
    parser.add_argument('--chain')
    args = parser.parse_args()
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    node_bins = map_nodes_to_bins(graph, data, LiftOver(args.chain) if args.chain else None)
    pairs = graph_distance_pairs(graph, 2000000, 5)
    pairs['pair_class'] = [pair_class(bool(graph.nodes[a].get('bait')), bool(graph.nodes[b].get('bait')))
                           for a, b in zip(pairs['node_a'], pairs['node_b'])]
    a = pairs['node_a'].map(node_bins)
    b = pairs['node_b'].map(node_bins)
    mapped = (a.notna() & b.notna()).to_numpy()
    offset = np.full(len(pairs), -1, dtype=np.int64)
    a_ids, b_ids = a[mapped].astype(np.int64).to_numpy(), b[mapped].astype(np.int64).to_numpy()
    same_chrom = data.bins.chrom_index_of_bins(a_ids) == data.bins.chrom_index_of_bins(b_ids)
    offset[np.nonzero(mapped)[0][same_chrom]] = np.abs(b_ids - a_ids)[same_chrom]
    pairs['window_offset'] = np.where(offset < 0, 'not mapped', np.where(offset == 0, '0', np.where(
        offset == 1, '1', '>=2')))
    rows = []
    for keys, group in pairs.groupby(['steps', 'pair_class']):
        rows.append(dict(steps=keys[0], pair_class=keys[1], **group['window_offset'].value_counts().to_dict()))
    for steps, group in pairs.groupby('steps'):
        rows.append(dict(steps=steps, pair_class='all', **group['window_offset'].value_counts().to_dict()))
    table = pd.DataFrame(rows).fillna(0)
    for column in ('0', '1', '>=2', 'not mapped'):
        if column not in table:
            table[column] = 0
    mapped_total = table['0'] + table['1'] + table['>=2']
    table['fraction_offset0'] = table['0'] / mapped_total
    table['fraction_offset1'] = table['1'] / mapped_total
    table['fraction_offset_le1'] = (table['0'] + table['1']) / mapped_total
    table.to_csv(args.output, sep='\t', index=False)
    print(table.to_string(), flush=True)


if __name__ == '__main__':
    main()
