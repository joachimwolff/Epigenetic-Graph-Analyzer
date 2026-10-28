"""Comparability of the input graphs: size, bait share, edge lengths and degrees of every graph used in the paper.

The manifest (TSV with columns name, path, assembly, calling, replicate_rule, source) lists the graphs. For a union
graph of the differential mode (edges with r_a / r_b), the two conditions are reported as well, an edge counting
for a condition when its strength there is positive. Per graph: nodes, baits, edges, cis edges, cis edges <= 2 Mb,
median and interquartile range of the cis edge length (node midpoints), mean degree of baits and of other ends,
largest connected component, and the number of node pairs at graph distance 2 within 2 Mb.

Usage: python graph_comparability.py <manifest.tsv> <output.tsv>
"""

import argparse

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.graph_nulls import graph_distance_pairs, graph_statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('manifest')
    parser.add_argument('output')
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep='\t', comment='#')
    rows = []
    for entry in manifest.itertuples(index=False):
        graph = read_graph(entry.path)
        variants = [('all', graph)]
        if graph.number_of_edges() and 'r_a' in next(iter(graph.edges(data=True)))[2]:
            names = str(entry.name).split('_vs_')
            for key, label in (('r_a', names[0]), ('r_b', names[-1])):
                keep = [(u, v) for u, v, d in graph.edges(data=True)
                        if np.isfinite(float(d.get(key, np.nan))) and float(d.get(key, 0)) > 0]
                sub = nx.Graph()
                sub.add_nodes_from((n, graph.nodes[n]) for e in keep for n in e)
                sub.add_edges_from(keep)
                variants.append((label, sub))
        for condition, g in variants:
            row = dict(name=entry.name, condition=condition, assembly=entry.assembly, calling=entry.calling,
                       replicate_rule=entry.replicate_rule, source=entry.source, path=entry.path)
            row.update(graph_statistics(g))
            row['pairs_step2_le_2mb'] = int((graph_distance_pairs(g, 2000000, 2)['steps'] == 2).sum())
            print(row, flush=True)
            rows.append(row)
        pd.DataFrame(rows).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
