"""Rewired interaction inputs of a differential run: the cis interactions of both conditions are joined into one
union graph, rewired by double-edge swaps that keep the degree of every fragment, the bait / other-end kind of both
ends, the chromosome and the length bin of every edge (graph_nulls.rewire_by_distance), and written back as one
WashU text file per condition. Every rewired edge keeps the scores of both conditions of the edge its bait had before,
so the numbers of gained, lost and shared interactions, their scores and their distances stay as in the real data;
only which fragments are connected changes. Only pairs called in at least one condition (score >= --minScore) are
rewired; a score below the threshold in the other condition (the reference for the absence of a contact) moves with
the pair, and pairs below the threshold in both conditions are dropped. Calibration null for the explanations beyond step 1
(42_calibration.sh).

Usage: python rewire_inputs.py <washU A> <washU B> <baitmap> <out A> <out B> [--seed 0] [--swapsPerEdge 10]
           [--logBinWidth 0.1] [--maxLength 2000000] [--minScore 5]
"""

import argparse
import gzip
import logging
import os
import sys

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer.graph_nulls import rewire_by_distance  # noqa: E402


def read_pairs(path):
    raw = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['a', 'b', 'score'], dtype=str)
    raw['score'] = raw.score.astype(float)
    return raw[raw.a.str.split(',').str[0] == raw.b.str.split(',').str[0]]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('washu_a', 'washu_b', 'baitmap', 'out_a', 'out_b'):
        parser.add_argument(name)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--swapsPerEdge', type=float, default=10)
    parser.add_argument('--logBinWidth', type=float, default=0.1)
    parser.add_argument('--maxLength', type=int, default=2000000)
    parser.add_argument('--minScore', type=float, default=5.0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    baits = pd.read_csv(args.baitmap, sep='\t', header=None, usecols=[0, 1, 2], dtype=str)
    bait_keys = set('chr' + baits[0].str.replace('^chr', '', regex=True) + ',' + baits[1] + ',' + baits[2])
    graph = nx.Graph()
    for column, path in (('score_a', args.washu_a), ('score_b', args.washu_b)):
        pairs = read_pairs(path)
        for a, b, score in zip(pairs.a, pairs.b, pairs.score):
            for key in (a, b):
                if key not in graph:
                    chrom, start, end = key.split(',')
                    graph.add_node(key, chrom=chrom, start=int(start), end=int(end), bait=key in bait_keys)
            if graph.has_edge(a, b):
                graph.edges[a, b][column] = max(score, graph.edges[a, b].get(column, score))
            else:
                graph.add_edge(a, b, **{column: score})
    graph.remove_edges_from([(u, v) for u, v, d in graph.edges(data=True)
                             if max(d.get('score_a', 0), d.get('score_b', 0)) < args.minScore])
    graph.remove_nodes_from([n for n in list(graph.nodes) if graph.degree(n) == 0])
    logging.info('union graph: %d fragments, %d edges', graph.number_of_nodes(), graph.number_of_edges())
    rewired, statistics = rewire_by_distance(graph, swaps_per_edge=args.swapsPerEdge, log_bin_width=args.logBinWidth,
                                             seed=args.seed, max_length=args.maxLength, keep_attributes=True)
    for column, path in (('score_a', args.out_a), ('score_b', args.out_b)):
        with gzip.open(path, 'wt') as out:
            for u, v, data in rewired.edges(data=True):
                if column not in data:
                    continue
                # bait first; between two baits or two other ends the upstream fragment
                if rewired.nodes[v]['bait'] and not rewired.nodes[u]['bait'] or (
                        rewired.nodes[u]['bait'] == rewired.nodes[v]['bait'] and
                        rewired.nodes[v]['start'] < rewired.nodes[u]['start']):
                    u, v = v, u
                out.write('{}\t{}\t{}\n'.format(u, v, data[column]))
    print('\t'.join('{}={}'.format(k, v) for k, v in statistics.items()))


if __name__ == '__main__':
    main()
