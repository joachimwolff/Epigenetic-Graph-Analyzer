"""Inputs for the ChAseR comparison: fragment pairs by graph distance and chromatin features per fragment.

For graph distance 1 to 5 (shortest path length, cis, within 2 Mb), up to 50,000 pairs are sampled as in
analysis/hop_distance.py and written as one ChAseR network (barebone format: chrom start end chrom start end) per
distance. The feature table holds the mean signal of every bigWig track over each fragment, log2(x + 1), with the
node names ChAseR expects (chrom:start-end).

Usage: python chaser_inputs.py <graph.graphml> <output dir> mark=track.bw [mark=track.bw ...]
"""

import os
import sys

import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig

from epigeneticgraphanalyzer.graph import read_graph

MAX_DISTANCE = 2000000
MAX_STEPS = 5
N_PAIRS = 50000


def main():
    graph_path, out = sys.argv[1:3]
    tracks = dict(arg.split('=', 1) for arg in sys.argv[3:])
    os.makedirs(out, exist_ok=True)
    graph = read_graph(graph_path)
    data = {n: (d['chrom'], int(d['start']), int(d['end'])) for n, d in graph.nodes(data=True)}
    middle = {n: (s + e) // 2 for n, (c, s, e) in data.items()}
    pairs = {k: [] for k in range(1, MAX_STEPS + 1)}
    for source in graph.nodes:
        for target, steps in nx.single_source_shortest_path_length(graph, source, cutoff=MAX_STEPS).items():
            if steps == 0 or target <= source or data[target][0] != data[source][0]:
                continue
            if abs(middle[target] - middle[source]) <= MAX_DISTANCE:
                pairs[steps].append((source, target))
    rng = np.random.default_rng(0)
    for steps, items in pairs.items():
        if len(items) > N_PAIRS:
            items = [items[i] for i in rng.choice(len(items), N_PAIRS, replace=False)]
        rows = [data[a] + data[b] for a, b in items]
        pd.DataFrame(rows).to_csv(os.path.join(out, 'network_step{}.txt'.format(steps)), sep='\t', header=False,
                                  index=False)
        print(steps, len(items))
    handles = {mark: pyBigWig.open(path) for mark, path in tracks.items()}
    features = {}
    for n, (chrom, start, end) in data.items():
        row = []
        for mark, handle in handles.items():
            if chrom not in handle.chroms():
                row.append(np.nan)
                continue
            value = handle.stats(chrom, start, min(end, handle.chroms()[chrom]), type='mean')[0]
            row.append(np.log2((value or 0.0) + 1.0))
        features['{}:{}-{}'.format(chrom, start, end)] = row
    table = pd.DataFrame.from_dict(features, orient='index', columns=list(handles))
    table.to_csv(os.path.join(out, 'features.txt'), sep='\t', index_label=False)


if __name__ == '__main__':
    main()
