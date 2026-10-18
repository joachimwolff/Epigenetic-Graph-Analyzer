"""Cell-type specificity of the multi-way support: interactions found only in one cell type, in both, tested in the
multi-way data of each cell type.

Two graphs on the same capture design (e.g. GM12878 and H9 hESC on the hg19 HindIII fragments) are split into
interactions present only in the first graph, only in the second, and in both. Each set is compared with matched
controls (same offset, same locus pools, same coverage classes and, with --compartments/--domains, the same A/B
compartments and domain membership) in each multi-way data set. If the support reflects cell-type-specific
proximity, interactions specific to a cell type are supported more in the multi-way data of that cell type than in
the other (a double dissociation); shared architecture supports all sets alike.

Usage: python celltype_specific_edges.py <graph1.graphml> <graph2.graphml> <name1> <name2> <output.tsv>
           --multiway NAME FILE.npz [CHAIN] [--multiway ...] [--domainsFor NAME COMPARTMENTS DOMAINS] [--pairs 50000]
           [--threads N]
"""

import argparse

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.indirect import pair_class
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import (LocusPools, block_bootstrap, domain_annotation, map_nodes_to_bins,
                                                pairs_to_bins, validate_pairs)

MAX_DISTANCE = 2000000


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph1')
    parser.add_argument('graph2')
    parser.add_argument('name1')
    parser.add_argument('name2')
    parser.add_argument('output')
    parser.add_argument('--multiway', nargs='+', action='append', required=True,
                        help='NAME FILE.npz [CHAIN from the graph assembly to the multi-way assembly]')
    parser.add_argument('--domainsFor', nargs=3, action='append', default=[],
                        metavar=('NAME', 'COMPARTMENTS', 'DOMAINS'))
    parser.add_argument('--pairs', type=int, default=50000)
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--threads', type=int, default=0,
                        help='Worker processes for the pair queries (default: available CPUs, at most 24)')
    args = parser.parse_args()
    graphs = [read_graph(args.graph1), read_graph(args.graph2)]
    keys = [{tuple(sorted((u, v))) for u, v in g.edges} for g in graphs]
    union = nx.compose(*graphs)
    middle = {n: (int(d['start']) + int(d['end'])) // 2 for n, d in union.nodes(data=True)}
    sets = {'{} only'.format(args.name1): keys[0] - keys[1], '{} only'.format(args.name2): keys[1] - keys[0],
            'both': keys[0] & keys[1]}
    domains = {name: (c, d) for name, c, d in args.domainsFor}
    rows = []
    for spec in args.multiway:
        name, path = spec[0], spec[1]
        chain = spec[2] if len(spec) > 2 else None
        data = MultiwayData.load(path).set_threads(args.threads)
        node_bins = map_nodes_to_bins(union, data, LiftOver(chain) if chain else None)
        pools = LocusPools(union, node_bins, data)
        compartment, tad_id = None, None
        if name in domains:
            compartment, tad_id = domain_annotation(data.bins, *domains[name])
        all_pairs = []
        for label, edges in sets.items():
            frame = pd.DataFrame([(u, v, union.nodes[u]['chrom'], abs(middle[u] - middle[v])) for u, v in edges
                                  if union.nodes[u]['chrom'] == union.nodes[v]['chrom']],
                                 columns=['node_a', 'node_b', 'chrom', 'distance'])
            frame = frame[frame['distance'] <= MAX_DISTANCE]
            frame['pair_class'] = [pair_class(bool(union.nodes[a].get('bait')), bool(union.nodes[b].get('bait')))
                                   for a, b in zip(frame['node_a'], frame['node_b'])]
            frame['set'] = label
            all_pairs.append(frame)
        pairs = pairs_to_bins(pd.concat(all_pairs, ignore_index=True), node_bins)
        excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
        max_offset = int(pairs['offset'].max())
        for k, (label, sub) in enumerate(pairs.groupby('set', sort=False)):
            if len(sub) > args.pairs:
                sub = sub.sample(n=args.pairs, random_state=k)
            sub = sub.reset_index(drop=True)
            table, summary = validate_pairs(sub, data, pools, excluded, label, max_offset=max_offset, seed=k,
                                            coverage_matched=args.coverageClasses, compartment=compartment,
                                            tad_id=tad_id)
            row = summary[(summary['metric'] == 'npmi') & (summary['pair_class'] == 'all') &
                          (summary['stratum'] == 'all')].iloc[0]
            low, high = block_bootstrap(table['npmi_concordance'].to_numpy(), table['chrom'].to_numpy(),
                                        n_boot=args.bootstrap, seed=k)
            rows.append(dict(multiway=name, interactions=label, n=int(row['n']),
                             median_distance=float(sub['distance'].median()),
                             concordance=float(row['fraction_above_control']), ci_low=low, ci_high=high,
                             domain_matched=name in domains))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
