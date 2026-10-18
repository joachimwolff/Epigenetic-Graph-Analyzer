"""Multi-way support of fragment pairs as a function of their distance in the graph (number of steps).

Pairs at graph distance 1 (direct contacts), 2, 3, 4 and 5 or more steps, and pairs in the same connected component
range are sampled within 2 Mb; each set is compared with control pairs at the same fragment offset drawn from the
same pools of graph loci and the same coverage classes of the multi-way data (the matched controls of the pair
validation, concordance of NPMI). With --compartments and --domains, controls also match the A/B compartments of both
ends and whether both ends lie in the same topological domain. The 95% interval of the concordance comes from a
block bootstrap over chromosomes, so that pairs sharing nodes and regions stay together.

Usage: python hop_distance.py <graph.graphml> <multiway.npz> <output.tsv> [--chain FILE] [--pairs 50000]
           [--coverageClasses 5] [--compartments BED] [--domains BED] [--bootstrap 1000] [--threads N]
"""

import argparse
import logging
import time

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
MAX_STEPS = 5


def bin_positions(bins, ids):
    """Midpoints of bins of a fragment index or of a fixed-size bin index."""
    ids = np.asarray(ids, dtype=np.int64)
    if hasattr(bins, 'bin_start'):
        return (bins.bin_start(ids) + bins.bin_end(ids)) // 2
    first = np.array([bins.offsets[c] for c in bins.chrom_names], dtype=np.int64)
    return (ids - first[bins.chrom_index_of_bins(ids)]) * bins.bin_size + bins.bin_size // 2


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('output')
    parser.add_argument('chain_file', nargs='?', default=None, help='(positional form of --chain, kept for old calls)')
    parser.add_argument('--chain', help='Chain from the graph assembly to the multi-way assembly')
    parser.add_argument('--pairs', type=int, default=50000, help='Pairs per graph distance')
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--compartments', help='A/B compartments (BED name A*/B* or eigenvector bedGraph), '
                                               'multi-way assembly')
    parser.add_argument('--domains', help='Topological domains (BED), multi-way assembly')
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--threads', type=int, default=0,
                        help='Worker processes for the pair queries (default: available CPUs, at most 24); the '
                             'results do not depend on the number')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
    clock = time.time()
    chain = args.chain or (args.chain_file if args.chain_file not in (None, '-') else None)
    liftover = LiftOver(chain) if chain else None
    n_pairs = args.pairs
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway).set_threads(args.threads)
    compartment, tad_id = None, None
    if args.compartments or args.domains:
        compartment, tad_id = domain_annotation(data.bins, args.compartments, args.domains)
        compartment = compartment if args.compartments else None
        tad_id = tad_id if args.domains else None
    node_bins = map_nodes_to_bins(graph, data, liftover)
    print('graph and multi-way data loaded, nodes mapped: {:.1f} s'.format(time.time() - clock), flush=True)
    clock = time.time()
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
    # after a liftover, a few pairs land on distant positions of the multi-way assembly; keep the distance limit there
    lifted = bin_positions(data.bins, pairs['bin_a']), bin_positions(data.bins, pairs['bin_b'])
    same_chrom = data.bins.chrom_index_of_bins(pairs['bin_a']) == data.bins.chrom_index_of_bins(pairs['bin_b'])
    keep = same_chrom & (np.abs(lifted[1] - lifted[0]) <= MAX_DISTANCE)
    if not keep.all():
        print('{} pairs beyond {} bp after mapping to the multi-way assembly removed'.format(int((~keep).sum()),
                                                                                          MAX_DISTANCE), flush=True)
        pairs = pairs[keep].reset_index(drop=True)
    pools = LocusPools(graph, node_bins, data)
    # every graph-connected pair within MAX_STEPS is excluded from the controls
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    max_offset = int(pairs['offset'].max())
    print('{} pairs within {} steps enumerated: {:.1f} s'.format(len(pairs), MAX_STEPS, time.time() - clock), flush=True)
    results = []
    for steps in range(1, MAX_STEPS + 1):
        subset = pairs[pairs['steps'] == steps]
        if len(subset) > n_pairs:
            subset = subset.sample(n=n_pairs, random_state=steps)
        subset = subset.reset_index(drop=True)
        table, summary = validate_pairs(subset, data, pools, excluded, 'steps_{}'.format(steps),
                                        max_offset=max_offset, seed=steps, coverage_matched=args.coverageClasses,
                                        compartment=compartment, tad_id=tad_id)
        clock = time.time()
        low, high = block_bootstrap(table['npmi_concordance'].to_numpy(), table['chrom'].to_numpy(),
                                    n_boot=args.bootstrap, seed=steps)
        print('step {}: block bootstrap {:.1f} s'.format(steps, time.time() - clock), flush=True)
        row = summary[(summary['metric'] == 'npmi') & (summary['pair_class'] == 'all') &
                      (summary['stratum'] == 'all')].iloc[0]
        results.append(dict(steps=steps, pairs_available=int((pairs['steps'] == steps).sum()), n=int(row['n']),
                            median_distance=float(subset['distance'].median()),
                            concordance=float(row['fraction_above_control']), ci_low=low, ci_high=high,
                            mean_npmi=float(row['test_mean']), mean_npmi_control=float(row['control_mean'])))
        print(results[-1], flush=True)
    pd.DataFrame(results).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
