"""Multi-way support of fragment pairs as a function of their distance in the graph (number of steps).

Pairs at graph distance 1 (direct contacts), 2, 3, 4 and 5 or more steps, and pairs in the same connected component
range are sampled within 2 Mb (--maxDistance; 0 keeps all pairs of one chromosome); each set is compared with control
pairs at the same fragment offset drawn from the same pools of graph loci and the same coverage classes of the
multi-way data (the matched controls of the pair validation, concordance of NPMI). With --compartments and
--domains, controls also match the A/B compartments of both ends and whether both ends lie in the same topological
domain. With --degreeClasses, controls also match the graph degree class of both ends (quantile classes among bait
and among other-end loci, graph_nulls.degree_bin_classes), so that pairs of high-degree hubs are compared with
pairs of equally connected loci. The 95% interval of the concordance comes from a block bootstrap over chromosomes,
so that pairs sharing nodes and regions stay together. control_fallback counts the pairs without a class-matched
control at their offset (matched on the offset only); then the concordance over the pairs with class-matched
controls only is reported as well (*_class_matched). --minCoverageQuantile keeps only graph loci with a multi-way
coverage at or above that quantile of all graph loci, for tests and controls. --minOffset drops pairs whose ends
are fewer bins apart (pairs in one bin are always dropped). With --rates, the absolute co-occurrence of tested and
control pairs is reported: the share found together in at least one group, the mean share of the groups holding
one end that also hold the other, and the mean number of groups holding both ends (unweighted).

Usage: python hop_distance.py <graph.graphml> <multiway.npz> <output.tsv> [--chain FILE] [--pairs 50000]
           [--coverageClasses 5] [--degreeClasses 0] [--compartments BED] [--domains BED] [--bootstrap 1000]
           [--maxDistance 2000000] [--minCoverageQuantile 0] [--minOffset 1] [--rates] [--pairTable FILE] [--threads N]
"""

import argparse
import logging

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.graph_nulls import degree_bin_classes, graph_distance_pairs
from epigeneticgraphanalyzer.indirect import pair_class
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import (LocusPools, block_bootstrap, cooccurrence_rates, domain_annotation,
                                                  map_nodes_to_bins, pairs_to_bins, validate_pairs)

MAX_DISTANCE = 2000000
MAX_STEPS = 5


def bin_positions(bins, ids):
    """Midpoints of bins of a fragment index or of a fixed-size bin index."""
    ids = np.asarray(ids, dtype=np.int64)
    if hasattr(bins, 'bin_start'):
        return (bins.bin_start(ids) + bins.bin_end(ids)) // 2
    first = np.array([bins.offsets[c] for c in bins.chrom_names], dtype=np.int64)
    return (ids - first[bins.chrom_index_of_bins(ids)]) * bins.bin_size + bins.bin_size // 2


def add_rates(rows, tables, data):
    """Mean absolute co-occurrence of the tested pairs and their controls added to the step rows."""
    for row in rows:
        table, controls = tables[row['steps']]
        test = cooccurrence_rates(data, table['bin_a'].to_numpy(), table['bin_b'].to_numpy())
        control = cooccurrence_rates(data, controls['bin_a'].to_numpy(), controls['bin_b'].to_numpy())
        for key in ('together', 'share', 'groups_both'):
            row[key] = float(np.mean(test[key]))
            row[key + '_control'] = float(np.mean(control[key]))
    return rows


def coverage_filtered(node_bins, data, quantile):
    """Node bins without the loci whose multi-way coverage is below the given quantile of the coverage of all graph
    loci (test pairs and control pools then use the better covered loci only)."""
    loci = np.unique(np.fromiter(node_bins.values(), dtype=np.int64))
    threshold = np.quantile(data.coverage_raw[loci], quantile)
    kept = {n: b for n, b in node_bins.items() if data.coverage_raw[b] >= threshold}
    print('coverage threshold {:g} (quantile {:g} of graph loci): {} of {} nodes kept'.format(
        threshold, quantile, len(kept), len(node_bins)), flush=True)
    return kept


def step_pairs(graph, node_bins, data, max_distance=MAX_DISTANCE, max_steps=MAX_STEPS):
    """Pairs at graph distance 1..max_steps within max_distance (None: all cis), mapped to multi-way bins."""
    pairs = graph_distance_pairs(graph, max_distance, max_steps)
    pairs['pair_class'] = [pair_class(bool(graph.nodes[a].get('bait')), bool(graph.nodes[b].get('bait')))
                           for a, b in zip(pairs['node_a'], pairs['node_b'])]
    pairs = pairs_to_bins(pairs, node_bins)
    # after a liftover, a few pairs land on distant positions of the multi-way assembly; keep the distance limit there
    same_chrom = data.bins.chrom_index_of_bins(pairs['bin_a']) == data.bins.chrom_index_of_bins(pairs['bin_b'])
    keep = same_chrom
    if max_distance is not None:
        lifted = bin_positions(data.bins, pairs['bin_a']), bin_positions(data.bins, pairs['bin_b'])
        keep = same_chrom & (np.abs(lifted[1] - lifted[0]) <= max_distance)
    if not keep.all():
        print('{} pairs beyond {} bp or on other chromosomes after mapping to the multi-way assembly removed'.format(
            int((~keep).sum()), max_distance), flush=True)
        pairs = pairs[keep].reset_index(drop=True)
    return pairs


def evaluate_steps(pairs, data, pools, excluded, n_pairs, coverage_classes=5, compartment=None, tad_id=None,
                   bin_class=None, n_boot=1000, steps_list=range(1, MAX_STEPS + 1), return_tables=False):
    """Concordance of NPMI against matched controls per graph distance (one row per number of steps)."""
    max_offset = int(pairs['offset'].max())
    results, tables = [], {}
    for steps in steps_list:
        subset = pairs[pairs['steps'] == steps]
        if len(subset) > n_pairs:
            subset = subset.sample(n=n_pairs, random_state=steps)
        subset = subset.reset_index(drop=True)
        table, summary, controls = validate_pairs(subset, data, pools, excluded, 'steps_{}'.format(steps),
                                                  max_offset=max_offset, seed=steps,
                                                  coverage_matched=coverage_classes, compartment=compartment,
                                                  tad_id=tad_id, bin_class=bin_class, with_expected=False,
                                                  return_controls=True)
        low, high = block_bootstrap(table['npmi_concordance'].to_numpy(), table['chrom'].to_numpy(),
                                    n_boot=n_boot, seed=steps)
        row = summary[(summary['metric'] == 'npmi') & (summary['pair_class'] == 'all') &
                      (summary['stratum'] == 'all')].iloc[0]
        results.append(dict(steps=steps, pairs_available=int((pairs['steps'] == steps).sum()), n=int(row['n']),
                            median_distance=float(subset['distance'].median()),
                            concordance=float(row['fraction_above_control']), ci_low=low, ci_high=high,
                            mean_npmi=float(row['test_mean']), mean_npmi_control=float(row['control_mean']),
                            control_fallback=int(table.attrs.get('control_fallback', 0))))
        if results[-1]['control_fallback']:
            # the same on the pairs with class-matched controls only
            matched = table['class_matched'].to_numpy()
            values = table['npmi_concordance'].to_numpy()[matched]
            low_m, high_m = block_bootstrap(values, table['chrom'].to_numpy()[matched], n_boot=n_boot, seed=steps)
            results[-1].update(n_class_matched=int(np.sum(~np.isnan(values))),
                               concordance_class_matched=float(np.nanmean(values)), ci_low_class_matched=low_m,
                               ci_high_class_matched=high_m,
                               median_distance_class_matched=float(np.median(table['distance'].to_numpy()[matched])))
        print(results[-1], flush=True)
        if return_tables:
            tables[steps] = (table, controls)
    return (results, tables) if return_tables else results


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('output')
    parser.add_argument('chain_file', nargs='?', default=None, help='(positional form of --chain, kept for old calls)')
    parser.add_argument('--chain', help='Chain from the graph assembly to the multi-way assembly')
    parser.add_argument('--pairs', type=int, default=50000, help='Pairs per graph distance')
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--degreeClasses', type=int, default=0,
                        help='Also match controls on this many graph degree classes of both ends (0: off)')
    parser.add_argument('--compartments', help='A/B compartments (BED name A*/B* or eigenvector bedGraph), '
                                               'multi-way assembly')
    parser.add_argument('--domains', help='Topological domains (BED), multi-way assembly')
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--maxDistance', type=int, default=MAX_DISTANCE,
                        help='Largest genomic distance of a pair (bp); 0 keeps all pairs of one chromosome')
    parser.add_argument('--minCoverageQuantile', type=float, default=0,
                        help='Use only graph loci whose multi-way coverage reaches this quantile of the coverage of '
                             'all graph loci (e.g. 0.2: top four quintiles)')
    parser.add_argument('--minOffset', type=int, default=1,
                        help='Test only pairs whose ends are at least this many bins apart (controls are matched on '
                             'the offset, so they follow); e.g. 2 drops pairs in adjacent GAM windows')
    parser.add_argument('--rates', action='store_true',
                        help='Also report the absolute co-occurrence of test and control pairs')
    parser.add_argument('--pairTable', help='Write the tested pairs with their NPMI, control mean and concordance '
                                            '(TSV, .gz allowed)')
    parser.add_argument('--threads', type=int, default=0,
                        help='Worker processes for the pair queries (default: available CPUs, at most 24); the '
                             'results do not depend on the number')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
    chain = args.chain or (args.chain_file if args.chain_file not in (None, '-') else None)
    liftover = LiftOver(chain) if chain else None
    max_distance = args.maxDistance if args.maxDistance > 0 else None
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway).set_threads(args.threads)
    compartment, tad_id = None, None
    if args.compartments or args.domains:
        compartment, tad_id = domain_annotation(data.bins, args.compartments, args.domains)
        compartment = compartment if args.compartments else None
        tad_id = tad_id if args.domains else None
    node_bins = map_nodes_to_bins(graph, data, liftover)
    if args.minCoverageQuantile > 0:
        node_bins = coverage_filtered(node_bins, data, args.minCoverageQuantile)
    pairs = step_pairs(graph, node_bins, data, max_distance)
    pools = LocusPools(graph, node_bins, data)
    # every graph-connected pair within MAX_STEPS is excluded from the controls
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    if args.minOffset > 1:
        close = (pairs['offset'] < args.minOffset).to_numpy()
        print('{} pairs with ends fewer than {} bins apart removed: {}'.format(
            int(close.sum()), args.minOffset, pairs[close].groupby('steps').size().to_dict()), flush=True)
        pairs = pairs[~close].reset_index(drop=True)
    bin_class = None
    if args.degreeClasses:
        bin_class = degree_bin_classes(graph, node_bins, data.bins.n_bins, args.degreeClasses)
        print('degree classes of graph loci: {}'.format(dict(zip(*np.unique(bin_class[bin_class >= 0],
                                                                            return_counts=True)))), flush=True)
    results, tables = evaluate_steps(pairs, data, pools, excluded, args.pairs, args.coverageClasses, compartment,
                                     tad_id, bin_class, args.bootstrap, return_tables=True)
    if args.rates:
        results = add_rates(results, tables, data)
    pd.DataFrame(results).to_csv(args.output, sep='\t', index=False)
    if args.pairTable:
        columns = ['steps', 'node_a', 'node_b', 'chrom', 'distance', 'pair_class', 'bin_a', 'bin_b', 'offset', 'npmi',
                   'npmi_control', 'npmi_concordance', 'class_matched']
        pd.concat([table[columns] for table, _ in tables.values()]).to_csv(args.pairTable, sep='\t', index=False)


if __name__ == '__main__':
    main()
