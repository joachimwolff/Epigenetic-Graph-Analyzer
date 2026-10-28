"""Balance of the matched controls of the step-2 pairs: standardised mean differences and eligible controls.

The step-2 pairs of hop_distance.py (same sample) are matched with controls under the coverage-matched design and,
with --compartments / --domains, the domain-matched design, and the design matched on coverage and degree classes
(hop_distance.py --degreeClasses). For tested pairs and their controls, the standardised
mean difference (SMD, mean test minus mean control over the pooled standard deviation) is reported for the bin
offset, the genomic distance, and per end (lower and higher value of the two ends, and their mean) for the
multi-way coverage, the node degree, the fragment (bin) length, the local interaction density (graph edges with an
end within +-100 kb) and, with --dnase, the mean DNase-seq signal over the locus. The number of eligible control
pairs (distinct candidates at the same offset and class) per tested pair is summarised as well.

Usage: python balance.py <graph.graphml> <multiway.npz> <output prefix> [--chain FILE] [--compartments BED]
           [--domains BED] [--dnase bigWig] [--dnaseChain FILE] [--pairs 50000] [--degreeClasses 3]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hop_distance import bin_positions, step_pairs  # noqa: E402

from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402
from epigeneticgraphanalyzer.graph import read_graph  # noqa: E402
from epigeneticgraphanalyzer.graph_nulls import degree_bin_classes, local_edge_density, per_bin  # noqa: E402
from epigeneticgraphanalyzer.multiway import MultiwayData  # noqa: E402
from epigeneticgraphanalyzer.validation import (LocusPools, domain_annotation, map_nodes_to_bins,  # noqa: E402
                                                  standardized_mean_difference, validate_pairs)

STEPS = 2


def dnase_per_bin(path, bins, ids, liftover=None):
    """Mean signal of a bigWig over the bins `ids` (lifted start and end; NaN if not lifted to one chromosome)."""
    import pyBigWig
    handle = pyBigWig.open(path)
    sizes = handle.chroms()
    values = {}
    chrom_names = np.array(bins.chrom_names)[bins.chrom_index_of_bins(ids)]
    starts, ends = np.asarray(bins.bin_start(ids)), np.asarray(bins.bin_end(ids))
    for bin_id, chrom, start, end in zip(ids, chrom_names, starts, ends):
        if liftover is not None:
            a, b = liftover.position(chrom, int(start)), liftover.position(chrom, int(max(start, end - 1)))
            if a is None or b is None or a[0] != b[0]:
                values[bin_id] = np.nan
                continue
            chrom, start, end = a[0], min(a[1], b[1]), max(a[1], b[1]) + 1
        if chrom not in sizes or end <= start:
            values[bin_id] = np.nan
            continue
        value = handle.stats(chrom, int(start), int(min(end, sizes[chrom])), type='mean')[0]
        values[bin_id] = np.nan if value is None else float(value)
    return values


def end_features(features, bin_a, bin_b):
    """Per feature, the lower and higher value of the two ends and their mean."""
    out = {}
    for name, per_bin_values in features.items():
        a, b = per_bin_values(bin_a), per_bin_values(bin_b)
        out[name + '_low'] = np.fmin(a, b)
        out[name + '_high'] = np.fmax(a, b)
        out[name + '_mean'] = 0.5 * (a + b)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('multiway')
    parser.add_argument('prefix')
    parser.add_argument('--chain')
    parser.add_argument('--compartments')
    parser.add_argument('--domains')
    parser.add_argument('--dnase')
    parser.add_argument('--dnaseChain', help='Chain from the multi-way assembly to the DNase assembly')
    parser.add_argument('--pairs', type=int, default=50000)
    parser.add_argument('--coverageClasses', type=int, default=5)
    parser.add_argument('--degreeClasses', type=int, default=3,
                        help='Also evaluate the design matched on coverage and this many degree classes (0: off)')
    args = parser.parse_args()
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    node_bins = map_nodes_to_bins(graph, data, LiftOver(args.chain) if args.chain else None)
    pairs = step_pairs(graph, node_bins, data)
    pools = LocusPools(graph, node_bins, data)
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    max_offset = int(pairs['offset'].max())
    subset = pairs[pairs['steps'] == STEPS]
    if len(subset) > args.pairs:
        subset = subset.sample(n=args.pairs, random_state=STEPS)
    subset = subset.reset_index(drop=True)
    n = data.bins.n_bins
    degree = per_bin(node_bins, dict(graph.degree()), n)
    density = per_bin(node_bins, local_edge_density(graph, 100000), n)
    length = np.asarray(data.bins.bin_end(np.arange(n))) - np.asarray(data.bins.bin_start(np.arange(n)))
    features = {'coverage': lambda ids: data.coverage_raw[ids].astype(float), 'degree': lambda ids: degree[ids],
                'fragment_length': lambda ids: length[ids].astype(float), 'density_100kb': lambda ids: density[ids]}
    designs = [('coverage', dict(coverage_matched=args.coverageClasses))]
    if args.degreeClasses:
        designs.append(('coverage_degree{}'.format(args.degreeClasses),
                        dict(coverage_matched=args.coverageClasses,
                             bin_class=degree_bin_classes(graph, node_bins, n, args.degreeClasses))))
    if args.compartments or args.domains:
        compartment, tad_id = domain_annotation(data.bins, args.compartments, args.domains)
        designs.append(('domains', dict(coverage_matched=args.coverageClasses,
                                        compartment=compartment if args.compartments else None,
                                        tad_id=tad_id if args.domains else None)))
    matched = {}
    for name, options in designs:
        table, _, controls = validate_pairs(subset, data, pools, excluded, name, max_offset=max_offset, seed=STEPS,
                                            with_expected=False, return_controls=True, **options)
        matched[name] = (table, controls)
    if args.dnase:
        ids = np.unique(np.concatenate([np.concatenate([t['bin_a'], t['bin_b'], c['bin_a'], c['bin_b']])
                                        for t, c in matched.values()]))
        values = dnase_per_bin(args.dnase, data.bins, ids, LiftOver(args.dnaseChain) if args.dnaseChain else None)
        dnase = np.full(n, np.nan)
        dnase[list(values)] = list(values.values())
        features['dnase'] = lambda ids: dnase[ids]
    smd_rows, eligible_rows = [], []
    for name, (table, controls) in matched.items():
        test = end_features(features, table['bin_a'].to_numpy(), table['bin_b'].to_numpy())
        control = end_features(features, controls['bin_a'].to_numpy(), controls['bin_b'].to_numpy())
        test['offset'] = (table['bin_b'] - table['bin_a']).to_numpy(dtype=float)
        control['offset'] = (controls['bin_b'] - controls['bin_a']).to_numpy(dtype=float)
        test['distance'] = np.abs(bin_positions(data.bins, table['bin_b']) - bin_positions(data.bins, table['bin_a']))
        control['distance'] = np.abs(bin_positions(data.bins, controls['bin_b']) -
                                     bin_positions(data.bins, controls['bin_a']))
        for feature in test:
            smd_rows.append(dict(design=name, feature=feature, smd=standardized_mean_difference(test[feature],
                                                                                                 control[feature]),
                                 mean_test=float(np.nanmean(test[feature])),
                                 mean_control=float(np.nanmean(control[feature])),
                                 median_test=float(np.nanmedian(test[feature])),
                                 median_control=float(np.nanmedian(control[feature])),
                                 missing_test=int(np.isnan(test[feature]).sum()),
                                 missing_control=int(np.isnan(control[feature]).sum())))
        eligible = controls.attrs['eligible'].to_numpy()
        eligible_rows.append(dict(design=name, pairs=len(table), controls=len(controls),
                                  class_fallback=int(table.attrs.get('control_fallback', 0)),
                                  eligible_min=int(eligible.min()), eligible_median=float(np.median(eligible)),
                                  eligible_q10=float(np.quantile(eligible, 0.1)), eligible_mean=float(eligible.mean()),
                                  fraction_below_10=float((eligible < 10).mean()),
                                  fraction_none=float((eligible == 0).mean()),
                                  concordance=float(np.nanmean(table['npmi_concordance']))))
        print(eligible_rows[-1], flush=True)
    pd.DataFrame(smd_rows).to_csv(args.prefix + '_smd.tsv', sep='\t', index=False)
    pd.DataFrame(eligible_rows).to_csv(args.prefix + '_eligible.tsv', sep='\t', index=False)
    print(pd.DataFrame(smd_rows)[['design', 'feature', 'smd']].to_string(), flush=True)


if __name__ == '__main__':
    main()
