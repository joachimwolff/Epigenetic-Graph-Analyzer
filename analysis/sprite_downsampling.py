"""Support of graph pairs as a function of the depth of the multi-way data.

Random subsets of the SPRITE clusters (or GAM tubes) are drawn; for each depth the concordance of direct and
indirect pairs with matched controls and the fraction of individually confirmed pairs (confirm_pairs) are
computed on the same tested pairs. A rise without plateau indicates that the depth of the multi-way data, not
the graph, limits the confirmation of single contacts.

Usage: python sprite_downsampling.py <graph.graphml> <multiway.npz> <output.tsv> [pairs per set=100000]
"""

import sys

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.indirect import direct_contacts, indirect_contacts
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import (LocusPools, confirm_pairs, map_nodes_to_bins, pairs_to_bins,
                                                validate_pairs)

FRACTIONS = (0.125, 0.25, 0.5, 1.0)


def main():
    graph_path, multiway_path, output = sys.argv[1:4]
    n_pairs = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
    graph = read_graph(graph_path)
    full = MultiwayData.load(multiway_path)
    rng = np.random.default_rng(0)
    order = rng.permutation(full.n_groups)
    rows = []
    for fraction in FRACTIONS:
        keep = np.sort(order[:int(round(fraction * full.n_groups))])
        data = MultiwayData(full.bins, full.incidence[:, keep], full.group_weight[keep], full.kind)
        # the same pairs at every depth: map nodes on the full data, drop bins that lost all coverage
        node_bins = {n: b for n, b in map_nodes_to_bins(graph, full).items() if data.coverage_raw[b] > 0}
        direct = pairs_to_bins(direct_contacts(graph), node_bins)
        indirect = pairs_to_bins(indirect_contacts(graph), node_bins)
        excluded = set(zip(direct['bin_a'], direct['bin_b'])) | set(zip(indirect['bin_a'], indirect['bin_b']))
        direct = direct.sample(n=min(n_pairs, len(direct)), random_state=0).reset_index(drop=True)
        indirect = indirect[indirect['pair_class'] == 'oe-oe']
        indirect = indirect.sample(n=min(n_pairs, len(indirect)), random_state=1).reset_index(drop=True)
        pools = LocusPools(graph, node_bins, data)
        max_offset = int(max(direct['offset'].max(), indirect['offset'].max()))
        for label, pairs, seed in (('direct', direct, 0), ('indirect_oe-oe', indirect, 1)):
            _, summary = validate_pairs(pairs, data, pools, excluded, label, max_offset=max_offset, seed=seed)
            npmi = summary[(summary['metric'] == 'npmi') & (summary['pair_class'] == 'all')].iloc[0]
            _, confirmed = confirm_pairs(pairs, data, pools, excluded, label, max_offset, seed=seed + 10)
            overall = confirmed.iloc[0]
            far = confirmed[confirmed['stratum'].isin(['500-1000kb', '1000-2000kb'])]
            rows.append(dict(fraction=fraction, groups=len(keep), set=label, n=len(pairs),
                             median_locus_coverage=float(np.median(data.coverage_raw[pools.other_array])),
                             concordance=npmi['fraction_above_control'],
                             observed=overall['observed'], control_observed=overall['control_observed'],
                             confirmed_p=overall['confirmed_p'], control_confirmed_p=overall['control_confirmed_p'],
                             confirmed_fdr=overall['confirmed_fdr'],
                             confirmed_p_500kb_plus=float(np.average(far['confirmed_p'], weights=far['n']))
                             if far['n'].sum() else np.nan))
            print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    main()
