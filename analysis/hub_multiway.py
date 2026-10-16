"""How many elements of each promoter hub occur together in single SPRITE clusters or GAM tubes.

A hub is a baited promoter fragment with all its direct partners in the graph. For every hub the script counts,
per multi-way group (SPRITE cluster, GAM tube), how many hub elements the group contains, and reports the largest
number in one group and the number of groups holding the promoter together with at least two partners (each such
group contains an indirect partner-partner contact together with its mediating promoter). As reference, the
whole hub is shifted along the chromosome by a random number of fragments (same shape and distances), ten times.

Usage: python hub_multiway.py <graph.graphml> <multiway.npz> <output prefix> [chain file for liftover]
"""

import sys

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import map_nodes_to_bins

N_SHIFTS = 10
SIZE_CLASSES = ((3, 3), (4, 5), (6, 10), (11, 1000))


def hub_counts(data, promoter_bin, partner_bins):
    """(largest number of hub elements in one group, groups with promoter + >= 2 partners).

    Counts only the groups that hold an element of the hub, so the cost does not grow with the number of groups.
    """
    promoter_groups = data.incidence[promoter_bin].indices
    groups, partners_in = np.unique(data.incidence[partner_bins].indices, return_counts=True)
    with_promoter = np.isin(groups, promoter_groups, assume_unique=True)
    largest = int((partners_in + with_promoter).max()) if len(groups) else 0
    largest = max(largest, int(len(promoter_groups) > 0))
    return largest, int(np.sum(with_promoter & (partners_in >= 2)))


def hub_counts_dense(data, promoter_bin, partner_bins):
    """Reference implementation of hub_counts over all groups (kept for the test)."""
    rows = data.incidence[np.concatenate([[promoter_bin], partner_bins])]
    per_group = np.asarray(rows.sum(axis=0)).ravel()
    with_promoter = np.zeros(data.n_groups, dtype=bool)
    with_promoter[data.incidence[promoter_bin].indices] = True
    partners_in = per_group - with_promoter
    return int(per_group.max()), int(np.sum(with_promoter & (partners_in >= 2)))


def main():
    graph_path, multiway_path, prefix = sys.argv[1:4]
    liftover = LiftOver(sys.argv[4]) if len(sys.argv) > 4 else None
    graph = read_graph(graph_path)
    data = MultiwayData.load(multiway_path)
    node_bins = map_nodes_to_bins(graph, data, liftover)
    chrom_of = data.bins._chrom_of_bin
    covered = data.coverage_raw > 0
    rng = np.random.default_rng(0)
    shift_range = 50 if not data.bins.fixed else max(10, 1500000 // data.bins.bin_size)
    rows = []
    for promoter, attributes in graph.nodes(data=True):
        if not attributes.get('bait') or promoter not in node_bins:
            continue
        partners = sorted({node_bins[n] for n in graph.neighbors(promoter) if n in node_bins} - {node_bins[promoter]})
        if len(partners) < 2:
            continue
        p = node_bins[promoter]
        partners = np.array(partners)
        largest, groups = hub_counts(data, p, partners)
        shifted_largest, shifted_groups = [], []
        tries = 0
        while len(shifted_largest) < N_SHIFTS and tries < 100:
            tries += 1
            shift = int(rng.integers(shift_range // 5, shift_range)) * int(rng.choice([-1, 1]))
            bins = np.concatenate([[p], partners]) + shift
            if bins.min() < 0 or bins.max() >= len(chrom_of) or np.any(chrom_of[bins] != chrom_of[p]) \
                    or not np.all(covered[bins]):
                continue
            a, b = hub_counts(data, bins[0], bins[1:])
            shifted_largest.append(a)
            shifted_groups.append(b)
        if not shifted_largest:
            continue
        rows.append(dict(promoter=promoter, genes=attributes.get('genes', ''), hub_size=1 + len(partners),
                         largest_in_one_group=largest, groups_promoter_plus2=groups,
                         shifted_largest_in_one_group=float(np.mean(shifted_largest)),
                         shifted_groups_promoter_plus2=float(np.mean(shifted_groups)),
                         shifted_any_promoter_plus2=float(np.mean(np.array(shifted_groups) > 0))))
    table = pd.DataFrame(rows)
    table.to_csv(prefix + '_hubs.tsv.gz', sep='\t', index=False)
    summary = []
    for low, high in SIZE_CLASSES:
        t = table[(table['hub_size'] >= low) & (table['hub_size'] <= high)]
        if len(t) == 0:
            continue
        summary.append(dict(
            hub_size='{}-{}'.format(low, high) if high < 1000 else '>={}'.format(low), hubs=len(t),
            hubs_with_promoter_plus2_in_one_group=float(np.mean(t['groups_promoter_plus2'] > 0)),
            shifted=float(t['shifted_any_promoter_plus2'].mean()),
            median_groups_promoter_plus2=float(t['groups_promoter_plus2'].median()),
            shifted_median_groups=float(t['shifted_groups_promoter_plus2'].median()),
            mean_largest_in_one_group=float(t['largest_in_one_group'].mean()),
            shifted_mean_largest=float(t['shifted_largest_in_one_group'].mean()),
            hubs_complete_in_one_group=float(np.mean(t['largest_in_one_group'] == t['hub_size']))))
    summary = pd.DataFrame(summary)
    summary.to_csv(prefix + '_hub_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
