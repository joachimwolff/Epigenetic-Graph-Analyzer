"""Aggregate multi-way co-occurrence around graph pairs at single-fragment resolution.

For every tested pair (a, b) the co-occurrence of (a + i, b + j), i, j in -W..W fragments, is divided by the
mean co-occurrence at its offset; the pile-up is the sum of observed over the sum of expected values of all
pairs per cell (positions on uncovered fragments are skipped). A focal maximum at (0, 0) shows that the
multi-way data see the contact of exactly these fragments and not only of their neighbourhood. Matched control
pairs (same offset, same locus pools) are piled up alike.

Usage: python sprite_pileup.py <graph.graphml> <multiway.npz> <output prefix> [window=10] [pairs per set=20000]
"""

import sys
from itertools import product

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.graph import read_graph
from epigeneticgraphanalyzer.indirect import direct_contacts, indirect_contacts
from epigeneticgraphanalyzer.multiway import MultiwayData
from epigeneticgraphanalyzer.validation import LocusPools, map_nodes_to_bins, pairs_to_bins, sample_matched_controls

OFFSET_STRATA = ((21, 50), (51, 150), (151, 100000))


def pileup(data, bin_a, bin_b, window, expected, coverage_normalised=False):
    """Pile-up of observed over expected co-occurrence; with `coverage_normalised` the expectation of a cell
    is also scaled by the coverage of both fragments relative to the mean coverage of covered fragments."""
    chrom_of = data.bins._chrom_of_bin
    covered = data.coverage_raw > 0
    relative = data.coverage_raw / data.coverage_raw[covered].mean()
    n_bins = len(chrom_of)
    size = 2 * window + 1
    observed = np.zeros((size, size))
    expect = np.zeros((size, size))
    for i in range(-window, window + 1):
        for j in range(-window, window + 1):
            a, b = bin_a + i, bin_b + j
            ok = (a >= 0) & (b < n_bins) & (b - a >= 1)
            a, b = a[ok], b[ok]
            ok = (chrom_of[a] == chrom_of[bin_a[ok]]) & (chrom_of[b] == chrom_of[bin_b[ok]]) & covered[a] & covered[b]
            a, b = a[ok], b[ok]
            offset = b - a
            inside = offset < len(expected)
            a, b, offset = a[inside], b[inside], offset[inside]
            observed[i + window, j + window] = data.cooccurrence(a, b).sum()
            scale = relative[a] * relative[b] if coverage_normalised else 1.0
            expect[i + window, j + window] = np.nansum(expected[offset] * scale)
    return observed / expect


def main():
    graph_path, multiway_path, prefix = sys.argv[1:4]
    window = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    n_pairs = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
    graph = read_graph(graph_path)
    data = MultiwayData.load(multiway_path)
    node_bins = map_nodes_to_bins(graph, data)
    direct = pairs_to_bins(direct_contacts(graph), node_bins)
    indirect = pairs_to_bins(indirect_contacts(graph), node_bins)
    excluded = set(zip(direct['bin_a'], direct['bin_b'])) | set(zip(indirect['bin_a'], indirect['bin_b']))
    pools = LocusPools(graph, node_bins, data)
    max_offset = int(max(direct['offset'].max(), indirect['offset'].max())) + 2 * window + 1
    expected = data.expected_by_offset(max_offset, 'cooc')
    sets = {'direct': direct, 'indirect_oe-oe': indirect[indirect['pair_class'] == 'oe-oe']}
    rows, matrices = [], {}
    rng = np.random.default_rng(0)
    for name, pairs in sets.items():
        for low, high in OFFSET_STRATA:
            subset = pairs[(pairs['offset'] >= low) & (pairs['offset'] <= high)]
            if len(subset) > n_pairs:
                subset = subset.iloc[rng.choice(len(subset), n_pairs, replace=False)]
            subset = subset.reset_index(drop=True)
            controls = sample_matched_controls(subset, pools, excluded, n_controls=1, seed=1)
            for (kind, a, b), normalisation in product(
                    (('test', subset['bin_a'], subset['bin_b']), ('control', controls['bin_a'], controls['bin_b'])),
                    ('distance', 'distance+coverage')):
                matrix = pileup(data, np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64), window, expected,
                                coverage_normalised=normalisation != 'distance')
                key = '{}|{}-{}|{}|{}'.format(name, low, high, kind, normalisation)
                matrices[key] = matrix
                centre = matrix[window, window]
                ring = matrix[window - 1:window + 2, window - 1:window + 2].copy()
                ring[1, 1] = np.nan
                corners = np.concatenate([matrix[:3, :3].ravel(), matrix[:3, -3:].ravel(),
                                          matrix[-3:, :3].ravel(), matrix[-3:, -3:].ravel()])
                rows.append(dict(set=name, offset_low=low, offset_high=high, kind=kind, normalisation=normalisation,
                                 n=len(a),
                                 centre=centre, neighbours=np.nanmean(ring), corners=np.mean(corners),
                                 centre_over_neighbours=centre / np.nanmean(ring),
                                 centre_over_corners=centre / np.mean(corners)))
                print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(prefix + '_pileup_summary.tsv', sep='\t', index=False)
    np.savez_compressed(prefix + '_pileup.npz', keys=np.array(list(matrices)), window=window,
                        matrices=np.stack([matrices[k] for k in matrices]))


if __name__ == '__main__':
    main()
