"""Multicore multi-way validation: identical results for 1 and 4 worker processes, and equal to the serial code."""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import parallel
from epigeneticgraphanalyzer.genome import IntervalIndex
from epigeneticgraphanalyzer.graph import write_graph
from epigeneticgraphanalyzer.indirect import direct_contacts, indirect_contacts, two_hop_triplets
from epigeneticgraphanalyzer.multiway import MultiwayData, read_sprite_clusters
from epigeneticgraphanalyzer.validation import (LocusPools, confirm_pairs, map_nodes_to_bins, pairs_to_bins,
                                                validate_communities, validate_pairs, validate_triplets)

from synthetic import CHROM_SIZES, simulate_sprite
from test_validation import hub_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope='module')
def fragment_data(tmp_path_factory):
    """SPRITE clusters on variable-length fragments (no co-occurrence band: every pair query is sparse)."""
    tmp = tmp_path_factory.mktemp('parallel')
    rng = np.random.default_rng(3)
    rows = []
    for chrom, size in CHROM_SIZES.items():
        position = 0
        while position < size - 20000:
            length = int(rng.integers(1000, 9000))
            rows.append((chrom, position, position + length))
            position += length
    bed = tmp / 'fragments.bed'
    bed.write_text(''.join('{}\t{}\t{}\n'.format(*x) for x in rows))
    graph, hubs = hub_graph(seed=4, n_hubs=30, n_background=60)
    for k, (_, _, nodes) in enumerate(hubs):
        for n in nodes:
            graph.nodes[n]['community'] = k
    clusters = str(tmp / 'clusters.txt')
    simulate_sprite(clusters, [(c, p) for c, p, _ in hubs], n_background=8000)
    data = read_sprite_clusters(clusters, CHROM_SIZES, None, index=IntervalIndex.from_bed(str(bed)))
    npz = str(tmp / 'fragments.npz')
    data.save(npz)
    graph_path = str(tmp / 'graph.graphml')
    write_graph(graph, graph_path)
    return graph, hubs, npz, graph_path, tmp


def _run_all(graph, hubs, npz, threads):
    data = MultiwayData.load(npz).set_threads(threads)
    data.batch_rows = 500   # several batches per query, so that the parallel path is used
    node_bins = map_nodes_to_bins(graph, data)
    direct = pairs_to_bins(direct_contacts(graph), node_bins)
    indirect = pairs_to_bins(indirect_contacts(graph), node_bins)
    excluded = set(zip(direct['bin_a'], direct['bin_b'])) | set(zip(indirect['bin_a'], indirect['bin_b']))
    pools = LocusPools(graph, node_bins, data)
    max_offset = int(max(direct['offset'].max(), indirect['offset'].max()))
    out = {}
    out['pairs'] = validate_pairs(indirect, data, pools, excluded, 'indirect', n_controls=10, max_offset=max_offset,
                                  coverage_matched=3)
    out['confirm'] = confirm_pairs(direct, data, pools, excluded, 'direct', max_offset, n_background=20, n_classes=3)
    out['triplets'] = validate_triplets(two_hop_triplets(graph), node_bins, data, pools, graph, n_controls=5)
    out['communities'] = validate_communities([set(h[2]) for h in hubs], node_bins, data, n_shifts=20)
    out['expected'] = data.expected_by_offset(max_offset, 'npmi')
    return out


def test_validation_identical_for_1_and_4_threads(fragment_data):
    graph, hubs, npz, _, _ = fragment_data
    one = _run_all(graph, hubs, npz, 1)
    four = _run_all(graph, hubs, npz, 4)
    for key in ('pairs', 'confirm', 'triplets', 'communities'):
        for a, b in zip(one[key], four[key]):
            pd.testing.assert_frame_equal(a, b, check_exact=True)
    np.testing.assert_array_equal(one['expected'], four['expected'])
    assert len(one['communities'][0]) == len(hubs)


def _sampled_expected_serial(data, max_offset, metric, samples_per_offset=20000, seed=0):
    """The single-stream implementation before the parallel version (one query per offset)."""
    rng = np.random.default_rng(seed)
    covered = np.nonzero(data.coverage_raw > 0)[0]
    chrom_of = data.bins._chrom_of_bin
    expected = np.full(max_offset + 1, np.nan)
    covered_mask = data.coverage_raw > 0
    for offset in range(1, max_offset + 1):
        left = covered[covered + offset < len(chrom_of)]
        left = left[covered_mask[left + offset] & (chrom_of[left] == chrom_of[left + offset])]
        if len(left) == 0:
            continue
        if len(left) > samples_per_offset:
            left = rng.choice(left, samples_per_offset, replace=False)
        values = data.cooccurrence(left, left + offset) if metric == 'cooc' else data.npmi(left, left + offset)
        expected[offset] = np.nanmean(values)
    return expected


def test_sampled_expected_equals_serial_implementation(fragment_data):
    _, _, npz, _, _ = fragment_data
    data = MultiwayData.load(npz)
    for metric in ('cooc', 'npmi'):
        serial = _sampled_expected_serial(data, 60, metric, samples_per_offset=300)
        data._expected.clear()
        data.set_threads(4)
        data.batch_rows = 500
        from epigeneticgraphanalyzer.multiway import _sampled_expected_impl
        np.testing.assert_array_equal(_sampled_expected_impl(data, 60, metric, samples_per_offset=300), serial)


def _script(name, *args):
    env = dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get('PYTHONPATH', ''))
    subprocess.run([sys.executable, os.path.join(ROOT, 'analysis', name)] + [str(a) for a in args], check=True,
                   env=env, stdout=subprocess.DEVNULL)


def test_hop_distance_and_hubs_identical_for_1_and_4_threads(fragment_data):
    _, _, npz, graph_path, tmp = fragment_data
    for threads in (1, 4):
        _script('hop_distance.py', graph_path, npz, tmp / 'hops_{}.tsv'.format(threads), '--pairs', 300,
                '--bootstrap', 50, '--threads', threads)
        _script('hub_multiway.py', graph_path, npz, tmp / 'hubs_{}'.format(threads), '--threads', threads)
    pd.testing.assert_frame_equal(pd.read_csv(tmp / 'hops_1.tsv', sep='\t'), pd.read_csv(tmp / 'hops_4.tsv', sep='\t'),
                                  check_exact=True)
    pd.testing.assert_frame_equal(pd.read_csv(tmp / 'hubs_1_hubs.tsv.gz', sep='\t'),
                                  pd.read_csv(tmp / 'hubs_4_hubs.tsv.gz', sep='\t'), check_exact=True)


def test_map_chunks_keeps_order_and_thread_bounds():
    parallel.share(offset=10)
    assert parallel.map_chunks(_add_offset, range(20), 4) == [x + 10 for x in range(20)]
    assert 1 <= parallel.default_threads() <= parallel.MAX_THREADS
    assert parallel.resolve_threads(3) == 3 and parallel.resolve_threads(0) == parallel.default_threads()


def _add_offset(x):
    return x + parallel.shared('offset')
