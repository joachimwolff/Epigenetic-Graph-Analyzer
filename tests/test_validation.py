import networkx as nx
import numpy as np
import pytest

from epigeneticgraphanalyzer.genome import BinIndex
from epigeneticgraphanalyzer.indirect import direct_contacts, indirect_contacts, two_hop_triplets
from epigeneticgraphanalyzer.multiway import MultiwayData, npmi_from_counts, read_gam_segregation, read_sprite_clusters
from epigeneticgraphanalyzer.validation import (LocusPools, confirm_pairs, coverage_classes, map_nodes_to_bins, pairs_to_bins,
                                                sample_matched_controls, validate_communities, validate_pairs,
                                                validate_triplets)

from synthetic import CHROM_SIZES, simulate_gam, simulate_sprite

BIN = 25000


def _add(graph, chrom, position, bait):
    name = '{}:{}-{}'.format(chrom, position, position + 999)
    graph.add_node(name, chrom=chrom, start=position, end=position + 999, bait=bait, genes='')
    return name


def hub_graph(seed=0, n_hubs=40, n_background=120):
    """Hubs a - v - b that are real (co-occur in the multi-way data) and background baits whose partners are not."""
    rng = np.random.default_rng(seed)
    graph = nx.Graph()
    hubs = []
    used = set()
    for h in range(n_hubs):
        chrom = 'chr1' if h % 2 == 0 else 'chr2'
        size = CHROM_SIZES[chrom]
        while True:
            v = int(rng.integers(30, size // BIN - 30))
            a = v - int(rng.integers(3, 15))
            b = v + int(rng.integers(3, 15))
            if not {(chrom, a), (chrom, v), (chrom, b)} & used:
                break
        used |= {(chrom, a), (chrom, v), (chrom, b)}
        positions = [x * BIN + 5000 for x in (a, v, b)]
        na, nv, nb = _add(graph, chrom, positions[0], False), _add(graph, chrom, positions[1], True), \
            _add(graph, chrom, positions[2], False)
        graph.add_edge(na, nv, weight=1.0)
        graph.add_edge(nv, nb, weight=1.0)
        hubs.append((chrom, positions, (na, nv, nb)))
    for k in range(n_background):
        chrom = 'chr1' if k % 2 == 0 else 'chr2'
        size = CHROM_SIZES[chrom]
        v = int(rng.integers(30, size // BIN - 30))
        nv = _add(graph, chrom, v * BIN + 7000, True)
        for _ in range(2):
            o = v + int(rng.choice([-1, 1]) * rng.integers(3, 15))
            no = _add(graph, chrom, o * BIN + 9000, False)
            if no != nv:
                graph.add_edge(nv, no, weight=1.0)
    return graph, hubs


@pytest.fixture(scope='module')
def sprite_data(tmp_path_factory):
    graph, hubs = hub_graph()
    path = str(tmp_path_factory.mktemp('sprite') / 'clusters.txt')
    simulate_sprite(path, [(c, p) for c, p, _ in hubs], bin_size=BIN)
    data = read_sprite_clusters(path, CHROM_SIZES, BIN)
    return graph, hubs, data


def test_bin_index():
    bins = BinIndex(CHROM_SIZES, BIN)
    assert bins.bin_of('chr1', 0) == 0
    assert bins.bin_of('chr2', 0) == CHROM_SIZES['chr1'] // BIN
    assert bins.bin_of('chr3', 0) == -1
    assert list(bins.bins_of(np.array(['chr1', 'chr2']), np.array([30000, 30000]))) == [1, 161]


def test_npmi():
    values = npmi_from_counts([0, 10, 5], [10, 10, 10], [10, 10, 0], 100)
    assert values[0] == -1
    assert values[1] == pytest.approx(1.0)
    assert np.isnan(values[2])


def test_sprite_parsing_and_cooccurrence(sprite_data, tmp_path):
    graph, hubs, data = sprite_data
    assert data.kind == 'sprite'
    chrom, positions, _ = hubs[0]
    bins = [data.bins.bin_of(chrom, p) for p in positions]
    triple = data.cooccurrence([bins[0]], [bins[1]], [bins[2]], weighted=False)
    assert triple[0] >= 30
    path = str(tmp_path / 'data.npz')
    data.save(path)
    loaded = MultiwayData.load(path)
    assert loaded.incidence.nnz == data.incidence.nnz
    assert np.allclose(loaded.group_weight, data.group_weight)


def test_expected_decays_with_distance(sprite_data):
    _, _, data = sprite_data
    expected = data.expected_by_offset(20)
    assert expected[1] > expected[10] > 0


def test_matched_controls_share_offset_and_pool(sprite_data):
    graph, _, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    pairs = pairs_to_bins(indirect_contacts(graph), node_bins)
    pools = LocusPools(graph, node_bins, data)
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    controls = sample_matched_controls(pairs, pools, excluded, n_controls=5)
    assert len(controls) > 0
    merged = controls.join(pairs[['offset', 'pair_class']], on='test_index')
    assert ((controls['bin_b'] - controls['bin_a']).to_numpy() == merged['offset'].to_numpy()).all()
    assert not any((a, b) in excluded for a, b in zip(controls['bin_a'], controls['bin_b']))
    oe = merged['pair_class'] == 'oe-oe'
    assert np.isin(controls.loc[oe, 'bin_a'], pools.other_array).all()
    assert np.isin(controls.loc[oe, 'bin_b'], pools.other_array).all()


def test_coverage_matched_controls(sprite_data):
    graph, _, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    pairs = pairs_to_bins(direct_contacts(graph), node_bins)
    pools = LocusPools(graph, node_bins, data)
    excluded = set(zip(pairs['bin_a'], pairs['bin_b']))
    classes = coverage_classes(pools, 3)
    assert set(np.unique(classes[np.union1d(pools.bait_array, pools.other_array)])) == {0, 1, 2}
    controls = sample_matched_controls(pairs, pools, excluded, n_controls=5, coverage_class=classes)
    merged = controls.join(pairs[['bin_a', 'bin_b']], on='test_index', rsuffix='_test')
    test_classes = np.sort(np.c_[classes[merged['bin_a_test']], classes[merged['bin_b_test']]], axis=1)
    control_classes = np.sort(np.c_[classes[merged['bin_a']], classes[merged['bin_b']]], axis=1)
    # most controls match the class pair of their test pair; the rest fell back to offset matching
    assert np.mean((test_classes == control_classes).all(axis=1)) > 0.5
    assert ((controls['bin_b'] - controls['bin_a']).to_numpy() ==
            (merged['bin_b_test'] - merged['bin_a_test']).to_numpy()).all()


def test_confirm_pairs_individually(sprite_data):
    graph, hubs, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    direct = pairs_to_bins(direct_contacts(graph), node_bins)
    indirect = pairs_to_bins(indirect_contacts(graph), node_bins)
    excluded = set(zip(direct['bin_a'], direct['bin_b'])) | set(zip(indirect['bin_a'], indirect['bin_b']))
    pools = LocusPools(graph, node_bins, data)
    table, summary = confirm_pairs(direct, data, pools, excluded, 'direct', int(direct['offset'].max()),
                                   n_background=20, n_classes=3)
    hub_nodes = {n for _, _, nodes in hubs for n in nodes}
    is_hub = table['node_a'].isin(hub_nodes) & table['node_b'].isin(hub_nodes)
    assert ((table['confirm_p'] >= 0) & (table['confirm_p'] <= 1)).all()
    assert (table.loc[table['groups_both'] == 0, 'confirm_p'] == 1).all()
    assert (table.loc[is_hub, 'confirm_p'] < 0.05).mean() > 0.8
    assert (table.loc[~is_hub, 'confirm_p'] < 0.05).mean() < 0.3
    overall = summary.iloc[0]
    assert overall['confirmed_p'] > overall['control_confirmed_p']


def test_real_hub_pairs_beat_controls_and_background_does_not(sprite_data):
    graph, hubs, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    pairs = pairs_to_bins(indirect_contacts(graph), node_bins)
    direct = pairs_to_bins(direct_contacts(graph), node_bins)
    excluded = set(zip(pairs['bin_a'], pairs['bin_b'])) | set(zip(direct['bin_a'], direct['bin_b']))
    pools = LocusPools(graph, node_bins, data)
    table, summary = validate_pairs(pairs, data, pools, excluded, 'indirect', n_controls=10, max_offset=40)
    hub_pairs = {(h[2][0], h[2][2]) for h in hubs}
    is_hub = np.array([(a, b) in hub_pairs or (b, a) in hub_pairs for a, b in zip(table['node_a'], table['node_b'])])
    assert is_hub.sum() >= 30
    hub = table[is_hub]
    background = table[~is_hub]
    assert (hub['cooc_oe'] > hub['cooc_oe_control']).mean() > 0.9
    assert abs((background['cooc_oe'] > background['cooc_oe_control']).mean() - 0.5) < 0.25
    row = summary[(summary['metric'] == 'cooc_oe') & (summary['pair_class'] == 'all') & (summary['stratum'] == 'all')]
    assert len(row) == 1


def test_triplet_cooperativity_of_hubs(sprite_data):
    graph, hubs, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    pools = LocusPools(graph, node_bins, data)
    triplets = two_hop_triplets(graph)
    table, summary = validate_triplets(triplets, node_bins, data, pools, graph, n_controls=5)
    hub_mediators = {h[2][1] for h in hubs}
    hub = table[table['mediator'].isin(hub_mediators)]
    assert (hub['c_uvw'] > hub['c_uvw_control']).mean() > 0.9
    assert (hub['log2_cooperativity'] > 0).mean() > 0.9
    assert set(summary['metric']) == {'c_uvw'}


def test_communities_against_shifted_sets(sprite_data):
    graph, hubs, data = sprite_data
    node_bins = map_nodes_to_bins(graph, data)
    communities = [set(h[2]) for h in hubs]
    table, summary = validate_communities(communities, node_bins, data, min_members=3, n_shifts=20)
    assert len(table) == len(hubs)
    assert (table['log2_enrichment_shift'] > 1).mean() > 0.9
    assert summary['wilcoxon_p'].iloc[0] < 1e-3


def test_gam_segregation_validation(tmp_path):
    graph, hubs = hub_graph(seed=1)
    path = str(tmp_path / 'gam.tsv')
    simulate_gam(path, [(c, p) for c, p, _ in hubs])
    data = read_gam_segregation([path], CHROM_SIZES)
    assert data.kind == 'gam'
    assert data.bins.bin_size == 25000
    node_bins = map_nodes_to_bins(graph, data)
    pairs = pairs_to_bins(indirect_contacts(graph), node_bins)
    pools = LocusPools(graph, node_bins, data)
    table, summary = validate_pairs(pairs, data, pools, set(zip(pairs['bin_a'], pairs['bin_b'])), 'indirect',
                                    max_offset=40)
    assert 'npmi_delta' in table
    hub_pairs = {(h[2][0], h[2][2]) for h in hubs}
    is_hub = np.array([(a, b) in hub_pairs for a, b in zip(table['node_a'], table['node_b'])])
    assert (table.loc[is_hub, 'npmi'] > table.loc[is_hub, 'npmi_control']).mean() > 0.9


def test_gam_replicates_concatenate(tmp_path):
    graph, hubs = hub_graph(seed=2, n_hubs=4, n_background=4)
    p1 = simulate_gam(str(tmp_path / 'r1.tsv'), [(c, p) for c, p, _ in hubs], n_tubes=50, seed=1)
    p2 = simulate_gam(str(tmp_path / 'r2.tsv'), [(c, p) for c, p, _ in hubs], n_tubes=70, seed=2)
    data = read_gam_segregation([p1, p2], CHROM_SIZES)
    assert data.n_groups == 120


def test_sprite_on_fragments(tmp_path):
    from epigeneticgraphanalyzer.genome import IntervalIndex
    # fragments of variable length tiling chr1 and chr2, with a gap on chr2
    rng = np.random.default_rng(3)
    chroms, starts, ends = [], [], []
    for chrom, size in CHROM_SIZES.items():
        position = 0
        while position < size - 20000:
            length = int(rng.integers(1000, 9000))
            if chrom == 'chr2' and 1000000 <= position < 1100000:
                position += length
                continue
            chroms.append(chrom)
            starts.append(position)
            ends.append(position + length)
            position += length
    bed = tmp_path / 'fragments.bed'
    bed.write_text(''.join('{}\t{}\t{}\n'.format(*x) for x in zip(chroms, starts, ends)))
    index = IntervalIndex.from_bed(str(bed))
    assert index.bin_of('chr2', 1050000) == -1
    graph, hubs = hub_graph(seed=4, n_hubs=20, n_background=40)
    clusters = str(tmp_path / 'clusters.txt')
    simulate_sprite(clusters, [(c, p) for c, p, _ in hubs], n_background=8000)
    data = read_sprite_clusters(clusters, CHROM_SIZES, None, index=index)
    path = str(tmp_path / 'fragments.npz')
    data.save(path)
    data = MultiwayData.load(path)
    assert not data.bins.fixed
    chrom, positions, nodes = hubs[0]
    fragments = [data.bins.bin_of(chrom, p) for p in positions]
    # reads are jittered by up to 2 kb, so part of a hub's reads falls into neighbouring fragments
    assert data.cooccurrence([fragments[0]], [fragments[1]], [fragments[2]], weighted=False)[0] >= 5
    node_bins = map_nodes_to_bins(graph, data)
    pairs = pairs_to_bins(indirect_contacts(graph), node_bins)
    pools = LocusPools(graph, node_bins, data)
    table, _ = validate_pairs(pairs, data, pools, set(zip(pairs['bin_a'], pairs['bin_b'])), 'indirect',
                              n_controls=5)
    hub_pairs = {(h[2][0], h[2][2]) for h in hubs}
    is_hub = np.array([(a, b) in hub_pairs for a, b in zip(table['node_a'], table['node_b'])])
    assert (table.loc[is_hub, 'cooc'] > table.loc[is_hub, 'cooc_control']).mean() > 0.9
