import networkx as nx
import numpy as np
import pytest

from epigeneticgraphanalyzer.validation import (assortativity_from_moments, assortativity_moments, block_bootstrap,
                                                block_resample_items, block_resamples, edge_assortativity)


def _reference_block_bootstrap(values, blocks, n_boot, seed):
    """The original implementation, kept to check that the refactoring does not change any interval."""
    labels, index = np.unique(blocks, return_inverse=True)
    members = [np.nonzero(index == k)[0] for k in range(len(labels))]
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_boot)
    for b in range(n_boot):
        chosen = rng.integers(0, len(labels), len(labels))
        estimates[b] = np.nanmean(values[np.concatenate([members[k] for k in chosen])])
    return float(np.nanpercentile(estimates, 2.5)), float(np.nanpercentile(estimates, 97.5))


def test_block_resamples_keep_blocks_together():
    blocks = np.array(['chr2', 'chr1', 'chr2', 'chr3', 'chr1', 'chr2'])
    labels, resamples = block_resamples(blocks, n_boot=20, seed=1)
    assert list(labels) == ['chr1', 'chr2', 'chr3']
    assert len(resamples) == 20 and all(len(r) == 3 for r in resamples)
    for chosen, items in zip(resamples, block_resample_items(blocks, n_boot=20, seed=1)):
        # every chosen block contributes all of its items, once per draw
        expected = sorted(i for k in chosen for i in np.nonzero(blocks == labels[k])[0])
        assert sorted(items) == expected


def test_block_bootstrap_unchanged():
    rng = np.random.default_rng(3)
    blocks = rng.choice(['chr{}'.format(i) for i in range(1, 12)], 400)
    values = rng.normal(size=400)
    assert block_bootstrap(values, blocks, n_boot=200, seed=5) == _reference_block_bootstrap(values, blocks, 200, 5)


def test_edge_assortativity_matches_networkx():
    rng = np.random.default_rng(0)
    graph = nx.gnm_random_graph(60, 200, seed=2)
    value = {n: float(rng.integers(0, 6)) for n in graph}
    nx.set_node_attributes(graph, value, 'x')
    u, v = zip(*graph.edges())
    ours = edge_assortativity([value[a] for a in u], [value[b] for b in v])
    assert ours == pytest.approx(nx.numeric_assortativity_coefficient(graph, 'x'), abs=1e-10)
    # moments of disjoint edge sets add up
    half = len(u) // 2
    m1 = assortativity_moments([value[a] for a in u[:half]], [value[b] for b in v[:half]])
    m2 = assortativity_moments([value[a] for a in u[half:]], [value[b] for b in v[half:]])
    assert assortativity_from_moments(*np.add(m1, m2)) == pytest.approx(ours, abs=1e-12)
    # edge direction does not matter
    assert edge_assortativity([value[b] for b in v], [value[a] for a in u]) == pytest.approx(ours, abs=1e-12)
