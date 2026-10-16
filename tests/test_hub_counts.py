import os
import sys
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from hub_multiway import hub_counts, hub_counts_dense  # noqa: E402


def test_hub_counts_matches_dense_reference():
    rng = np.random.default_rng(0)
    incidence = sp.csr_matrix((rng.random((60, 400)) < 0.05).astype(np.float64))
    data = SimpleNamespace(incidence=incidence, n_groups=400)
    for _ in range(200):
        bins = rng.choice(60, size=int(rng.integers(3, 12)), replace=False)
        assert hub_counts(data, bins[0], np.sort(bins[1:])) == hub_counts_dense(data, bins[0], np.sort(bins[1:]))


def test_hub_counts_promoter_only():
    incidence = sp.csr_matrix(np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64))
    data = SimpleNamespace(incidence=incidence, n_groups=3)
    assert hub_counts(data, 0, np.array([1, 2])) == hub_counts_dense(data, 0, np.array([1, 2])) == (1, 0)
