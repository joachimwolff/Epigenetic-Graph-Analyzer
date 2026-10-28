"""Domain-matched controls and the block bootstrap."""

import numpy as np

from epigeneticgraphanalyzer.genome import BinIndex
from epigeneticgraphanalyzer.validation import _pair_code, block_bootstrap, domain_annotation


def test_domain_annotation_and_pair_codes(tmp_path):
    bins = BinIndex({'chr1': 100000}, 10000)
    compartments = tmp_path / 'compartments.bed'
    compartments.write_text('chr1\t0\t50000\tA1\nchr1\t50000\t100000\tB2\n')
    domains = tmp_path / 'domains.txt'
    # Rao-style domain list: header, chromosome without 'chr', a nested inner domain
    domains.write_text('chr1\tx1\tx2\tchr2\ty1\ty2\n1\t0\t60000\t1\t0\t60000\n1\t10000\t30000\t1\t10000\t30000\n')
    compartment, tad_id = domain_annotation(bins, str(compartments), str(domains))
    assert list(compartment) == [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    # bins 1 and 2 share the inner domain, bin 0 and bin 4 only the outer one, bins 7-9 none
    assert tad_id[1] == tad_id[2] != tad_id[0]
    assert tad_id[0] == tad_id[4] and tad_id[8] == -1
    classes = np.zeros(10, dtype=np.int64)
    codes = _pair_code(classes, tad_id, np.array([1, 1, 7]), np.array([2, 5, 8]))
    assert codes[0] % 2 == 1 and codes[1] % 2 == 0 and codes[2] % 2 == 0


def test_block_bootstrap_keeps_blocks_together():
    rng = np.random.default_rng(0)
    blocks = np.repeat(np.arange(20), 50)
    # strong block effect: items of one block share their value
    values = np.repeat(rng.normal(0.5, 0.1, 20), 50)
    low, high = block_bootstrap(values, blocks, n_boot=500)
    naive = 1.96 * values.std() / np.sqrt(len(values))
    assert high - low > 4 * naive          # wider than the interval that treats items as independent
    assert low < values.mean() < high
