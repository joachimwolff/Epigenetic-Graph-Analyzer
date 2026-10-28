"""LiftOver.hits and the node liftover report."""

import os
import sys

from epigeneticgraphanalyzer.genome import LiftOver

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'analysis'))
from liftover_report import node_liftover  # noqa: E402


def chain_file(tmp_path):
    # chr1:0-10000 -> chr1:100-10100; chr1:20000-30000 -> chr2:0-10000 (gap 10000-20000 unmapped)
    path = tmp_path / 'test.chain'
    path.write_text('chain 1000 chr1 50000 + 0 10000 chr1 50000 + 100 10100 1\n10000\n\n'
                    'chain 900 chr1 50000 + 20000 30000 chr2 50000 + 0 10000 2\n10000\n\n')
    return str(path)


def test_hits_and_position(tmp_path):
    liftover = LiftOver(chain_file(tmp_path))
    assert liftover.hits('chr1', 50) == [('chr1', 150)]
    assert liftover.position('chr1', 50) == ('chr1', 150)
    assert liftover.hits('chr1', 15000) == [] and liftover.position('chr1', 15000) is None
    assert liftover.hits('chr1', 25000) == [('chr2', 5000)]


def test_node_liftover_flags(tmp_path):
    liftover = LiftOver(chain_file(tmp_path))
    inside = node_liftover(liftover, 'chr1', 1000, 2000)
    assert not any(inside.values())
    spanning = node_liftover(liftover, 'chr1', 9000, 21000)     # midpoint in the gap, ends on two chromosomes
    assert spanning['midpoint_unmapped'] and spanning['end_other_chrom']
    moved = node_liftover(liftover, 'chr1', 21000, 22000)
    assert moved['midpoint_other_chrom'] and not moved['length_changed_10pct']
