"""Liftover of graph nodes to the assembly of the multi-way data: unmapped, multi-mapped and resized nodes.

The pair validation lifts the midpoint of every node with the chain file and uses the first hit
(genome.LiftOver.position); nodes whose midpoint is not lifted are dropped. This report lifts midpoint, start and
end of every node and counts nodes whose midpoint is unmapped, whose midpoint has more than one hit, whose midpoint
lands on another chromosome, whose start or end is unmapped or lands on another chromosome than the midpoint, and
whose lifted length (|lifted end - lifted start|) differs by more than 10% from the original length.

Usage: python liftover_report.py <chain> <output.tsv> <name>=<graph.graphml> [<name>=<graph.graphml> ...]
"""

import argparse

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.graph import read_graph


def node_liftover(liftover, chrom, start, end):
    """Liftover status of one node: dict of flags and the lifted length."""
    middle = (start + end) // 2
    hits = liftover.hits(chrom, middle)
    flags = dict(midpoint_unmapped=not hits, midpoint_multi=len(hits) > 1,
                 midpoint_other_chrom=bool(hits) and hits[0][0] != chrom)
    first = liftover.hits(chrom, start)
    last = liftover.hits(chrom, max(start, end - 1))
    flags['end_unmapped'] = not first or not last
    target = hits[0][0] if hits else None
    flags['end_other_chrom'] = bool(first and last) and (first[0][0] != last[0][0] or
                                                           (target is not None and first[0][0] != target))
    flags['length_changed_10pct'] = False
    if first and last and not flags['end_other_chrom']:
        length = max(1, end - 1 - start)
        lifted = abs(last[0][1] - first[0][1])
        flags['length_changed_10pct'] = abs(lifted - length) > 0.1 * length
    return flags


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('chain')
    parser.add_argument('output')
    parser.add_argument('graphs', nargs='+')
    args = parser.parse_args()
    liftover = LiftOver(args.chain)
    rows = []
    for item in args.graphs:
        name, path = item.split('=', 1)
        graph = read_graph(path)
        flags = pd.DataFrame([node_liftover(liftover, d['chrom'], int(d['start']), int(d['end']))
                              for _, d in graph.nodes(data=True)])
        bait = np.array([bool(d.get('bait')) for _, d in graph.nodes(data=True)])
        row = dict(graph=name, chain=args.chain, nodes=len(flags), baits=int(bait.sum()))
        for column in flags.columns:
            row[column] = int(flags[column].sum())
            row[column + '_fraction'] = float(flags[column].mean())
            row[column + '_baits'] = int(flags[column][bait].sum())
        any_issue = flags.any(axis=1)
        row['any_issue'] = int(any_issue.sum())
        row['any_issue_fraction'] = float(any_issue.mean())
        print(row, flush=True)
        rows.append(row)
    pd.DataFrame(rows).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
