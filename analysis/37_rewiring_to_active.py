"""Gained and lost contacts to elements that are active in both conditions, as explanations of expression changes.

Only genes that are not explained otherwise are analysed: differentially expressed genes without a concordant level
at the promoter (marks, contacts) or at step 1 (marks) (--earlyLevels; unchanged states such as 'poised in both' do not
count). For these genes, a rewiring path is a shortest path promoter - ... - target (target at step 1 or
step 2) of which at least one interaction changed (gained or lost) without a call in the other condition (support
'stable' or 'not testable'; 'threshold-sensitive' interactions, called in a replicate of the other condition or just
below the threshold there, are excluded; --support), whose target is active in both conditions (an
active mark at or above the --activeQuantile of all nodes in both conditions) and whose marks do not change
(|log2 change| < --maxChange for every mark). A path is concordant when its changed interactions are gained for a gene
that increases or lost for a gene that decreases.

Candidates: the unexplained genes with a concordant rewiring path at step 2, one row per gene.

Outputs: <out>_paths.tsv (every rewiring path of the unexplained genes) and <out>_candidates.tsv.

Usage: python 37_rewiring_to_active.py <cmp prefix> <out prefix> --conditions A B --active H3K27ac [H3K4me3 ...]
           [--marks H3K27ac H3K4me1 ...] [--activeQuantile 0.9] [--maxChange 1]
"""

import argparse

import networkx as nx
import numpy as np
import pandas as pd

EARLY = ['promoter_marks', 'promoter_contacts', 'step1_marks']


def log2_change(attrs, mark, a, b):
    return float(np.log2((float(attrs.get('{}_{}'.format(mark, b), 0) or 0) + 1) /
                         (float(attrs.get('{}_{}'.format(mark, a), 0) or 0) + 1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('prefix')
    parser.add_argument('output')
    parser.add_argument('--conditions', nargs=2, required=True)
    parser.add_argument('--active', nargs='+', required=True, help='Marks that define an active element')
    parser.add_argument('--marks', nargs='+', help='All marks (default: the active marks)')
    parser.add_argument('--activeQuantile', type=float, default=0.9)
    parser.add_argument('--maxChange', type=float, default=1.0)
    parser.add_argument('--support', nargs='+', default=['stable', 'not testable'])
    parser.add_argument('--earlyLevels', nargs='+', default=EARLY,
                        help='Levels that explain a gene before step 2; only genes with none of them are analysed')
    args = parser.parse_args()
    a, b = args.conditions
    marks = args.marks or args.active
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    levels['explains'] = levels['explains_declared_mean'].astype(str) == 'True'
    early = levels[levels.level.isin(args.earlyLevels)].groupby('node').explains.any()
    genes = levels.drop_duplicates('node').set_index('node')[['genes', 'rna_log2FC']]
    genes = genes[~genes.index.map(early).fillna(False).astype(bool)]
    print('{} differentially expressed genes without an explanation at the promoter or step 1'.format(len(genes)),
          flush=True)
    nodes = pd.DataFrame.from_dict(dict(graph.nodes(data=True)), orient='index')
    thresholds = {m: np.nanquantile(np.concatenate([nodes['{}_{}'.format(m, c)].astype(float).values
                                                    for c in (a, b)]), args.activeQuantile) for m in args.active}
    print('active thresholds', {m: round(v, 2) for m, v in thresholds.items()}, flush=True)

    def active(attrs):
        return any(min(float(attrs.get('{}_{}'.format(m, c), 0) or 0) for c in (a, b)) >= t
                   for m, t in thresholds.items())

    def unchanged(attrs):
        return all(abs(log2_change(attrs, m, a, b)) < args.maxChange for m in marks)

    rows = []
    for node, info in genes.iterrows():
        if node not in graph:
            continue
        sign = np.sign(info.rna_log2FC)
        distance = nx.single_source_shortest_path_length(graph, node, cutoff=2)
        for target, k in distance.items():
            if k == 0 or not active(graph.nodes[target]) or not unchanged(graph.nodes[target]):
                continue
            for path in nx.all_shortest_paths(graph, node, target):
                edges = [graph.edges[u, v] for u, v in zip(path[:-1], path[1:])]
                changed = [e for e in edges if e.get('status') in ('gained', 'lost')]
                if not changed or any(e.get('support') not in args.support for e in changed):
                    continue
                statuses = {e['status'] for e in changed}
                if len(statuses) > 1:
                    continue
                status = statuses.pop()
                attrs = graph.nodes[target]
                rows.append(dict(node=node, genes=info.genes, rna_log2FC=info.rna_log2FC, step=k, target=target,
                                 target_genes=attrs.get('genes') or '', mediator=path[1] if k == 2 else '',
                                 status=status, concordant=(status == 'gained') == (sign > 0),
                                 distance_kb=round(abs((int(attrs['start']) + int(attrs['end'])) / 2 -
                                                       (int(graph.nodes[node]['start']) +
                                                        int(graph.nodes[node]['end'])) / 2) / 1000, 1),
                                 **{'{}_{}'.format(m, c): round(float(attrs.get('{}_{}'.format(m, c), 0) or 0), 2)
                                    for m in marks for c in (a, b)}))
    paths = pd.DataFrame(rows, columns=['node', 'genes', 'rna_log2FC', 'step', 'target', 'target_genes', 'mediator',
                                        'status', 'concordant', 'distance_kb'] +
                         ['{}_{}'.format(m, c) for m in marks for c in (a, b)])
    paths.to_csv(args.output + '_paths.tsv', sep='\t', index=False)
    concordant = paths[(paths.step == 2) & paths.concordant]
    candidates = concordant.groupby(['node', 'genes', 'rna_log2FC']).agg(
        paths=('target', 'size'), targets=('target', 'nunique'),
        target_genes=('target_genes', lambda x: ','.join(sorted({g for v in x for g in str(v).split(',') if g}))),
        mediators=('mediator', 'nunique'), status=('status', 'first'),
        min_distance_kb=('distance_kb', 'min')).reset_index().sort_values(['targets', 'rna_log2FC'],
                                                                         ascending=[False, False])
    candidates.to_csv(args.output + '_candidates.tsv', sep='\t', index=False)
    print('{} of them with a concordant rewiring path at step 2'.format(len(candidates)))
    print(candidates.head(25).to_string(index=False))


if __name__ == '__main__':
    main()
