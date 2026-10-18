"""Support classes of the gained and lost contacts and of the explanations of egaDifferential runs (29_threshold_support.sh).

Usage: support_summary.py <dir> <comparison> [<comparison> ...] [--gene comparison=GENE ...]

Writes <dir>/support_summary.tsv (per comparison: gained and lost edges per support class with fractions, the
differential genes, the explained genes, the genes whose first explaining level rests only on threshold-sensitive
contacts and those whose every explaining level does; declared roles, fragment mean), <dir>/support_levels.tsv (per
comparison and level: genes explained there and how many of them rest only on threshold-sensitive contacts, by the
class of the carrying contacts) and <dir>/support_genes.tsv (every level of the named genes with the support of its
carrying contacts and the contacts of the strongest node with the promoter: status, scores, support, reason).
"""

import argparse
import os

import pandas as pd

KEY = 'declared_mean'


def comparison_rows(directory, name):
    prefix = os.path.join(directory, name)
    support = pd.read_csv(prefix + '_support.tsv', sep='\t')
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t')
    row = dict(comparison=name, threshold=support['threshold'].iloc[0], floor=support['floor'].iloc[0])
    for status in ('gained', 'lost', 'gained+lost'):
        part = support[support['status'] == status].set_index('support')
        row[status.replace('+', '_') + '_edges'] = int(part['edges'].sum())
        for cls in ('stable', 'threshold-sensitive', 'not testable'):
            label = '{}_{}'.format(status.replace('+', '_'), cls.replace(' ', '_').replace('-', '_'))
            row[label] = int(part.loc[cls, 'edges'])
            row[label + '_fraction'] = float(part.loc[cls, 'fraction'])
    explained = classes['first_level_' + KEY] != 'none'
    row['differential_genes'] = int(len(classes))
    row['explained_genes'] = int(explained.sum())
    row['first_level_with_contacts'] = int((classes['first_level_support_' + KEY].fillna('') != '').sum())
    for cls in ('stable', 'threshold-sensitive', 'not testable'):
        row['first_level_support_' + cls.replace(' ', '_').replace('-', '_')] = int(
            (classes['first_level_support_' + KEY].fillna('') == cls).sum())
    row['first_level_threshold_sensitive_only'] = int(classes['first_level_threshold_sensitive_' + KEY].sum())
    row['all_levels_threshold_sensitive_only'] = int(classes['only_threshold_sensitive_' + KEY].sum())
    row['any_level_threshold_sensitive_only'] = int(
        (classes['threshold_sensitive_levels_' + KEY].fillna('') != '').sum())
    by_level = []
    hit = levels[levels['explains_' + KEY].astype(bool)]
    for level, part in hit.groupby('level', sort=False):
        support_class = part['support_' + KEY].fillna('')
        by_level.append(dict(comparison=name, rank=int(part['rank'].iloc[0]), level=level, genes_explained=len(part),
                             with_contacts=int((support_class != '').sum()),
                             stable=int((support_class == 'stable').sum()),
                             threshold_sensitive=int((support_class == 'threshold-sensitive').sum()),
                             not_testable=int((support_class == 'not testable').sum()),
                             threshold_sensitive_only=int(part['threshold_sensitive_only_' + KEY].sum()),
                             edges_stable=int(part['edges_stable_' + KEY].sum()),
                             edges_threshold_sensitive=int(part['edges_threshold_sensitive_' + KEY].sum()),
                             edges_not_testable=int(part['edges_not_testable_' + KEY].sum())))
    return row, sorted(by_level, key=lambda r: r['rank'])


def gene_rows(directory, name, gene):
    prefix = os.path.join(directory, name)
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t')
    edges = pd.read_csv(prefix + '_edges.tsv.gz', sep='\t')
    names = edges.columns[edges.columns.str.startswith('score_')].tolist()
    part = levels[levels['genes'].astype(str).str.split(',').map(lambda g: gene in g)]
    rows = []
    for _, r in part.sort_values(['node', 'rank']).iterrows():
        row = dict(comparison=name, gene=gene, node=r['node'], rank=r['rank'], level=r['level'],
                   provenance_chain=r.get('provenance_chain', ''), explains=bool(r['explains_' + KEY]),
                   support=r.get('support_' + KEY, ''), threshold_sensitive_only=bool(r['threshold_sensitive_only_' + KEY]),
                   edges_stable=r.get('edges_stable_' + KEY), edges_threshold_sensitive=r.get('edges_threshold_sensitive_' + KEY),
                   edges_not_testable=r.get('edges_not_testable_' + KEY), support_detail=r.get('support_detail_' + KEY, ''),
                   strongest_node=r.get('strongest_node_declared', ''), detail=r['detail_' + KEY])
        strongest = r.get('strongest_node_declared')
        if isinstance(strongest, str) and strongest:
            pair = edges[((edges['node_a'] == r['node']) & (edges['node_b'] == strongest)) |
                         ((edges['node_b'] == r['node']) & (edges['node_a'] == strongest))]
            if len(pair):
                p = pair.iloc[0]
                row.update(strongest_contact_status=p['status'], strongest_contact_support=p['support'],
                           strongest_contact_reason=p['support_reason'],
                           strongest_contact_other_max_score=p['other_max_score'],
                           **{'strongest_contact_' + c: p[c] for c in names})
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('directory')
    parser.add_argument('comparisons', nargs='+')
    parser.add_argument('--gene', nargs='*', default=[], help='comparison=GENE')
    args = parser.parse_args()
    rows, levels, genes = [], [], []
    for name in args.comparisons:
        row, by_level = comparison_rows(args.directory, name)
        rows.append(row)
        levels += by_level
    for item in args.gene:
        name, gene = item.split('=', 1)
        genes += gene_rows(args.directory, name, gene)
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(args.directory, 'support_summary.tsv'), sep='\t', index=False)
    pd.DataFrame(levels).to_csv(os.path.join(args.directory, 'support_levels.tsv'), sep='\t', index=False)
    pd.DataFrame(genes).to_csv(os.path.join(args.directory, 'support_genes.tsv'), sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_columns', 60, 'display.max_colwidth', 120):
        print(summary.T.to_string())
        print(pd.DataFrame(levels).to_string(index=False))
        if genes:
            print(pd.DataFrame(genes).drop(columns=['detail']).to_string(index=False))


if __name__ == '__main__':
    main()
