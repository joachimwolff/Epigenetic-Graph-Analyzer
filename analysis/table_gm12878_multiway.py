"""Supplementary table: multi-way support of the explaining paths of genes explained at step 2 (naive B vs GM12878).

Rows: genes without an explanation at the promoter (marks, contacts) or step 1 (marks) and explained at step 2 with a
per-gene P below --maxP; columns: expression change, step-2 P and q, explaining paths present in GM12878 and present
only in naive B cells, and per multi-way data set the number of present paths that occur together above their
shifted controls (P < 0.05). The output holds the table rows; the longtable frame and caption are in the supplement.

Usage: python table_gm12878_multiway.py <cmp prefix> <multiway_paths prefix> <output.tex> [--maxP 0.01]
"""

import argparse

import pandas as pd

DATA_SETS = ['SPRITE', 'HiPore-C', 'Pore-C']


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('prefix')
    parser.add_argument('paths')
    parser.add_argument('output')
    parser.add_argument('--maxP', type=float, default=0.01)
    args = parser.parse_args()
    triplets = pd.read_csv(args.paths + '_triplets.tsv', sep='\t')
    triplets['significant'] = triplets.p_empirical < 0.05
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    step2 = levels[levels.level == 'step2'].set_index('node')
    tested = triplets.dropna(subset=['control_mean'])
    lines = []
    explaining = tested[tested.set == 'explaining']
    rows = []
    for promoter, sub in explaining.groupby('promoter'):
        info = step2.loc[promoter]
        if not info['p_declared'] < args.maxP:
            continue
        present = sub[sub.present]
        counts = {d.replace('-', ''): int(present[present.data_set == d].significant.sum()) for d in DATA_SETS}
        rows.append(dict(genes=info['genes'], fc=info['rna_log2FC'], p=info['p_declared'], q=info['q_declared'],
                         present=present.drop_duplicates(['mediator', 'target']).shape[0],
                         absent=sub[~sub.present].drop_duplicates(['mediator', 'target']).shape[0], **counts))
    table = pd.DataFrame(rows).sort_values('p')
    for r in table.itertuples():
        genes = ', '.join(g for g in str(r.genes).split(',') if not g.startswith(('RP11-', 'RP1-', 'RP4-', 'RP5-',
                                                                                      'CTD-', 'CTB-', 'AC0', 'AL')))
        lines.append(r'{} & {:+.2f} & {:.1e} & {:.2g} & {} & {} & {} & {} & {} \\'.format(
            genes.replace('_', r'\_') or r.genes, r.fc, r.p, r.q, r.present, r.absent, r.SPRITE, r.HiPoreC, r.PoreC))
    open(args.output, 'w').write('\n'.join(lines) + '\n')
    print('{} genes with step-2 P < {}'.format(len(table), args.maxP))


if __name__ == '__main__':
    main()
