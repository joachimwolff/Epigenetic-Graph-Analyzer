"""Genes explained only beyond the promoter and its direct contacts, over several comparisons, one row per gene.

A gene is listed when it has no explanation at promoter_marks, promoter_contacts or step1_marks, is explained at
step 2 or step 3 (first such level) with q < --q beyond the mirrored loci, and its expression change does not rest on
zero counts (baseMean >= --baseMean, lfcSE < --maxSE). Columns: comparison, gene, DESeq2 log2 fold change, adjusted P
and baseMean, step, explaining marks, concordant nodes of the nodes at the step, per-gene P and q, the direct contacts
(shared / gained / lost), and with a multi-way table of 36_gm12878_multiway_paths.py for the comparison the number of
explaining paths present in the condition of the multi-way data and of these the paths with P < 0.05 in each data set.

Usage: python table_explained_beyond.py <output prefix> --comparison LABEL=<cmp prefix>,<deseq2.tsv>[,<multiway prefix>]
           ... [--q 0.05] [--baseMean 50] [--maxSE 1]
Outputs: <output prefix>.tsv and <output prefix>_rows.tex (longtable rows).
"""

import argparse
import collections

import networkx as nx
import pandas as pd

EARLY = ['promoter_marks', 'promoter_contacts', 'step1_marks']
DATA_SETS = ['SPRITE', 'HiPore-C', 'Pore-C']


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('output')
    parser.add_argument('--comparison', nargs='+', required=True)
    parser.add_argument('--q', type=float, default=0.05)
    parser.add_argument('--baseMean', type=float, default=50)
    parser.add_argument('--maxSE', type=float, default=1)
    args = parser.parse_args()
    rows = []
    for item in args.comparison:
        label, paths = item.split('=', 1)
        prefix, deseq2, *multiway = paths.split(',')
        levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
        levels['explains'] = levels['explains_declared_mean'].astype(str) == 'True'
        early = set(levels[levels.explains & levels.level.isin(EARLY)].node)
        expression = pd.read_csv(deseq2, sep='\t').dropna(subset=['padj']).drop_duplicates('gene_name')
        expression = expression.set_index('gene_name')
        graph = nx.read_graphml(prefix + '_union.graphml')
        triplets = pd.read_csv(multiway[0] + '_triplets.tsv', sep='\t') if multiway else None
        for node, sub in levels[levels.level.isin(['step2', 'step3'])].groupby('node', sort=False):
            if node in early:
                continue
            hit = sub[sub.explains & (sub.q_declared < args.q) &
                      (sub.beyond_mirror_declared.astype(str) == 'True')].sort_values('level')
            if not len(hit):
                continue
            r = hit.iloc[0]
            genes = [g for g in str(r.genes).split(',') if g in expression.index]
            if not genes:
                continue
            gene = min(genes, key=lambda g: abs(expression.loc[g, 'log2FoldChange'] - r.rna_log2FC))
            e = expression.loc[gene]
            if not (e.baseMean >= args.baseMean and e.lfcSE < args.maxSE):
                continue
            status = collections.Counter(graph.edges[node, m].get('status') for m in graph.neighbors(node))
            row = dict(comparison=label, gene=gene, node=node, log2FC=e.log2FoldChange, padj=e.padj,
                       baseMean=e.baseMean, step=r.level.replace('step', ''),
                       marks=str(r.concordant_marks_declared), concordant=int(r.n_concordant_nodes_declared),
                       tested=int(r.nodes_tested), P=r.p_declared, q=r.q_declared,
                       direct='{}/{}/{}'.format(status['shared'], status['gained'], status['lost']))
            if triplets is not None:
                t = triplets[(triplets.promoter == node) & (triplets.set == 'explaining') & triplets.present]
                row['present_paths'] = t.drop_duplicates(['mediator', 'target']).shape[0]
                for d in DATA_SETS:
                    row[d.replace('-', '')] = int(t[(t.data_set == d) & (t.p_empirical < 0.05)].shape[0])
            rows.append(row)
    # one row per gene: the promoter fragment with the lowest q
    table = pd.DataFrame(rows).sort_values(['comparison', 'q', 'P']).drop_duplicates(['comparison', 'gene'])
    table.to_csv(args.output + '.tsv', sep='\t', index=False)
    with open(args.output + '_rows.tex', 'w') as out:
        for r in table.itertuples():
            multi = '{} & {} & {} & {}'.format(int(r.present_paths), int(r.SPRITE), int(r.HiPoreC), int(r.PoreC)) \
                if 'present_paths' in table and pd.notna(getattr(r, 'present_paths', float('nan'))) else '-- & -- & -- & --'
            out.write('{} & {} & {:+.2f} & {:.1e} & {} & {} & {}/{} & {:.1e} & {:.2g} & {} & {} \\\\\n'.format(
                r.comparison, r.gene.replace('_', r'\_'), r.log2FC, r.padj, r.step, r.marks.replace(',', ', '),
                r.concordant, r.tested, r.P, r.q, r.direct, multi))
    print(table.groupby('comparison').size().to_string())


if __name__ == '__main__':
    main()
