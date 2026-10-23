"""Genes explained only through the graph, as candidates for a report example.

A candidate gene is differentially expressed (padj < --padj, |log2FC| >= --lfc, baseMean >= --baseMean, lfcSE <
--maxSE so that the change does not rest on zero counts in one condition), its promoter fragment is explained at no
promoter level and not at step 1, and step 2 or step 3 explains it with q < --q beyond the mirrored loci. For each
candidate the table gives the gene of the fragment that carries the expression change, its DESeq2 statistics, the
first explaining step and its statistics, the number of direct contacts and their status, and the strongly changed
genes whose bodies overlap the candidate gene (which could carry its reads or its marks).

Usage: python screen_report_examples.py <cmp prefix> <deseq2.tsv> <genes.gtf(.gz)> <output.tsv> [--baseMean 50]
           [--padj 0.05] [--lfc 1] [--maxSE 1] [--q 0.05] [--geneColumn gene_name]
"""

import argparse
import gzip

import numpy as np
import pandas as pd

EARLY = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state']


def gene_bodies(path):
    rows = []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] != 'gene':
                continue
            name = f[8].split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in f[8] else None
            if name:
                rows.append((name, f[0], int(f[3]), int(f[4])))
    return pd.DataFrame(rows, columns=['gene', 'chrom', 'start', 'end']).drop_duplicates('gene')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('prefix', 'deseq2', 'gtf', 'output'):
        parser.add_argument(name)
    parser.add_argument('--baseMean', type=float, default=50)
    parser.add_argument('--padj', type=float, default=0.05)
    parser.add_argument('--lfc', type=float, default=1)
    parser.add_argument('--maxSE', type=float, default=1)
    parser.add_argument('--q', type=float, default=0.05)
    parser.add_argument('--geneColumn', default='gene_name')
    args = parser.parse_args()
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    levels['explains'] = levels['explains_declared_mean'].astype(str) == 'True'
    wide = levels.pivot_table(index='node', columns='level', values='explains', aggfunc='first').fillna(False)
    deseq = pd.read_csv(args.deseq2, sep='\t').dropna(subset=['padj'])
    deseq = deseq.drop_duplicates(args.geneColumn).set_index(args.geneColumn)
    bodies = gene_bodies(args.gtf).set_index('gene')
    changed = deseq[(deseq['padj'] < args.padj) & (deseq['log2FoldChange'].abs() >= args.lfc)]
    rows = []
    for node, flags in wide.iterrows():
        if flags[EARLY].any():
            continue
        step = next((s for s in ('step2', 'step3') if flags.get(s, False)), None)
        if step is None:
            continue
        info = levels[(levels.node == node) & (levels.level == step)].iloc[0]
        if not (info['q_declared'] < args.q and str(info['beyond_mirror_declared']) == 'True'):
            continue
        # the gene of the fragment whose change was assigned to it
        genes = [g for g in str(info['genes']).split(',') if g in deseq.index]
        gene = min(genes, key=lambda g: abs(deseq.loc[g, 'log2FoldChange'] - info['rna_log2FC']), default=None)
        if gene is None:
            continue
        d = deseq.loc[gene]
        if not (d['baseMean'] >= args.baseMean and d['padj'] < args.padj and abs(d['log2FoldChange']) >= args.lfc
                and d['lfcSE'] < args.maxSE):
            continue
        overlap = []
        if gene in bodies.index:
            b = bodies.loc[gene]
            near = bodies[(bodies.chrom == b.chrom) & (bodies.start <= b.end) & (bodies.end >= b.start)]
            overlap = [g for g in near.index if g != gene and g in changed.index]
        contacts = levels[(levels.node == node) & (levels.level == 'promoter_contacts')].iloc[0]
        step1 = levels[(levels.node == node) & (levels.level == 'step1_marks')].iloc[0]
        rows.append(dict(node=node, gene=gene, baseMean=round(d['baseMean'], 1),
                         log2FC=round(d['log2FoldChange'], 2), lfcSE=round(d['lfcSE'], 2), padj=d['padj'],
                         first_step=step, p_step=info['p_declared'], q_step=info['q_declared'],
                         concordant_nodes=info['n_concordant_nodes_declared'], nodes_tested=info['nodes_tested'],
                         strongest=info['detail_declared_mean'], p_step1=step1['p_declared'],
                         direct=contacts['nodes_tested'], gained=contacts['contacts_gained'],
                         lost=contacts['contacts_lost'], shared=contacts['contacts_shared'],
                         overlapping_changed=','.join(overlap)))
    table = pd.DataFrame(rows)
    if len(table):
        table = table.sort_values(['q_step', 'padj'])
    table.to_csv(args.output, sep='\t', index=False)
    print('{} candidates'.format(len(table)))


if __name__ == '__main__':
    main()
