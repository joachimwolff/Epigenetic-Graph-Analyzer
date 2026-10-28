"""Gene expression of AML patients from RNA-seq coverage bigWigs (Assi et al. 2019, GSE108266).

Expression of a gene is the summed coverage over the union of its exons (GENCODE) divided by the total coverage
of the sample, per million. The table holds one column per sample, the mean log2 expression per group, the
log2 fold change of group B over group A and a Welch t-test across patients with BH adjustment.

Usage: python assi_expression.py <genes.gtf.gz> <output.tsv> A=<a1.bw,a2.bw,...> B=<b1.bw,...>
A sample given as plus.bw+minus.bw sums the coverage of both strands.
"""

import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats

from epigeneticgraphanalyzer.genome import open_text

GENE_TYPES = {'protein_coding', 'lincRNA', 'antisense'}


def exons(gtf):
    merged = defaultdict(list)
    names, chroms = {}, {}
    with open_text(gtf) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'exon':
                continue
            attributes = dict(item.strip().split(' ', 1) for item in fields[8].split(';') if ' ' in item.strip())
            if attributes.get('gene_type', '').strip('"') not in GENE_TYPES:
                continue
            gene = attributes['gene_id'].strip('"').split('.')[0]
            names[gene] = attributes['gene_name'].strip('"')
            chroms[gene] = fields[0]
            merged[gene].append((int(fields[3]) - 1, int(fields[4])))
    for gene, intervals in merged.items():
        intervals.sort()
        union = [list(intervals[0])]
        for start, end in intervals[1:]:
            if start <= union[-1][1]:
                union[-1][1] = max(union[-1][1], end)
            else:
                union.append([start, end])
        merged[gene] = union
    return merged, names, chroms


def coverage(path, merged, chroms):
    """Raw summed coverage per gene and the total coverage of the track."""
    handle = pyBigWig.open(path)
    sizes = handle.chroms()
    total = sum(float(handle.stats(c, 0, n, type='sum', exact=True)[0] or 0) for c, n in sizes.items())
    values = {}
    for gene, intervals in merged.items():
        chrom = chroms[gene]
        if chrom not in sizes:
            continue
        values[gene] = sum(float(handle.stats(chrom, s, min(e, sizes[chrom]), type='sum', exact=True)[0] or 0)
                           for s, e in intervals if s < min(e, sizes[chrom]))
    handle.close()
    return pd.Series(values), total


def main():
    gtf, output = sys.argv[1:3]
    groups = dict(arg.split('=', 1) for arg in sys.argv[3:])
    merged, names, chroms = exons(gtf)
    columns = {}
    for group, paths in groups.items():
        for sample in paths.split(','):
            parts = sample.split('+')
            tracks = [coverage(p, merged, chroms) for p in parts]
            values = sum(t[0] for t in tracks) / sum(t[1] for t in tracks) * 1e6
            columns['{}:{}'.format(group, parts[0].rsplit('/', 1)[-1].split('.')[0])] = values
    table = pd.DataFrame(columns).fillna(0)
    log_values = np.log2(table + 1)
    a = log_values[[c for c in table if c.startswith('A:')]]
    b = log_values[[c for c in table if c.startswith('B:')]]
    table['mean_log2_A'] = a.mean(axis=1)
    table['mean_log2_B'] = b.mean(axis=1)
    table['log2FoldChange'] = table['mean_log2_B'] - table['mean_log2_A']
    expressed = (table[[c for c in columns]] >= 1).any(axis=1)
    pvalue = stats.ttest_ind(b, a, axis=1, equal_var=False).pvalue
    table['pvalue'] = np.where(np.isfinite(pvalue), pvalue, 1.0)
    table = table[expressed]
    order = np.argsort(table['pvalue'].to_numpy())
    ranked = table['pvalue'].to_numpy()[order] * len(table) / np.arange(1, len(table) + 1)
    padj = np.empty(len(table))
    padj[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    table['padj'] = padj
    table.insert(0, 'gene_name', [names[g] for g in table.index])
    table.index.name = 'gene_id'
    table.reset_index().to_csv(output, sep='\t', index=False)
    print('{} expressed genes, {} with padj <= 0.05 and |log2FC| >= 1'.format(
        len(table), int(((table['padj'] <= 0.05) & (table['log2FoldChange'].abs() >= 1)).sum())))


if __name__ == '__main__':
    main()
