"""Differential expression from STAR gene counts (ReadsPerGene.out.tab of rnaseq_coverage.sh), one run per replicate.

samples.tsv (of rnaseq_coverage.sh) gives the condition of every row; row k of the table is the run
<counts dir>/<condition>.<k>.ReadsPerGene.out.tab. The unstranded counts (column 2) of the gene IDs are named with the
gene_name of a GENCODE GTF and summed per name, genes with fewer than 10 counts in total are removed, and PyDESeq2
tests condition B against condition A (log2FoldChange is B relative to A).

Usage: python star_counts_deseq2.py <samples.tsv> <counts dir> <genes.gtf.gz> <condition A> <condition B> <output.tsv>
"""

import gzip
import os
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def gene_names(gtf):
    names = {}
    with gzip.open(gtf, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] == 'gene':
                names[f[8].split('gene_id "', 1)[1].split('"', 1)[0]] = \
                    f[8].split('gene_name "', 1)[1].split('"', 1)[0]
    return names


def main():
    samples, directory, gtf, a, b, output = sys.argv[1:7]
    table = pd.read_csv(samples, sep='\t')
    counts, conditions = {}, {}
    for row, condition in enumerate(table['condition'], start=1):
        if condition not in (a, b):
            continue
        name = '{}.{}'.format(condition, row)
        tab = pd.read_csv(os.path.join(directory, name + '.ReadsPerGene.out.tab'), sep='\t', header=None,
                          index_col=0)
        counts[name] = tab.loc[~tab.index.str.startswith('N_'), 1]
        conditions[name] = condition
    counts = pd.DataFrame(counts).fillna(0).astype(int)
    names = gene_names(gtf)
    counts = counts.groupby(counts.index.map(lambda g: names.get(g, g))).sum()
    counts = counts[counts.sum(axis=1) >= 10]
    metadata = pd.DataFrame({'condition': [conditions[c] for c in counts.columns]}, index=counts.columns)
    dds = DeseqDataSet(counts=counts.T, metadata=metadata, design='~condition', quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', b, a], quiet=True)
    stats.summary()
    result = stats.results_df.reset_index().rename(columns={'index': 'gene_name'})
    result.columns = ['gene_name'] + list(result.columns[1:])
    result.to_csv(output, sep='\t', index=False)
    significant = result[(result['padj'] <= 0.05) & (result['log2FoldChange'].abs() >= 1)]
    print('{} genes tested, {} with padj <= 0.05 and |log2FC| >= 1 ({} up in {})'.format(
        len(result), len(significant), int((significant['log2FoldChange'] > 0).sum()), b))


if __name__ == '__main__':
    main()
