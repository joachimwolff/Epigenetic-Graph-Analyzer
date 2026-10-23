"""Differential expression between two conditions from ENCODE RSEM gene quantifications.

files.tsv lists one quantification per row (columns condition, filename, ...); expected counts are rounded, gene IDs
(ENSG, version removed) are named with the gene_name of a GENCODE GTF and summed per name, genes with fewer than 10
counts in total are removed, and PyDESeq2 tests condition B against condition A (log2FoldChange is B relative to A).

Usage: python encode_deseq2.py <files.tsv> <genes.gtf.gz> <condition A> <condition B> <output.tsv>
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
            if f[2] != 'gene':
                continue
            gene_id = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split('.')[0]
            names[gene_id] = f[8].split('gene_name "', 1)[1].split('"', 1)[0]
    return names


def main():
    files, gtf, a, b, output = sys.argv[1:6]
    samples = pd.read_csv(files, sep='\t')
    samples = samples[samples['condition'].isin([a, b])]
    counts = {}
    for row in samples.itertuples():
        table = pd.read_csv(os.path.join(os.path.dirname(files), row.filename), sep='\t',
                            usecols=['gene_id', 'expected_count'])
        table = table[table['gene_id'].str.startswith('ENSG')]
        counts[row.file_accession] = table.set_index(table['gene_id'].str.split('.').str[0])['expected_count']
    counts = pd.DataFrame(counts).fillna(0).round().astype(int)
    names = gene_names(gtf)
    counts = counts[counts.index.isin(names)]
    counts = counts.groupby(counts.index.map(names)).sum()
    counts = counts[counts.sum(axis=1) >= 10]
    metadata = pd.DataFrame({'condition': list(samples['condition'])}, index=list(samples['file_accession']))
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
