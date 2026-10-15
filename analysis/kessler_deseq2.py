"""Differential expression of E10.5 PA2 versus Md cranial neural crest (Minoux et al. 2017, GSE89434).

Usage: python kessler_deseq2.py <gene_counts.csv.gz> <output.tsv>

Columns that occur twice in the count table hold identical counts of the same
library and are used once. log2FoldChange is PA2 relative to Md.
"""

import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def main():
    counts_path, output = sys.argv[1], sys.argv[2]
    table = pd.read_csv(counts_path, compression='infer')
    table = table.loc[:, ~table.columns.duplicated()]
    samples = [c for c in table.columns if c.startswith('RNAseq_E10.5-Md_') or c.startswith('RNAseq_E10.5-PA2_')]
    genes = table.groupby('Associated.Gene.Name')[samples].sum()
    genes = genes[genes.sum(axis=1) >= 10]
    metadata = pd.DataFrame({'condition': ['PA2' if '-PA2_' in s else 'Md' for s in samples]}, index=samples)
    dds = DeseqDataSet(counts=genes.T.astype(int), metadata=metadata, design='~condition', quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', 'PA2', 'Md'], quiet=True)
    stats.summary()
    result = stats.results_df.reset_index().rename(columns={'index': 'gene_name', 'Associated.Gene.Name': 'gene_name'})
    result.columns = ['gene_name'] + list(result.columns[1:])
    result.to_csv(output, sep='\t', index=False)
    significant = result[(result['padj'] <= 0.05) & (result['log2FoldChange'].abs() >= 1)]
    print('{} genes tested, {} with padj <= 0.05 and |log2FC| >= 1 ({} up in PA2)'.format(
        len(result), len(significant), int((significant['log2FoldChange'] > 0).sum())))


if __name__ == '__main__':
    main()
