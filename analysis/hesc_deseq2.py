"""Differential expression of H9 hESC-derived neuroectoderm (hNEC) versus hESC (Freire-Pritchett et al. 2017,
GSE86821 raw counts, two replicates each).

Usage: python hesc_deseq2.py <GSE86821_RNA-seq_rawcounts_hESChNPC.txt.gz> <output.tsv>

log2FoldChange is hNEC relative to hESC; genes are keyed by Ensembl id (column gene_id) with the symbol in
gene_name.
"""

import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def main():
    counts_path, output = sys.argv[1], sys.argv[2]
    table = pd.read_csv(counts_path, sep='\t', compression='infer')
    samples = [c for c in table.columns if c.startswith('RNA_Seq_')]
    table = table[table['ID'].astype(str).str.startswith('ENSG')]
    genes = table.groupby('ID')[samples].sum().round().astype(int)
    genes = genes[genes.sum(axis=1) >= 10]
    names = table.drop_duplicates('ID').set_index('ID')['Probe']
    metadata = pd.DataFrame({'condition': ['hNEC' if '_hNEC_' in s else 'hESC' for s in samples]}, index=samples)
    dds = DeseqDataSet(counts=genes.T, metadata=metadata, design='~condition', quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', 'hNEC', 'hESC'], quiet=True)
    stats.summary()
    result = stats.results_df
    result.insert(0, 'gene_name', names.reindex(result.index).to_numpy())
    result.index.name = 'gene_id'
    result.reset_index().to_csv(output, sep='\t', index=False)
    significant = result[(result['padj'] <= 0.05) & (result['log2FoldChange'].abs() >= 1)]
    print('{} genes tested, {} with padj <= 0.05 and |log2FC| >= 1 ({} up in hNEC)'.format(
        len(result), len(significant), int((significant['log2FoldChange'] > 0).sum())))


if __name__ == '__main__':
    main()
