"""Differential expression over keratinocyte differentiation (Rubin et al. 2017, GSE84661 RSEM expected counts, two
donors at day 0, 3 and 6): one PyDESeq2 model ~ donor + day, contrasts day 3 vs 0, day 6 vs 3 and day 6 vs 0.

Usage: python rubin_deseq2.py <GSE84661_GV19_TSS_RNAseq_RSEM_counts.bed.gz> <out_dir>

Outputs: <out_dir>/deseq2_D3_vs_D0.tsv, deseq2_D6_vs_D3.tsv, deseq2_D6_vs_D0.tsv (gene_name first, gene_id,
baseMean, log2FoldChange = later day relative to earlier day, padj), normalised_counts.tsv (DESeq2 size-factor
normalised counts per sample) and deseq2_summary.tsv. Expected counts are rounded; Ensembl ids lose the version.
"""

import os
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

SAMPLES = ['B1_D0', 'B2_D0', 'B1_D3', 'B2_D3', 'B1_D6', 'B2_D6']
CONTRASTS = [('D3', 'D0'), ('D6', 'D3'), ('D6', 'D0')]


def main():
    counts_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    table = pd.read_csv(counts_path, sep='\t', comment='#', header=None,
                        names=['chrom', 'tss', 'tss2', 'gene_id', 'gene_name', 'strand', 'tts'] + SAMPLES)
    table['gene_id'] = table['gene_id'].str.replace(r'\.\d+$', '', regex=True)
    names = table.drop_duplicates('gene_id').set_index('gene_id')['gene_name']
    genes = table.groupby('gene_id')[SAMPLES].sum().round().astype(int)
    genes = genes[genes.sum(axis=1) >= 10]
    metadata = pd.DataFrame({'donor': [s.split('_')[0] for s in SAMPLES], 'day': [s.split('_')[1] for s in SAMPLES]},
                            index=SAMPLES)
    dds = DeseqDataSet(counts=genes.T, metadata=metadata, design='~donor + day', quiet=True)
    dds.deseq2()
    normalised = pd.DataFrame(dds.layers['normed_counts'].T, index=genes.index, columns=SAMPLES)
    normalised.insert(0, 'gene_name', names.reindex(normalised.index).to_numpy())
    normalised.index.name = 'gene_id'
    normalised.reset_index().to_csv(os.path.join(out_dir, 'normalised_counts.tsv'), sep='\t', index=False)
    rows = []
    for later, earlier in CONTRASTS:
        stats = DeseqStats(dds, contrast=['day', later, earlier], quiet=True)
        stats.summary()
        result = stats.results_df.copy()
        result.insert(0, 'gene_id', result.index)
        result.insert(0, 'gene_name', names.reindex(result.index).to_numpy())
        name = 'deseq2_{}_vs_{}.tsv'.format(later, earlier)
        result.to_csv(os.path.join(out_dir, name), sep='\t', index=False)
        significant = result[(result['padj'] <= 0.05) & (result['log2FoldChange'].abs() >= 1)]
        rows.append({'contrast': '{}_vs_{}'.format(later, earlier), 'tested': len(result),
                     'differential': len(significant), 'up': int((significant['log2FoldChange'] > 0).sum()),
                     'down': int((significant['log2FoldChange'] < 0).sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(out_dir, 'deseq2_summary.tsv'), sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
