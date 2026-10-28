"""Gene-level PyDESeq2 of Kasumi-1 RUNX1-ETO knockdown (siRE) versus mismatch control (siMM) from salmon output.

Transcript read estimates (NumReads) are summed per gene (Ensembl id without version, from the GENCODE transcript
names) and rounded. log2FoldChange is siRE relative to siMM.

Usage: python kasumi_deseq2.py <samples.tsv> <salmon output dir> <output.tsv> <gencode transcripts.fa.gz>
The transcript-to-gene map is read from the GENCODE FASTA headers (salmon --gencode keeps only the transcript id).
"""

import gzip
import os
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def transcript_genes(fasta):
    gene_of, names = {}, {}
    with gzip.open(fasta, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                fields = line[1:].strip().split('|')
                gene = fields[1].split('.')[0]
                gene_of[fields[0]] = gene
                names[gene] = fields[5]
    return gene_of, names


def main():
    samples_path, quant_dir, output, fasta = sys.argv[1:5]
    samples = pd.read_csv(samples_path, sep='\t')
    gene_of, names = transcript_genes(fasta)
    counts = {}
    for run in samples['run']:
        quant = pd.read_csv(os.path.join(quant_dir, run, 'quant.sf'), sep='\t')
        gene = quant['Name'].map(gene_of)
        counts[run] = quant.groupby(gene)['NumReads'].sum()
    genes = pd.DataFrame(counts).fillna(0).round().astype(int)
    genes = genes[genes.sum(axis=1) >= 10]
    metadata = pd.DataFrame({'condition': samples['condition'].to_numpy()}, index=samples['run'].to_numpy())
    dds = DeseqDataSet(counts=genes.T, metadata=metadata, design='~condition', quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', 'siRE', 'siMM'], quiet=True)
    stats.summary()
    result = stats.results_df
    result.insert(0, 'gene_name', [names.get(g, g) for g in result.index])
    result.index.name = 'gene_id'
    result.reset_index().to_csv(output, sep='\t', index=False)
    significant = result[(result['padj'] <= 0.05) & (result['log2FoldChange'].abs() >= 1)]
    print('{} genes tested, {} with padj <= 0.05 and |log2FC| >= 1 ({} up after RUNX1-ETO knockdown)'.format(
        len(result), len(significant), int((significant['log2FoldChange'] > 0).sum())))


if __name__ == '__main__':
    main()
