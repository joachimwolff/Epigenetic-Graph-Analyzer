"""Expression table of unchanged genes posing as differentially expressed: the genes that pass --padj and
--minLog2FC in the real table are set to unchanged (padj 1, log2 fold change 0), and genes that do not change
(padj >= --unchangedPadj and |log2 fold change| <= --unchangedLog2FC) receive real fold changes, drawn at random
without replacement, with padj 1e-10; all unchanged genes are used when there are fewer of them than changed genes.
The run of the explanation levels on this table measures how often EGA offers an explanation for a gene whose
expression does not change (calibration of the per-gene explanations, 42_calibration.sh).

Usage: python pseudo_de.py <deseq2.tsv> <out.tsv> [--seed 0] [--padj 0.05] [--minLog2FC 1] [--unchangedPadj 0.2]
           [--unchangedLog2FC 0.5]
"""

import argparse

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('deseq2')
    parser.add_argument('output')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--padj', type=float, default=0.05)
    parser.add_argument('--minLog2FC', type=float, default=1.0)
    parser.add_argument('--unchangedPadj', type=float, default=0.2)
    parser.add_argument('--unchangedLog2FC', type=float, default=0.5)
    args = parser.parse_args()
    table = pd.read_csv(args.deseq2, sep='\t')
    changed = (table.padj < args.padj) & (table.log2FoldChange.abs() >= args.minLog2FC)
    unchanged = np.flatnonzero(((table.padj >= args.unchangedPadj) &
                                (table.log2FoldChange.abs() <= args.unchangedLog2FC)).to_numpy())
    rng = np.random.default_rng(args.seed)
    fold_changes = table.loc[changed, 'log2FoldChange'].to_numpy()
    size = min(len(unchanged), len(fold_changes))
    chosen = rng.choice(unchanged, size=size, replace=False)
    fold_changes = rng.choice(fold_changes, size=size, replace=False)
    out = table.copy()
    out.loc[changed, ['log2FoldChange', 'padj']] = [0.0, 1.0]
    out.loc[out.index[chosen], 'log2FoldChange'] = rng.permutation(fold_changes)
    out.loc[out.index[chosen], 'padj'] = 1e-10
    out.to_csv(args.output, sep='\t', index=False)
    print('{} genes with real fold changes, {} unchanged genes chosen from {}'.format(
        int(changed.sum()), len(chosen), len(unchanged)))


if __name__ == '__main__':
    main()
