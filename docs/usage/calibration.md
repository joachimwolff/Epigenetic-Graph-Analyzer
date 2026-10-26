# Calibration with unchanged genes

A level offers an explanation when a mark changes concordantly with the expression change of a gene. Where the
chromatin of two conditions differs at many loci, such concordance also arises for genes whose expression does not
change. The calibration measures this rate for a comparison: the comparison is repeated with unchanged genes that
receive the fold changes of the differentially expressed genes, with all other inputs and settings unchanged. A level
offers explanations beyond chance where the rate of the differentially expressed genes exceeds that of the unchanged
genes.

In the comparisons of the manuscript, the marks at the promoter and two steps away exceeded the rate of unchanged
genes in all three comparisons; the third step stayed close to it.

## 1. Expression table of unchanged genes

```bash
for seed in 0 1 2; do
  python analysis/pseudo_de.py deseq2_B_vs_A.tsv calibration/unchanged_$seed/expression.tsv --seed $seed
done
```

The genes passing `--padj` (0.05) and `--minLog2FC` (1) are set to unchanged; genes with an adjusted $P$ of at least
`--unchangedPadj` (0.2) and an absolute log2 fold change of at most `--unchangedLog2FC` (0.5) receive the real fold
changes, drawn at random without replacement.

## 2. The same run on each table

Run `egaDifferential` with exactly the options of the real run and only `--expression` and `-o` changed:

```bash
for seed in 0 1 2; do
  egaDifferential <options of the real run> \
      --expression calibration/unchanged_$seed/expression.tsv -o calibration/unchanged_$seed/A_vs_B
done
```

## 3. Compare

```bash
python analysis/43_calibration_summary.py differential/A_vs_B calibration A_vs_B 'unchanged_*'
```

The table gives, for the real run and every unchanged run, the promoter fragments without an explanation at the
promoter or step 1 and the share of them with a candidate at step 2 or 3 (any, with $P<0.05$ beyond the mirrored
loci, with $q<0.05$). The share of genes per level is in `<prefix>_explanation_levels.tsv` of each run (column
`fraction_any`).

Report the rate of unchanged genes next to the candidates: a candidate explanation at a level whose rate is close to
that of unchanged genes is weak evidence on its own.
