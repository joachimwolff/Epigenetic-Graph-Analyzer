"""DNA methylation of one condition as an annotation of the element classes (blind test, analysis/21_blind_roles.sh).

Roadmap has WGBS for E008 (H9 hESC) only, not for E009 (hNEC), so methylation cannot enter the comparison as a
change; it is read as a state of condition A. Per fragment: the mean fractional CpG methylation over the covered
CpGs (pyBigWig mean over covered bases; fragments without a covered CpG are NA). Outputs (TSV):

  <out>_classes.tsv    methylation of the fragments of every element class in condition A (declared and learned
                       roles, fragment mean presence): fragments, with a covered CpG, median, quartiles, fraction
                       >= 0.6 (methylated) and <= 0.2 (unmethylated)
  <out>_promoters.tsv  the baited promoters of genes with an expression change: Spearman rho of the promoter
                       methylation in A with the expression change B / A (all, and per element class of the
                       promoter in A), and the mean expression change of unmethylated (<= 0.2), intermediate and
                       methylated (>= 0.6) promoters with chromosome block-bootstrap intervals

Known biology: promoter CpG methylation is repressive; a promoter methylated in hESC is expected to stay silent or
be activated only with demethylation, so its expression change is not predicted by the state alone.

Usage: python blind_methylation.py <egaDifferential prefix> <methylation bigWig> <condition A> <output prefix>
"""

import sys

import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats

sys.path.insert(0, __file__.rsplit('/', 2)[0])
from epigeneticgraphanalyzer.validation import block_bootstrap  # noqa: E402


def fragment_methylation(nodes, path):
    handle = pyBigWig.open(path)
    chroms = handle.chroms()
    values = []
    for node in nodes:
        chrom, rest = node.split(':')
        start, end = (int(x) for x in rest.split('-'))
        if chrom not in chroms:
            values.append(np.nan)
            continue
        value = handle.stats(chrom, max(0, start), min(end, chroms[chrom]), type='mean')[0]
        values.append(np.nan if value is None else float(value))
    handle.close()
    return np.array(values, dtype=float)


def main():
    prefix, bigwig, condition, out = sys.argv[1:5]
    elements = pd.read_csv(prefix + '_elements.tsv', sep='\t')
    a = elements[elements['condition'].astype(str) == condition].reset_index(drop=True)
    a['methylation'] = fragment_methylation(a['node'], bigwig)
    rows = []
    for roles in ('declared', 'learned'):
        column = 'class_{}_mean'.format(roles)
        for cls, part in a.groupby(column):
            m = part['methylation'].dropna()
            rows.append(dict(roles=roles, element_class=cls, fragments=len(part), with_cpg=len(m),
                             median=m.median(), q25=m.quantile(0.25), q75=m.quantile(0.75),
                             fraction_methylated=(m >= 0.6).mean() if len(m) else np.nan,
                             fraction_unmethylated=(m <= 0.2).mean() if len(m) else np.nan))
    pd.DataFrame(rows).to_csv(out + '_classes.tsv', sep='\t', index=False)
    features = pd.read_csv(prefix + '_gene_features.tsv', sep='\t', usecols=['node', 'rna_log2FC'])
    genes = features.merge(a[['node', 'methylation', 'class_declared_mean']], on='node').dropna()
    chroms = genes['node'].str.split(':').str[0].to_numpy()
    rows = []
    for cls, part in [('all', genes)] + list(genes.groupby('class_declared_mean')):
        if len(part) < 30:
            continue
        rho, p = stats.spearmanr(part['methylation'], part['rna_log2FC'])
        rows.append(dict(kind='spearman', promoter_class=cls, group='', n_genes=len(part), rho=rho, p=p))
    for group, mask in (('unmethylated <= 0.2', genes['methylation'] <= 0.2),
                        ('intermediate', (genes['methylation'] > 0.2) & (genes['methylation'] < 0.6)),
                        ('methylated >= 0.6', genes['methylation'] >= 0.6)):
        values = genes.loc[mask, 'rna_log2FC'].to_numpy()
        low, high = block_bootstrap(values, chroms[mask.to_numpy()], 1000, 0) if len(values) else (np.nan, np.nan)
        rows.append(dict(kind='group', promoter_class='all', group=group, n_genes=int(mask.sum()),
                         mean_rna_log2FC=values.mean() if len(values) else np.nan, ci_low=low, ci_high=high,
                         fraction_up=(values >= 1).mean() if len(values) else np.nan,
                         fraction_down=(values <= -1).mean() if len(values) else np.nan))
    table = pd.DataFrame(rows)
    table.to_csv(out + '_promoters.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 200):
        print(pd.DataFrame(pd.read_csv(out + '_classes.tsv', sep='\t')).round(3).to_string(index=False))
        print(table.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
