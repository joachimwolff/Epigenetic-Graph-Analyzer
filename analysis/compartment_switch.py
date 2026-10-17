"""Prototype: condition-specific compartment switches from chromatin tracks as an explanation level.

For two conditions with the same marks the binned log enrichment of every track (100 kb) is z-scored jointly over
both conditions, the first principal component is fitted on the stacked bins of both conditions (one set of
loadings, oriented by the active-role marks) and applied to each; scores are smoothed over 3 bins, centred on the
per-chromosome median of the mean of both conditions and divided by the pooled per-chromosome standard deviation.
Two normalisations: 'joint' (as above) and 'quantile' (main): per chromosome both scores are mapped onto their mean
quantile distribution, which removes differences in location and spread between the two data sets (ChIP efficiency,
sequencing depth), so a switch is a change of rank along the chromosome.
A bin switches (B->A or A->B) when the sign changes and |score_B - score_A| >= threshold (default 0.5 SD).

Genes (TSS of the GTF, median TSS per gene name) are assigned to bins; for differential genes (padj < alpha) a switch
is concordant when B->A goes with up-regulation in condition B (A->B with down-regulation). Reference: differential
genes in non-switching bins, matched by baseMean decile, scored by the sign of their delta (concordant when the sign
of delta and log2FC agree), all of them and only those with |delta| >= threshold without a sign change, and the Spearman correlation of delta with log2FC over all genes.
Control comparison without biology expected: E003 (H1) vs E008 (H9), two hESC lines (switch rate under noise).

Consistency check without Hi-C: per condition with an accessibility track (H9 DNase E008 for hESC, E11.5 facial
prominence ATAC as proxy for Md; H1 E003 and H9 E008 DNase for the control) the chromatin score is compared with the
binned accessibility level (log enrichment, 3-bin smoothing, centred per chromosome): Pearson, sign agreement, AUC
(accessibility_consistency.tsv). hNEC (E009) and PA2 have no accessibility data.

Output per comparison and normalisation (<results>/compartments_domains/switch/<comparison>_<normalisation>_genes.tsv), readable as one explanation
level: gene, node (promoter fragment from the hops_n tables where available), chrom, tss, compartment_A,
compartment_B, score_A, score_B, delta, switch, rna_log2FC, padj, differential, concordant, in_cascade.
Plus <comparison>_<normalisation>_bins.tsv.gz and switch_summary.tsv.

Usage: python compartment_switch.py [threshold] [alpha]
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402
from epigeneticgraphanalyzer.compartments import (bin_bigwig, chromatin_pc1, compartment_switches,  # noqa: E402
                                                  compare_to_reference, log_enrichment, make_bins,
                                                  quantile_normalise_pair,
                                                  smooth_per_chrom)
from epigeneticgraphanalyzer.io import read_tss_from_gtf  # noqa: E402

ROAD = os.path.join(cc.RAW, 'roadmap')
FOLD = os.path.join(cc.RESULTS, 'minoux_fold')


def roadmap(eid, marks):
    return {m: os.path.join(ROAD, '{}-{}.fc.signal.bigwig'.format(eid, m)) for m in marks}


MARKS3 = ('H3K4me1', 'H3K4me3', 'H3K27me3')
COMPARISONS = [
    # name, names (A, B), tracks A, tracks B, assembly, gtf, expression table (log2FC of B over A), gene types
    ('hESC_vs_hNEC', ('hESC', 'hNEC'), roadmap('E008', MARKS3), roadmap('E009', MARKS3), 'hg19',
     os.path.join(cc.RAW, 'annotation', 'gencode.v19.annotation.gtf.gz'),
     os.path.join(cc.RESULTS, 'differential', 'deseq2_hNEC_vs_hESC.tsv'), ('protein_coding', 'lincRNA', 'antisense')),
    ('Md_vs_PA2', ('Md', 'PA2'),
     {m: os.path.join(FOLD, 'Md_{}.fold.bw'.format(m)) for m in ('H3K27ac', 'H3K4me2', 'H3K27me3')},
     {m: os.path.join(FOLD, 'PA2_{}.fold.bw'.format(m)) for m in ('H3K27ac', 'H3K4me2', 'H3K27me3')}, 'mm10',
     os.path.join(cc.RAW, 'kessler', 'mm10.ncbiRefSeq.gtf'),
     os.path.join(cc.RESULTS, 'kessler', 'deseq2_PA2_vs_Md.tsv'), None),
    ('E003_vs_E008_control', ('H1', 'H9'), roadmap('E003', MARKS3), roadmap('E008', MARKS3), 'hg19',
     os.path.join(cc.RAW, 'annotation', 'gencode.v19.annotation.gtf.gz'), None, None),
]


ACC = os.path.join(cc.RAW, 'accessibility')
# accessibility of a condition (consistency check where no Hi-C exists); hNEC (E009) and PA2 have none
ACCESSIBILITY = {
    ('hESC_vs_hNEC', 'A'): ('H9 DNase (Roadmap E008)', os.path.join(ACC, 'H9_DNase-seq.hg19.roadmap_E008.bigWig')),
    ('Md_vs_PA2', 'A'): ('E11.5 facial prominence ATAC (ENCODE, proxy for Md)',
                         os.path.join(ACC, 'mouseE11.5FacialProminence_ATAC-seq.mm10.ENCFF788BFR.bigWig')),
    ('E003_vs_E008_control', 'A'): ('H1 DNase (Roadmap E003)', os.path.join(ACC, 'H1_DNase-seq.hg19.roadmap_E003.bigWig')),
    ('E003_vs_E008_control', 'B'): ('H9 DNase (Roadmap E008)', os.path.join(ACC, 'H9_DNase-seq.hg19.roadmap_E008.bigWig')),
}


def bins_for(assembly, tracks, bin_size=100000):
    if assembly in ('hg19', 'hg38'):
        return cc.bins_for(assembly, bin_size)
    import pyBigWig
    sizes = pyBigWig.open(next(iter(tracks.values()))).chroms()
    keep = ['chr{}'.format(c) for c in list(range(1, 20)) + ['X']]
    return make_bins({c: sizes[c] for c in keep if c in sizes}, bin_size)


def joint_scores(bins, tracks_a, tracks_b):
    marks = list(tracks_a)
    fa = np.column_stack([log_enrichment(bin_bigwig(tracks_a[m], bins)) for m in marks])
    fb = np.column_stack([log_enrichment(bin_bigwig(tracks_b[m], bins)) for m in marks])
    stacked = np.vstack([fa, fb])
    active = np.array([cc.ROLES.get(m, 'other') in ('active', 'enhancer', 'transcription') for m in marks])
    score, loadings, explained = chromatin_pc1(stacked, active=active)
    n = len(bins)
    a, b = score[:n], score[n:]
    chroms = bins.chrom.values
    a, b = smooth_per_chrom(a, chroms, 1), smooth_per_chrom(b, chroms, 1)
    for chrom in pd.unique(chroms):
        k = chroms == chrom
        mean = (a[k] + b[k]) / 2
        ok = np.isfinite(mean)
        if ok.sum() < 3:
            a[k], b[k] = np.nan, np.nan
            continue
        centre = np.median(mean[ok])
        sd = np.nanstd(np.concatenate([a[k], b[k]]))
        a[k] = (a[k] - centre) / sd
        b[k] = (b[k] - centre) / sd
    return a, b, dict(zip(marks, loadings)), explained


def gene_tss(gtf, gene_types):
    tss = read_tss_from_gtf(gtf, gene_types=gene_types)
    out = tss.groupby('gene_name').agg(chrom=('chrom', 'first'), tss=('tss', 'median')).reset_index()
    out['tss'] = out.tss.round().astype(np.int64)
    return out


def matched_concordance(sw, pool, rng, n_draws=200):
    """Concordance (sign of delta vs sign of log2FC) of genes drawn from `pool` with the baseMean deciles of `sw`."""
    values = []
    for _ in range(n_draws):
        picks = []
        for d, n in sw.decile.value_counts().items():
            cand = pool[pool.decile == d]
            if len(cand):
                picks.append(cand.sample(n=n, replace=len(cand) < n, random_state=rng.integers(1 << 31)))
        if picks:
            p = pd.concat(picks)
            values.append(np.mean(np.sign(p.delta) == np.sign(p.rna_log2FC)))
    return np.array(values)


def evaluate(name, na, nb, bins, a, b, threshold, normalisation, loadings, explained, genes_base, node_of, alpha,
             out_dir, rng):
    delta = b - a
    switch = compartment_switches(a, b, threshold)
    tag = '{}_{}'.format(name, normalisation)
    bins.assign(score_A=a, score_B=b, delta=delta, switch=switch).to_csv(
        os.path.join(out_dir, '{}_bins.tsv.gz'.format(tag)), sep='\t', index=False, float_format='%.4g')
    valid = switch != 'NA'
    is_sw = np.isin(switch, ['B->A', 'A->B'])
    row = {'comparison': name, 'normalisation': normalisation, 'condition_A': na, 'condition_B': nb,
           'threshold_sd': threshold, 'loadings': ';'.join('{}:{:+.2f}'.format(k, v) for k, v in loadings.items()),
           'explained_pc1': explained, 'bins': int(valid.sum()),
           'B_to_A_bins': int(np.sum(switch == 'B->A')), 'A_to_B_bins': int(np.sum(switch == 'A->B')),
           'switch_fraction': float(np.mean(is_sw[valid])),
           'score_correlation_A_B': float(pd.Series(a).corr(pd.Series(b)))}
    for t in (0.25, 0.5, 1.0):
        s = compartment_switches(a, b, t)
        row['switch_fraction_sd{}'.format(t)] = float(np.mean(np.isin(s[valid], ['B->A', 'A->B'])))
    if genes_base is None:
        return row
    i = genes_base.i.values
    genes = genes_base.assign(score_A=a[i], score_B=b[i], delta=delta[i], switch=switch[i])
    genes['compartment_A'] = np.where(genes.score_A > 0, 'A', 'B')
    genes['compartment_B'] = np.where(genes.score_B > 0, 'A', 'B')
    genes.loc[~np.isfinite(genes.score_A), ['compartment_A', 'compartment_B']] = 'NA'
    genes['rna_log2FC'] = genes.log2FoldChange
    genes['differential'] = genes.padj < alpha
    gsw = genes.switch.isin(['B->A', 'A->B'])
    conc = np.where(gsw & genes.differential,
                    np.where(genes.switch == 'B->A', genes.rna_log2FC > 0, genes.rna_log2FC < 0), np.nan)
    genes['concordant'] = pd.Series(conc, index=genes.index).map({1.0: True, 0.0: False})
    genes['node'] = genes.gene_name.map(node_of) if node_of is not None else np.nan
    genes['in_cascade'] = genes.node.notna()
    genes = genes.rename(columns={'gene_name': 'gene'})
    cols = ['gene'] + (['gene_id'] if 'gene_id' in genes else []) + [
            'node', 'chrom', 'tss', 'compartment_A', 'compartment_B', 'score_A', 'score_B', 'delta', 'switch', 'rna_log2FC', 'padj', 'baseMean', 'differential', 'concordant', 'in_cascade']
    genes[cols].sort_values(['chrom', 'tss']).to_csv(os.path.join(out_dir, '{}_genes.tsv'.format(tag)),
                                                    sep='\t', index=False, float_format='%.4g')
    g = genes[np.isfinite(genes.delta) & np.isfinite(genes.rna_log2FC)].copy()
    g['decile'] = pd.qcut(np.log10(g.baseMean + 1), 10, labels=False, duplicates='drop')
    g_sw = g.switch.isin(['B->A', 'A->B'])
    sw = g[g_sw & g.differential]
    k_conc = int(sw.concordant.sum())
    row.update(genes=len(g), differential_genes=int(g.differential.sum()), genes_in_switch_bins=int(g_sw.sum()),
               differential_in_switch_bins=len(sw), concordant=k_conc,
               concordance=k_conc / len(sw) if len(sw) else np.nan,
               binom_p=stats.binomtest(k_conc, len(sw), 0.5, alternative='greater').pvalue if len(sw) else np.nan)
    # references of the same expression level (baseMean decile): differential genes outside switch bins, all of
    # them and only those with |delta| >= threshold (a large change without a change of sign)
    for label, pool in (('nonswitch', g[~g_sw & g.differential & (g.delta != 0)]),
                        ('nonswitch_large_delta', g[~g_sw & g.differential & (g.delta.abs() >= threshold)])):
        values = matched_concordance(sw, pool, rng) if len(sw) else np.array([])
        row['matched_{}_concordance'.format(label)] = float(values.mean()) if len(values) else np.nan
        row['matched_{}_n_pool'.format(label)] = len(pool)
        row['matched_{}_p'.format(label)] = float((1 + np.sum(values >= row['concordance'])) / (1 + len(values))) \
            if len(values) and np.isfinite(row['concordance']) else np.nan
    row['differential_rate_switch'] = float(g[g_sw].differential.mean()) if g_sw.any() else np.nan
    row['differential_rate_nonswitch'] = float(g[~g_sw].differential.mean())
    rho = stats.spearmanr(g.delta, g.rna_log2FC)
    row.update(spearman_delta_log2FC_all=rho[0], spearman_p_all=rho[1])
    d = g[g.differential]
    rho = stats.spearmanr(d.delta, d.rna_log2FC)
    row.update(spearman_delta_log2FC_differential=rho[0], spearman_p_differential=rho[1])
    row['cascade_genes_in_switch_bins'] = int((g_sw & g.in_cascade).sum())
    return row


def main():
    threshold = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    out_dir = os.path.join(cc.OUT, 'switch')
    os.makedirs(out_dir, exist_ok=True)
    summary, consistency = [], []
    rng = np.random.default_rng(0)
    for name, (na, nb), ta, tb, assembly, gtf, expression, gene_types in COMPARISONS:
        print(name, flush=True)
        bins = bins_for(assembly, ta)
        a, b, loadings, explained = joint_scores(bins, ta, tb)
        genes_base, node_of = None, None
        if expression is not None:
            de = pd.read_csv(expression, sep='\t')
            genes_base = de.merge(gene_tss(gtf, gene_types), on='gene_name', how='inner')
            size = int((bins.end - bins.start).max())
            key = pd.DataFrame({'chrom': bins.chrom.values, 'b': bins.start.values // size, 'i': np.arange(len(bins))})
            genes_base['b'] = (genes_base.tss // size).astype(np.int64)
            genes_base = genes_base.merge(key, on=['chrom', 'b'], how='inner')
            path = os.path.join(cc.HOPS, '{}_gene_features.tsv'.format(name))
            if os.path.exists(path):
                lv = pd.read_csv(path, sep='\t', usecols=['node', 'genes'])
                lv = lv.assign(gene=lv.genes.str.split(',')).explode('gene')
                node_of = lv.drop_duplicates('gene').set_index('gene').node
        for side, condition, score in (('A', na, a), ('B', nb, b)):
            label, path = ACCESSIBILITY.get((name, side), (None, None))
            if path is None or not os.path.exists(path):
                consistency.append({'comparison': name, 'condition': condition, 'accessibility': 'none available'})
                continue
            level = smooth_per_chrom(log_enrichment(bin_bigwig(path, bins)), bins.chrom.values, 1)
            centred = level.copy()
            for chrom in pd.unique(bins.chrom.values):
                k = (bins.chrom == chrom).values & np.isfinite(level)
                centred[k] = level[k] - np.median(level[k])
            table = compare_to_reference(score, centred, bins.chrom.values, threshold='median')
            row = table[table.chrom == 'genome'].iloc[0].to_dict()
            row.update(comparison=name, condition=condition, accessibility=label,
                       median_of_chroms_pearson=float(table[table.chrom == 'median_of_chroms'].pearson.iloc[0]))
            consistency.append(row)
            print('consistency', name, condition, label, {k: row[k] for k in ('pearson', 'sign_agreement', 'auc')},
                  flush=True)
        qa, qb = quantile_normalise_pair(a, b, bins.chrom.values)
        for normalisation, (x, y) in (('quantile', (qa, qb)), ('joint', (a, b))):
            row = evaluate(name, na, nb, bins, x, y, threshold, normalisation, loadings, explained, genes_base,
                           node_of, alpha, out_dir, rng)
            summary.append(row)
            print({k: row[k] for k in row if k != 'loadings'}, flush=True)
    pd.DataFrame(summary).to_csv(os.path.join(out_dir, 'switch_summary.tsv'), sep='\t', index=False,
                                 float_format='%.4g')
    pd.DataFrame(consistency).to_csv(os.path.join(out_dir, 'accessibility_consistency.tsv'), sep='\t', index=False,
                                     float_format='%.4g')


if __name__ == '__main__':
    main()
