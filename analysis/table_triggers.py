"""Supplementary Table S31: motifs enriched at the explaining nodes per gene group, and the candidate trigger factors
of selected genes (egaTriggers, analysis/20_triggers.sh).

Part A: per comparison and gene group the motifs with q < 0.05 against both backgrounds (GC- and length-matched
random marked nodes and the non-explaining nodes of the same genes), at most five per group, with the fraction of
foreground and background windows carrying the motif, the log2 odds ratio against random marked nodes and both q.
Groups without such a motif are listed with their sizes. Part B: the five best ranked candidates of the genes shown
in the figures (figure_triggers.py defaults).

Usage: python table_triggers.py <triggers dir> <output dir>
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_triggers import COMPARISONS, combined, default_genes  # noqa: E402

LABELS = {'hESC_vs_hNEC': 'hESC vs hNEC', 'Md_vs_PA2': 'Md vs PA2', 'siMM_vs_siRE': 'Kasumi-1 siMM vs siRE',
          't821_vs_FLT3ITD': 't(8;21) vs FLT3-ITD'}


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('#', r'\#')


def pval(p):
    if pd.isna(p):
        return '--'
    if p >= 0.01:
        return '{:.2f}'.format(p)
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))


def num(x, digits=2):
    return '--' if pd.isna(x) else '{:.{}f}'.format(x, digits).replace('-', '$-$')


def longtable(columns, header, rows, caption, label):
    lines = [r'{\footnotesize', r'\begin{longtable}{@{}' + columns + '@{}}', r'\caption{' + caption + r'}\label{'
             + label + r'}\\', r'\toprule', header + r' \\', r'\midrule', r'\endfirsthead', r'\toprule',
             header + r' \\', r'\midrule', r'\endhead']
    return lines + rows + [r'\bottomrule', r'\end{longtable}}']


def part_a(directory, top=5, alpha=0.05):
    rows = []
    for comparison in COMPARISONS:
        path = os.path.join(directory, comparison + '_motif_enrichment.tsv')
        if not os.path.exists(path) or os.path.getsize(path) < 2:
            continue
        enrichment = pd.read_csv(path, sep='\t')
        sizes = enrichment[enrichment['background'] == 'random_marked'].drop_duplicates('group').set_index('group')
        table = combined(enrichment)
        if rows:
            rows.append(r'\midrule')
        for group in sizes.index:
            hits = table[(table['group'] == group) & (table['q'] < alpha)].sort_values('q').head(top)
            size = '{} / {:,}'.format(int(sizes.at[group, 'n_genes']), int(sizes.at[group, 'fg_windows']))
            if hits.empty:
                rows.append(' & '.join([LABELS[comparison], tex(group), size, 'none', '', '', '', '', '']) + r' \\')
                continue
            for i, r in enumerate(hits.itertuples()):
                rows.append(' & '.join([LABELS[comparison] if i == 0 else '', tex(group) if i == 0 else '',
                                        size if i == 0 else '', tex(r.motif), num(r.fg_fraction), num(r.bg_fraction),
                                        num(r.log2_odds), pval(r.q_random), pval(r.q_same)]) + r' \\')
    header = (r'Comparison & Group & Genes / windows & Motif & Fg & Bg & log2 OR & $q_\mathrm{random}$ & '
              r'$q_\mathrm{same}$')
    caption = ('Motifs enriched at the explaining nodes of groups of differentially expressed genes (JASPAR 2024 '
               'CORE vertebrates, hit at $P<10^{-4}$ in the 1\\,kb window of highest active-mark signal per '
               'fragment). Groups: genes up or down explained by promoter marks (promoter), by marks at step $k$ '
               '(step$k$), by a poised contacting fragment (step1\\_state), with a promoter bivalent in both '
               'conditions (poised\\_both: promoter fragments; \\_steps: their explaining nodes) and promoter state '
               'transitions involving a bivalent state. Genes / windows: genes of the group and scanned foreground '
               'windows. Fg / Bg: fraction of foreground and of GC- and length-matched random marked windows with '
               'a hit; log2 OR against random marked windows; $q$ (one-sided Fisher test, BH across motifs) against '
               'random marked windows and against the non-explaining nodes at the same steps of the same genes '
               '(promoter groups: promoters of the other differential genes). Listed: $q<0.05$ against both, at '
               'most five per group.')
    return longtable('llrlrrrrr', header, rows, caption, 'tab:s31')


def part_b(directory, top=5):
    rows = []
    for comparison, gene in default_genes(directory):
        path = os.path.join(directory, comparison + '_gene_triggers.tsv')
        if not os.path.exists(path):
            continue
        t = pd.read_csv(path, sep='\t')
        t = t[t['gene'].astype(str).str.split(',').apply(lambda gs: gene in gs)].sort_values('rank')
        if t.empty:
            rows.append(' & '.join([LABELS[comparison], tex(gene), '', 'no candidate'] + [''] * 6) + r' \\')
            continue
        head = t.iloc[0]
        t = t[t['gene_node'] == head['gene_node']]  # several promoter fragments: the first listed
        for i, r in enumerate(t.drop_duplicates('tf').head(top).itertuples()):
            what = '{} ({:+.1f}, {})'.format(gene, head['gene_log2FC'], head['explained_at'])
            rows.append(' & '.join([LABELS[comparison] if i == 0 else '', tex(what) if i == 0 else '',
                                    str(r.rank), tex(r.tf) + ('' if r.evidence == 'named factor' else '$^f$'),
                                    tex(r.motif), num(r.tf_log2FC, 1), tex(r.relation.split('-')[0]),
                                    '{}/{}'.format(r.n_nodes_with_motif, r.n_nodes_scanned), pval(r.gene_motif_p),
                                    tex(r.steps_with_motif)]) + r' \\')
    header = 'Comparison & Gene (log2FC, explained at) & Rank & Factor & Motif & log2FC & Relation & Nodes & P & Steps'
    caption = ('Candidate trigger factors of selected genes: factors whose motif occurs at the promoter fragment or '
               'the explaining nodes of the gene and whose expression changes (padj $<0.05$; t(8;21) vs FLT3-ITD: '
               '$P<0.05$). Relation: activator (factor changes with the gene) or repressor (against it). Nodes: '
               'scanned nodes of the gene with the motif. $^f$ family member of the factor named by the motif. '
               'P: binomial P of the nodes with the motif against the motif frequency in marked windows. Ranked by '
               'the mean rank of P, motif score, factor change, graph distance to the factor and group enrichment; '
               'candidates, not decisions.')
    return longtable('lllllrlrrl', header, rows, caption, 'tab:s31b')


def main(argv):
    directory, output = argv[1], argv[2]
    os.makedirs(output, exist_ok=True)
    lines = part_a(directory) + [''] + part_b(directory)
    with open(os.path.join(output, 'table_s31.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv)
