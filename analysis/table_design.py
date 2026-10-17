"""Supplementary Table S34: experiment design (egaDesign, analysis/25_design.sh).

Part A: the three best ranked candidates of the genes shown in the figures (figure_design.py defaults). Part B: per
comparison and level the explained genes and those whose level is removed by flattening a single node. Part C: K562
CRISPRi validation of the single-condition ranking (fraction of significant pairs per design rank bin and the
likelihood-ratio tests against log distance and element DNase).

Usage: python table_design.py <design dir> <output dir>
"""

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_design import COMPARISONS, DEFAULT_GENES, LABELS, LEVEL_LABELS, LEVELS, gene_rows  # noqa: E402
from table_triggers import longtable, num, pval, tex  # noqa: E402


def stat_change(r):
    if pd.isna(r['statistic_before']):
        return 'state ' + ('lost' if str(r['level_lost_flatten']).lower() == 'true' else 'kept')
    return '{} $\\to$ {}'.format(num(r['statistic_before'], 1), num(r['statistic_flatten'], 1))


def part_a(directory, genes=DEFAULT_GENES, top=3):
    rows = []
    for comparison, gene in genes:
        path = os.path.join(directory, comparison + '_design.tsv')
        if not os.path.exists(path):
            continue
        t = gene_rows(pd.read_csv(path, sep='\t'), gene).head(top)
        if rows:
            rows.append(r'\midrule')
        for i, r in enumerate(t.to_dict('records')):
            what = '{} ({:+.1f})'.format(r['gene'], r['rna_log2FC'])
            level = LEVEL_LABELS.get(r['level'], r['level'])
            p = '{} / {}'.format(pval(r['p_before_reduced']), pval(r['p_after_reduced'])) \
                if not pd.isna(r['p_before_reduced']) else '--'
            rows.append(' & '.join([LABELS[comparison] if i == 0 else '', tex(what) if i == 0 else '',
                                    str(int(r['rank'])), tex('{}:{:,}-{:,}'.format(r['chrom'], r['start'], r['end'])),
                                    tex(r['element_class']), tex(r['position']), tex(level), stat_change(r),
                                    num(r['effect_fraction']), p, str(int(r['evidence_points'])),
                                    str(int(r['n_genes_strong'])), tex(r['expected_direction'])]) + r' \\')
    header = (r'Comparison & Gene (log2FC) & Rank & Node & Element & Position & Level & Statistic & Removed & '
              r'$P$ before / after & Ev. & Genes & Expected')
    caption = ('Best ranked perturbation candidates of selected genes (egaDesign). Statistic: the level statistic '
               '(largest concordant log2 change over the marks, contact-weighted at step $k$) before and after the '
               'mark change of the node is set to 0; state: the step-1 state level is lost when the node is the only '
               'contacting fragment with the poised transition. Removed: fraction of the statistic removed. $P$: '
               'per-gene permutation $P$ of the level before and after (1000 permutations). Ev.: evidence points '
               '(own concordant mark change, level $P<0.05$, beyond mirror, 3D support, contact in both conditions, '
               'non-promoter element, not pleiotropic). Genes: genes for which the node removes at least half of a '
               'level. Expected: direction of the gene\'s change when the element is silenced (active mark) or '
               'deleted (repressive mark). Candidates, not decisions.')
    return longtable('lllp{2.6cm}p{2.6cm}p{2.2cm}lllllll', header, rows, caption, 'tab:s34')


def part_b(directory):
    rows = []
    for comparison in COMPARISONS:
        path = os.path.join(directory, comparison + '_design_summary.json')
        if not os.path.exists(path):
            continue
        s = json.load(open(path))
        first = True
        for level in LEVELS:
            if level not in s['by_level']:
                continue
            v = s['by_level'][level]
            rows.append(' & '.join([LABELS[comparison] if first else '',
                                    '{:,}'.format(s['genes']) if first else '', tex(LEVEL_LABELS[level]),
                                    '{:,}'.format(v['genes']), '{:,}'.format(v['genes_single_node_removes_level']),
                                    num(v['genes_single_node_removes_level'] / v['genes'])]) + r' \\')
            first = False
    header = r'Comparison & Explained genes & Level & Genes & One node removes the level & Fraction'
    caption = ('Explanations that depend on a single node: per comparison and explaining level the genes whose '
               'level no longer explains after flattening one node (its mark change set to 0; own contacts: the '
               'edge change; step-1 state: the only fragment with the transition).')
    return longtable('lrlrrr', header, rows, caption, 'tab:s34b')


def part_c(directory):
    prefix = os.path.join(directory, 'validation', 'K562')
    if not os.path.exists(prefix + '_bins.tsv'):
        return []
    bins = pd.read_csv(prefix + '_bins.tsv', sep='\t')
    rows = []
    for name, part in bins.groupby('set', sort=False):
        first = True
        for r in part.to_dict('records'):
            rows.append(' & '.join([tex(name) if first else '', tex(r['rank_bin']), '{:,}'.format(r['pairs']),
                                    '{:,}'.format(r['positives']),
                                    '{} ({}--{})'.format(num(r['fraction']), num(r['ci_low']), num(r['ci_high']))])
                        + r' \\')
            first = False
        rows.append(r'\midrule')
    tests = pd.read_csv(prefix + '_added_value.tsv', sep='\t')
    for r in tests.to_dict('records'):
        if pd.isna(r.get('lr_p')):
            continue
        rows.append(' & '.join(['LR test, {}'.format(tex(r['set'])), tex(r['feature']),
                                '{:,}'.format(int(r['n'])), '{:,}'.format(int(r['positives'])),
                                'z {} , $P$ {}'.format(num(r['z']), pval(r['lr_p']))]) + r' \\')
    summary = json.load(open(prefix + '_summary.json'))
    if 'within_genes' in summary:
        rows.append(r'\midrule')
        rows.append(' & '.join(['within gene ({} genes)'.format(summary['within_genes']), 'AUROC design / distance',
                                '', '', '{} / {}, $P$ {}'.format(num(summary['within_mean_auroc_design']),
                                                                 num(summary['within_mean_auroc_distance']),
                                                                 pval(summary['within_wilcoxon_design_vs_distance_p']))
                                ]) + r' \\')
    header = r'Set & Design rank / feature & Pairs & Significant & Fraction (95\% CI) / test'
    caption = ('K562 CRISPRi validation of the single-condition design ranking: the nodes within three steps of a '
               'gene on the K562 pcHi-C graph ranked by their contribution (path weight times DNase signal) to the '
               'gene\'s DNase activity; ENCODE-rE2G K562 training pairs mapped as in the CRISPRi benchmark (TSS '
               'fragment a bait, element not on the TSS fragment). Not reached: element farther than three steps or '
               'not in the graph. LR test: logistic likelihood-ratio test of the feature added to log distance and '
               'element DNase (given\\_direct: also a direct edge); share: contribution share; top3 / top10: design '
               'rank at most 3 / 10. Within gene: mean AUROC over genes with significant and non-significant reached '
               'elements, Wilcoxon signed-rank $P$.')
    return longtable('llrrl', header, rows, caption, 'tab:s34c')


def main(argv):
    directory, output = argv[1], argv[2]
    os.makedirs(output, exist_ok=True)
    lines = part_a(directory) + [''] + part_b(directory) + [''] + part_c(directory)
    with open(os.path.join(output, 'table_s34.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv)
