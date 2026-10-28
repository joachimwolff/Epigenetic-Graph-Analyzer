"""Supplementary Table S36: TAD boundaries crossed by gained, lost and shared contacts (analysis/28_tad_crossing.sh).

Part A: per comparison and boundary set, the fraction of contacts crossing a boundary per condition and per status,
and the promoters / boundaries with a significant enrichment. Part B: the focus promoter (Hoxa2 for Md / PA2, EPHA7
for hESC / hNEC) and the boundary crossed by its gained contacts. Part C: the five best promoters per comparison
(gained contacts, boundaries at or above the median strength). Part D: the five best boundaries. Part E: egaBoundaries
rules (no Hi-C) against the hicFindTADs boundaries of the same condition and their call at chr6:51.175 Mb.

Written as a longtable (LaTeX package longtable), caption below the last row.

Usage: python table_tad_crossing.py <tad_crossing dir> <output dir>
"""

import os
import sys

import numpy as np
import pandas as pd

SETTING = {'all': 'all', 's50': '>= median'}
TADS = {'PA2': 'PA2', 'Md': 'Md', 'H1': 'H1'}
RULE = {'dip': 'capture dip', 'ctcf': 'CTCF', 'ctcf_pattern': 'CTCF pattern', 'loop_anchor': 'loop anchor',
        'dip_or_ctcf': 'dip OR CTCF', 'dip_and_ctcf': 'dip AND CTCF', 'dip_or_ctcf_pattern': 'dip OR pattern',
        'dip_and_ctcf_pattern': 'dip AND pattern'}
NCOL = 11


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('>=', r'$\geq$')


def f(x, d=2):
    return '--' if x is None or pd.isna(x) else '{:.{}f}'.format(float(x), d)


def p(x):
    if x is None or pd.isna(x):
        return '--'
    x = float(x)
    if x == 0:
        return '0'
    if x >= 0.01:
        return '{:.2f}'.format(x)
    m, e = '{:.1e}'.format(x).split('e')
    return r'${}\times10^{{{}}}$'.format(m, int(e))


def i(x):
    return '--' if x is None or pd.isna(x) else '{:d}'.format(int(x))


def row(cells):
    cells = list(cells) + [''] * (NCOL - len(cells))
    return ' & '.join(str(c) for c in cells) + r' \\'


def heading(text):
    return r'\multicolumn{%d}{p{0.9\textwidth}}{\textbf{%s}}\\' % (NCOL, text)


def label(setting):
    tads, kind = setting.split('_tads_')
    return '{}, {}'.format(TADS.get(tads, tads), SETTING.get(kind, kind))


def comparison_names(comparison):
    return comparison.split('_vs_')


def part_a(s):
    lines = [heading('A. Contacts crossing a boundary (fraction) and significant enrichments (q < 0.05)'),
             row(['comparison', 'boundaries', 'cond. A', 'cond. B', 'gained', 'lost', 'shared',
                  'prom. tested', 'prom. gained', 'prom. lost', 'bound. gained / lost']), r'\midrule']
    for r in s.itertuples():
        a, b = comparison_names(r.comparison)
        lines.append(row([tex(r.comparison.replace('_vs_', ' / ')), tex(label(r.setting)),
                          f(getattr(r, 'crossing_' + a)), f(getattr(r, 'crossing_' + b)), f(r.crossing_gained),
                          f(r.crossing_lost), f(r.crossing_shared), i(r.promoters_tested_gained),
                          i(r.promoters_q05_gained), i(r.promoters_q05_lost),
                          '{} / {}'.format(i(r.boundaries_q05_gained), i(r.boundaries_q05_lost))]))
    return lines


def part_b(s):
    lines = [r'\midrule', heading('B. Focus promoter (Hoxa2: Md / PA2, EPHA7: hESC / hNEC) and the boundary crossed by most of its gained contacts (strength rank)'),
             row(['boundaries', 'gene', 'gained c / w', 'lost c / w', 'OR', 'P', 'q', 'rank',
                  'boundary (rank)', 'gained obs / exp, q']), r'\midrule']
    for r in s.itertuples():
        if pd.isna(getattr(r, 'focus_gene', np.nan)):
            continue
        boundary = getattr(r, 'focus_boundary', '')
        if isinstance(boundary, str) and ':' in boundary:
            chrom, pos = boundary.split(':')
            boundary = '{}:{:.3f} ({})'.format(chrom, int(pos) / 1e6,
                                                  i(getattr(r, 'focus_boundary_strength_rank', np.nan)))
        else:
            boundary = '--'
        lines.append(row([tex(label(r.setting)), r'\textit{%s}' % r.focus_gene,
                          '{} / {}'.format(i(r.focus_gained_crossing), i(r.focus_gained_within)),
                          '{} / {}'.format(i(r.focus_lost_crossing), i(r.focus_lost_within)),
                          f(r.focus_gained_odds_ratio, 1), p(r.focus_gained_p), p(r.focus_gained_q),
                          i(r.focus_gained_rank), boundary,
                          '{} / {}, {}'.format(i(getattr(r, 'focus_boundary_gained', np.nan)),
                                               f(getattr(r, 'focus_boundary_gained_expected', np.nan), 1),
                                               p(getattr(r, 'focus_boundary_gained_q', np.nan)))]))
    return lines


def main_setting(comparison):
    return 'PA2_tads_s50' if comparison == 'Md_vs_PA2' else 'H1_tads_s50'


def part_c(top):
    lines = [r'\midrule', heading('C. Five best promoters per comparison: gained contacts enriched among the '
                                  'boundary-crossing contacts (boundaries at or above the median strength)'),
             row(['comparison', 'rank', 'genes', 'fragment', 'gained cross / within', 'lost cross / within',
                  'P', 'q']), r'\midrule']
    for comparison in sorted(top.comparison.unique()):
        t = top[(top.comparison == comparison) & (top.setting == main_setting(comparison)) &
                (top.direction == 'gained')].head(5)
        for k, r in enumerate(t.itertuples()):
            genes = ','.join(str(r.genes).split(',')[:2]) if isinstance(r.genes, str) else ''
            lines.append(row([tex(comparison.replace('_vs_', ' / ')) if k == 0 else '', i(r.gained_rank),
                              r'\textit{%s}' % tex(genes), tex(r.node),
                              '{} / {}'.format(i(r.gained_crossing), i(r.gained_within)),
                              '{} / {}'.format(i(r.lost_crossing), i(r.lost_within)),
                              p(r.gained_p), p(r.gained_q)]))
    return lines


def part_d(top):
    lines = [r'\midrule', heading('D. Five best boundaries per comparison: gained contacts crossing the boundary '
                                  'against the expectation from contacts of the same length'),
             row(['comparison', 'rank', 'boundary', 'strength rank', 'crossing', 'gained obs / exp',
                  'lost obs / exp', 'P', 'q']), r'\midrule']
    for comparison in sorted(top.comparison.unique()):
        t = top[(top.comparison == comparison) & (top.setting == main_setting(comparison)) &
                (top.direction == 'gained')].head(5)
        for k, r in enumerate(t.itertuples()):
            lines.append(row([tex(comparison.replace('_vs_', ' / ')) if k == 0 else '', i(r.gained_rank),
                              '{}:{:.3f} Mb'.format(r.chrom, r.pos / 1e6), i(r.strength_rank), i(r.crossing),
                              '{} / {}'.format(i(r.gained_crossing), f(r.gained_expected, 1)),
                              '{} / {}'.format(i(r.lost_crossing), f(r.lost_expected, 1)),
                              p(r.gained_p), p(r.gained_q)]))
    return lines


def part_e(ega):
    lines = [r'\midrule', heading('E. egaBoundaries rules (capture graph, Kessler CTCF, CTCF motif; no Hi-C) against '
                                  'the hicFindTADs boundaries of the same condition (tolerance 1 bin)'),
             row(['condition', 'rule', 'calls', 'P', 'R', 'F1', 'P covered', 'R covered', 'call at 51.175 Mb']),
             r'\midrule']
    sub = ega[(ega.rule != 'focus_bin') & (ega.tolerance_bins == 1)]
    for cond in sorted(sub.condition.unique()):
        for k, r in enumerate(sub[sub.condition == cond].itertuples()):
            lines.append(row([cond if k == 0 else '', tex(RULE.get(r.rule, r.rule)), i(r.calls), f(r.precision),
                              f(r.recall), f(r.f1), f(r.precision_covered), f(r.recall_covered),
                              'yes' if str(r.calls_focus_boundary).lower() == 'true' else 'no']))
    return lines


def main():
    base, out = sys.argv[1], sys.argv[2]
    summary = os.path.join(base, 'summary')
    s = pd.read_csv(os.path.join(summary, 'settings.tsv'), sep='\t')
    top_p = pd.read_csv(os.path.join(summary, 'top_promoters.tsv'), sep='\t')
    top_b = pd.read_csv(os.path.join(summary, 'top_boundaries.tsv'), sep='\t')
    ega_path = os.path.join(summary, 'egaboundaries.tsv')
    parts = [part_a(s), part_b(s), part_c(top_p), part_d(top_b)]
    if os.path.exists(ega_path):
        parts.append(part_e(pd.read_csv(ega_path, sep='\t')))
    blocks = []
    for part in parts:
        # one longtable per part (own column widths); the leading \midrule of a part is dropped
        part = part[1:] if part and part[0] == r'\midrule' else part
        blocks.append([r'\begin{longtable}{' + 'l' * NCOL + '}', r'\toprule'] + part + [r'\bottomrule'])
    caption = (r'\caption{TAD boundaries crossed by gained, lost and shared capture Hi-C contacts. A contact crosses '
               r'a boundary when the boundary midpoint lies between the midpoints of its two fragments. Boundaries: the '
               r'hicFindTADs calls of the named sample, every boundary (all) or those at or above the median strength '
               r'(PA2, Md: Kessler et al. 2023, 25 kb; H1: 4DN Micro-C, 50 kb; c / w: crossing / within its domain) '
               r'($\geq$ median); boundary '
               r'strength is the negative hicFindTADs TAD-separation score (strength rank 1: the strongest). '
               r'Promoter enrichment: one-sided Fisher test of gained (lost) contacts among the boundary-crossing '
               r'against the non-crossing contacts of a promoter fragment, Benjamini-Hochberg over the tested '
               r'promoters (OR with 0.5 added to every cell). Boundary enrichment: gained (lost) contacts crossing a '
               r'boundary against the sum over these contacts of the genome-wide fraction of gained (lost) union '
               r'contacts of the same length (20 length quantiles), exact upper tail of the sum of binomials, '
               r'Benjamini-Hochberg over the boundaries. Hi-C: Kessler et al. 2023 (GSE211901), balanced 25 kb, '
               r'hicFindTADs minDepth 75 kb, maxDepth 250 kb, FDR 0.05; H1: 4DN Micro-C, hicFindTADs 50 kb, '
               r'lifted from GRCh38 to hg19. egaBoundaries rules as in Table S35, without cohesin tracks (every CTCF '
               r'summit treated as cohesin-bound).}')
    lines = []
    for k, block in enumerate(blocks):
        lines += block
        if k == len(blocks) - 1:
            lines.append(caption + r'\label{tab:s36}\\')
        lines.append(r'\end{longtable}')
        if k < len(blocks) - 1:
            lines.append(r'\addtocounter{table}{-1}')     # uncaptioned parts do not take a table number
    body = '\n'.join([r'\begingroup', r'\scriptsize', r'\setlength{\tabcolsep}{3pt}'] + lines + [r'\endgroup'])
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, 'table_s36.tex')
    with open(path, 'w') as h:
        h.write(body + '\n')
    print(path)


if __name__ == '__main__':
    main()
