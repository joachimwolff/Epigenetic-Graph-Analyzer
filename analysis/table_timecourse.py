"""Supplementary Table S33: keratinocyte differentiation time course (Rubin et al. 2017; analysis/22_timecourse.sh).

Part A: graphs per day, differential genes per pair and genes explained per level (declared roles, fragment mean).
Part B: order of the contact, H3K27ac and expression changes over the two intervals (day 0-3, day 3-6) against the
within-gene permutation null and the layer-marginal null; sensitivity rows for the thresholds 0.5 and 1.5 and for
H3K27ac of donor 2 only. Part C: H3K27ac change in the siRNA knockdowns at day 3 at the explaining nodes of the genes
whose candidate trigger is the knocked-down factor (or at explaining nodes carrying the factor's motif), against the
other explaining nodes.

Usage: python table_timecourse.py <timecourse dir> <output dir>
"""

import os
import re
import sys

import pandas as pd

NCOL = 10
PAIRS = ['D0_vs_D3', 'D3_vs_D6', 'D0_vs_D6']
LEVELS = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state', 'step2', 'step3']
CATEGORY = {'forward': 'contact $\\leq$ mark $\\leq$ expr.', 'simultaneous': 'same interval',
            'reversed': 'expr.\\ $\\leq$ mark $\\leq$ contact', 'other': 'other',
            'before': 'first layer earlier', 'same': 'same interval', 'after': 'first layer later'}
LAYERS = {'contact_mark_expression': 'contact, H3K27ac, expression', 'contact_expression': 'contact, expression',
          'mark_expression': 'H3K27ac, expression', 'contact_mark': 'contact, H3K27ac'}


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')


def pval(p):
    if pd.isna(p):
        return '--'
    if p >= 0.01:
        return '{:.2f}'.format(p)
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))


def num(x, digits=2):
    return '--' if pd.isna(x) else '{:.{}f}'.format(x, digits).replace('-', '$-$')


def row(cells):
    cells = [str(c) for c in cells] + [''] * (NCOL - len(cells))
    return ' & '.join(cells) + r' \\'


def heading(text):
    return r'\multicolumn{%d}{l}{\textbf{%s}}\\' % (NCOL, text)


def graph_sizes(base):
    rows = []
    for day in ('D0', 'D3', 'D6'):
        text = open(os.path.join(base, 'graphs', day + '.log')).read()
        match = re.search(r'graph: (\d+) nodes \((\d+) baits\), (\d+) edges', text)
        rows.append((day, *match.groups()) if match else (day, '--', '--', '--'))
    return rows


def part_a(base):
    lines = [heading(r'A. Graphs per day (contacts with reads in both donors, score $\geq$ 5 in one) and '
                     r'genes explained per level'),
             row(['day', 'nodes', 'baits', 'edges']), r'\midrule']
    for day, nodes, baits, edges in graph_sizes(base):
        lines.append(row([day.replace('D', 'day '), nodes, baits, edges]))
    lines += [r'\midrule', row(['pair', 'DE genes'] + [tex(lv.replace('_', ' ')) for lv in LEVELS]), r'\midrule']
    for pair in PAIRS:
        path = os.path.join(base, 'differential', pair + '_explanation_levels.tsv')
        if not os.path.exists(path):
            continue
        table = pd.read_csv(path, sep='\t')
        table = table[(table['roles'] == 'declared') & (table['presence'] == 'mean')].set_index('level')
        cells = ['{} ({})'.format(int(table.at[lv, 'genes_any']), int(table.at[lv, 'genes_first']))
                 if lv in table.index else '--' for lv in LEVELS]
        lines.append(row([pair.replace('_vs_', ' vs ').replace('D', 'day '),
                          int(table['differential_genes'].iloc[0])] + cells))
    return lines


def order_rows(sub, label):
    lines = []
    for i, r in enumerate(sub.itertuples()):
        lines.append(row([label if i == 0 else '', CATEGORY[r.category], int(r.units) if i == 0 else '',
                          '{} ({}--{})'.format(int(r.observed), int(r.boot_low), int(r.boot_high)), num(r.fraction),
                          '{} ({}--{})'.format(num(r.null_mean, 1), int(r.null_low), int(r.null_high)),
                          num(r.ratio), pval(r.p), num(r.marginal_ratio), pval(r.marginal_p)]))
    return lines


def part_b(base):
    ordering = pd.read_csv(os.path.join(base, 'timecourse_ordering.tsv'), sep='\t')
    lines = [r'\midrule', heading('B. Order of the changes (event interval day 0--3 or day 3--6; units = gene and '
                                  'explaining node)'),
             row(['units: layers', 'order', 'units', r'observed (95\% CI)', 'fraction',
                  r'within-gene null (95\%)', 'ratio', '$P$', 'marginal ratio', 'marginal $P$']), r'\midrule']
    main = ordering[ordering['chromatin_threshold'] == 1.0]
    selection = [('all', 'contact_mark_expression'), ('promoter', 'contact_mark_expression'),
                 ('step1', 'contact_mark_expression'), ('step2', 'contact_mark_expression'),
                 ('step3', 'contact_mark_expression'), ('all', 'contact_expression'), ('all', 'mark_expression'),
                 ('all', 'contact_mark')]
    for subset, layers in selection:
        lines += order_rows(main[(main['subset'] == subset) & (main['layers'] == layers)],
                            '{}: {}'.format(subset, LAYERS[layers]))
    for threshold in sorted(set(ordering['chromatin_threshold']) - {1.0}):
        sub = ordering[(ordering['chromatin_threshold'] == threshold) & (ordering['subset'] == 'all') &
                       (ordering['layers'] == 'contact_mark_expression')]
        lines += order_rows(sub, 'all, threshold {:g}'.format(threshold))
    donor2 = os.path.join(base, 'timecourse_donor2_ordering.tsv')
    if os.path.exists(donor2):
        d2 = pd.read_csv(donor2, sep='\t')
        d2 = d2[d2['chromatin_threshold'] == 1.0]
        for layers in ('contact_mark_expression', 'mark_expression'):
            lines += order_rows(d2[(d2['subset'] == 'all') & (d2['layers'] == layers)],
                                'H3K27ac donor 2 only: {}'.format(LAYERS[layers]))
    return lines


def part_c(base):
    path = os.path.join(base, 'timecourse_knockdown_test.tsv')
    if not os.path.exists(path):
        return []
    test = pd.read_csv(path, sep='\t')
    lines = [r'\midrule', heading('C. H3K27ac siRNA knockdown vs control (day 3) at explaining nodes (day 0 vs 3 run)'),
             row(['factor: grouping (genes)', 'cand.\\ genes', r'nodes cand.\ / other',
                  r'median log2 cand.\ / other', r'median diff.\ (95\% CI)', r'Cliff $\delta$ (95\% CI)',
                  'rep1 / rep2 diff.', '$P$']), r'\midrule']
    main = test[(test['measure'] == 'mean') & (test['evidence'].isin(['any', 'motif at node']))]
    for r in main.itertuples():
        reps = test[(test['factor'] == r.factor) & (test['evidence'] == r.evidence) & (test['genes'] == r.genes) &
                    (test['measure'] != 'mean')].sort_values('measure')
        rep_text = ' / '.join(num(x, 3) for x in reps['median_difference']) if len(reps) else '--'
        motif = '' if r.motif_in_jaspar in (True, 'True') else ' (no JASPAR motif)'
        grouping = 'candidate trigger' if r.evidence == 'any' else 'motif at node'
        lines.append(row(['{}{}: {} ({})'.format(r.factor, motif, grouping, r.genes), int(r.genes_candidate),
                          '{} / {}'.format(int(r.nodes_candidate), int(r.nodes_control)),
                          '{} / {}'.format(num(r.median_candidate, 3), num(r.median_control, 3)),
                          '{} ({}--{})'.format(num(r.median_difference, 3), num(r.md_low, 3), num(r.md_high, 3)),
                          '{} ({}--{})'.format(num(r.cliffs_delta, 3), num(r.delta_low, 3), num(r.delta_high, 3)),
                          rep_text, pval(r.p_two_sided)]))
    return lines


def main():
    base, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    lines = [r'\begin{tabular}{' + 'l' * NCOL + '}', r'\toprule'] + part_a(base) + part_b(base) + part_c(base) + \
        [r'\bottomrule', r'\end{tabular}']
    caption = (r'\caption{Keratinocyte differentiation time course (Rubin et al.\ 2017, two donors, days 0, 3 and 6). '
               r'A: graphs per day and genes explained per level in the pairwise runs (any level; first level in '
               r'parentheses). Capture Hi-C contacts are the deposited CHiCAGO calls (FDR $<$ 0.01 at any day); '
               r'bait-to-non-bait contacts were deposited only where the other end overlaps an H3K27ac peak. '
               r'B: event interval of a layer = first interval with $|\log_2|$ change $\geq 1$ (expression: '
               r'DESeq2 padj $\leq 0.05$); with two intervals an order is resolved only as earlier versus same '
               r'interval. Within-gene null: event intervals shuffled among the layers of each gene; marginal null: '
               r'intervals of each layer shuffled across units (10,000 permutations each); CI: bootstrap over genes. '
               r'C: two-sided Mann-Whitney $P$; intervals from a bootstrap over genes.}')
    body = '\n'.join([r'\begin{table*}', r'\centering', r'\scriptsize', r'\resizebox{\textwidth}{!}{%'] + lines +
                     ['}', caption, r'\label{tab:s33}', r'\end{table*}'])
    with open(os.path.join(out, 'table_s33.tex'), 'w') as handle:
        handle.write(body + '\n')
    print(body)


if __name__ == '__main__':
    main()
