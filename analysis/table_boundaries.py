"""Supplementary Table S35: rule-based domain boundaries from capture contacts, CTCF / cohesin and motif orientation
(analysis/27_boundaries.sh).

Part A: every rule (a priori thresholds) against the Hi-C insulation boundaries of the same cell type, all bins:
precision, recall, F1 (tolerance 1 bin) and AUPRC per resolution and cell type. Part B: the same, covered bins only.
Part C: adding the capture dip to CTCF: change of precision, recall, F1 and AUPRC with 95 % chromosome
block-bootstrap intervals. Part D: other Hi-C references (and the agreement of the Hi-C insulation boundaries with
them). Part E: F1 over the threshold grid. Part F: inputs, loop orientation and domain agreement.

Written as a longtable (LaTeX package longtable) because the table spans several pages.

Usage: python table_boundaries.py <boundaries dir> <output dir>
"""

import os
import sys

import pandas as pd

RULES = ['dip', 'ctcf', 'ctcf_pattern', 'loop_anchor', 'dip_or_ctcf', 'dip_and_ctcf', 'dip_or_ctcf_pattern',
         'dip_and_ctcf_pattern']
LABEL = {'dip': 'dip', 'ctcf': 'CTCF', 'ctcf_pattern': 'CTCF pattern', 'loop_anchor': 'loop anchor',
         'dip_or_ctcf': 'dip OR CTCF', 'dip_and_ctcf': 'dip AND CTCF', 'dip_or_ctcf_pattern': 'dip OR pattern',
         'dip_and_ctcf_pattern': 'dip AND pattern', 'hic_insulation_minima': 'Hi-C insulation boundaries'}
CELLS = ['GM12878', 'K562', 'H1']
NCOL = 14


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')


def f(x, d=2):
    return '--' if pd.isna(x) else '{:.{}f}'.format(x, d)


def row(cells):
    cells = list(cells) + [''] * (NCOL - len(cells))
    return ' & '.join(str(c) for c in cells) + r' \\'


def heading(text):
    return r'\multicolumn{%d}{p{0.9\textwidth}}{\textbf{%s}}\\' % (NCOL, text)


def part_metrics(ev, region, title):
    sub = ev[(ev.reference == 'hic_insulation') & (ev.region == region)]
    cells = ' & '.join(['', ''] + [r'\multicolumn{4}{c}{%s}' % c for c in CELLS]) + r' \\'
    lines = [r'\midrule', heading(title), cells,
             row(['bin', 'rule'] + ['P', 'R', 'F1', 'AUPRC'] * 3), r'\midrule']
    for res in (10, 25, 50):
        for i, rule in enumerate(RULES):
            vals = []
            for c in CELLS:
                g = sub[(sub.resolution_kb == res) & (sub.rule == rule) & (sub.cell == c)]
                vals += [f(g.precision_tol1.mean()), f(g.recall_tol1.mean()), f(g.f1_tol1.mean()),
                         f(g.auprc.mean(), 3)]
            lines.append(row(['{} kb'.format(res) if i == 0 else '', tex(LABEL[rule])] + vals))
        ncalls = []
        for c in CELLS:
            g = sub[(sub.resolution_kb == res) & (sub.cell == c)]
            ncalls += ['', '', 'n ref {}'.format(int(g.n_reference_tol1.iloc[0])), f(g.positive_rate.iloc[0], 3)]
        lines.append(row(['', 'reference / positive rate'] + ncalls))
    return lines


def part_delta(delta):
    lines = [r'\midrule', heading('C. Adding the capture dip to CTCF (reference Hi-C insulation, tolerance 1 bin; '
                                  '95\\% chromosome bootstrap interval)'),
             row(['bin', 'cell', 'region', r'$\Delta$P', 'CI', r'$\Delta$R', 'CI', r'$\Delta$F1', 'CI',
                  r'$\Delta$AUPRC', 'CI', 'calls']), r'\midrule']
    d = delta[delta.reference == 'hic_insulation']
    ci = '[{},{}]'
    for (big, small), g in d.groupby(['with', 'without'], sort=False):
        lines.append(r'\multicolumn{%d}{l}{\textit{%s vs %s}}\\' % (NCOL, tex(LABEL[big]), tex(LABEL[small])))
        for _, r in g.iterrows():
            lines.append(row(['{} kb'.format(int(r.resolution_kb)), r.cell, r.region,
                              f(r.delta_precision), ci.format(f(r.precision_lo), f(r.precision_hi)),
                              f(r.delta_recall), ci.format(f(r.recall_lo), f(r.recall_hi)),
                              f(r.delta_f1), ci.format(f(r.f1_lo), f(r.f1_hi)),
                              f(r.delta_auprc, 3), ci.format(f(r.auprc_lo, 3), f(r.auprc_hi, 3)),
                              '{}/{}'.format(int(r.n_called_with), int(r.n_called_without))]))
    return lines


def part_refs(ev):
    lines = [r'\midrule', heading('D. Other Hi-C references (all bins)'),
             row(['reference', 'cell', 'bin', 'rule', 'calls', 'P', 'R', 'F1', 'F1 tol. 2', 'AUPRC',
                  'positive rate']), r'\midrule']
    sub = ev[ev.region == 'all']
    for ref in ('hicFindTADs_50kb', 'rao2014_arrowhead', '4dn_H1_microc'):
        for c in CELLS:
            for res in (10, 25, 50):
                for rule in ('hic_insulation_minima', 'dip', 'ctcf_pattern', 'dip_and_ctcf_pattern'):
                    g = sub[(sub.reference == ref) & (sub.cell == c) & (sub.resolution_kb == res) & (sub.rule == rule)]
                    if g.empty:
                        continue
                    g = g.iloc[0]
                    lines.append(row([tex(ref), c, '{} kb'.format(res), tex(LABEL[rule]), int(g.n_called_tol1),
                                      f(g.precision_tol1), f(g.recall_tol1), f(g.f1_tol1), f(g.f1_tol2),
                                      f(g.auprc, 3), f(g.positive_rate, 3)]))
    return lines


def part_grid(grid):
    pts = grid[['dip_prominence', 'ctcf_min']].drop_duplicates().sort_values(['dip_prominence', 'ctcf_min'])
    lines = [r'\midrule', heading('E. F1 (tolerance 1, all bins, reference Hi-C insulation) over the threshold grid '
                                  '(dip prominence / CTCF enrichment; a priori 0.3 / 10)'),
             row(['bin', 'cell', 'rule'] + ['{:g} / {:g}'.format(a, b) for a, b in pts.values]), r'\midrule']
    sub = grid[(grid.region == 'all') & (grid.reference == 'hic_insulation')]
    for res in (10, 25, 50):
        for c in CELLS:
            for i, rule in enumerate(('dip', 'ctcf', 'ctcf_pattern', 'dip_and_ctcf', 'dip_and_ctcf_pattern',
                                      'dip_or_ctcf')):
                vals = []
                for a, b in pts.values:
                    g = sub[(sub.resolution_kb == res) & (sub.cell == c) & (sub.rule == rule) &
                            (sub.dip_prominence == a) & (sub.ctcf_min == b)]
                    vals.append(f(g.f1_tol1.mean()))
                lines.append(row(['{} kb'.format(res) if i == 0 and c == CELLS[0] else '', c if i == 0 else '',
                                  tex(LABEL[rule])] + vals))
    return lines


def part_inputs(tables):
    cov = pd.read_csv(os.path.join(tables, 'coverage.tsv'), sep='\t')
    qc = pd.read_csv(os.path.join(tables, 'sites_qc.tsv'), sep='\t')
    lo = pd.read_csv(os.path.join(tables, 'loop_orientation.tsv'), sep='\t')
    dom = pd.read_csv(os.path.join(tables, 'eval_rules_domains.tsv'), sep='\t')
    lines = [r'\midrule', heading('F. Inputs: coverage, CTCF sites, orientation of edges with cohesin-bound CTCF at '
                                  'both ends, domain agreement of dip AND pattern (50 kb, Hi-C insulation domains)'),
             row(['set', 'covered', 'CTCF sites', 'cohesin-b.', 'motif', 'cis edges', 'both ends', 'conv.',
                  'tandem', 'div.', 'VI (rnd)', 'Rand (rnd)']), r'\midrule']
    for s in qc.set:
        q = qc[qc.set == s].iloc[0]
        cv = cov[(cov.cell == s) & (cov.resolution_kb == 50)].covered_fraction
        o = lo[lo.set == s].set_index('class')
        d = dom[(dom.cell == s) & (dom.rule == 'dip_and_ctcf_pattern') & (dom.resolution_kb == 50) &
                (dom.reference == 'hic_insulation')]
        lines.append(row([tex(s), f(cv.iloc[0]) if len(cv) else '--', int(q.n_sites), int(q.n_cohesin_bound),
                          f(q.motif_fraction_cohesin_bound), int(q.n_edges_cis), int(o.n_edges_both_ends.iloc[0]),
                          f(o.fraction['convergent']),
                          f(o.fraction['tandem_forward'] + o.fraction['tandem_reverse']), f(o.fraction['divergent']),
                          '{} ({})'.format(f(d.vi_norm.iloc[0], 3), f(d.vi_norm_random.iloc[0], 3)) if len(d) else '--',
                          '{} ({})'.format(f(d.local_rand_1mb.iloc[0], 3), f(d.local_rand_1mb_random.iloc[0], 3))
                          if len(d) else '--']))
    lines.append(row(['', '', '', '', '', '', 'independence:', '0.25', '0.50', '0.25']))
    return lines


def main():
    base, out = sys.argv[1], sys.argv[2]
    tables = os.path.join(base, 'tables')
    ev = pd.read_csv(os.path.join(tables, 'eval_rules.tsv'), sep='\t')
    grid = pd.read_csv(os.path.join(tables, 'eval_rules_grid.tsv'), sep='\t')
    delta = pd.read_csv(os.path.join(tables, 'eval_rules_delta.tsv'), sep='\t')
    global NCOL
    parts = [(14, lambda: part_metrics(ev, 'all', 'A. Rules vs Hi-C insulation boundaries of the same cell type '
                                       '(all bins; P, R, F1 at tolerance 1 bin)')),
             (14, lambda: part_metrics(ev, 'covered', 'B. As A, covered bins only (>= 5 expected spanning edges)')),
             (12, lambda: part_delta(delta)), (11, lambda: part_refs(ev)), (12, lambda: part_grid(grid)),
             (12, lambda: part_inputs(tables))]
    blocks = []
    for ncol, make in parts:
        NCOL = ncol
        # one longtable per part (own column widths); each part starts with a rule and its heading
        blocks.append([r'\begin{longtable}{' + 'l' * ncol + '}', r'\toprule'] + make()[1:] + [r'\bottomrule'])
    caption = (r'\caption{Domain boundaries called by rules from three evidence sources, evaluated against Hi-C. '
               r'Dip: local minimum of the capture insulation (log2 observed / expected number of graph edges '
               r'spanning the bin, expectation from the same edges mirrored at their bait) with prominence '
               r'$\geq$ 0.3 and value $\leq$ 0. CTCF: a cohesin-bound CTCF summit (enrichment $\geq$ 10) in the bin '
               r'(dip AND CTCF: within 1 bin of the dip). CTCF pattern: reverse CTCF motif (JASPAR MA0139.2) sites '
               r'in the bin or the bin to the left and '
               r'forward sites in the bin or the bin to the right. Loop anchor: an end of a graph edge with '
               r'convergent cohesin-bound CTCF within 1 bin. One call per run of consecutive flagged bins. Thresholds '
               r'were fixed before the evaluation; no weights were fitted. AUPRC: per-bin score of the rule (AND / OR: '
               r'minimum / maximum of the percentile ranks of the two scores). VI: normalised variation of '
               r'information between called and Hi-C domains; local Rand: agreement on same-domain status of bin '
               r'pairs up to 1 Mb apart; random: the same number of boundaries placed at random.}')
    # longtable (LaTeX package longtable): the table runs over several pages; caption below the last row
    lines = []
    for k, block in enumerate(blocks):
        lines += block
        if k == len(blocks) - 1:
            lines.append(caption + r'\label{tab:s35}\\')
        lines.append(r'\end{longtable}')
        if k < len(blocks) - 1:
            lines.append(r'\addtocounter{table}{-1}')     # uncaptioned parts do not take a table number
    body = '\n'.join([r'\begingroup', r'\scriptsize', r'\setlength{\tabcolsep}{3pt}'] + lines + [r'\endgroup'])
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'table_s35.tex'), 'w') as h:
        h.write(body + '\n')
    print(os.path.join(out, 'table_s35.tex'))


if __name__ == '__main__':
    main()
