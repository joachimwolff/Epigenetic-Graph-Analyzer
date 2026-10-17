"""Supplementary Table S35: domain boundaries from capture contacts, CTCF / cohesin and motif orientation
(analysis/27_boundaries.sh).

Part A: cross-cell boundary calling against Hi-C insulation boundaries (all bins), per resolution and method:
precision, recall and F1 (tolerance 1 bin) and AUPRC, mean over the two training cell types per test cell type;
positive rate = fraction of bins holding a reference boundary. Part B: the same restricted to covered bins
(>= 5 expected spanning edges). Part C: gain of the contacts over CTCF (cross-cell and chromosome split, all bins):
delta AUPRC and delta F1 with 95 % chromosome block-bootstrap intervals. Part D: combined and CTCF + orientation
against the other Hi-C references (cross-cell, 50 kb, mean over training cells). Part E: coverage, CTCF sites and
loop orientation. Part F: keratinocyte time course.

Usage: python table_boundaries.py <boundaries dir> <output dir>
"""

import os
import sys

import pandas as pd

METHODS = ['density', 'capture', 'ctcf_cohesin', 'ctcf_orientation', 'ctcf_orientation_density', 'combined']
LABEL = {'density': 'bait density', 'capture': 'density + capture insulation', 'ctcf_cohesin': 'CTCF + cohesin',
         'ctcf_orientation': 'CTCF + cohesin + orientation', 'ctcf_orientation_density': '... + bait density',
         'combined': 'combined', 'capture_minima': 'capture minima (no training)'}
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
    return r'\multicolumn{%d}{l}{\textbf{%s}}\\' % (NCOL, text)


def part_metrics(ev, region, title):
    sub = ev[(ev.design == 'cross_cell') & (ev.reference == 'hic_insulation') & (ev.region == region)]
    uns = ev[(ev.design == 'unsupervised') & (ev.reference == 'hic_insulation') & (ev.region == region)]
    lines = [r'\midrule', heading(title),
             row(['', ''] + [r'\multicolumn{4}{c}{%s}' % c for c in CELLS]),
             row(['bin', 'method'] + ['P', 'R', 'F1', 'AUPRC'] * 3), r'\midrule']
    for res in (10, 25, 50):
        for i, method in enumerate(METHODS + ['capture_minima']):
            src = uns if method == 'capture_minima' else sub
            vals = []
            for c in CELLS:
                g = src[(src.resolution_kb == res) & (src.method == method) & (src.test == c)]
                vals += [f(g.precision_tol1.mean()), f(g.recall_tol1.mean()), f(g.f1_tol1.mean()),
                         f(g.auprc.mean(), 3)]
            lines.append(row(['{} kb'.format(res) if i == 0 else '', tex(LABEL[method])] + vals))
        pr = [sub[(sub.resolution_kb == res) & (sub.test == c)].positive_rate.mean() for c in CELLS]
        lines.append(row(['', 'positive rate'] + sum([['', '', '', f(p, 3)] for p in pr], [])))
    return lines


def part_delta(delta):
    lines = [r'\midrule', heading('C. Contacts beyond CTCF: delta AUPRC / delta F1 (tolerance 1), reference Hi-C '
                                  'insulation, all bins and covered bins, 95\\% chromosome bootstrap interval'),
             row(['comparison', 'bin', 'design', 'train', 'test', r'$\Delta$AUPRC all', '95\\% CI',
                  r'$\Delta$F1 all', '95\\% CI', r'$\Delta$AUPRC covered', '95\\% CI', r'$\Delta$F1 covered',
                  '95\\% CI']), r'\midrule']
    d = delta[delta.reference == 'hic_insulation']
    comps = [('combined', 'ctcf_orientation'), ('combined', 'ctcf_orientation_density'),
             ('ctcf_orientation_density', 'ctcf_orientation'), ('capture', 'density')]
    for big, small in comps:
        first = True
        for res in (10, 25, 50):
            sel = d[(d.full == big) & (d.reduced == small) & (d.resolution_kb == res)]
            for (design, train, test), g in sel.groupby(['design', 'train', 'test'], sort=False):
                a = g[g.region == 'all'].iloc[0]
                c = g[g.region == 'covered'].iloc[0]
                lines.append(row([tex('{} vs {}'.format(LABEL[big], LABEL[small])) if first else '',
                                  '{} kb'.format(res), 'split' if design == 'chromosome_split' else 'cross',
                                  train, test,
                                  f(a.delta_auprc, 3), '{}, {}'.format(f(a.delta_auprc_lo, 3), f(a.delta_auprc_hi, 3)),
                                  f(a.delta_f1_tol1, 3), '{}, {}'.format(f(a.delta_f1_lo, 3), f(a.delta_f1_hi, 3)),
                                  f(c.delta_auprc, 3), '{}, {}'.format(f(c.delta_auprc_lo, 3), f(c.delta_auprc_hi, 3)),
                                  f(c.delta_f1_tol1, 3), '{}, {}'.format(f(c.delta_f1_lo, 3), f(c.delta_f1_hi, 3))]))
                first = False
    return lines


def part_refs(ev):
    lines = [r'\midrule', heading('D. Other Hi-C references (cross-cell, mean over training cells, all bins)'),
             row(['reference', 'test', 'bin', 'method', 'P', 'R', 'F1', 'F1 tol. 2', 'AUPRC', 'positive rate']),
             r'\midrule']
    sub = ev[(ev.design == 'cross_cell') & (ev.region == 'all')]
    for ref in ('hic_insulation', 'hicFindTADs_50kb', 'rao2014_arrowhead', '4dn_H1_microc'):
        for c in CELLS:
            for res in (25, 50):
                for method in ('ctcf_orientation', 'combined'):
                    g = sub[(sub.reference == ref) & (sub.test == c) & (sub.resolution_kb == res) &
                            (sub.method == method)]
                    if g.empty:
                        continue
                    lines.append(row([tex(ref), c, '{} kb'.format(res), tex(LABEL[method]),
                                      f(g.precision_tol1.mean()), f(g.recall_tol1.mean()), f(g.f1_tol1.mean()),
                                      f(g.f1_tol2.mean()), f(g.auprc.mean(), 3), f(g.positive_rate.mean(), 3)]))
    return lines


def part_inputs(tables):
    cov = pd.read_csv(os.path.join(tables, 'coverage.tsv'), sep='\t')
    qc = pd.read_csv(os.path.join(tables, 'sites_qc.tsv'), sep='\t')
    lo = pd.read_csv(os.path.join(tables, 'loop_orientation.tsv'), sep='\t')
    dom = pd.read_csv(os.path.join(tables, 'eval_domains.tsv'), sep='\t')
    lines = [r'\midrule', heading('E. Inputs: coverage, CTCF sites, orientation of edges with cohesin-bound CTCF '
                                  'at both ends, domain agreement (combined, cross-cell, 50 kb, mean over training '
                                  'cells)'),
             row(['set', 'covered bins', 'CTCF sites', 'cohesin-bound', 'with motif', 'cis edges',
                  'edges CTCF both ends', 'convergent', 'tandem', 'divergent', 'VI (random)', 'local Rand (random)']),
             r'\midrule']
    for s in qc.set:
        q = qc[qc.set == s].iloc[0]
        cv = cov[(cov.cell == s) & (cov.resolution_kb == 50)].covered_fraction
        o = lo[lo.set == s].set_index('class')
        d = dom[(dom.design == 'cross_cell') & (dom.test == s) & (dom.method == 'combined') &
                (dom.resolution_kb == 50) & (dom.reference == 'hic_insulation')]
        lines.append(row([tex(s), f(cv.iloc[0]) if len(cv) else '--', int(q.n_sites), int(q.n_cohesin_bound),
                          f(q.motif_fraction_cohesin_bound), int(q.n_edges_cis), int(o.n_edges_both_ends.iloc[0]),
                          f(o.fraction['convergent']),
                          f(o.fraction['tandem_forward'] + o.fraction['tandem_reverse']), f(o.fraction['divergent']),
                          '{} ({})'.format(f(d.vi_norm.mean(), 3), f(d.vi_norm_random.mean(), 3)) if len(d) else '--',
                          '{} ({})'.format(f(d.local_rand_1mb.mean(), 3), f(d.local_rand_1mb_random.mean(), 3))
                          if len(d) else '--']))
    lines.append(row(['', '', '', '', '', '', 'independence:', '0.25', '0.50', '0.25']))
    return lines


def part_timecourse(tables):
    ev = pd.read_csv(os.path.join(tables, 'timecourse_eval.tsv'), sep='\t')
    corr = pd.read_csv(os.path.join(tables, 'timecourse_change_correlation.tsv'), sep='\t')
    ch = pd.read_csv(os.path.join(tables, 'timecourse_changes.tsv'), sep='\t')
    ex = pd.read_csv(os.path.join(tables, 'timecourse_expression.tsv'), sep='\t')
    lines = [r'\midrule', heading('F. Keratinocyte days 0, 3, 6 (hg19, 50 kb, chr2--22): calls against the Hi-C '
                                  'insulation boundaries of the same day (all bins)'),
             row(['day', 'method', 'P', 'R', 'F1', 'AUPRC', 'positive rate']), r'\midrule']
    for day in ('D0', 'D3', 'D6'):
        for i, m in enumerate(['capture_minima', 'ctcf_orientation_encode', 'combined_encode', 'ctcf_orientation_split',
                               'combined_split']):
            g = ev[(ev.day == day) & (ev.method == m) & (ev.region == 'all')]
            if g.empty:
                continue
            g = g.iloc[0]
            lines.append(row([day if i == 0 else '', tex(m.replace('_encode', ' (ENCODE models)')
                                                         .replace('_split', ' (chromosome split)')),
                              f(g.precision_tol1), f(g.recall_tol1), f(g.f1_tol1), f(g.auprc, 3),
                              f(g.positive_rate, 3)]))
    lines += [r'\midrule', row(['pair', 'track', 'bins', r'Spearman vs $\Delta$Hi-C insulation', '$P$']), r'\midrule']
    for _, r in corr[corr.region == 'covered_both'].iterrows():
        lines.append(row([r.pair.replace('_', '--'), tex(r.track), int(r.n_bins), f(r.spearman_vs_hic_change, 3),
                          '{:.1e}'.format(r.p)]))
    lines += [r'\midrule', row(['pair', 'caller', 'change', 'Hi-C changes', 'caller changes', 'matched', 'shifted',
                                '$P$']), r'\midrule']
    for _, r in ch.iterrows():
        lines.append(row([r.pair.replace('_', '--'), tex(r.caller), r.change, int(r.n_hic_changes),
                          int(r.n_caller_changes), f(r.fraction_hic_changes_matched), f(r.expected_shifted),
                          f(r.p_shift, 3)]))
    lines += [r'\midrule', row(['pair', 'change', 'boundaries', 'genes near', 'genes stable',
                                r'median $|$LFC$|$ near', 'stable', 'DE near', 'DE stable', '$P$']), r'\midrule']
    for _, r in ex.iterrows():
        lines.append(row([r.pair.replace('_', '--'), r.change, int(r.n_changed_boundaries), int(r.n_genes_changed),
                          int(r.n_genes_stable), f(r.median_abs_lfc_changed), f(r.median_abs_lfc_stable),
                          f(r.de_fraction_changed), f(r.de_fraction_stable),
                          '{:.2g}'.format(r.p_mannwhitney_greater)]))
    return lines


def main():
    base, out = sys.argv[1], sys.argv[2]
    tables = os.path.join(base, 'tables')
    ev = pd.read_csv(os.path.join(tables, 'eval_boundaries.tsv'), sep='\t')
    delta = pd.read_csv(os.path.join(tables, 'eval_delta.tsv'), sep='\t')
    lines = [r'\begin{tabular}{' + 'l' * NCOL + '}', r'\toprule']
    lines += part_metrics(ev, 'all', 'A. Cross-cell boundary calls vs Hi-C insulation boundaries (all bins; P, R, '
                                     'F1 at tolerance 1 bin; mean over the two training cell types)')[1:]
    lines += part_metrics(ev, 'covered', 'B. As A, covered bins only (>= 5 expected spanning edges)')
    lines += part_delta(delta) + part_refs(ev) + part_inputs(tables)
    if os.path.exists(os.path.join(tables, 'timecourse_eval.tsv')):
        lines += part_timecourse(tables)
    lines += [r'\bottomrule', r'\end{tabular}']
    caption = (r'\caption{Domain boundaries from capture contacts, CTCF and cohesin binding and CTCF motif '
               r'orientation, evaluated against Hi-C. Methods are logistic regressions on per-bin features z-scored '
               r'within each cell type, trained on the Hi-C insulation boundaries of one cell type and applied to '
               r'another (cross) or to the other half of the chromosomes (split). Bait density: number of baits and '
               r'expected spanning edges per bin, without the side to which the edges point. Capture insulation: '
               r'log2 observed / expected number of edges spanning the bin, expectation from the same edges '
               r'mirrored at their bait. Orientation: forward / reverse CTCF motif (JASPAR MA0139.2) at '
               r'cohesin-bound CTCF summits. VI: normalised variation of information between called and Hi-C '
               r'domains; local Rand: agreement on same-domain status of bin pairs up to 1 Mb apart; random: the '
               r'same number of boundaries placed at random.}')
    body = '\n'.join([r'\begin{table*}', r'\centering', r'\scriptsize', r'\resizebox{\textwidth}{!}{%'] + lines +
                     ['}', caption, r'\label{tab:s35}', r'\end{table*}'])
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'table_s35.tex'), 'w') as h:
        h.write(body + '\n')
    print(os.path.join(out, 'table_s35.tex'))


if __name__ == '__main__':
    main()
