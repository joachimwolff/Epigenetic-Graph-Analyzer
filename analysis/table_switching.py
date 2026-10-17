"""Supplementary Table S29: enhancer switching and promoter competition per comparison (egaSwitching,
analysis/18_enhancer_switching.sh).

table_s29.tex: per comparison and event set (switch events with a non-promoter element, with any fragment, hub
switches within two steps, competition of promoters sharing elements) the number of events, the number with the
expression of both promoters, the concordant fraction (lost promoter down and gained promoter up; for competition the
promoter gaining more of the shared elements up and the other down), the discordant fraction, the means of the two
null models with their one-sided permutation P values, and the concordant fraction of promoter pairs that keep their
contacts with the same element (no switch).
table_s29_events.tex: the three highest-ranked switch events with a non-promoter element per comparison.

Usage: python table_switching.py <enhancer_switching dir> <output dir>
"""

import os
import sys

import pandas as pd

COMPARISONS = [('hESC_vs_hNEC', 'hESC / hNEC'), ('Md_vs_PA2', 'Md / PA2'), ('siMM_vs_siRE', 'Kasumi-1 siMM / siRE'),
               ('t821_vs_FLT3ITD', 'AML t(8;21) / FLT3-ITD')]
SETS = [('non-promoter element', 'switch, element'), ('any fragment', 'switch, any fragment'),
        ('hub switch (two steps)', 'hub switch'), ('competition', 'competition')]


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('#', r'\#')


def pct(x):
    return '--' if pd.isna(x) else '{:.1f}'.format(100 * x)


def pval(p):
    if pd.isna(p):
        return '--'
    if p >= 0.01:
        return '{:.2f}'.format(p)
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))


def integer(x):
    return '--' if pd.isna(x) else '{:,}'.format(int(x))


def num(x, digits=1):
    return '--' if pd.isna(x) else '{:+.{}f}'.format(x, digits).replace('-', '$-$')


def summary_table(src, out):
    rows = []
    for key, label in COMPARISONS:
        path = os.path.join(src, key + '_switch_summary.tsv')
        if not os.path.exists(path):
            continue
        summary = pd.read_csv(path, sep='\t').set_index('set')
        first = True
        for name, set_label in SETS:
            if name not in summary.index:
                continue
            r = summary.loc[name]
            rows.append(' & '.join([
                label if first else '', set_label, integer(r.get('n_events')), integer(r.get('n_tested')),
                pct(r.get('concordant')), pct(r.get('discordant')), pct(r.get('concordant_null_permutation')),
                pval(r.get('concordant_p_permutation')), pct(r.get('concordant_null_rotation')),
                pval(r.get('concordant_p_rotation')), pct(r.get('stable_pairs_concordant'))]) + r' \\')
            first = False
        rows.append(r'\addlinespace')
    if rows and rows[-1] == r'\addlinespace':
        rows.pop()
    header = (r'Comparison & Set & Events & Tested & Conc.\ (\%) & Disc.\ (\%) & Perm.\ null (\%) & $P_\mathrm{perm}$ '
              r'& Rot.\ null (\%) & $P_\mathrm{rot}$ & No switch (\%)')
    caption = (r'Enhancer switching and promoter competition. Switch: a fragment loses its contact with promoter A '
               r'(status lost or $\log_2$ change $\leq -1$) and gains a contact with promoter B (gained or '
               r'$\geq +1$); element: non-promoter fragments, any fragment: including promoters. Hub switch: '
               r'a non-promoter element leaves the hubs of the promoters it reaches within two steps only in the first '
               r'condition and joins those it reaches only in the second (one row per element, mean expression change '
               r'of each set). Competition: promoter pairs that share non-promoter elements with a different net '
               r'change of the shared contacts. Tested: events with the expression of both promoters. Concordant: A '
               r'down and B up (competition: the promoter gaining more of the shared elements up, the other down); '
               r'discordant: the reverse. Null models, 10,000 draws each: permutation of the expression changes among '
               r'the promoters of the events within strata of expression level and contact distance, and rotation of '
               r'the expression changes along each chromosome; one-sided permutation $P$. No switch: concordant '
               r'fraction of ordered promoter pairs that keep their contacts with the same element.')
    lines = [r'\begin{table}[!h]', r'\centering\scriptsize\setlength{\tabcolsep}{3pt}',
             r'\begin{adjustbox}{max width=\textwidth}', r'\begin{tabular}{@{}llrrrrrrrrr@{}}', r'\toprule',
             header + r' \\', r'\midrule'] + rows + [
             r'\bottomrule', r'\end{tabular}', r'\end{adjustbox}', r'\caption{' + caption + '}', r'\label{tab:s29}',
             r'\end{table}']
    with open(os.path.join(out, 'table_s29.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def short_genes(text, limit=2):
    genes = [g for g in str(text).split(',') if g and g != 'nan']
    return ','.join(genes[:limit]) + (',...' if len(genes) > limit else '')


def events_table(src, out, top=3):
    rows = []
    for key, label in COMPARISONS:
        path = os.path.join(src, key + '_switch_events.tsv')
        if not os.path.exists(path):
            continue
        events = pd.read_csv(path, sep='\t')
        events = events[(events['element_type'] == 'non-promoter element') & events['expression_available']]
        events = events.drop_duplicates('element').head(top)
        first = True
        for _, e in events.iterrows():
            rows.append(' & '.join([
                label if first else '', tex(e['element']),
                r'\textit{' + tex(short_genes(e['genes_lost'])) + '}', num(e['rna_log2FC_lost']),
                '{:.0f}'.format(e['distance_lost'] / 1000), num(e['log2FC_lost']),
                r'\textit{' + tex(short_genes(e['genes_gained'])) + '}', num(e['rna_log2FC_gained']),
                '{:.0f}'.format(e['distance_gained'] / 1000), num(e['log2FC_gained']),
                tex(e['element_chromatin_changes']) if isinstance(e['element_chromatin_changes'], str) else '--',
                '{:.2f}'.format(e['evidence'])]) + r' \\')
            first = False
    header = (r'Comparison & Element E & A & RNA & kb & Contact & B & RNA & kb & Contact & Chromatin at E & Evidence')
    caption = (r'Highest-ranked switch events with a non-promoter element (one per element). A: promoter losing the '
               r'contact, B: promoter gaining it; RNA: $\log_2$ expression change; kb: distance to E; contact: '
               r'$\log_2$ change of the contact; chromatin at E: marks changing by at least 2-fold. Evidence: contact '
               r'change, expression change in the expected direction and chromatin change at E (egaSwitching).')
    lines = [r'\begin{table}[!h]', r'\centering\scriptsize\setlength{\tabcolsep}{3pt}',
             r'\begin{adjustbox}{max width=\textwidth}', r'\begin{tabular}{@{}llllrrlrrrlr@{}}', r'\toprule',
             header + r' \\', r'\midrule'] + rows + [
             r'\bottomrule', r'\end{tabular}', r'\end{adjustbox}', r'\caption{' + caption + '}',
             r'\label{tab:s29events}', r'\end{table}']
    with open(os.path.join(out, 'table_s29_events.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    summary_table(src, out)
    events_table(src, out)


if __name__ == '__main__':
    main()
