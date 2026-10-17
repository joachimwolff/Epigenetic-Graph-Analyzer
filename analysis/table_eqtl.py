"""Supplementary Table S30: eQTL variant-gene pairs against distance-matched controls at 1, 2 and 3 graph steps
(analysis/18_eqtl.sh), and a tab-separated summary of the same numbers.

Per graph and eQTL data set, observed and mirrored graph: eQTL pairs and controls, the fraction with a direct edge,
Mantel-Haenszel odds ratios over distance bins (chromosome block-bootstrap 95% intervals) for a direct edge, a
2-step path among pairs without a direct edge and a 3-step path among pairs without a 1- or 2-step path, and the
likelihood-ratio P of the 2- and 3-step terms in the logistic model (log distance, direct, 2-step, 3-step) on all
pairs and on pairs without a direct edge.

Usage: python table_eqtl.py <eqtl_validation dir> <output dir>
"""

import os
import sys

import pandas as pd

GRAPHS = {'GM_chicago': 'GM12878 (CHiCAGO)', 'GM_mifsud': 'GM12878 (Mifsud)', 'VC_GM12878': 'GM12878 (virtual)',
          'K562': 'K562'}
DATA = {'geuvadis_susie': 'GEUVADIS LCL SuSiE', 'gtex_lcl': 'GTEx LCL', 'gtex_wb': 'GTEx blood',
        'eqtlgen': 'eQTLGen blood'}


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')


def num(x, digits=2):
    if pd.isna(x):
        return '--'
    return '{:.{}f}'.format(x, digits)


def pval(p):
    if pd.isna(p):
        return '--'
    if p == 0:
        return '$<10^{-300}$'
    if p >= 0.01:
        return '{:.2f}'.format(p)
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))


def odds_cell(r):
    if r is None:
        return '--'
    return '{} ({}--{})'.format(num(r['odds_ratio']), num(r['or_low']), num(r['or_high']))


def collect(root):
    rows = []
    for graph in GRAPHS:
        for data in DATA:
            prefix = os.path.join(root, 'runs', '{}__{}'.format(graph, data))
            if not os.path.exists(prefix + '_models.tsv'):
                continue
            odds = pd.read_csv(prefix + '_odds.tsv', sep='\t')
            models = pd.read_csv(prefix + '_models.tsv', sep='\t')
            mapping = pd.read_csv(prefix + '_mapping.tsv', sep='\t').set_index('step')['n']
            for version in ('observed', 'mirrored'):
                o = odds[(odds['graph'] == version) & (odds['control'] == 'control')]

                def pick(contrast, bin_label='all_mh'):
                    r = o[(o['contrast'] == contrast) & (o['distance_bin'] == bin_label)]
                    return r.iloc[0] if len(r) else None

                def lrt(model, term):
                    r = models[(models['graph'] == version) & (models['model'] == model) & (models['term'] == term)]
                    return r.iloc[0]['lr_p'] if len(r) else float('nan')

                def coef(model, term, column):
                    r = models[(models['graph'] == version) & (models['model'] == model) & (models['term'] == term)]
                    return r.iloc[0][column] if len(r) else float('nan')
                direct = pick('direct')
                two = pick('two_step_without_direct')
                three = pick('three_step_without_1_2')
                rows.append({
                    'graph': graph, 'data': data, 'version': version,
                    'eqtl_pairs': int(direct['eqtl_pairs']), 'control_pairs': int(direct['control_pairs']),
                    'excluded_variant_in_promoter_fragment': int(mapping.get('variant_in_promoter_fragment', 0)),
                    'excluded_promoter_not_baited': int(mapping.get('gene_promoter_not_baited', 0)),
                    'eqtl_direct': direct['eqtl_fraction'], 'control_direct': direct['control_fraction'],
                    'or_direct': direct['odds_ratio'], 'or_direct_low': direct['or_low'],
                    'or_direct_high': direct['or_high'],
                    'eqtl_two_without_direct': two['eqtl_fraction'],
                    'control_two_without_direct': two['control_fraction'],
                    'or_two_without_direct': two['odds_ratio'], 'or_two_low': two['or_low'],
                    'or_two_high': two['or_high'],
                    'eqtl_three_without_1_2': three['eqtl_fraction'],
                    'control_three_without_1_2': three['control_fraction'],
                    'or_three_without_1_2': three['odds_ratio'], 'or_three_low': three['or_low'],
                    'or_three_high': three['or_high'],
                    'model_or_two': coef('all_pairs', 'two_step', 'odds_ratio'),
                    'model_or_two_low': coef('all_pairs', 'two_step', 'or_low'),
                    'model_or_two_high': coef('all_pairs', 'two_step', 'or_high'),
                    'model_or_three': coef('all_pairs', 'three_step', 'odds_ratio'),
                    'model_or_three_low': coef('all_pairs', 'three_step', 'or_low'),
                    'model_or_three_high': coef('all_pairs', 'three_step', 'or_high'),
                    'lrt_p_direct': lrt('all_pairs', 'direct'),
                    'lrt_p_two': lrt('all_pairs', 'two_step'), 'lrt_p_three': lrt('all_pairs', 'three_step'),
                    'lrt_p_two_adjusted': lrt('all_pairs_adjusted', 'two_step'),
                    'lrt_p_three_adjusted': lrt('all_pairs_adjusted', 'three_step'),
                    'lrt_p_two_without_direct': lrt('without_direct', 'two_step'),
                    'lrt_p_three_without_direct': lrt('without_direct', 'three_step'),
                    'lrt_p_two_conditional': lrt('conditional_all_pairs', 'two_step'),
                    'lrt_p_three_conditional': lrt('conditional_all_pairs', 'three_step'),
                    'lrt_p_two_conditional_without_direct': lrt('conditional_without_direct', 'two_step'),
                    'lrt_p_three_conditional_without_direct': lrt('conditional_without_direct', 'three_step'),
                })
    return pd.DataFrame(rows)


def main():
    root, out = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    summary = collect(root)
    summary.to_csv(os.path.join(root, 'eqtl_summary.tsv'), sep='\t', index=False)
    lines = [r'\begin{table}[!h]', r'\centering\scriptsize', r'\begin{adjustbox}{max width=\textwidth}',
             r'\begin{tabular}{@{}lllrrlllllll@{}}', r'\toprule',
             r'Graph & eQTL set & Graph version & eQTL pairs & Controls & Direct, eQTL / control (\%) & '
             r'OR direct & OR 2-step, no direct & OR 3-step, no 1--2 step & $P_\mathrm{LR}$ 2-step & '
             r'$P_\mathrm{LR}$ 3-step & $P_\mathrm{LR}$ 2-step, no direct \\', r'\midrule']
    previous = None
    for _, r in summary.iterrows():
        if previous is not None and r['graph'] != previous:
            lines.append(r'\midrule')
        previous = r['graph']
        lines.append(' & '.join([
            tex(GRAPHS[r['graph']]), tex(DATA[r['data']]), r['version'], '{:,}'.format(r['eqtl_pairs']),
            '{:,}'.format(r['control_pairs']),
            '{} / {}'.format(num(100 * r['eqtl_direct'], 1), num(100 * r['control_direct'], 1)),
            odds_cell({'odds_ratio': r['or_direct'], 'or_low': r['or_direct_low'], 'or_high': r['or_direct_high']}),
            odds_cell({'odds_ratio': r['or_two_without_direct'], 'or_low': r['or_two_low'],
                       'or_high': r['or_two_high']}),
            odds_cell({'odds_ratio': r['or_three_without_1_2'], 'or_low': r['or_three_low'],
                       'or_high': r['or_three_high']}),
            pval(r['lrt_p_two']), pval(r['lrt_p_three']), pval(r['lrt_p_two_without_direct'])]) + r' \\')
    caption = (r'eQTL variant-gene pairs against distance-matched non-eQTL pairs in the promoter interaction graphs. '
               r'Controls per eQTL pair: the same variant fragment with baited promoters of other eGenes, and the same '
               r'gene with other fragments carrying a significant variant, within a factor $10^{0.1}$ of the eQTL '
               r'distance (up to 1~Mb) and never a significant pair of the data set. OR: Mantel-Haenszel odds ratio over '
               r'distance bins with chromosome block-bootstrap 95\% interval (1,000 resamples). $P_\mathrm{LR}$: '
               r'likelihood-ratio test of the 2- or 3-step term in the logistic model of eQTL status on log distance, '
               r'direct edge, 2-step path and 3-step path, on all pairs or on pairs without a direct edge. Mirrored: '
               r'every edge reflected around its bait, keeping bait degree and edge length.')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{adjustbox}', r'\caption{' + caption + '}',
              r'\label{tab:s30}', r'\end{table}']
    with open(os.path.join(out, 'table_s30.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
