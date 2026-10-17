"""Figures of the eQTL validation (analysis/18_eqtl.sh), one PDF per panel, no panel letters.

  eqtl_fraction_<graph>__<data>_<step>.pdf   fraction of pairs connected within 1, 2 or 3 steps per distance bin,
                                             eQTL pairs against matched controls, observed and mirrored graph,
                                             chromosome block-bootstrap 95% intervals
  eqtl_odds_<contrast>.pdf                   odds ratio eQTL vs controls (Mantel-Haenszel over distance bins) per
                                             graph and eQTL data set, observed and mirrored graph

Usage: python figure_eqtl.py <eqtl_validation dir> <output dir>
"""

import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

EQTL = '#1f5a96'
CONTROL = '#8c8c8c'
GRAPHS = {'GM_chicago': 'GM12878 pcHi-C (CHiCAGO)', 'GM_mifsud': 'GM12878 pcHi-C (Mifsud)',
          'VC_GM12878': 'GM12878 virtual capture', 'K562': 'K562 pcHi-C'}
DATA = {'geuvadis_susie': 'GEUVADIS LCL, SuSiE', 'gtex_lcl': 'GTEx LCL', 'gtex_wb': 'GTEx whole blood',
        'eqtlgen': 'eQTLGen blood'}
STEPS = {'step1': '1 step', 'within2': 'within 2 steps', 'within3': 'within 3 steps'}
CONTRASTS = {'direct': 'direct edge', 'two_step_without_direct': '2-step path, pairs without direct edge',
             'three_step_without_1_2': '3-step path, pairs without 1- or 2-step path',
             'two_step_any': '2-step path (any)', 'three_step_any': '3-step path (any)'}


def fraction_panel(fractions, step, path):
    bins = [b for b in fractions['distance_bin'].unique() if b != 'all']
    x = np.arange(len(bins))
    fig, ax = plt.subplots(figsize=(2.8, 2.6))
    for graph, marker, fill in (('observed', 'o', True), ('mirrored', 's', False)):
        for group, color, offset in (('eqtl', EQTL, -0.08), ('control', CONTROL, 0.08)):
            t = fractions[(fractions['graph'] == graph) & (fractions['group'] == group)].set_index('distance_bin')
            t = t.reindex(bins)
            y = t[step].to_numpy() * 100
            err = np.vstack([y - t[step + '_low'].to_numpy() * 100, t[step + '_high'].to_numpy() * 100 - y])
            ax.errorbar(x + offset + (0.02 if graph == 'mirrored' else 0), y, yerr=err, color=color, marker=marker,
                        markersize=3.5, lw=0.9 if fill else 0.7, ls='-' if fill else ':', capsize=1.5,
                        markerfacecolor=color if fill else 'white',
                        label='{} ({} graph)'.format('eQTL pairs' if group == 'eqtl' else 'controls', graph))
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace('kb', '') for b in bins], rotation=40, ha='right')
    ax.set_xlabel('Distance variant - TSS (kb)')
    ax.set_ylabel('Pairs connected {} (%)'.format(STEPS[step]))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=5.5, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2,
              handlelength=2.2, columnspacing=1.0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def odds_panel(odds, contrast, path):
    t = odds[(odds['contrast'] == contrast) & (odds['control'] == 'control') & (odds['distance_bin'] == 'all_mh')]
    rows = [(g, d) for g in GRAPHS for d in DATA if ((t['graph_name'] == g) & (t['data'] == d)).any()]
    fig, ax = plt.subplots(figsize=(4.4, 0.22 * len(rows) + 1.0))
    for i, (g, d) in enumerate(rows):
        for graph, color, offset, fill in (('observed', EQTL, -0.15, True), ('mirrored', CONTROL, 0.15, False)):
            r = t[(t['graph_name'] == g) & (t['data'] == d) & (t['graph'] == graph)]
            if r.empty:
                continue
            r = r.iloc[0]
            y = len(rows) - 1 - i + offset
            ax.plot([r['or_low'], r['or_high']], [y, y], color=color, lw=0.9)
            ax.plot(r['odds_ratio'], y, 'o', color=color, markersize=3.5,
                    markerfacecolor=color if fill else 'white', label='{} graph'.format(graph) if i == 0 else None)
    ax.axvline(1, color='black', lw=0.5, ls='--')
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.25, 0.5, 1, 2, 4]))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(['0.25', '0.5', '1', '2', '4']))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(['{} / {}'.format(GRAPHS[g], DATA[d]) for g, d in rows][::-1], fontsize=5.5)
    ax.set_xlabel('Odds ratio, eQTL vs control pairs\n({})'.format(CONTRASTS[contrast]))
    ax.legend(frameon=False, fontsize=5.5, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    root, out = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    all_odds = []
    for path in sorted(glob.glob(os.path.join(root, 'runs', '*_fractions.tsv'))):
        name = os.path.basename(path)[:-len('_fractions.tsv')]
        graph, data = name.split('__')
        fractions = pd.read_csv(path, sep='\t')
        for step in STEPS:
            fraction_panel(fractions, step, os.path.join(out, 'eqtl_fraction_{}_{}.pdf'.format(name, step)))
        odds = pd.read_csv(path.replace('_fractions.tsv', '_odds.tsv'), sep='\t')
        all_odds.append(odds.assign(graph_name=graph, data=data))
    odds = pd.concat(all_odds, ignore_index=True)
    for contrast in CONTRASTS:
        odds_panel(odds, contrast, os.path.join(out, 'eqtl_odds_{}.pdf'.format(contrast)))


if __name__ == '__main__':
    main()
