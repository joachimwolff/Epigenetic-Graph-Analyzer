"""Explanation cascade of the differentially expressed genes, one PDF per comparison (3.6 x 2.4 in plus legend), no panel letters.

Reads <prefix>_explanation_cascade.tsv of egaDifferential (analysis/16_hops_n.sh) and draws one horizontal bar per
level. Gene (step 0) and step 1 are broken down by explanation type: marks, own interactions (gene only), both
(gene only), a chromatin state transition from or to the bivalent promoter / poised enhancer state concordant with
the expression change ('bivalent or poised state resolved'), and a promoter or contacting fragment poised in both
conditions ('poised in both conditions'); the indirect steps (2, 3, ...) and the unexplained genes follow. Bars give the number of
baited promoter fragments; the hatched part of a mark bar at a step are the genes beyond the mirrored loci with a
per-gene permutation P < 0.05. The number of genes and, where a step curve was given (--proximitySupport), the
probability that graph neighbours at that step are close in the nucleus are written next to each bar. Comparisons
whose marks allow no state calls (no repressive mark) say so.

Usage: python figure_explanation_cascade.py <results dir> <figure dir> [name=prefix ...]
       (default: the four comparisons of analysis/16_hops_n.sh under <results dir>/hops_n)
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

COMPARISONS = (('hESC_vs_hNEC', 'hops_n/hESC_vs_hNEC'), ('Md_vs_PA2', 'hops_n/Md_vs_PA2'),
               ('Kasumi1_siMM_vs_siRE', 'hops_n/siMM_vs_siRE'), ('AML_t821_vs_FLT3ITD', 'hops_n/t821_vs_FLT3ITD'))
# categorical slots 1-5 of the reference palette, stacked in slot order
TYPES = (('marks', '#2a78d6', 'marks'), ('interactions', '#eb6834', 'own interactions'),
         ('both', '#1baf7a', 'marks and interactions'), ('state', '#eda100', 'bivalent or poised state resolved'),
         ('poised_unresolved', '#e87ba4', 'poised in both conditions'))
ROWS = {'gene': {'marks': 'gene_marks_only', 'interactions': 'gene_interactions_only', 'both': 'gene_both',
                 'state': 'gene_state', 'poised_unresolved': 'gene_poised_unresolved'},
        'step1': {'marks': 'step1_marks', 'state': 'step1_state', 'poised_unresolved': 'step1_poised_unresolved'}}
UNEXPLAINED = '#e3e2dc'
INK = '#1f1f1e'
MUTED = '#6b6a66'


def label(level):
    if level == 'gene':
        return 'gene (step 0)'
    if level == 'unexplained':
        return 'unexplained'
    k = int(level[4:])
    return '{} step{}'.format(k, '' if k == 1 else 's')


def value(table, row, column, default=0):
    if row not in table.index or column not in table.columns:
        return default
    v = table.loc[row, column]
    return default if pd.isna(v) else v


def draw(table, path):
    table = table.set_index('level')
    steps = [lv for lv in table.index if lv.startswith('step') and lv[4:].isdigit()]
    levels = ['gene'] + steps + ['unexplained']
    not_callable = 'note' in table.columns and (table['note'] == 'states not callable').any()
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    largest = max(table.loc[levels, 'genes'].max(), 1)
    used = set()
    for y, level in enumerate(levels):
        n = int(table.loc[level, 'genes'])
        if level in ROWS:
            left = 0
            for kind, colour, _ in TYPES:
                if kind not in ROWS[level]:
                    continue
                width = int(value(table, ROWS[level][kind], 'genes'))
                if width:
                    used.add(kind)
                ax.barh(y, width, left=left, color=colour, height=0.62, lw=0)
                if kind == 'marks' and level != 'gene':
                    dark = int(value(table, level, 'beyond_mirror_p_lt_0.05'))
                    ax.barh(y, dark, left=left, color='none', height=0.62, lw=0, hatch='//////', edgecolor='white')
                left += width
        elif level == 'unexplained':
            ax.barh(y, n, color=UNEXPLAINED, height=0.62, lw=0.4, edgecolor='#8c8c88')
        else:
            dark = int(value(table, level, 'beyond_mirror_p_lt_0.05'))
            ax.barh(y, n, color=TYPES[0][1], height=0.62, lw=0)
            ax.barh(y, dark, color='none', height=0.62, lw=0, hatch='//////', edgecolor='white')
        text = '{:,}'.format(n)
        if level.startswith('step'):
            dark = value(table, level, 'beyond_mirror_p_lt_0.05', float('nan'))
            if not pd.isna(dark):
                text += ' ({:,})'.format(int(dark))
            support = value(table, level, 'proximity_support', float('nan'))
            interval = value(table, level, 'proximity_ci', '')
            if not pd.isna(support):
                text += ', 3D {:.2f}'.format(support)
                if isinstance(interval, str) and interval:
                    text += ' ({})'.format(interval)
        ax.text(n + 0.02 * largest, y, text, va='center', ha='left', fontsize=6, color=INK)
    ax.set_yticks(range(len(levels)))
    ax.set_yticklabels([label(lv) for lv in levels], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, largest * 1.8)
    ax.set_xlabel('differentially expressed promoter fragments', fontsize=7)
    ax.tick_params(axis='x', labelsize=6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    kinds = [t for t in TYPES if t[0] in used or (t[0] in ('marks',))]
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in kinds]
    names = [name for _, _, name in kinds]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=TYPES[0][1], hatch='//////', edgecolor='white', lw=0))
    names.append('marks beyond mirror, P < 0.05')
    ax.legend(handles, names, fontsize=5.5, frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.2), ncol=2)
    if not_callable:
        ax.text(0.99, 0.99, 'states not callable (no repressive mark)', transform=ax.transAxes, ha='right',
                va='top', fontsize=5.5, color=MUTED)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def main():
    results, figures = sys.argv[1:3]
    comparisons = [tuple(x.split('=', 1)) for x in sys.argv[3:]] or COMPARISONS
    os.makedirs(figures, exist_ok=True)
    for name, prefix in comparisons:
        source = os.path.join(results, prefix + '_explanation_cascade.tsv')
        if not os.path.exists(source):
            print('missing', source)
            continue
        target = os.path.join(figures, 'explanation_cascade_{}.pdf'.format(name))
        draw(pd.read_csv(source, sep='\t', keep_default_na=False, na_values=['']), target)
        print(target)


if __name__ == '__main__':
    main()
