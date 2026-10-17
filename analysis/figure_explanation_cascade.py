"""Explanation cascade of the differentially expressed genes, one PDF per comparison (3.4 x 2.6 in).

Reads <prefix>_explanation_cascade.tsv of egaDifferential (analysis/16_hops_n.sh) and draws one horizontal bar per
level: explained at the gene (by its marks, its own interactions or both), first explained 1, 2, ... steps away
in the graph, and unexplained. Bars give the number of baited promoter fragments; the darker part of a step bar
are the genes beyond the mirrored loci with a per-gene permutation P < 0.05; the number of genes and, where a
step curve was given (--proximitySupport), the probability that graph neighbours at that step are close in the
nucleus are written next to each bar.

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
GENE_COLOURS = {'gene_marks_only': '#4d4d4a', 'gene_both': '#7a7a76', 'gene_interactions_only': '#b0b0aa'}
STEP_LIGHT = '#9ecae1'
STEP_DARK = '#2171b5'
UNEXPLAINED = '#e3e2dc'
INK = '#1f1f1e'


def label(level):
    if level == 'gene':
        return 'gene'
    if level == 'unexplained':
        return 'unexplained'
    k = int(level[4:])
    return '{} step{}'.format(k, '' if k == 1 else 's')


def draw(table, path):
    table = table.set_index('level')
    levels = ['gene'] + [lv for lv in table.index if lv.startswith('step')] + ['unexplained']
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    largest = max(table.loc[levels, 'genes'].max(), 1)
    for y, level in enumerate(levels):
        n = int(table.loc[level, 'genes'])
        if level == 'gene':
            left = 0
            for part, colour in GENE_COLOURS.items():
                width = int(table.loc[part, 'genes']) if part in table.index else 0
                ax.barh(y, width, left=left, color=colour, height=0.62, lw=0)
                left += width
        elif level == 'unexplained':
            ax.barh(y, n, color=UNEXPLAINED, height=0.62, lw=0.4, edgecolor='#8c8c88')
        else:
            dark = table.loc[level].get('beyond_mirror_p_lt_0.05', 0)
            dark = 0 if pd.isna(dark) else int(dark)
            ax.barh(y, n, color=STEP_LIGHT, height=0.62, lw=0)
            ax.barh(y, dark, color=STEP_DARK, height=0.62, lw=0)
        text = '{:,}'.format(n)
        if level.startswith('step'):
            dark = table.loc[level].get('beyond_mirror_p_lt_0.05', float('nan'))
            if not pd.isna(dark):
                text += ' ({:,})'.format(int(dark))
            support = table.loc[level].get('proximity_support', float('nan'))
            interval = table.loc[level].get('proximity_ci', '')
            if not pd.isna(support):
                text += ', 3D {:.2f}'.format(support)
                if isinstance(interval, str) and interval:
                    text += ' ({})'.format(interval)
        ax.text(n + 0.02 * largest, y, text, va='center', ha='left', fontsize=6, color=INK)
    ax.set_yticks(range(len(levels)))
    ax.set_yticklabels([label(lv) for lv in levels], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, largest * 1.75)
    ax.set_xlabel('differentially expressed promoter fragments', fontsize=7)
    ax.tick_params(axis='x', labelsize=6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in GENE_COLOURS.values()] + \
        [plt.Rectangle((0, 0), 1, 1, color=STEP_DARK)]
    ax.legend(handles, ['gene: marks', 'gene: marks and interactions', 'gene: interactions',
                        'beyond mirror, P < 0.05'], fontsize=5.5, frameon=False, loc='lower right')
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
