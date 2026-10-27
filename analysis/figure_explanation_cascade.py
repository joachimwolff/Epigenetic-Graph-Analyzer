"""Explanation cascade of the differentially expressed genes, one PDF per comparison (3.6 x 2.4 in plus legend), no panel letters.

Reads <prefix>_explanation_cascade.tsv of egaDifferential (analysis/16_hops_n.sh) and draws one horizontal bar per
level. Gene (step 0) and step 1 are broken down by explanation type: marks, own interactions (gene only), both
(gene only), a chromatin state transition from or to the bivalent promoter / poised enhancer state concordant with
the expression change ('bivalent or poised state resolved'), and a promoter or contacting fragment poised in both
conditions ('poised in both conditions'); the indirect steps (2, 3, ...) and the unexplained genes follow. Bars give the number of
baited promoter fragments; the hatched part of a mark bar at a step are the genes beyond the mirrored loci with a
per-gene permutation P < 0.05. The number of genes and, where a step curve was given (--proximitySupport), the
probability that graph neighbours at that step are close in the nucleus are written next to each bar. Comparisons
whose marks allow no state calls (no repressive mark) say so. With the roles learned at every level
(<prefix>_explanation_cascade_learned.tsv) a second PDF per comparison, explanation_cascade_<name>_learned.pdf; with the
states under the maximum 1 kb window statistic (<prefix>_explanation_cascade_max1kb.tsv) explanation_cascade_<name>_max1kb.pdf
(and _learned_max1kb).

Genes explained at each level regardless of earlier levels (<prefix>_explanation_levels.tsv of egaDifferential):
explanation_levels_<name>.pdf (declared roles) and explanation_levels_<name>_learned.pdf, one bar per level in rank
order (promoter marks, own contacts, promoter state, step 1 marks, step 1 states, step 2, step 3), the state levels
and the genes explained at no level once per presence statistic (fragment mean, maximum 1 kb window). State bars are
split into a concordant transition and poised in both conditions; the hatched part of a step bar are the genes beyond
the mirrored loci with a per-gene permutation P < 0.05 at that step; the 3D support of each step is written next to it.

The title names the comparison (the name with '_' read as an arrow, hESC_hNEC as hESC -> hNEC).

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
TYPES = (('marks', '#2a78d6', 'mark change'), ('interactions', '#eb6834', 'interaction change'),
         ('both', '#1baf7a', 'mark and interaction change'), ('state', '#eda100', 'change to or from bivalent'),
         ('poised_unresolved', '#e87ba4', 'bivalent in both conditions'))
ROWS = {'gene': {'marks': 'gene_marks_only', 'interactions': 'gene_interactions_only', 'both': 'gene_both',
                 'state': 'gene_state', 'poised_unresolved': 'gene_poised_unresolved'},
        'step1': {'marks': 'step1_marks', 'state': 'step1_state', 'poised_unresolved': 'step1_poised_unresolved'}}
# panel titles of the comparisons of the manuscript; other names are shown with '_' read as an arrow
TITLES = {'hESC_hNEC': 'hESC \u2192 hNEC (neural induction)', 'Md_PA2': 'Md \u2192 PA2 (mouse pharyngeal arches)'}
UNEXPLAINED = '#7f7f7a'
INK = '#1f1f1e'
MUTED = '#6b6a66'


def label(level):
    if level == 'gene':
        return 'promoter'
    if level == 'unexplained':
        return 'no explanation'
    return 'step {}'.format(int(level[4:]))


def value(table, row, column, default=0):
    if row not in table.index or column not in table.columns:
        return default
    v = table.loc[row, column]
    return default if pd.isna(v) else v


def draw(table, path, subject=None):
    table = table.set_index('level')
    steps = [lv for lv in table.index if lv.startswith('step') and lv[4:].isdigit()]
    levels = ['gene'] + steps + ['unexplained']
    not_callable = 'note' in table.columns and (table['note'] == 'states not callable').any()
    fig, ax = plt.subplots(figsize=(2.9, 2.6))
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
            ax.barh(y, n, color=UNEXPLAINED, height=0.62, lw=0)
        else:
            dark = int(value(table, level, 'beyond_mirror_p_lt_0.05'))
            ax.barh(y, n, color=TYPES[0][1], height=0.62, lw=0)
            ax.barh(y, dark, color='none', height=0.62, lw=0, hatch='//////', edgecolor='white')
        text = '{:,}'.format(n)
        if level.startswith('step'):
            dark = value(table, level, 'beyond_mirror_p_lt_0.05', float('nan'))
            if not pd.isna(dark):
                text += ' ({:,})'.format(int(dark))
        ax.text(n + 0.02 * largest, y, text, va='center', ha='left', fontsize=8.5, color=INK)
    ax.set_yticks(range(len(levels)))
    ax.set_yticklabels([label(lv) for lv in levels], fontsize=9.5)
    ax.invert_yaxis()
    ax.set_xlim(0, largest * 1.8)
    ax.set_xlabel('promoter fragments of differentially expressed genes', fontsize=9.5)
    if subject:
        ax.set_title(subject, loc='left', fontsize=10.5, color=INK)
    ax.tick_params(axis='x', labelsize=8.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    kinds = [t for t in TYPES if t[0] in used or (t[0] in ('marks',))]
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in kinds]
    names = [name for _, _, name in kinds]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=TYPES[0][1], hatch='//////', edgecolor='white', lw=0))
    names.append('mark change beyond mirrored loci, P < 0.05 (in parentheses)')
    ax.legend(handles, names, fontsize=8, frameon=False, loc='upper center', bbox_to_anchor=(0.4, -0.24), ncol=2)
    if not_callable:
        ax.text(0.99, 0.99, 'states not callable (no repressive mark)', transform=ax.transAxes, ha='right',
                va='top', fontsize=8, color=MUTED)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


LEVEL_ROWS = (('promoter_marks', 'mean', 'promoter marks'), ('promoter_contacts', 'mean', 'own contacts'),
              ('promoter_state', 'mean', 'promoter state, mean'), ('promoter_state', 'max1kb', 'promoter state, max 1 kb'),
              ('step1_marks', 'mean', 'step 1 marks'), ('step1_state', 'mean', 'step 1 states, mean'),
              ('step1_state', 'max1kb', 'step 1 states, max 1 kb'))


def draw_levels(table, roles, path):
    """Genes explained at each level regardless of earlier levels (one role set), state levels and 'no level' under
    both presence statistics."""
    table = table[table['roles'] == roles].set_index(['presence', 'level'])
    steps = sorted({lv for _, lv in table.index if lv.startswith('step') and lv[4:].isdigit()}, key=lambda x: int(x[4:]))
    rows = list(LEVEL_ROWS) + [(lv, 'mean', '{} steps'.format(lv[4:])) for lv in steps if lv != 'step1'] + \
        [('none', 'mean', 'no level, mean'), ('none', 'max1kb', 'no level, max 1 kb')]
    rows = [r for r in rows if (r[1], r[0]) in table.index]
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    largest = max(max(int(table.loc[(p, lv), 'genes_any']) for lv, p, _ in rows), 1)
    used = set()
    for y, (level, presence, _) in enumerate(rows):
        r = table.loc[(presence, level)]
        n = int(r['genes_any'])
        if level in ('promoter_state', 'step1_state'):
            left = 0
            for kind, column in (('state', 'any_transition'), ('poised_unresolved', 'any_poised_in_both')):
                width = int(value_of(r, column))
                colour = dict((t[0], t[1]) for t in TYPES)[kind]
                if width:
                    used.add(kind)
                ax.barh(y, width, left=left, color=colour, height=0.62, lw=0)
                left += width
        elif level == 'none':
            ax.barh(y, n, color=UNEXPLAINED, height=0.62, lw=0)
        else:
            colour = TYPES[1][1] if level == 'promoter_contacts' else TYPES[0][1]
            used.add('interactions' if level == 'promoter_contacts' else 'marks')
            ax.barh(y, n, color=colour, height=0.62, lw=0)
        text = '{:,}'.format(n)
        if level.startswith('step') and level[4:].isdigit() or level == 'step1_marks':
            dark = value_of(r, 'any_beyond_mirror_p_lt_0.05', float('nan'))
            if not pd.isna(dark):
                ax.barh(y, int(dark), color='none', height=0.62, lw=0, hatch='//////', edgecolor='white')
                text += ' ({:,})'.format(int(dark))
        if level.startswith('step') and level != 'none':
            support = value_of(r, 'proximity_support', float('nan'))
            interval = value_of(r, 'proximity_ci', '')
            if not pd.isna(support):
                text += ', 3D {:.2f}'.format(support)
                if isinstance(interval, str) and interval and interval != 'NA':
                    text += ' ({})'.format(interval)
        ax.text(n + 0.02 * largest, y, text, va='center', ha='left', fontsize=6, color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for _, _, label in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, largest * 1.8)
    ax.set_xlabel('differentially expressed promoter fragments explained at the level', fontsize=7)
    ax.tick_params(axis='x', labelsize=6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    names = {'marks': 'marks', 'interactions': 'own interactions', 'state': 'bivalent or poised state resolved',
             'poised_unresolved': 'poised in both conditions'}
    kinds = [t for t in TYPES if t[0] in used or t[0] == 'marks']
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in kinds]
    labels = [names[k] for k, _, _ in kinds]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=TYPES[0][1], hatch='//////', edgecolor='white', lw=0))
    labels.append('marks beyond mirror, P < 0.05')
    ax.legend(handles, labels, fontsize=5.5, frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.16), ncol=2)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def value_of(row, column, default=0):
    if column not in row.index:
        return default
    v = row[column]
    return default if pd.isna(v) else v


def main():
    results, figures = sys.argv[1:3]
    comparisons = [tuple(x.split('=', 1)) for x in sys.argv[3:]] or COMPARISONS
    os.makedirs(figures, exist_ok=True)
    for name, prefix in comparisons:
        # first-level cascade: declared roles; roles learned separately at every level (egaDifferential --roles
        # both); states under the fragment mean and under the maximum 1 kb window
        for suffix in ('', '_learned', '_max1kb', '_learned_max1kb'):
            source = os.path.join(results, prefix + '_explanation_cascade' + suffix + '.tsv')
            if not os.path.exists(source):
                print('missing', source)
                continue
            target = os.path.join(figures, 'explanation_cascade_{}{}.pdf'.format(name, suffix))
            # the comparison in the axis label, as 'hESC / hNEC' for the name hESC_hNEC
            draw(pd.read_csv(source, sep='\t', keep_default_na=False, na_values=['']), target,
                 subject=TITLES.get(name, name.replace('_', ' \u2192 ')))
            print(target)
        # genes explained at each level regardless of earlier levels
        source = os.path.join(results, prefix + '_explanation_levels.tsv')
        if not os.path.exists(source):
            print('missing', source)
            continue
        table = pd.read_csv(source, sep='\t', keep_default_na=False, na_values=[''])
        for roles, suffix in (('declared', ''), ('learned', '_learned')):
            if (table['roles'] == roles).any():
                target = os.path.join(figures, 'explanation_levels_{}{}.pdf'.format(name, suffix))
                draw_levels(table, roles, target)
                print(target)


if __name__ == '__main__':
    main()
