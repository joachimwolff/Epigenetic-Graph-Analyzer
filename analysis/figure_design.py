"""Figures of the experiment design (egaDesign, analysis/25_design.sh), one PDF per panel.

  candidates_<comparison>_<gene>.pdf  the eight best ranked candidates of one gene: statistic of the explaining
                                      level before (open circle) and after flattening the node (filled circle),
                                      coloured by level, the explanation threshold dashed; state candidates as a
                                      label (level lost or kept).
  single_node_dependence.pdf          per comparison and level the fraction of the explained genes whose level is
                                      removed by flattening a single node.
  pleiotropy.pdf                      per comparison the number of strong candidate nodes (removing >= half of a
                                      level statistic) by the number of genes they explain.
  validation_rank_bins.pdf            K562 CRISPRi: fraction of element-gene pairs with a significant effect by the
                                      design rank of the element (Wilson 95% CI), all scorable pairs.
  validation_rank_bins_distance.pdf   the same within three distance strata.
  validation_within_gene.pdf          per gene with positive and negative reached elements: within-gene AUROC of
                                      the design score against that of distance.

Usage: python figure_design.py <design dir> <figure dir> [comparison:gene ...]
Default genes: Hoxa3 and Hoxa2 (Md_vs_PA2), STAMBPL1 and NR2F1 (hESC_vs_hNEC).
"""

import json
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COMPARISONS = ['hESC_vs_hNEC', 'Md_vs_PA2', 'siMM_vs_siRE', 't821_vs_FLT3ITD']
LABELS = {'hESC_vs_hNEC': 'hESC vs hNEC', 'Md_vs_PA2': 'Md vs PA2', 'siMM_vs_siRE': 'Kasumi-1 siMM vs siRE',
          't821_vs_FLT3ITD': 't(8;21) vs FLT3-ITD'}
LEVELS = ['promoter_contacts', 'step1_marks', 'step1_state', 'step2', 'step3']
LEVEL_LABELS = {'promoter_contacts': 'own contacts', 'step1_marks': 'step 1 marks', 'step1_state': 'step 1 state',
                'step2': 'step 2', 'step3': 'step 3', 'promoter_marks+promoter_state': 'own promoter',
                'promoter_marks': 'own promoter', 'promoter_state': 'own promoter'}
# categorical slots in fixed order (one per level / comparison)
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']
LEVEL_COLOURS = dict(zip(LEVELS, SERIES))
COMPARISON_COLOURS = dict(zip(COMPARISONS, SERIES))
INK, MUTED = '#2b2b2a', '#8a8a85'
DEFAULT_GENES = [('Md_vs_PA2', 'Hoxa3'), ('Md_vs_PA2', 'Hoxa2'), ('hESC_vs_hNEC', 'STAMBPL1'),
                 ('hESC_vs_hNEC', 'NR2F1')]
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'axes.edgecolor': MUTED, 'xtick.color': INK, 'ytick.color': INK, 'axes.labelcolor': INK})


def gene_rows(design, gene):
    hit = design['gene'].astype(str).str.split(',').apply(lambda names: gene in names)
    rows = design[hit]
    if rows.empty:
        return rows
    first = rows.sort_values('rank').iloc[0]['gene_node']
    return rows[rows['gene_node'] == first].sort_values('rank')


def short_label(r):
    import re
    cls = str(r['element_class']).replace(' in both', '')
    both = re.match(r'(.+) in \S+, (.+) in \S+$', cls)
    if both:
        cls = '{} -> {}'.format(both.group(1), both.group(2))
    position = str(r['position']).rsplit(' ', 1)[0]
    position = position[:-3] if position.endswith(' of') else position
    return '{}. {}, {}'.format(int(r['rank']), cls, position)


def candidate_panel(design, comparison, gene, path, top=8):
    rows = gene_rows(design, gene).head(top)
    if rows.empty:
        return False
    fig, ax = plt.subplots(figsize=(4.6, 0.32 * len(rows) + 0.9))
    y = np.arange(len(rows))[::-1]
    xmax = 1.0
    for yi, (_, r) in zip(y, rows.iterrows()):
        level = r['level'].split(',')[0]
        colour = LEVEL_COLOURS.get(level, MUTED)
        if level.startswith('promoter_marks') or level.startswith('promoter_state'):
            colour = MUTED
        before, after = r['statistic_before'], r['statistic_flatten']
        if pd.isna(before):
            lost = str(r['level_lost_flatten']).lower() == 'true'
            ax.text(0.02, yi, 'state level {}'.format('lost' if lost else 'kept'), va='center', fontsize=7,
                    color=INK, transform=ax.get_yaxis_transform())
            ax.plot([0], [yi], marker='s', ms=6, color=colour, transform=ax.get_yaxis_transform(), clip_on=False)
            continue
        xmax = max(xmax, before, after)
        ax.plot([after, before], [yi, yi], color=colour, lw=2, solid_capstyle='round')
        ax.plot(before, yi, 'o', ms=8, mfc='white', mec=colour, mew=2)
        ax.plot(after, yi, 'o', ms=8, color=colour, mec='white', mew=1)
    ax.axvline(1.0, color=MUTED, lw=1, ls='--')
    ax.set_yticks(y)
    ax.set_yticklabels([short_label(r) for _, r in rows.iterrows()], fontsize=7)
    ax.set_xlim(min(0, rows['statistic_flatten'].min(skipna=True) if rows['statistic_flatten'].notna().any()
                    else 0) - 0.1, xmax * 1.08)
    ax.set_xlabel('level statistic (open: observed, filled: node flattened)')
    first = rows.iloc[0]
    ax.set_title('{} ({}), log2FC {:+.1f}'.format(gene, LABELS[comparison], first['rna_log2FC']), fontsize=8,
                 loc='left', color=INK)
    present = [lv for lv in LEVELS if any(rows['level'].str.startswith(lv))]
    handles = [plt.Line2D([], [], marker='o', ls='', color=LEVEL_COLOURS[lv], label=LEVEL_LABELS[lv])
               for lv in present]
    if any(rows['level'].str.startswith('promoter_marks') | rows['level'].str.startswith('promoter_state')):
        handles.append(plt.Line2D([], [], marker='o', ls='', color=MUTED, label='own promoter'))
    height = 0.32 * len(rows) + 0.9
    ax.legend(handles=handles, frameon=False, fontsize=7, loc='upper center',
              bbox_to_anchor=(0.5, -0.55 / height), ncol=min(5, len(handles)))
    ax.grid(axis='x', color='#e6e5e0', lw=0.6)
    ax.set_axisbelow(True)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return True


def dependence_panel(directory, path):
    fig, ax = plt.subplots(figsize=(4.6, 2.4))
    width = 0.8 / len(COMPARISONS)
    for i, comparison in enumerate(COMPARISONS):
        summary_path = os.path.join(directory, comparison + '_design_summary.json')
        if not os.path.exists(summary_path):
            continue
        by_level = json.load(open(summary_path))['by_level']
        x, values, counts = [], [], []
        for j, level in enumerate(LEVELS):
            if level not in by_level or by_level[level]['genes'] == 0:
                continue
            x.append(j + (i - (len(COMPARISONS) - 1) / 2) * width)
            values.append(by_level[level]['genes_single_node_removes_level'] / by_level[level]['genes'])
            counts.append(by_level[level]['genes'])
        ax.bar(x, values, width=width - 0.02, color=COMPARISON_COLOURS[comparison], label=LABELS[comparison],
               edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(len(LEVELS)))
    ax.set_xticklabels([LEVEL_LABELS[lv] for lv in LEVELS])
    ax.set_ylabel('explained genes whose level is\nremoved by one node (fraction)')
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc='upper right')
    ax.grid(axis='y', color='#e6e5e0', lw=0.6)
    ax.set_axisbelow(True)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def pleiotropy_panel(directory, path, strong=0.5):
    fig, ax = plt.subplots(figsize=(4.6, 2.4))
    bins = ['1', '2', '3', '4', '5+']
    width = 0.8 / len(COMPARISONS)
    for i, comparison in enumerate(COMPARISONS):
        candidates_path = os.path.join(directory, comparison + '_design_candidates.tsv.gz')
        if not os.path.exists(candidates_path):
            continue
        c = pd.read_csv(candidates_path, sep='\t', usecols=['node', 'genes', 'n_genes_strong', 'own_promoter'])
        nodes = c[(c['n_genes_strong'] > 0) & ~c['own_promoter'].astype(bool)].drop_duplicates('node')
        counts = nodes['n_genes_strong'].clip(upper=5).value_counts()
        values = [counts.get(k, 0) for k in range(1, 6)]
        total = max(1, sum(values))
        x = np.arange(5) + (i - (len(COMPARISONS) - 1) / 2) * width
        ax.bar(x, np.array(values) / total, width=width - 0.02, color=COMPARISON_COLOURS[comparison],
               label='{} ({:,} nodes)'.format(LABELS[comparison], total), edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(5))
    ax.set_xticklabels(bins)
    ax.set_xlabel('genes for which the node removes >= half of a level')
    ax.set_ylabel('strong candidate nodes (fraction)')
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis='y', color='#e6e5e0', lw=0.6)
    ax.set_axisbelow(True)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def rank_bin_panel(bins, sets, path, labels=None):
    fig, ax = plt.subplots(figsize=(4.6, 2.4))
    order = ['1', '2-3', '4-10', '11-50', '>50', 'not reached']
    width = 0.8 / len(sets)
    for i, name in enumerate(sets):
        part = bins[bins['set'] == name].set_index('rank_bin').reindex(order)
        x = np.arange(len(order)) + (i - (len(sets) - 1) / 2) * width
        frac = part['fraction'].to_numpy(dtype=float)
        low = np.clip(frac - part['ci_low'].to_numpy(dtype=float), 0, None)
        high = np.clip(part['ci_high'].to_numpy(dtype=float) - frac, 0, None)
        colour = SERIES[i]
        ax.errorbar(x, frac, yerr=[low, high], fmt='o', ms=6, color=colour, ecolor=colour, elinewidth=1.2,
                    capsize=0, label=(labels or {}).get(name, name))
        for xi, n in zip(x, part['pairs']):
            if not pd.isna(n) and len(sets) == 1:
                ax.text(xi, -0.02, 'n={:,}'.format(int(n)), ha='center', va='top', fontsize=6, color=MUTED,
                        transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.tick_params(axis='x', pad=10 if len(sets) == 1 else 3)
    ax.set_xlabel('design rank of the element among the gene\'s neighbourhood nodes')
    ax.set_ylabel('pairs with a significant\nCRISPRi effect (fraction)')
    ax.set_ylim(bottom=0)
    if len(sets) > 1:
        ax.legend(frameon=False, fontsize=7)
    ax.grid(axis='y', color='#e6e5e0', lw=0.6)
    ax.set_axisbelow(True)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def within_gene_panel(within, path):
    fig, ax = plt.subplots(figsize=(2.8, 2.8))
    ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls='--')
    ax.scatter(within['auroc_distance'], within['auroc_design'], s=18 + 4 * within['pairs'].clip(upper=20),
               color=SERIES[0], edgecolor='white', linewidth=1, zorder=3)
    ax.set_xlabel('within-gene AUROC, distance')
    ax.set_ylabel('within-gene AUROC, design score')
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_aspect('equal')
    ax.text(0.03, 0.97, '{} genes; mean {:.2f} vs {:.2f}'.format(len(within), within['auroc_design'].mean(),
                                                                within['auroc_distance'].mean()),
            transform=ax.transAxes, va='top', fontsize=7, color=INK)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def main(argv):
    directory, output = argv[1], argv[2]
    os.makedirs(output, exist_ok=True)
    genes = [tuple(a.split(':', 1)) for a in argv[3:]] or DEFAULT_GENES
    for comparison, gene in genes:
        path = os.path.join(directory, comparison + '_design.tsv')
        if not os.path.exists(path):
            continue
        design = pd.read_csv(path, sep='\t')
        candidate_panel(design, comparison, gene, os.path.join(output, 'candidates_{}_{}.pdf'.format(comparison,
                                                                                                    gene)))
    dependence_panel(directory, os.path.join(output, 'single_node_dependence.pdf'))
    pleiotropy_panel(directory, os.path.join(output, 'pleiotropy.pdf'))
    validation = os.path.join(directory, 'validation', 'K562')
    if os.path.exists(validation + '_bins.tsv'):
        bins = pd.read_csv(validation + '_bins.tsv', sep='\t')
        rank_bin_panel(bins, ['scorable'], os.path.join(output, 'validation_rank_bins.pdf'))
        sets = [s for s in bins['set'].unique() if s.startswith('scorable, distance')]
        rank_bin_panel(bins, sets, os.path.join(output, 'validation_rank_bins_distance.pdf'),
                       {s: s.replace('scorable, distance ', '') for s in sets})
    if os.path.exists(validation + '_within_gene.tsv'):
        within = pd.read_csv(validation + '_within_gene.tsv', sep='\t')
        if len(within):
            within_gene_panel(within, os.path.join(output, 'validation_within_gene.pdf'))


if __name__ == '__main__':
    main(sys.argv)
