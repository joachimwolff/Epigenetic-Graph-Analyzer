"""Figure 2: evaluation of the interaction graph with ligation-free multi-way data.

(a) SPRITE clusters and GAM tubes as groups of loci; the loci of a hub or community occur together in several groups.
(b) Controls: matched pairs at the same genomic offset, and hubs or communities shifted along the chromosome, which
keeps the spacing of their loci. (c) Probability that a pair at graph distance k co-occurs more than matched pairs
without graph connection (analysis/hop_distance.py) for three graph / assay combinations. (d) Fraction of promoter
hubs with at least one SPRITE cluster holding the promoter and two of its partners, for the hubs and for the same
hubs shifted along the chromosome (analysis/hub_multiway.py).

Usage: python figure_3d_support.py <results dir> <output.pdf>
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

AQUA, INK, MUTED, GRID = '#1baf7a', '#1f1f1e', '#6b6a64', '#d9d8d2'
SERIES = (('GM12878 graph, GM12878 SPRITE', 'hops/GM_chicago_sprite_GM12878_hg19_HindIII_fragments.tsv', '#1f5f8b', 'o'),
          ('H9 hESC graph, H1-hESC GAM', 'hops/hESC_gam_H1_hg38_30kb.tsv', '#c0504d', 's'),
          ('H1-hESC graph, H1-hESC GAM', 'hops/H1_gam_H1_hg38_30kb.tsv', '#d99a4e', '^'))


def label(ax, letter, title=None):
    ax.annotate(letter, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(-14, 6), textcoords='offset points',
                fontweight='bold', fontsize=10, color=INK, va='bottom')
    if title:
        ax.annotate(title, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, 6), textcoords='offset points',
                    fontsize=8, color=INK, va='bottom')


def panel_groups(ax):
    rng = np.random.default_rng(3)
    loci = (2, 5, 8)
    for k in range(5):
        y = 4 - k
        members = set(rng.choice(np.arange(12), size=rng.integers(2, 5), replace=False))
        if k in (0, 2):
            members |= set(loci)
        if k == 3:
            members |= {2, 5}
        for m in sorted(members):
            ax.add_patch(mpatches.FancyBboxPatch((m - 0.35, y - 0.3), 0.7, 0.6,
                                                 boxstyle='round,pad=0.02,rounding_size=0.1',
                                                 color=AQUA if m in loci else GRID, lw=0))
    for m in loci:
        ax.plot([m, m], [4.55, 4.85], color=INK, lw=0.8)
    ax.text(5, 4.95, 'loci of a hub or community', ha='center', va='bottom', fontsize=6.5, color=INK)
    ax.set_xlim(-0.6, 11.6)
    ax.set_ylim(-0.6, 5.5)
    ax.set_yticks(range(5))
    ax.set_yticklabels(['group {}'.format(5 - i) for i in range(5)], fontsize=6.5)
    ax.set_xticks([])
    ax.set_xlabel('genomic loci', fontsize=7)
    for side in ('left', 'top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(left=False)
    label(ax, 'a', 'SPRITE clusters / GAM tubes')


def panel_controls(ax):
    loci = np.array([2, 5, 8])
    ax.scatter(loci, [1.0] * 3, s=40, color=AQUA, zorder=3)
    for shift, y in ((3, 0.6), (6, 0.2)):
        shifted = (loci + shift) % 12
        ax.scatter(shifted, [y] * 3, s=40, facecolor='white', edgecolor=MUTED, zorder=3)
    ax.hlines([1.0, 0.6, 0.2], -0.5, 11.5, color=GRID, lw=0.8, zorder=0)
    ax.text(12.0, 1.0, 'graph loci', va='center', fontsize=6.5, color=INK)
    ax.text(12.0, 0.6, 'shifted', va='center', fontsize=6.5, color=MUTED)
    ax.text(12.0, 0.2, 'shifted', va='center', fontsize=6.5, color=MUTED)
    # matched pair at the same offset
    ax.plot([2, 5], [1.45, 1.45], color=AQUA, lw=1.5)
    ax.plot([7, 10], [1.45, 1.45], color=MUTED, lw=1.5, ls=(0, (3, 2)))
    ax.scatter([2, 5], [1.45, 1.45], s=18, color=AQUA, zorder=3)
    ax.scatter([7, 10], [1.45, 1.45], s=18, facecolor='white', edgecolor=MUTED, zorder=3)
    ax.text(12.0, 1.45, 'pair and matched pair\n(same offset, no edge)', va='center', fontsize=6.5, color=INK)
    ax.set_xlim(-0.6, 18.5)
    ax.set_ylim(0.0, 1.75)
    ax.axis('off')
    label(ax, 'b', 'Controls')


def panel_steps(ax, results):
    ax.axhline(0.5, color='#999999', lw=0.8, ls='--')
    for name, path, color, marker in SERIES:
        table = pd.read_csv(os.path.join(results, path), sep='\t')
        ax.plot(table['steps'], table['concordance'], marker=marker, color=color, lw=1.5, ms=5, label=name)
    ax.set_xticks(range(1, 6))
    ax.set_xlabel('graph distance (steps)')
    ax.set_ylabel('P(graph pair > matched pair)')
    ax.set_ylim(0.44, 0.69)
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    label(ax, 'c')


def panel_hubs(ax, results):
    hubs = pd.read_csv(os.path.join(results, 'hubs/GM_chicago_sprite_hub_summary.tsv'), sep='\t')
    names = ['3' if s == '3-3' else s.replace('>=', '≥') for s in hubs['hub_size']]
    x = np.arange(len(hubs))
    width = 0.38
    ax.bar(x - width / 2, hubs['hubs_with_promoter_plus2_in_one_group'] * 100, width, color='#1f5f8b',
           label='promoter hubs')
    ax.bar(x + width / 2, hubs['shifted'] * 100, width, color='#b8c7d6', label='shifted hubs')
    for i, n in enumerate(hubs['hubs']):
        ax.text(i, 102, 'n={:,}'.format(n), ha='center', fontsize=6.5, color='#555555')
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_xlabel('hub size (promoter and partners)')
    ax.set_ylabel('hubs with promoter + 2 partners\nin one SPRITE cluster (%)')
    ax.set_ylim(0, 110)
    ax.legend(frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(0.0, 0.93))
    label(ax, 'd')


def main():
    results, output = sys.argv[1:3]
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(7.0, 5.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.75, 1], width_ratios=[1, 1.1], hspace=0.45, wspace=0.35)
    panel_groups(fig.add_subplot(grid[0, 0]))
    panel_controls(fig.add_subplot(grid[0, 1]))
    panel_steps(fig.add_subplot(grid[1, 0]), results)
    panel_hubs(fig.add_subplot(grid[1, 1]), results)
    fig.savefig(output, bbox_inches='tight', pad_inches=0.03)
    print(output)


if __name__ == '__main__':
    main()
