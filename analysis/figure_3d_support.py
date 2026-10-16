"""Figure 2: evaluation of the interaction graph with ligation-free multi-way data.

(a) SPRITE: crosslinked complexes are barcoded by split-pool rounds, so all fragments of one complex share a barcode
and are read out together as one cluster, without ligation. (b) GAM: thin slices of single nuclei; loci that are close
in three dimensions are found in the same slice (tube) more often than distant loci. (c) How the multi-way groups test
the graph: an indirect contact (E1-E2, two steps in the graph) is observed directly when both loci occur in one group,
and is compared with a matched pair at the same genomic offset; a community is compared with its loci shifted along
the chromosome. (d) Probability that a pair at graph distance k co-occurs more than matched pairs without graph
connection (analysis/hop_distance.py) for three graph / assay combinations. (e) Fraction of promoter hubs with at
least one SPRITE cluster holding the promoter and two of its partners, for the hubs and for the same hubs shifted
along the chromosome (analysis/hub_multiway.py).

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

from figure_concept import (AQUA, BLUE, DASH, EDGES, GRID, INK, LAYOUT, MUTED, ORANGE, YELLOW, other_end,  # noqa: E402
                            promoter)

SERIES = (('GM12878 graph, GM12878 SPRITE', 'hops/GM_chicago_sprite_GM12878_hg19_HindIII_fragments.tsv', '#1f5f8b', 'o'),
          ('H9 hESC graph, H1-hESC GAM', 'hops/hESC_gam_H1_hg38_30kb.tsv', '#c0504d', 's'),
          ('H1-hESC graph, H1-hESC GAM', 'hops/H1_gam_H1_hg38_30kb.tsv', '#d99a4e', '^'))
BARCODE = ('#7b5ea7', '#2a9d8f', '#e9c46a')
ASSEMBLY = {'P': (0.0, 0.0), 'E1': (-0.32, 0.22), 'E2': (0.3, 0.26), 'E3': (0.05, -0.3)}


def label(ax, letter, title=None):
    ax.annotate(letter, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(-12, 4), textcoords='offset points',
                fontweight='bold', fontsize=10, color=INK, va='bottom')
    if title:
        ax.annotate(title, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(2, 4), textcoords='offset points',
                    fontsize=8, color=INK, va='bottom')


def assembly(ax, cx, cy, scale=1.0, names=True, links=True):
    points = {n: (cx + scale * x, cy + scale * y) for n, (x, y) in ASSEMBLY.items()}
    if links:
        for name in ('E1', 'E2', 'E3'):
            ax.plot(*zip(points['P'], points[name]), color=MUTED, lw=0.6, zorder=2)
    promoter(ax, *points['P'], 'P', size=45 * scale, label=False)
    for name in ('E1', 'E2', 'E3'):
        other_end(ax, *points[name], name, size=38 * scale, label=False)
    if names:
        for name, (dx, dy) in (('P', (-0.3, -0.2)), ('E1', (-0.26, 0.0)), ('E2', (0.1, 0.05)), ('E3', (0.12, -0.08))):
            ax.text(points[name][0] + dx * scale, points[name][1] + dy * scale, name, fontsize=6, color=INK, zorder=5)
    return points


def barcode(ax, x, y, width=0.12):
    for i, color in enumerate(BARCODE):
        ax.add_patch(mpatches.Rectangle((x + i * width, y - 0.05), width, 0.1, color=color, lw=0, zorder=5))


def panel_sprite(ax):
    # crosslinked complex in a nucleus
    ax.add_patch(mpatches.Ellipse((0.75, 0.95), 1.35, 1.45, facecolor='#f4f3ef', edgecolor=MUTED, lw=0.8))
    ax.add_patch(mpatches.Circle((0.75, 0.95), 0.48, facecolor=AQUA, alpha=0.15, lw=0))
    points = assembly(ax, 0.75, 0.95)
    for a, b in (('E1', 'E2'), ('E2', 'E3'), ('E3', 'E1')):
        (x1, y1), (x2, y2) = points[a], points[b]
        ax.plot([x1, x2], [y1, y2], color=MUTED, lw=0.5, ls=(0, (1, 1.5)), zorder=1)
    ax.text(0.75, 0.08, 'crosslinked complex', ha='center', fontsize=6.5, color=INK)
    ax.annotate('', xy=(1.95, 0.95), xytext=(1.5, 0.95), arrowprops=dict(arrowstyle='->', color=INK, lw=1.0))
    ax.text(1.62, 0.85, 'split-pool\nbarcoding', ha='center', va='top', fontsize=6.5, color=INK)
    # every fragment of the complex carries the same barcode
    for i, name in enumerate(('P', 'E1', 'E2', 'E3')):
        y = 1.55 - 0.28 * i
        if name == 'P':
            ax.add_patch(mpatches.Rectangle((2.05, y - 0.07), 0.32, 0.14, color=INK, lw=0))
        else:
            ax.add_patch(mpatches.Rectangle((2.05, y - 0.07), 0.32, 0.14, facecolor='white', edgecolor=INK, lw=1.0))
        barcode(ax, 2.4, y)
        ax.text(2.82, y, name, va='center', fontsize=6.5, color=INK)
    ax.text(2.05, 0.38, 'one barcode = one cluster', ha='left', fontsize=6.5, color=INK)
    ax.text(0.0, -0.25, 'E1 and E2 are read out together, without ligation', fontsize=6.5, color=ORANGE)
    ax.set_xlim(0.0, 3.4)
    ax.set_ylim(-0.35, 1.75)
    ax.set_aspect('equal')
    ax.axis('off')
    label(ax, 'a', 'SPRITE: split-pool barcoding of complexes')


def panel_gam(ax):
    # three nuclei, sliced; the slice of the first one cuts through the assembly
    for k, (cx, inside) in enumerate(((0.65, True), (1.75, False), (2.85, False))):
        nucleus = mpatches.Ellipse((cx, 0.95), 0.95, 1.3, facecolor='#f4f3ef', edgecolor=MUTED, lw=0.8)
        ax.add_patch(nucleus)
        slice_y = 0.95 if inside else (0.55 if k == 1 else 1.35)
        band = mpatches.Rectangle((cx - 0.5, slice_y - 0.13), 1.0, 0.26, facecolor=YELLOW, alpha=0.3, lw=0)
        ax.add_patch(band)
        band.set_clip_path(nucleus)
        assembly(ax, cx, 0.95, scale=0.75, names=(k == 0))
        ax.text(cx, 0.08, 'slice: P, E1, E2, E3' if inside else ('slice: E3' if k == 1 else 'slice: E1, E2'),
                ha='center', fontsize=6, color=INK if inside else MUTED)
    ax.text(0.0, -0.25, 'close loci share a slice more often than distant loci', fontsize=6.5, color=ORANGE)
    ax.set_xlim(0.0, 3.4)
    ax.set_ylim(-0.35, 1.75)
    ax.set_aspect('equal')
    ax.axis('off')
    label(ax, 'b', 'GAM: thin nuclear slices')


def panel_test(ax):
    # left: graph with grey edges, the indirect contact E1-E2 and community 1
    ox, oy = 0.0, 0.25
    pos = {n: (x + ox, y + oy) for n, (x, y) in LAYOUT.items()}
    ax.add_patch(mpatches.FancyBboxPatch((ox + 0.3, oy + 0.0), 2.2, 1.4, boxstyle='round,pad=0.02,rounding_size=0.25',
                                         facecolor=AQUA, alpha=0.13, lw=0))
    ax.text(ox + 0.35, oy - 0.2, 'community 1', fontsize=6.5, color=AQUA)
    for u, v in EDGES:
        ax.plot(*zip(pos[u], pos[v]), color=GRID, lw=1.3, zorder=2)
    t = np.linspace(0, 1, 60)[:, None]
    (x1, y1), (x2, y2) = pos['E1'], pos['E2']
    control = np.array([(x1 + x2) / 2 - 0.25 * (y2 - y1) / np.hypot(x2 - x1, y2 - y1),
                        (y1 + y2) / 2 + 0.25 * (x2 - x1) / np.hypot(x2 - x1, y2 - y1)])
    curve = (1 - t) ** 2 * np.array([x1, y1]) + 2 * (1 - t) * t * control + t ** 2 * np.array([x2, y2])
    ax.plot(curve[:, 0], curve[:, 1], color=ORANGE, lw=1.6, ls=DASH, zorder=3)
    for name, (x, y) in pos.items():
        if name.startswith('P'):
            promoter(ax, x, y, name, size=45, label=False)
        else:
            other_end(ax, x, y, name, size=38, label=False,
                      color='#fde3d6' if name in ('E1', 'E2') else 'white')
        dy = 0.11 if name in ('E2', 'P3') else -0.25
        ax.text(x, y + dy, name, fontsize=6, color=INK, zorder=5, ha='center')
    ax.text(ox + 1.1, oy + 1.62, 'graph: E1-E2 is an indirect contact', ha='center', fontsize=6.5, color=ORANGE)
    ax.annotate('', xy=(4.65, 0.85), xytext=(4.05, 0.85), arrowprops=dict(arrowstyle='->', color=INK, lw=1.0))
    ax.text(4.35, 0.97, 'test', ha='center', fontsize=6.5, color=INK)
    # right: multi-way groups (clusters or tubes) over loci in genomic order and a matched control pair
    loci = ['E1', 'P1', 'E2', 'E3', 'P2', 'E4', 'E5', 'P3', '', 'c1', 'c2']
    groups = [{'E1', 'P1', 'E2', 'P2'}, {'E1', 'E2', 'E3'}, {'P1', 'E3', 'P2', 'E4'},
              {'E1', 'P1', 'E2', 'E3', 'P2'}, {'E4', 'P3', 'E5', 'c1'}, {'E2', 'c2', 'E5'}]
    x0, y0, w, h = 5.35, 0.05, 0.3, 0.24
    for j, name in enumerate(loci):
        if name:
            ax.text(x0 + j * w + w / 2, y0 + len(groups) * h + 0.05, name, ha='center', va='bottom', fontsize=6,
                    color=ORANGE if name in ('E1', 'E2') else (MUTED if name.startswith('c') else INK))
    for i, members in enumerate(groups):
        y = y0 + (len(groups) - 1 - i) * h
        ax.text(x0 - 0.08, y + h / 2, 'group {}'.format(i + 1), ha='right', va='center', fontsize=6, color=MUTED)
        both = {'E1', 'E2'} <= members
        community = {'E1', 'P1', 'E2', 'E3', 'P2'} <= members
        if community:
            ax.add_patch(mpatches.Rectangle((x0 - 0.02, y - 0.01), 5 * w + 0.04, h + 0.02, facecolor='none',
                                            edgecolor=AQUA, lw=1.4, zorder=4))
        for j, name in enumerate(loci):
            if name in members:
                color = ORANGE if (both and name in ('E1', 'E2')) else (MUTED if name.startswith('c') else '#9fa3a8')
                ax.add_patch(mpatches.FancyBboxPatch((x0 + j * w + 0.04, y + 0.03), w - 0.08, h - 0.06,
                                                     boxstyle='round,pad=0.0,rounding_size=0.04', color=color, lw=0,
                                                     zorder=3))
    right = x0 + len(loci) * w + 0.12
    ax.text(right, 1.5, 'E1 with E2: 3 groups', fontsize=6.5, color=ORANGE, va='center')
    ax.text(right, 1.16, 'matched pair c1 with c2\n(same offset, no edge): 0 groups', fontsize=6.5, color=MUTED,
            va='center')
    ax.text(right, 0.62, 'community 1 complete\nin one group', fontsize=6.5, color=AQUA, va='center')
    ax.text(right, 0.18, 'compared with shifted\ncommunities', fontsize=6.5, color=INK, va='center')
    ax.set_xlim(-0.05, 11.6)
    ax.set_ylim(-0.1, 2.0)
    ax.set_aspect('equal')
    ax.axis('off')
    label(ax, 'c', 'Multi-way groups observe indirect contacts and communities directly')


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
    label(ax, 'd')


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
    label(ax, 'e')


def main():
    results, output = sys.argv[1:3]
    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(7.2, 7.0))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.05, 1.35], width_ratios=[1, 1.1], hspace=0.32, wspace=0.3)
    panel_sprite(fig.add_subplot(grid[0, 0]))
    panel_gam(fig.add_subplot(grid[0, 1]))
    panel_test(fig.add_subplot(grid[1, :]))
    for column, panel in ((0, panel_steps), (1, panel_hubs)):
        ax = fig.add_subplot(grid[2, column])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        panel(ax, results)
    fig.savefig(output, bbox_inches='tight', pad_inches=0.03)
    print(output)


if __name__ == '__main__':
    main()
