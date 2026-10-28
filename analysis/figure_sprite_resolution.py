"""Supplementary figure: SPRITE support of single interactions at restriction fragment resolution.

(a) Multi-way support of direct interactions and of indirect contacts between two other ends against the fraction of
SPRITE clusters used (analysis/sprite_downsampling.py). (b) Fraction of contacts exceeding the 95th percentile of their
matched controls, overall and for direct interactions beyond 500 kb. (c) Aggregate co-occurrence around the two
fragments of direct interactions (+-10 fragments) over matched pairs, log2, by offset in fragments, normalised for
distance (top) and for distance and fragment coverage (bottom) (analysis/sprite_pileup.py).

Usage: python figure_sprite_resolution.py <downsampling.tsv> <pileup.npz> <output.pdf>
"""

import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OFFSETS = (('21-50', '21 to 50 fragments'), ('51-150', '51 to 150 fragments'),
           ('151-100000', 'more than 150 fragments'))
SETS = (('direct', 'direct interactions', '#1f5f8b'), ('indirect_oe-oe', 'indirect, other-other', '#d99a4e'))


def main():
    downsampling, pileup, output = sys.argv[1:4]
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(7.0, 5.6))
    grid = fig.add_gridspec(2, 5, height_ratios=[1, 1.05], hspace=0.55, wspace=0.6)
    t = pd.read_csv(downsampling, sep='\t')
    ax = fig.add_subplot(grid[0, 0:2])
    for name, label, color in SETS:
        sub = t[t['set'] == name]
        ax.plot(sub['groups'] / 1e6, sub['concordance'], marker='o', color=color, lw=1.5, ms=4, label=label)
    ax.axhline(0.5, color='#999999', lw=0.8, ls='--')
    ax.set_xscale('log', base=2)
    ax.set_xticks(t['groups'].unique() / 1e6)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in t['groups'].unique() / 1e6])
    ax.set_xlabel('SPRITE clusters (millions)')
    ax.set_ylabel('P(graph pair > matched pair)')
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title('a', loc='left', fontweight='bold')
    ax = fig.add_subplot(grid[0, 2:4])
    for name, label, color in SETS:
        sub = t[t['set'] == name]
        ax.plot(sub['groups'] / 1e6, 100 * sub['confirmed_p'], marker='o', color=color, lw=1.5, ms=4, label=label)
        ax.plot(sub['groups'] / 1e6, 100 * sub['control_confirmed_p'], marker='o', color=color, lw=1, ms=3, ls=':',
                label=label + ', controls')
    sub = t[t['set'] == 'direct']
    ax.plot(sub['groups'] / 1e6, 100 * sub['confirmed_p_500kb_plus'], marker='s', color='#1f5f8b', lw=1.5, ms=4,
            ls='--', label='direct interactions > 500 kb')
    ax.set_xscale('log', base=2)
    ax.set_xticks(t['groups'].unique() / 1e6)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in t['groups'].unique() / 1e6])
    ax.set_xlabel('SPRITE clusters (millions)')
    ax.set_ylabel('above 95th percentile of controls (%)')
    ax.legend(frameon=False, fontsize=6, loc='upper left')
    ax.set_title('b', loc='left', fontweight='bold')
    data = np.load(pileup)
    matrices = dict(zip(data['keys'], data['matrices']))
    window = int(data['window'])
    images = []
    for row, normalisation in enumerate(('distance', 'distance+coverage')):
        for col, (offset, title) in enumerate(OFFSETS):
            test = matrices['direct|{}|test|{}'.format(offset, normalisation)]
            control = matrices['direct|{}|control|{}'.format(offset, normalisation)]
            ratio = np.log2((test + 1e-9) / (control + 1e-9))
            images.append((row, col, ratio, title, normalisation))
    sub = grid[1, :].subgridspec(2, 3, hspace=0.35, wspace=0.15)
    limit = max(np.abs(r).max() for _, _, r, _, _ in images)
    for row, col, ratio, title, normalisation in images:
        ax = fig.add_subplot(sub[row, col])
        image = ax.imshow(ratio, cmap='RdBu_r', vmin=-limit, vmax=limit, origin='upper',
                          extent=(-window - 0.5, window + 0.5, window + 0.5, -window - 0.5))
        ax.set_xticks([-window, 0, window])
        ax.set_yticks([-window, 0, window])
        ax.tick_params(labelsize=6, length=2)
        if row == 0:
            ax.set_title(title, fontsize=7)
        else:
            ax.set_xlabel('fragments from end 2', fontsize=6.5)
        if col == 0:
            ax.set_ylabel(('distance' if normalisation == 'distance' else 'distance and\ncoverage')
                          + '\nfrom end 1', fontsize=6.5)
        if row == 0 and col == 0:
            ax.text(-0.55, 1.25, 'c', transform=ax.transAxes, fontweight='bold', fontsize=10)
    bar = fig.colorbar(image, ax=fig.axes[-6:], shrink=0.8, pad=0.02)
    bar.set_label('$\\log_2$ co-occurrence over matched pairs', fontsize=6.5)
    bar.ax.tick_params(labelsize=6)
    fig.savefig(output, bbox_inches='tight')
    print(output)


if __name__ == '__main__':
    main()
