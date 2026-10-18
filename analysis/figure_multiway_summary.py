"""Multi-way and functional support of the interaction graph in one figure (one PDF per panel, no panel letters).

  <prefix>_test.pdf      schematic of the pair test (figure_3d_support.panel_test without communities)
  <prefix>_a.pdf         probability that pairs at graph distance k co-occur more than matched pairs in the ligation-free
                         assays (GM12878 SPRITE, H1-hESC GAM); solid, controls matched for offset and coverage, with 95%
                         chromosome bootstrap band; dotted, controls also matched for compartment and domain
  <prefix>_b.pdf         the same for the concatemer assays (HiPore-C of GM12878 and K562, Pore-C of GM12878)
  <prefix>_c.pdf         hubs of three fragments with the promoter and both partners in one multi-way group, against
                         the same hubs shifted along the chromosome, in SPRITE, HiPore-C and Pore-C
  <prefix>_d.pdf         eQTL odds ratios (Mantel-Haenszel, 95% block-bootstrap interval) for a direct edge, a two-step
                         path without a direct edge and a three-step path without a shorter one, GM12878 and K562
                         graphs, observed and reflected edges

Usage: python figure_multiway_summary.py <results dir> <output prefix>
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from figure_3d_support import panel_test  # noqa: E402

plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

LIGATION_FREE = (('GM12878 graph, SPRITE', 'reviewer_controls/GM_chicago_sprite_hops', '#1f5f8b', 'o'),
                 ('H9 hESC graph, H1-hESC GAM', 'reviewer_controls/hESC_gam_hops', '#c0504d', 's'),
                 ('H1-hESC graph, H1-hESC GAM', 'reviewer_controls/H1_gam_hops', '#d99a4e', '^'))
CONCATEMER = (('GM12878 graph, HiPore-C', 'porec/GM_chicago_hiporec_hops.tsv', '#1f5f8b', 'o'),
              ('K562 graph, HiPore-C', 'porec/K562_hiporec_hops.tsv', '#6a4c93', 'D'),
              ('GM12878 graph, Pore-C', 'porec/GM_chicago_porec_hops.tsv', '#2a9d8f', 'v'))
HUBS = (('SPRITE\nGM12878', 'hubs/GM_chicago_sprite_hub_summary.tsv'),
        ('HiPore-C\nGM12878', 'porec/GM_chicago_hiporec_hubs_hub_summary.tsv'),
        ('HiPore-C\nK562', 'porec/K562_hiporec_hubs_hub_summary.tsv'),
        ('Pore-C\nGM12878', 'porec/GM_chicago_porec_hubs_hub_summary.tsv'))
EQTL_GRAPHS = (('GM_chicago', 'GM12878 graph'), ('K562', 'K562 graph'))
EQTL_SETS = (('geuvadis_susie', 'GEUVADIS LCL', 'o'), ('gtex_lcl', 'GTEx LCL', 's'), ('gtex_wb', 'GTEx blood', '^'),
             ('eqtlgen', 'eQTLGen', 'D'))
EQTL_STEPS = (('or_direct', 'or_direct_low', 'or_direct_high', 'direct'),
              ('or_two_without_direct', 'or_two_low', 'or_two_high', '2 steps'),
              ('or_three_without_1_2', 'or_three_low', 'or_three_high', '3 steps'))
OBSERVED, REFLECTED = '#1f5a96', '#a6a6a6'


def steps_axis(ax, ylim):
    ax.axhline(0.5, color='#999999', lw=0.7, ls='--')
    ax.set_xticks(range(1, 6))
    ax.set_xlabel('graph distance (steps)')
    ax.set_ylabel('P(graph pair > matched pair)')
    ax.set_ylim(*ylim)
    ax.legend(frameon=False, fontsize=6, loc='upper right')


def panel_ligation_free(ax, results):
    for name, stem, color, marker in LIGATION_FREE:
        table = pd.read_csv(os.path.join(results, stem + '_coverage.tsv'), sep='\t')
        strict = pd.read_csv(os.path.join(results, stem + '_domains.tsv'), sep='\t')
        ax.fill_between(table['steps'], table['ci_low'], table['ci_high'], color=color, alpha=0.15, lw=0)
        ax.plot(table['steps'], table['concordance'], marker=marker, color=color, lw=1.3, ms=4, label=name)
        ax.plot(strict['steps'], strict['concordance'], color=color, lw=1.0, ls=(0, (1.5, 1.5)))
    ax.plot([], [], color='#555555', lw=1.0, ls=(0, (1.5, 1.5)), label='+ compartment, domain')
    steps_axis(ax, (0.44, 0.70))


def panel_concatemer(ax, results):
    for name, path, color, marker in CONCATEMER:
        table = pd.read_csv(os.path.join(results, path), sep='\t')
        ax.fill_between(table['steps'], table['ci_low'], table['ci_high'], color=color, alpha=0.15, lw=0)
        ax.plot(table['steps'], table['concordance'], marker=marker, color=color, lw=1.3, ms=4, label=name)
    steps_axis(ax, (0.44, 0.86))


def panel_hubs(ax, results):
    observed, shifted, counts = [], [], []
    for _, path in HUBS:
        table = pd.read_csv(os.path.join(results, path), sep='\t')
        row = table[table['hub_size'] == '3-3'].iloc[0]
        observed.append(row['hubs_with_promoter_plus2_in_one_group'] * 100)
        shifted.append(row['shifted'] * 100)
        counts.append(int(row['hubs']))
    x = np.arange(len(HUBS))
    width = 0.38
    ax.bar(x - width / 2, observed, width, color='#1f5f8b', label='promoter hubs')
    ax.bar(x + width / 2, shifted, width, color='#b8c7d6', label='shifted hubs')
    for i, n in enumerate(counts):
        ax.text(i, max(observed[i], shifted[i]) + 2, 'n={:,}'.format(n), ha='center', fontsize=5.5, color='#555555')
    ax.set_xticks(x)
    ax.set_xticklabels([name for name, _ in HUBS])
    ax.set_ylabel('hubs of three fragments with\npromoter + 2 partners in one group (%)')
    ax.set_ylim(0, 60)
    ax.legend(frameon=False, fontsize=6, loc='upper right')


def panel_eqtl(fig, results):
    table = pd.read_csv(os.path.join(results, 'eqtl_validation/eqtl_summary.tsv'), sep='\t')
    axes = fig.subplots(1, 2, sharey=True)
    for ax, (graph, title) in zip(axes, EQTL_GRAPHS):
        ax.axhline(1.0, color='#999999', lw=0.7, ls='--')
        for k, (data, _, marker) in enumerate(EQTL_SETS):
            for version, color, shift in (('observed', OBSERVED, -0.09), ('mirrored', REFLECTED, 0.09)):
                row = table[(table['graph'] == graph) & (table['data'] == data) & (table['version'] == version)]
                row = row.iloc[0]
                for s, (value, low, high, _) in enumerate(EQTL_STEPS):
                    xpos = s + shift + (k - 1.5) * 0.045
                    ax.errorbar(xpos, row[value], yerr=[[row[value] - row[low]], [row[high] - row[value]]],
                                fmt=marker, color=color, ms=3, lw=0.7, capsize=0)
        ax.set_yscale('log')
        ax.set_yticks([0.5, 1, 2])
        ax.set_yticklabels(['0.5', '1', '2'])
        ax.minorticks_off()
        ax.set_xticks(range(3))
        ax.set_xticklabels([label for *_, label in EQTL_STEPS])
        ax.set_xlim(-0.45, 2.45)
        ax.set_title(title, fontsize=7)
    axes[0].set_ylabel('eQTL odds ratio')
    handles = [axes[1].plot([], [], marker, color=OBSERVED, ms=3, label=name)[0] for data, name, marker in EQTL_SETS]
    handles.append(axes[1].plot([], [], 's', color=REFLECTED, ms=3, label='reflected edges')[0])
    fig.legend(handles=handles, frameon=False, fontsize=5.5, loc='lower center', ncol=3, handletextpad=0.2,
               columnspacing=0.8, bbox_to_anchor=(0.5, -0.12))


def save(fig, path):
    fig.savefig(path, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    print(path)


def main():
    results, prefix = sys.argv[1:3]
    fig, ax = plt.subplots(figsize=(7.2, 1.6))
    panel_test(ax, communities=False)
    save(fig, prefix + '_test.pdf')
    for letter, panel in (('a', panel_ligation_free), ('b', panel_concatemer), ('c', panel_hubs)):
        fig, ax = plt.subplots(figsize=(3.3, 2.3))
        panel(ax, results)
        save(fig, '{}_{}.pdf'.format(prefix, letter))
    fig = plt.figure(figsize=(3.3, 2.3))
    panel_eqtl(fig, results)
    save(fig, prefix + '_d.pdf')


if __name__ == '__main__':
    main()
