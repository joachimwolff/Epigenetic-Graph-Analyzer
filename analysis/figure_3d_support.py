"""Figure: multi-way support by graph distance and promoter hubs in single SPRITE clusters.

(a) Probability that a pair at graph distance k co-occurs more than matched pairs without graph connection
(analysis/hop_distance.py) for three graph / assay combinations. (b) Fraction of promoter hubs with at least one
SPRITE cluster holding the promoter and two of its partners, for the hubs and for the same hubs shifted along the
chromosome (analysis/hub_multiway.py).

Usage: python figure_3d_support.py <results dir> <output.pdf>
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

SERIES = (('GM12878 graph, GM12878 SPRITE', 'hops/GM_chicago_sprite_GM12878_hg19_HindIII_fragments.tsv', '#1f5f8b', 'o'),
          ('H9 hESC graph, H1-hESC GAM', 'hops/hESC_gam_H1_hg38_30kb.tsv', '#c0504d', 's'),
          ('H1-hESC graph, H1-hESC GAM', 'hops/H1_gam_H1_hg38_30kb.tsv', '#d99a4e', '^'))


def main():
    results, output = sys.argv[1:3]
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={'width_ratios': [1, 1.1]})
    ax = axes[0]
    ax.axhline(0.5, color='#999999', lw=0.8, ls='--')
    for label, path, color, marker in SERIES:
        table = pd.read_csv(os.path.join(results, path), sep='\t')
        ax.plot(table['steps'], table['concordance'], marker=marker, color=color, lw=1.5, ms=5, label=label)
    ax.set_xticks(range(1, 6))
    ax.set_xlabel('graph distance (steps)')
    ax.set_ylabel('P(graph pair > matched pair)')
    ax.set_ylim(0.44, 0.69)
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    ax.set_title('a', loc='left', fontweight='bold')

    hubs = pd.read_csv(os.path.join(results, 'hubs/GM_chicago_sprite_hub_summary.tsv'), sep='\t')
    labels = [s.replace('-', '-').replace('>=', '>') for s in hubs['hub_size']]
    labels = ['3' if s == '3-3' else s.replace('>=', '≥') for s in hubs['hub_size']]
    ax = axes[1]
    x = range(len(hubs))
    width = 0.38
    ax.bar([i - width / 2 for i in x], hubs['hubs_with_promoter_plus2_in_one_group'] * 100, width, color='#1f5f8b',
           label='promoter hubs')
    ax.bar([i + width / 2 for i in x], hubs['shifted'] * 100, width, color='#b8c7d6', label='shifted hubs')
    for i, n in enumerate(hubs['hubs']):
        ax.text(i, 102, 'n={:,}'.format(n), ha='center', fontsize=6.5, color='#555555')
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_xlabel('hub size (promoter and partners)')
    ax.set_ylabel('hubs with promoter + 2 partners\nin one SPRITE cluster (%)')
    ax.set_ylim(0, 110)
    ax.legend(frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(0.0, 0.93))
    ax.set_title('b', loc='left', fontweight='bold')
    fig.tight_layout()
    fig.savefig(output)
    print(output)


if __name__ == '__main__':
    main()
