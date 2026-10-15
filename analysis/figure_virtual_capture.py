"""Figure: virtual promoter capture from Hi-C and Micro-C.

(a) Precision and recall of virtual calls against promoter capture calls of the same cell type (exact bin and one
bin tolerance; analysis/virtual_capture_benchmark.py), with two capture experiments of GM12878 as reference.
(b) Multi-way support by graph distance for graphs from virtual calls and from capture calls of the same cell type
(analysis/hop_distance.py).

Usage: python figure_virtual_capture.py <results dir> <output.pdf>
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

BENCHMARKS = (('GM12878, two\ncapture sets', 'GM12878_cairns_vs_mifsud20'),
              ('GM12878\nHi-C', 'GM12878'), ('K562\nHi-C', 'K562'), ('H1-hESC\nMicro-C', 'H1_microc_matched'))
CURVES = (('GM12878, SPRITE', 'virtual_capture/graphs/VC_GM12878_sprite_hops.tsv',
           'hops/GM_chicago_sprite_GM12878_hg19_HindIII_fragments.tsv', '#1f5f8b'),
          ('H1-hESC, GAM', 'virtual_capture/graphs/VC_H1_microc_gam_hops.tsv', 'hops/H1_gam_H1_hg38_30kb.tsv',
           '#d99a4e'))


def main():
    results, output = sys.argv[1:3]
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={'width_ratios': [1.3, 1]})
    ax = axes[0]
    width = 0.2
    for i, (label, name) in enumerate(BENCHMARKS):
        table = pd.read_csv(os.path.join(results, 'virtual_capture', name + '_benchmark.tsv'), sep='\t')
        table = table[table['stratum'] == 'all'].set_index('tolerance_bins')
        values = [table.loc[0, 'precision'], table.loc[1, 'precision'], table.loc[0, 'recall'], table.loc[1, 'recall']]
        colors = ['#1f5f8b', '#8fb3cf', '#c0504d', '#e6a8a6']
        for k, (value, color) in enumerate(zip(values, colors)):
            ax.bar(i + (k - 1.5) * width, value * 100, width, color=color,
                   label=['precision', 'precision, 1 bin tolerance', 'recall', 'recall, 1 bin tolerance'][k]
                   if i == 0 else None)
    ax.set_xticks(range(len(BENCHMARKS)))
    ax.set_xticklabels([b[0] for b in BENCHMARKS], fontsize=7)
    ax.set_ylabel('agreement with capture calls (%)')
    ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    ax.set_title('a', loc='left', fontweight='bold')
    ax = axes[1]
    ax.axhline(0.5, color='#999999', lw=0.8, ls='--')
    for label, virtual, capture, color in CURVES:
        v = pd.read_csv(os.path.join(results, virtual), sep='\t')
        c = pd.read_csv(os.path.join(results, capture), sep='\t')
        ax.plot(v['steps'], v['concordance'], marker='o', color=color, lw=1.5, ms=4.5, label=label + ', virtual')
        ax.plot(c['steps'], c['concordance'], marker='s', color=color, lw=1.2, ms=4, ls=':', label=label + ', capture')
    ax.set_xticks(range(1, 6))
    ax.set_xlabel('graph distance (steps)')
    ax.set_ylabel('P(graph pair > matched pair)')
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title('b', loc='left', fontweight='bold')
    fig.tight_layout()
    fig.savefig(output)
    print(output)


if __name__ == '__main__':
    main()
