"""Multi-way support of the candidate explaining paths of one gene, two PDFs.

<out>_paths_matrix.pdf  every explaining path of the gene (promoter - mediator - element, from
                        36_gm12878_multiway_paths.py) by multi-way data set: colour, log2((observed + 1) / (control
                        mean + 1)); number, groups holding all three fragments; black frame, P < 0.05 against the
                        controls shifted by whole restriction fragments; rows 'mediator > element', in two blocks:
                        paths present in the condition of the multi-way data (both interactions called there) and
                        paths present only in the other condition
<out>_paths_matrix_supported.pdf  the same for the paths with P < 0.05 in at least one data set
<out>_paths_groups.pdf  for the paths with P < 0.05 in at least one data set (up to --maxPaths), every multi-way group
                        that holds the three fragments, drawn along the chromosome: a line from the first to the last
                        fragment of the group within the window, ticks at its fragments, the promoter, mediator and
                        element fragments marked

Usage: python figure_gene_path_multiway.py <cmp prefix> <multiway_paths prefix> <promoter node> <out prefix>
           --data LABEL=FILE.npz ... [--name ASUN] [--maxPaths 4] [--padding 50000]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_gene_multiway import ChromGroups, intersect  # noqa: E402

from epigeneticgraphanalyzer.multiway import MultiwayData  # noqa: E402

plt.rcParams.update({'font.size': 7, 'pdf.fonttype': 42})
INK, MUTED = '#222222', '#7f7f7a'
ROLE_COLOURS = {'promoter': '#222222', 'mediator': '#7f7f7a', 'element': '#b07f2a'}
DATA_COLOURS = ['#2a6fb0', '#b2182b', '#4d9221', '#7b3294']
NORM = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
CMAP = plt.get_cmap('RdBu_r')


def middle(attrs):
    return (int(attrs['start']) + int(attrs['end'])) / 2


def node_label(graph, node, max_genes=1):
    genes = [g for g in str(graph.nodes[node].get('genes') or '').split(',') if g and not g.startswith(('RP', 'CT', 'AC'))]
    mb = '{:.3f} Mb'.format(middle(graph.nodes[node]) / 1e6)
    return mb + (' ' + ','.join(genes[:max_genes]) if genes else '')


def panel_matrix(graph, paths, data_sets, path, condition='GM12878', other='naive B'):
    paths = paths.copy()
    paths['mediator_pos'] = paths.mediator.map(lambda n: middle(graph.nodes[n]))
    paths['target_pos'] = paths.target.map(lambda n: middle(graph.nodes[n]))
    order = paths.drop_duplicates(['mediator', 'target']).sort_values(['present', 'mediator_pos', 'target_pos'],
                                                                       ascending=[False, True, True])
    blocks = [('present in ' + condition, order[order.present]), ('present only in ' + other, order[~order.present])]
    layout, y = [], 0
    for title, block in blocks:
        if not len(block):
            continue
        layout.append(('title', title, y))
        y += 1
        for m, t in zip(block.mediator, block.target):
            layout.append(('row', (m, t), y))
            y += 1
        y += 0.4
    fig, ax = plt.subplots(figsize=(0.36 * len(data_sets) + 1.6, 0.14 * y + 0.9))
    labels, ticks = [], []
    for kind, item, yy in layout:
        if kind == 'title':
            ax.text(-0.1, yy + 0.6, item, ha='right', va='center', fontsize=7.5, fontweight='bold')
            continue
        m, t = item
        ticks.append(yy + 0.5)
        labels.append('{} > {}'.format(node_label(graph, m), node_label(graph, t)))
        for j, d in enumerate(data_sets):
            hit = paths[(paths.mediator == m) & (paths.target == t) & (paths.data_set == d)]
            if not len(hit) or pd.isna(hit.control_mean.iloc[0]):
                ax.add_patch(Rectangle((j, yy), 1, 1, facecolor='#eeeeee', edgecolor='white', lw=0.5))
                continue
            r = hit.iloc[0]
            ax.add_patch(Rectangle((j, yy), 1, 1, facecolor=CMAP(NORM(r.log2_ratio)), edgecolor='white', lw=0.5))
            if r.p_empirical < 0.05:
                ax.add_patch(Rectangle((j + 0.06, yy + 0.06), 0.88, 0.88, fill=False, edgecolor=INK, lw=1.0))
            ax.text(j + 0.5, yy + 0.5, '{:d}'.format(int(r.observed)), ha='center', va='center', fontsize=6.5,
                    color=INK)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xticks(np.arange(len(data_sets)) + 0.5)
    ax.set_xticklabels(data_sets, fontsize=7, rotation=30, ha='right')
    ax.set_xlim(0, len(data_sets))
    ax.set_ylim(y, 0)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    height = fig.get_size_inches()[1]
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), ax=ax, orientation='horizontal',
                       fraction=0.12 / height, pad=0.55 / height, shrink=0.9, aspect=14)
    bar.set_label('log$_2$ (observed + 1) / (control + 1)', fontsize=6)
    bar.ax.tick_params(labelsize=6)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return [item for kind, item, _ in layout if kind == 'row']


def groups_of(data, chrom, spans):
    groups = ChromGroups(data, chrom)
    units = [groups.overlapping(s + 1, e - 1) for s, e in spans]
    holding = intersect(*[groups.groups(u) for u in units])
    first = data.bins.offsets[chrom]
    block = data.incidence[first:first + groups.n].tocsc()
    members = [block.indices[block.indptr[g]:block.indptr[g + 1]] for g in holding]
    return groups, units, members


def panel_groups(graph, chrom, selected, data, path, padding):
    if not selected:
        return
    fig, axes = plt.subplots(len(selected), 1, figsize=(4.2, 0.9 + 0.75 * len(selected)), squeeze=False)
    for ax, (m, t, promoter) in zip(axes[:, 0], selected):
        spans = [(int(graph.nodes[n]['start']), int(graph.nodes[n]['end'])) for n in (promoter, m, t)]
        lo, hi = min(s for s, _ in spans) - padding, max(e for _, e in spans) + padding
        y = 0
        for k, (label, d) in enumerate(data.items()):
            groups, units, members = groups_of(d, chrom, spans)
            for frags in members:
                pos = np.array([(groups.starts[f] + groups.ends[f]) / 2 for f in frags])
                pos = pos[(pos >= lo) & (pos <= hi)]
                ax.plot([pos.min(), pos.max()], [y, y], color=DATA_COLOURS[k], lw=0.8, zorder=1)
                ax.scatter(pos, [y] * len(pos), s=4, color=DATA_COLOURS[k], zorder=2)
                y += 1
            if members:
                ax.text(hi, y - 0.5 * len(members) - 0.5, '{} {}'.format(len(members), label), fontsize=6,
                        color=DATA_COLOURS[k], ha='left', va='center')
        for k, ((s, e), role) in enumerate(zip(spans, ('promoter', 'mediator', 'element'))):
            ax.axvspan(s, e, color=ROLE_COLOURS[role], alpha=0.25, lw=0, zorder=0)
            # labels staggered in height, so that close fragments do not overlap
            ax.text((s + e) / 2, -0.8 - 0.9 * (k % 2), role, fontsize=6, ha='center', va='top',
                    color=ROLE_COLOURS[role])
        ax.set_xlim(lo, hi)
        ax.set_ylim(-2.6, max(y, 1))
        ax.set_yticks([])
        ax.set_title('via {} to {}'.format(node_label(graph, m), node_label(graph, t)), fontsize=7, loc='left')
        ticks = ax.get_xticks()
        ax.set_xticks(ticks)
        ax.set_xticklabels(['{:.2f}'.format(v / 1e6) for v in ticks], fontsize=6)
        ax.set_xlim(lo, hi)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
    axes[-1, 0].set_xlabel('{} (Mb)'.format(chrom), fontsize=7)
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('prefix', 'paths', 'promoter', 'output'):
        parser.add_argument(name)
    parser.add_argument('--data', nargs='+', required=True, help='LABEL=FILE.npz, in the column order')
    parser.add_argument('--name', default=None)
    parser.add_argument('--maxPaths', type=int, default=4)
    parser.add_argument('--padding', type=int, default=50000)
    args = parser.parse_args()
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    triplets = pd.read_csv(args.paths + '_triplets.tsv', sep='\t')
    paths = triplets[(triplets.promoter == args.promoter) & (triplets.set == 'explaining')]
    labels = [item.split('=', 1)[0] for item in args.data]
    rows = panel_matrix(graph, paths, labels, args.output + '_paths_matrix.pdf')
    keep = paths.groupby(['mediator', 'target']).p_empirical.transform('min') < 0.05
    panel_matrix(graph, paths[keep], labels, args.output + '_paths_matrix_supported.pdf')
    supported = paths[paths.p_empirical < 0.05].groupby(['mediator', 'target']).agg(
        n=('data_set', 'nunique'), groups=('observed', 'sum'), present=('present', 'first')).reset_index()
    supported = supported.sort_values(['n', 'present', 'groups'], ascending=False).head(args.maxPaths)
    data = {label: MultiwayData.load(path) for label, path in (item.split('=', 1) for item in args.data)}
    chrom = graph.nodes[args.promoter]['chrom']
    selected = [(m, t, args.promoter) for m, t in zip(supported.mediator, supported.target)]
    panel_groups(graph, chrom, selected, data, args.output + '_paths_groups.pdf', args.padding)
    print('{}: {} paths, {} with P < 0.05 in at least one data set; shown {}'.format(
        args.name or args.promoter, len(rows), int((paths.groupby(['mediator', 'target']).p_empirical.min() < 0.05)
                                                   .sum()), [node_label(graph, t) for _, t, _ in selected]))


if __name__ == '__main__':
    main()
