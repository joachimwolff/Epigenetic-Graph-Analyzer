"""Paths from a promoter to its concordant step-2 nodes, drawn in three columns (promoter, mediating fragments,
step-2 nodes), one PDF.

Mediating fragments within --merge bp of each other are drawn as one node, labelled with their range and number.
The union graph of a differential comparison (egaDifferential) is read; the step-2 nodes whose active marks
increase at least two-fold (log2 change of (signal + 1) >= 1 for any of the given marks) are kept, together with
every fragment on a shortest path from the promoter to them. Nodes are coloured by the log2 change of the first
mark, edges by status (gained, lost, shared) with a width by the strength in condition B. Mediating fragments are
ordered by position.

Usage: python figure_hoxa3_paths.py <union.graphml> <promoter node prefix> <output.pdf>
           [--marks H3K27ac H3K4me2] [--conditions Md PA2] [--merge 50000] [--title TEXT]
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

STATUS_STYLE = {'gained': ('#b2182b', '--'), 'lost': ('#2166ac', '--'), 'shared': ('#8c8c8c', '-')}


def log2_change(data, mark, a, b):
    return float(np.log2((data.get('{}_{}'.format(mark, b), 0.0) + 1) / (data.get('{}_{}'.format(mark, a), 0.0) + 1)))


def node_label(node, data):
    genes = (data.get('genes') or '').split(',')[0]
    if genes:
        return genes
    start, end = int(data['start']), int(data['end'])
    return '{:.3f} Mb'.format((start + end) / 2e6)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('promoter')
    parser.add_argument('output')
    parser.add_argument('--marks', nargs='+', default=['H3K27ac', 'H3K4me2'])
    parser.add_argument('--conditions', nargs=2, default=['Md', 'PA2'])
    parser.add_argument('--merge', type=int, default=50000)
    parser.add_argument('--title', default='')
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.graph)
    source = [n for n in graph if n.startswith(args.promoter)][0]
    distance = nx.single_source_shortest_path_length(graph, source, cutoff=2)
    targets = [n for n, k in distance.items() if k == 2
               and max(log2_change(graph.nodes[n], m, a, b) for m in args.marks) >= 1]
    mediators = sorted({p[1] for t in targets for p in nx.all_shortest_paths(graph, source, t)},
                       key=lambda n: int(graph.nodes[n]['start']))
    targets = sorted(targets, key=lambda n: int(graph.nodes[n]['start']))
    # mediating fragments within --merge bp form one group
    groups, current = [], []
    for n in mediators:
        mid = (int(graph.nodes[n]['start']) + int(graph.nodes[n]['end'])) // 2
        if current and mid - current[-1][1] > args.merge:
            groups.append(current)
            current = []
        current.append((n, mid))
    if current:
        groups.append(current)
    groups = [[n for n, _ in grp] for grp in groups]
    norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    cmap = plt.get_cmap('RdBu_r')

    def mean_change(nodes):
        return float(np.mean([log2_change(graph.nodes[n], args.marks[0], a, b) for n in nodes]))

    def ys(k):
        return np.linspace(1, 0, k) if k > 1 else [0.5]

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    src_pos = (0.0, 0.5)
    grp_pos = {i: (1.0, y) for i, y in enumerate(ys(len(groups)))}
    tgt_pos = {t: (2.0, y) for t, y in zip(targets, ys(len(targets)))}

    def draw_edges(edges, p0, p1):
        status = [graph.edges[u, v].get('status', 'shared') for u, v in edges]
        for kind in ('shared', 'lost', 'gained'):
            k = status.count(kind)
            if not k:
                continue
            colour, style = STATUS_STYLE[kind]
            offset = {'shared': 0.0, 'gained': 0.012, 'lost': -0.012}[kind]
            ax.plot([p0[0], p1[0]], [p0[1] + offset, p1[1] + offset], color=colour, linestyle=style,
                    linewidth=0.5 + 0.5 * k, alpha=0.8, zorder=1)

    for i, grp in enumerate(groups):
        draw_edges([(source, n) for n in grp if graph.has_edge(source, n)], src_pos, grp_pos[i])
        for t in targets:
            draw_edges([(n, t) for n in grp if graph.has_edge(n, t)], grp_pos[i], tgt_pos[t])
    ax.scatter([src_pos[0]], [src_pos[1]], s=170, marker='s', c=[cmap(norm(mean_change([source])))],
               edgecolors='black', linewidths=0.6, zorder=3)
    ax.text(src_pos[0] - 0.09, src_pos[1], node_label(source, graph.nodes[source]), ha='right', va='center',
            fontsize=9, fontweight='bold')
    for i, grp in enumerate(groups):
        x, y = grp_pos[i]
        ax.scatter([x], [y], s=60 + 25 * len(grp), marker='o', c=[cmap(norm(mean_change(grp)))],
                   edgecolors='black', linewidths=0.6, zorder=3)
        lo = min(int(graph.nodes[n]['start']) for n in grp) / 1e6
        hi = max(int(graph.nodes[n]['end']) for n in grp) / 1e6
        text = '{:.2f}-{:.2f} Mb'.format(lo, hi) if hi - lo > 0.005 else '{:.2f} Mb'.format(lo)
        ax.text(x, y - 0.075, '{}, {} fragment{}'.format(text, len(grp), 's' if len(grp) > 1 else ''),
                ha='center', va='top', fontsize=6.5, color='#333333', zorder=4,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.9))
    for t, (x, y) in tgt_pos.items():
        bait = int(graph.nodes[t].get('bait', 0) or 0) == 1
        ax.scatter([x], [y], s=170, marker='s' if bait else 'o', c=[cmap(norm(mean_change([t])))],
                   edgecolors='black', linewidths=0.6, zorder=3)
        ax.text(x + 0.09, y, node_label(t, graph.nodes[t]), ha='left', va='center', fontsize=9)
    for x, text in ((0.0, 'promoter\n'), (1.0, 'step 1\nmediating fragments'), (2.0, 'step 2\nconcordant promoters')):
        ax.text(x, 1.10, text, ha='center', va='bottom', fontsize=8.5)
    ax.set_xlim(-0.75, 2.75)
    ax.set_ylim(-0.16, 1.22)
    ax.axis('off')
    if args.title:
        ax.set_title(args.title, fontsize=9, loc='left')
    handles = [Line2D([0], [0], color=c, linestyle=st, label=k + (' in ' + b if k != 'shared' else ''))
               for k, (c, st) in STATUS_STYLE.items()]
    handles += [Line2D([0], [0], marker='s', color='w', markerfacecolor='white', markeredgecolor='black',
                       markersize=7, label='baited promoter'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',
                       markersize=6, label='other end')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.13), ncol=5, fontsize=7, frameon=False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.025, pad=0.01)
    bar.set_label('{} log$_2$ change {}/{}'.format(args.marks[0], b, a), fontsize=8)
    bar.ax.tick_params(labelsize=7)
    fig.savefig(args.output, bbox_inches='tight')
    print('{}: {} mediating fragments, {} step-2 nodes: {}'.format(
        args.output, len(mediators), len(targets), ', '.join(node_label(t, graph.nodes[t]) for t in targets)))
    print('groups: ' + '; '.join('{} fragments'.format(len(g)) for g in groups))


if __name__ == '__main__':
    main()
