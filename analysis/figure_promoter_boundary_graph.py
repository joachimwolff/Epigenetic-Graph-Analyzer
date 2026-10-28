"""Interactions of one promoter on both sides of a domain boundary, drawn as a graph, one PDF.

The promoter of a differential union graph (egaDifferential) is drawn in the middle; its partners are grouped by
genomic position (fragments within --merge bp form one node) and placed on the side of the boundary they lie on,
the far side (across the boundary from the promoter) on the left, the near side on the right, ordered by position.
Edges are coloured by status (gained, lost, shared) with a width by their number; node size gives the number of
fragments, and the boundary is drawn as a vertical line between the two sides.

Usage: python figure_promoter_boundary_graph.py <union.graphml> <promoter node prefix> <boundary position> <output.pdf>
           [--conditions Md PA2] [--merge 50000] [--far-label TEXT] [--near-label TEXT]
"""

import argparse
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D

STATUS_STYLE = {'gained': ('#b2182b', '--'), 'lost': ('#2166ac', '--'), 'shared': ('#8c8c8c', '-')}


def middle(data):
    return (int(data['start']) + int(data['end'])) // 2


def group(nodes, graph, merge):
    groups, current = [], []
    for n in sorted(nodes, key=lambda n: middle(graph.nodes[n])):
        if current and middle(graph.nodes[n]) - middle(graph.nodes[current[-1]]) > merge:
            groups.append(current)
            current = []
        current.append(n)
    if current:
        groups.append(current)
    return groups


def label(grp, graph):
    lo = min(int(graph.nodes[n]['start']) for n in grp) / 1e6
    hi = max(int(graph.nodes[n]['end']) for n in grp) / 1e6
    return '{:.2f}-{:.2f} Mb'.format(lo, hi) if hi - lo > 0.005 else '{:.2f} Mb'.format(lo)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('promoter')
    parser.add_argument('boundary', type=int)
    parser.add_argument('output')
    parser.add_argument('--conditions', nargs=2, default=['Md', 'PA2'])
    parser.add_argument('--merge', type=int, default=50000)
    parser.add_argument('--far-label', default='across the boundary')
    parser.add_argument('--near-label', default='same side as the promoter')
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.graph)
    source = [n for n in graph if n.startswith(args.promoter)][0]
    promoter_side = middle(graph.nodes[source]) > args.boundary
    partners = list(graph.neighbors(source))
    far = [n for n in partners if (middle(graph.nodes[n]) > args.boundary) != promoter_side]
    near = [n for n in partners if (middle(graph.nodes[n]) > args.boundary) == promoter_side]
    far_groups, near_groups = group(far, graph, args.merge), group(near, graph, args.merge)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    src_pos = (0.0, 0.5)

    def positions(groups, x):
        ys = np.linspace(1, 0, len(groups)) if len(groups) > 1 else [0.5]
        return [(x, y) for y in ys]

    far_pos, near_pos = positions(far_groups, -1.1), positions(near_groups, 1.0)
    ax.axvline(-0.5, color='black', linestyle=':', linewidth=1.2, zorder=0)
    ax.text(-0.5, -0.06, 'boundary\n{:.3f} Mb'.format(args.boundary / 1e6), ha='center', va='top', fontsize=7.5)
    for groups, pos in ((far_groups, far_pos), (near_groups, near_pos)):
        for grp, (x, y) in zip(groups, pos):
            counts = Counter(graph.edges[source, n].get('status', 'shared') for n in grp)
            for kind, offset in (('shared', 0.0), ('gained', 0.013), ('lost', -0.013)):
                if counts[kind]:
                    colour, style = STATUS_STYLE[kind]
                    ax.plot([src_pos[0], x], [src_pos[1] + offset, y + offset], color=colour, linestyle=style,
                            linewidth=0.6 + 0.35 * counts[kind], alpha=0.85, zorder=1)
            ax.scatter([x], [y], s=40 + 18 * len(grp), c='white', edgecolors='black', linewidths=0.7, zorder=3)
            text = '{}, {} fragment{}'.format(label(grp, graph), len(grp), 's' if len(grp) > 1 else '')
            ax.text(x + (-0.08 if x < 0 else 0.08), y, text, ha='right' if x < 0 else 'left', va='center',
                    fontsize=6.5, color='#333333')
    ax.scatter([src_pos[0]], [src_pos[1]], s=200, marker='s', c='#f4a582', edgecolors='black', linewidths=0.8,
               zorder=4)
    name = (graph.nodes[source].get('genes') or '').split(',')[0]
    ax.text(src_pos[0], src_pos[1] - 0.07, name, ha='center', va='top', fontsize=9, fontweight='bold')
    ax.text(-1.1, 1.12, args.far_label, ha='center', va='bottom', fontsize=8.5)
    ax.text(1.0, 1.12, args.near_label, ha='center', va='bottom', fontsize=8.5)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-0.2, 1.22)
    ax.axis('off')
    handles = [Line2D([0], [0], color=c, linestyle=st, label=k + (' in ' + b if k != 'shared' else ''))
               for k, (c, st) in STATUS_STYLE.items()]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=7, frameon=False)
    fig.savefig(args.output, bbox_inches='tight')
    count = Counter((('far' if n in far else 'near'), graph.edges[source, n].get('status')) for n in partners)
    print('{}: {}'.format(args.output, dict(count)))


if __name__ == '__main__':
    main()
