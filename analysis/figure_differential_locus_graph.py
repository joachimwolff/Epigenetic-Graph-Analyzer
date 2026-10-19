"""Differential graph around a promoter, laid out along the genome, one PDF.

Nodes within --steps of the promoter in the union graph of a differential comparison (egaDifferential) are placed
in rows by their graph distance to the promoter (promoter, step 1, step 2, ...) and, within a row, evenly spaced in
the order of their genomic position (left to right), so that neighbouring fragments do not overlap.
Shared interactions are coloured by their log2 fold change between the conditions, interactions called in one
condition only are dashed (gained red, lost blue); line width follows the larger strength of the two conditions.
Baited promoters are squares and labelled with their first gene, other ends are circles. Optional vertical lines
mark the position of domain boundaries (--boundaries, in bp) within each row.

Usage: python figure_differential_locus_graph.py <union.graphml> <promoter node prefix> <output.pdf>
           [--steps 2] [--conditions Md PA2] [--boundaries 51175000 ...] [--width 7.0] [--height 3.2]
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


def middle(data):
    return (int(data['start']) + int(data['end'])) / 2e6


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('promoter')
    parser.add_argument('output')
    parser.add_argument('--steps', type=int, default=2)
    parser.add_argument('--conditions', nargs=2, default=['Md', 'PA2'])
    parser.add_argument('--boundaries', type=int, nargs='*', default=[])
    parser.add_argument('--width', type=float, default=7.0)
    parser.add_argument('--height', type=float, default=3.2)
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.graph)
    source = [n for n in graph if n.startswith(args.promoter)][0]
    distance = nx.single_source_shortest_path_length(graph, source, cutoff=args.steps)
    sub = graph.subgraph(distance)
    # rows from the top: promoter, step 1, step 2, ...; a small vertical spread separates neighbouring fragments
    pos = {}
    for k in range(args.steps + 1):
        nodes = sorted((n for n, d in distance.items() if d == k), key=lambda n: middle(graph.nodes[n]))
        xs = np.linspace(0, 1, len(nodes)) if len(nodes) > 1 else [0.5]
        for x, n in zip(xs, nodes):
            pos[n] = (x, -k)
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = plt.get_cmap('RdBu_r')
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    for boundary in args.boundaries:
        for k in range(1, args.steps + 1):
            row = [n for n, d in distance.items() if d == k]
            left = [pos[n][0] for n in row if middle(graph.nodes[n]) < boundary / 1e6]
            right = [pos[n][0] for n in row if middle(graph.nodes[n]) >= boundary / 1e6]
            if left and right:
                x = (max(left) + min(right)) / 2
                ax.plot([x, x], [-k - 0.3, -k + 0.3], color='black', linestyle=':', linewidth=1.0, zorder=0)
                ax.text(x, -k + 0.32, 'boundary', ha='center', va='bottom', fontsize=6)
    order = {'shared': 0, 'lost': 1, 'gained': 2}
    edges = sorted(sub.edges(data=True), key=lambda e: order.get(e[2].get('status', 'shared'), 0))
    for u, v, data in edges:
        status = data.get('status', 'shared')
        strength = max(float(data.get('r_a', 0) or 0), float(data.get('r_b', 0) or 0))
        width = 0.3 + 1.0 * min(strength, 3.0) / 3.0
        if status == 'gained':
            colour, style = '#b2182b', (0, (3, 2))
        elif status == 'lost':
            colour, style = '#2166ac', (0, (3, 2))
        else:
            colour, style = cmap(norm(float(data.get('log2FC', 0) or 0))), '-'
        (x0, y0), (x1, y1) = pos[u], pos[v]
        rad = 0.4 if abs(y0 - y1) < 0.3 else 0.0
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-', connectionstyle='arc3,rad={}'.format(rad),
                                     color=colour, linestyle=style, linewidth=width, alpha=0.7, zorder=1))
    for n, (x, y) in pos.items():
        bait = int(graph.nodes[n].get('bait', 0) or 0) == 1
        size = 70 if n == source else (28 if bait else 9)
        ax.scatter([x], [y], s=size, marker='s' if bait else 'o', c='black' if bait else 'white',
                   edgecolors='black', linewidths=0.5, zorder=3)
    labelled = [n for n in pos if int(graph.nodes[n].get('bait', 0) or 0) == 1 and graph.nodes[n].get('genes')]
    for n in sorted(labelled, key=lambda n: pos[n][0]):
        x, y = pos[n]
        name = graph.nodes[n]['genes'].split(',')[0]
        ax.text(x, y + (0.16 if n == source else -0.12), name, ha='center' if n == source else 'right',
                va='bottom' if n == source else 'top',
                fontsize=8 if n == source else 6.5, rotation=0 if n == source else 50, rotation_mode='anchor',
                fontweight='bold' if n == source else 'normal', zorder=4)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-args.steps - 0.9, 0.45)
    ax.set_yticks([-k for k in range(args.steps + 1)])
    ax.set_yticklabels(['promoter'] + ['step {}'.format(k) for k in range(1, args.steps + 1)], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['{:.2f} Mb'.format(min(middle(graph.nodes[n]) for n in pos)),
                        '{:.2f} Mb'.format(max(middle(graph.nodes[n]) for n in pos))], fontsize=7)
    ax.set_xlabel('fragments in genomic order ({})'.format(graph.nodes[source]['chrom']), fontsize=8)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='y', length=0)
    handles = [Line2D([0], [0], color='#b2182b', linestyle=(0, (3, 2)), label='gained in ' + b),
               Line2D([0], [0], color='#2166ac', linestyle=(0, (3, 2)), label='lost in ' + b),
               Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=6, label='baited promoter'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',
                      markersize=5, label='other end')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 1.18), ncol=4, fontsize=7, frameon=False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.02, pad=0.01)
    bar.set_label('shared: log$_2$ change {}/{}'.format(b, a), fontsize=7)
    bar.ax.tick_params(labelsize=6)
    fig.savefig(args.output, bbox_inches='tight')
    counts = {s: sum(1 for _, _, d in sub.edges(data=True) if d.get('status') == s) for s in order}
    print('{}: {} nodes, edges {}'.format(args.output, len(pos), counts))


if __name__ == '__main__':
    main()
