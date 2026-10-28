"""Differential graph of a locus in a ForceAtlas2 layout, one PDF.

The nodes within --steps of a promoter in the union graph of a differential comparison (egaDifferential) are laid
out with ForceAtlas2 (networkx, distributed attraction against hubs, node sizes against overlap, fixed seed), weighted by the larger strength of each
interaction. Shared interactions are coloured by their log2 fold change between the conditions, interactions called
in one condition only are dashed (gained red, lost blue). Baited promoters are squares labelled with their first
genes (at most two, LOC entries skipped), so that a gene with start sites on several fragments is told apart;
other ends are small circles.

Usage: python figure_differential_graph_fa2.py <union.graphml> <promoter node prefix> <output.pdf>
           [--steps 2] [--conditions Md PA2] [--seed 1] [--label-space 3.0] [--width 6.5] [--height 4.8]
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


def repel_labels(fig, ax, labels, pos, iterations=300):
    """Move overlapping text labels apart in display space and connect moved labels to their node."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_data = ax.transData.inverted()
    offsets = {n: np.array([4.0, 3.0]) for n, _ in labels}
    anchor = {n: ax.transData.transform(pos[n]) for n, _ in labels}
    for _ in range(iterations):
        moved = False
        for n, text in labels:
            text.set_position(to_data.transform(anchor[n] + offsets[n]))
        boxes = {n: text.get_window_extent(renderer).expanded(1.05, 1.15) for n, text in labels}
        names = [n for n, _ in labels]
        for i, u in enumerate(names):
            for v in names[i + 1:]:
                if boxes[u].overlaps(boxes[v]):
                    shift = np.array([0.0, 2.0]) if boxes[u].y0 >= boxes[v].y0 else np.array([0.0, -2.0])
                    offsets[u] += shift
                    offsets[v] -= shift
                    moved = True
        if not moved:
            break
    for n, text in labels:
        text.set_position(to_data.transform(anchor[n] + offsets[n]))
        if abs(offsets[n][1] - 3.0) > 6:
            end = to_data.transform(anchor[n] + offsets[n])
            ax.plot([pos[n][0], end[0]], [pos[n][1], end[1]], color='#666666', linewidth=0.4, zorder=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('promoter')
    parser.add_argument('output')
    parser.add_argument('--steps', type=int, default=2)
    parser.add_argument('--conditions', nargs=2, default=['Md', 'PA2'])
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--label-space', type=float, default=3.0,
                        help='Size of baited promoters in the layout, which keeps room for their labels')
    parser.add_argument('--width', type=float, default=6.5)
    parser.add_argument('--height', type=float, default=4.8)
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.graph)
    source = [n for n in graph if n.startswith(args.promoter)][0]
    nodes = nx.single_source_shortest_path_length(graph, source, cutoff=args.steps)
    sub = nx.Graph(graph.subgraph(nodes))
    for u, v, data in sub.edges(data=True):
        data['layout_weight'] = max(float(data.get('r_a', 0) or 0), float(data.get('r_b', 0) or 0), 0.1)
    pos = nx.forceatlas2_layout(sub, max_iter=3000, scaling_ratio=12.0, gravity=0.3, linlog=False,
                                distributed_action=True, weight='layout_weight', seed=args.seed,
                                node_size={n: (args.label_space if int(graph.nodes[n].get('bait', 0) or 0) else 0.2) for n in sub})
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = plt.get_cmap('RdBu_r')
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    order = {'shared': 0, 'lost': 1, 'gained': 2}
    for u, v, data in sorted(sub.edges(data=True), key=lambda e: order.get(e[2].get('status', 'shared'), 0)):
        status = data.get('status', 'shared')
        width = 0.3 + 1.0 * min(data['layout_weight'], 3.0) / 3.0
        if status == 'gained':
            colour, style = '#b2182b', (0, (3, 2))
        elif status == 'lost':
            colour, style = '#2166ac', (0, (3, 2))
        else:
            colour, style = cmap(norm(float(data.get('log2FC', 0) or 0))), '-'
        (x0, y0), (x1, y1) = pos[u], pos[v]
        ax.plot([x0, x1], [y0, y1], color=colour, linestyle=style, linewidth=width, alpha=0.7, zorder=1)
    for n, (x, y) in pos.items():
        bait = int(graph.nodes[n].get('bait', 0) or 0) == 1
        ax.scatter([x], [y], s=42 if bait else 12, marker='s' if bait else 'o', c='black' if bait else 'white',
                   edgecolors='black', linewidths=0.5, zorder=3)
    ax.set_aspect('equal')
    ax.axis('off')
    labels = []
    for n, (x, y) in pos.items():
        if int(graph.nodes[n].get('bait', 0) or 0) == 1 and graph.nodes[n].get('genes'):
            name = '/'.join([g for g in graph.nodes[n]['genes'].split(',') if not g.startswith('LOC')][:2])
            labels.append((n, ax.text(x, y, name, fontsize=7, ha='left', va='bottom', zorder=4,
                                      bbox=dict(boxstyle='round,pad=0.1', facecolor='white', edgecolor='none',
                                                alpha=0.85))))
    repel_labels(fig, ax, labels, pos)
    handles = [Line2D([0], [0], color='#b2182b', linestyle=(0, (3, 2)), label='gained in ' + b),
               Line2D([0], [0], color='#2166ac', linestyle=(0, (3, 2)), label='lost in ' + b),
               Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=6, label='baited promoter'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',
                      markersize=5, label='other end')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=7, frameon=False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', fraction=0.035,
                       pad=0.12, shrink=0.45)
    bar.set_label('shared: log$_2$ change {}/{}'.format(b, a), fontsize=7)
    bar.ax.tick_params(labelsize=6)
    fig.savefig(args.output, bbox_inches='tight')
    counts = {s: sum(1 for _, _, d in sub.edges(data=True) if d.get('status') == s) for s in order}
    print('{}: {} nodes, edges {}'.format(args.output, sub.number_of_nodes(), counts))


if __name__ == '__main__':
    main()
