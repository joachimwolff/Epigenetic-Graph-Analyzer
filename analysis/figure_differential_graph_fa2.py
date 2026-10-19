"""Differential graph of a locus in a ForceAtlas2 layout, one PDF.

The nodes within --steps of a promoter in the union graph of a differential comparison (egaDifferential) are laid
out with ForceAtlas2 (networkx, distributed attraction against hubs, node sizes against overlap, fixed seed), weighted by the larger strength of each
interaction. Shared interactions are coloured by their log2 fold change between the conditions, interactions called
in one condition only are dashed (gained red, lost blue). Baited promoters are squares labelled with their first
gene, other ends small circles.

Usage: python figure_differential_graph_fa2.py <union.graphml> <promoter node prefix> <output.pdf>
           [--steps 2] [--conditions Md PA2] [--seed 1] [--width 6.5] [--height 4.8]
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('graph')
    parser.add_argument('promoter')
    parser.add_argument('output')
    parser.add_argument('--steps', type=int, default=2)
    parser.add_argument('--conditions', nargs=2, default=['Md', 'PA2'])
    parser.add_argument('--seed', type=int, default=1)
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
                                node_size={n: (0.6 if int(graph.nodes[n].get('bait', 0) or 0) else 0.2) for n in sub})
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
        if bait and graph.nodes[n].get('genes'):
            ax.annotate(graph.nodes[n]['genes'].split(',')[0], (x, y), xytext=(4, 3), textcoords='offset points',
                        fontsize=7, zorder=4,
                        bbox=dict(boxstyle='round,pad=0.1', facecolor='white', edgecolor='none', alpha=0.8))
    ax.set_aspect('equal')
    ax.axis('off')
    handles = [Line2D([0], [0], color='#b2182b', linestyle=(0, (3, 2)), label='gained in ' + b),
               Line2D([0], [0], color='#2166ac', linestyle=(0, (3, 2)), label='lost in ' + b),
               Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=6, label='baited promoter'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markeredgecolor='black',
                      markersize=5, label='other end')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=4, fontsize=7, frameon=False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.01, shrink=0.7)
    bar.set_label('shared: log$_2$ change {}/{}'.format(b, a), fontsize=7)
    bar.ax.tick_params(labelsize=6)
    fig.savefig(args.output, bbox_inches='tight')
    counts = {s: sum(1 for _, _, d in sub.edges(data=True) if d.get('status') == s) for s in order}
    print('{}: {} nodes, edges {}'.format(args.output, sub.number_of_nodes(), counts))


if __name__ == '__main__':
    main()
