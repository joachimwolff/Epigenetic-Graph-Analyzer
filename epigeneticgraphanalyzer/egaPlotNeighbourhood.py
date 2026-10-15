"""egaPlotNeighbourhood: draw the one- and two-hop neighbourhood of a gene."""

import logging

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ._cli import base_parser, ensure_parent, setup_logging  # noqa: E402
from .graph import neighbourhood_layers, read_graph  # noqa: E402

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Plot the neighbourhood of a gene: nodes on a genomic axis, arcs for direct contacts, '
                         'dashed arcs for indirect contacts of the gene promoter.')
    parser.add_argument('--graph', '-g', required=True)
    parser.add_argument('--gene', required=True, help='Gene name as annotated on the bait node')
    parser.add_argument('--steps', type=int, default=2, choices=(1, 2))
    parser.add_argument('--colorBy', help='Numeric node attribute used for node colour (e.g. H3K27ac)')
    parser.add_argument('--output', '-o', required=True, help='Image file (pdf, png, svg)')
    return parser.parse_args(args)


def find_gene_node(graph, gene):
    for node, genes in graph.nodes(data='genes'):
        if genes and gene in str(genes).split(','):
            return node
    raise ValueError('Gene {} is not annotated on any node'.format(gene))


def _arc(ax, x1, x2, height, **kwargs):
    xs = np.linspace(x1, x2, 50)
    centre = (x1 + x2) / 2
    half = abs(x2 - x1) / 2 or 1.0
    ys = height * np.sqrt(np.clip(1 - ((xs - centre) / half) ** 2, 0, 1))
    ax.plot(xs, ys, **kwargs)


def plot_neighbourhood(graph, source, steps, color_by, output):
    layers = neighbourhood_layers(graph, source, steps)
    nodes = {source} | set().union(*layers)
    sub = graph.subgraph(nodes)
    mid = {n: (graph.nodes[n]['start'] + graph.nodes[n]['end']) / 2 / 1e6 for n in nodes}
    fig, ax = plt.subplots(figsize=(8, 3.2))
    max_span = max(abs(mid[u] - mid[v]) for u, v in sub.edges) if sub.number_of_edges() else 1.0
    for u, v in sub.edges:
        on_source = source in (u, v)
        _arc(ax, mid[u], mid[v], abs(mid[u] - mid[v]) / max_span,
             color='#1f4e79' if on_source else '#9a9a9a', lw=1.2 if on_source else 0.6)
    if steps >= 2:
        for m in layers[1]:
            _arc(ax, mid[source], mid[m], -abs(mid[source] - mid[m]) / max_span, color='#b8541c', lw=0.8, ls='--')
    values = None
    if color_by:
        values = np.array([float(graph.nodes[n].get(color_by, 0.0)) for n in nodes])
    xs = np.array([mid[n] for n in nodes])
    markers = ['s' if graph.nodes[n].get('bait') else 'o' for n in nodes]
    for marker in ('s', 'o'):
        mask = np.array([m == marker for m in markers])
        if not mask.any():
            continue
        scatter = ax.scatter(xs[mask], np.zeros(mask.sum()), marker=marker, s=28, zorder=3,
                             c=values[mask] if values is not None else '#333333', cmap='viridis',
                             edgecolors='black', linewidths=0.4)
    if values is not None:
        fig.colorbar(scatter, ax=ax, label=color_by, fraction=0.04)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_yticks([])
    ax.set_xlabel('{} (Mb)'.format(graph.nodes[source]['chrom']))
    ax.set_title('{}: {} direct, {} indirect'.format(graph.nodes[source].get('genes', source), len(layers[0]),
                                                     len(layers[1]) if steps >= 2 else 0), fontsize=10)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(ensure_parent(output), dpi=300)
    plt.close(fig)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    plot_neighbourhood(graph, find_gene_node(graph, args.gene), args.steps, args.colorBy, args.output)


if __name__ == '__main__':
    main()
