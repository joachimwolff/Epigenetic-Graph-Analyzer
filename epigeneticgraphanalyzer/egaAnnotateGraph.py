"""egaAnnotateGraph: attach chromatin signal and expression values to graph nodes."""

import logging

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .annotation import annotate_signals, gene_values_to_nodes, parse_track_arguments
from .graph import node_table, read_graph, write_graph
from .io import read_expression

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Annotate graph nodes with the mean bigWig signal over each fragment and baits with '
                         'gene expression values.')
    parser.add_argument('--graph', '-g', required=True, help='GraphML input')
    parser.add_argument('--tracks', '-t', nargs='*', default=[],
                        help='name=file.bw[,replicate.bw ...]; replicates are averaged')
    parser.add_argument('--flank', type=int, default=0, help='Extend fragments by this many bp on both sides')
    parser.add_argument('--expression', help='Expression table (gene identifier column + value columns)')
    parser.add_argument('--expressionColumns', nargs='*', default=None,
                        help='Value columns of the expression table to attach (default: all numeric)')
    parser.add_argument('--geneKey', choices=('gene_ids', 'genes'), default='gene_ids',
                        help='Node attribute used to match the expression table')
    parser.add_argument('--outputPrefix', '-o', required=True, help='Writes <prefix>.graphml and _nodes.tsv')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    tracks = parse_track_arguments(args.tracks)
    annotate_signals(graph, tracks, flank=args.flank)
    columns = []
    if args.expression:
        table = read_expression(args.expression)
        columns = args.expressionColumns or [c for c in table.columns if table[c].dtype.kind in 'fi']
        for column in columns:
            values = gene_values_to_nodes(graph, table[column], key=args.geneKey)
            for node, value in values.items():
                graph.nodes[node]['expression_' + column] = value
    write_graph(graph, ensure_parent(args.outputPrefix + '.graphml'))
    nodes = node_table(graph)
    for name in list(tracks) + ['expression_' + c for c in columns]:
        nodes[name] = [graph.nodes[n].get(name) for n in nodes['node']]
    write_table(nodes, args.outputPrefix + '_nodes.tsv')


if __name__ == '__main__':
    main()
