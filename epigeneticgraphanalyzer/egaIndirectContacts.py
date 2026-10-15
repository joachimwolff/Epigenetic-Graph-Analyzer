"""egaIndirectContacts: indirect contacts and two-hop triplets of an interaction graph."""

import logging

from ._cli import base_parser, setup_logging, write_table
from .graph import read_graph
from .indirect import direct_contacts, indirect_contacts, two_hop_triplets

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Infer indirect contacts (pairs sharing a neighbour without a direct significant '
                         'interaction) and list the two-hop paths (triplets) of a graph.')
    parser.add_argument('--graph', '-g', required=True, help='GraphML of egaCreateGraph')
    parser.add_argument('--maxDistance', type=int, default=2000000, help='Maximum distance of a pair (bp)')
    parser.add_argument('--minDistance', type=int, default=0, help='Minimum distance of a pair (bp)')
    parser.add_argument('--maxMediatorDegree', type=int, default=None,
                        help='Ignore mediators with a larger degree (limits hub explosion)')
    parser.add_argument('--triplets', action='store_true', help='Also write the open two-hop triplets')
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>_indirect.tsv, _direct.tsv and optionally _triplets.tsv')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    write_table(indirect_contacts(graph, max_distance=args.maxDistance, min_distance=args.minDistance,
                                  max_mediator_degree=args.maxMediatorDegree), args.outputPrefix + '_indirect.tsv')
    write_table(direct_contacts(graph, max_distance=args.maxDistance), args.outputPrefix + '_direct.tsv')
    if args.triplets:
        write_table(two_hop_triplets(graph, max_distance=args.maxDistance,
                                     max_mediator_degree=args.maxMediatorDegree),
                    args.outputPrefix + '_triplets.tsv')


if __name__ == '__main__':
    main()
