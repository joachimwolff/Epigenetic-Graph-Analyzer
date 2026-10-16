"""egaValidateMultiway: test direct contacts, indirect contacts, triplets and communities against SPRITE or GAM."""

import logging
from collections import defaultdict

import pandas as pd

from ._cli import base_parser, setup_logging, write_table
from .genome import LiftOver
from .graph import read_graph
from .indirect import direct_contacts, indirect_contacts, two_hop_triplets
from .multiway import MultiwayData
from .validation import (LocusPools, confirm_pairs, domain_annotation, map_nodes_to_bins, pairs_to_bins, validate_communities,
                         validate_pairs, validate_triplets)

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Validate graph-derived direct and indirect contacts, two-hop triplets and communities '
                         'with ligation-free multi-way contact data against distance-matched controls.')
    parser.add_argument('--graph', '-g', required=True, help='GraphML of egaCreateGraph')
    parser.add_argument('--multiway', '-m', required=True, help='.npz of egaBuildMultiway')
    parser.add_argument('--chain', help='UCSC chain file to lift graph coordinates to the multi-way assembly')
    parser.add_argument('--maxDistance', type=int, default=2000000)
    parser.add_argument('--maxMediatorDegree', type=int, default=None)
    parser.add_argument('--controls', type=int, default=10, help='Matched controls per tested pair')
    parser.add_argument('--compartments', help='A/B compartments in the multi-way assembly (BED with names A*/B* or '
                                               'eigenvector bedGraph); controls match the compartments of both ends')
    parser.add_argument('--domains', help='Topological domains in the multi-way assembly (BED); controls match whether '
                                          'both ends lie in the same domain')
    parser.add_argument('--coverageClasses', type=int, default=0,
                        help='Also match pair controls on the coverage class (quantiles over graph loci) of both '
                             'ends; 0 matches on offset and locus pools only')
    parser.add_argument('--confirmPairs', action='store_true',
                        help='Test every direct and indirect pair individually against background pairs of the same '
                             'distance and coverage (writes <prefix>_confirm_*.tsv)')
    parser.add_argument('--maxPairs', type=int, default=500000,
                        help='Random subset size of direct and of indirect pairs (all pairs if fewer)')
    parser.add_argument('--maxTriplets', type=int, default=500000, help='Random subset size of triplets')
    parser.add_argument('--communityMinMembers', type=int, default=3,
                        help='A group supports a community if it holds at least this many of its bins')
    parser.add_argument('--skipPairs', action='store_true')
    parser.add_argument('--skipTriplets', action='store_true')
    parser.add_argument('--skipCommunities', action='store_true')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--outputPrefix', '-o', required=True)
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    liftover = LiftOver(args.chain) if args.chain else None
    node_bins = map_nodes_to_bins(graph, data, liftover)
    log.info('%d of %d nodes mapped to covered bins', len(node_bins), graph.number_of_nodes())

    direct = pairs_to_bins(direct_contacts(graph, args.maxDistance), node_bins)
    indirect = pairs_to_bins(indirect_contacts(graph, max_distance=args.maxDistance,
                                               max_mediator_degree=args.maxMediatorDegree), node_bins)
    excluded = set(zip(direct['bin_a'], direct['bin_b'])) | set(zip(indirect['bin_a'], indirect['bin_b']))
    log.info('%d direct and %d indirect pairs on distinct bins', len(direct), len(indirect))
    if len(direct) > args.maxPairs:
        direct = direct.sample(n=args.maxPairs, random_state=args.seed).reset_index(drop=True)
    if len(indirect) > args.maxPairs:
        indirect = indirect.sample(n=args.maxPairs, random_state=args.seed).reset_index(drop=True)
    pools = LocusPools(graph, node_bins, data)
    compartment, tad_id = None, None
    if args.compartments or args.domains:
        compartment, tad_id = domain_annotation(data.bins, args.compartments, args.domains)
        compartment = compartment if args.compartments else None
        tad_id = tad_id if args.domains else None
    max_offset = int(max(direct['offset'].max() if len(direct) else 1,
                         indirect['offset'].max() if len(indirect) else 1))

    summaries = []
    if not args.skipPairs:
        direct_table, summary = validate_pairs(direct, data, pools, excluded, 'direct', n_controls=args.controls,
                                               max_offset=max_offset, seed=args.seed,
                                               coverage_matched=args.coverageClasses, compartment=compartment,
                                               tad_id=tad_id)
        summaries.append(summary)
        write_table(direct_table, args.outputPrefix + '_direct_pairs.tsv.gz')
        indirect_table, summary = validate_pairs(indirect, data, pools, excluded, 'indirect',
                                                 n_controls=args.controls, max_offset=max_offset, seed=args.seed + 1,
                                                 coverage_matched=args.coverageClasses, compartment=compartment,
                                                 tad_id=tad_id)
        summaries.append(summary)
        write_table(indirect_table, args.outputPrefix + '_indirect_pairs.tsv.gz')

    if args.confirmPairs:
        confirmations = []
        for label, pairs, seed in (('direct', direct, args.seed + 10), ('indirect', indirect, args.seed + 11)):
            table, summary = confirm_pairs(pairs, data, pools, excluded, label, max_offset,
                                           n_classes=max(args.coverageClasses, 5), seed=seed)
            write_table(table, '{}_confirm_{}.tsv.gz'.format(args.outputPrefix, label))
            confirmations.append(summary)
        confirmations = pd.concat([c for c in confirmations if len(c)], ignore_index=True)
        confirmations.insert(0, 'assay', data.kind)
        write_table(confirmations, args.outputPrefix + '_confirm_summary.tsv')

    if not args.skipTriplets:
        triplets = two_hop_triplets(graph, max_distance=args.maxDistance,
                                    max_mediator_degree=args.maxMediatorDegree)
        triplet_table, summary = validate_triplets(triplets, node_bins, data, pools, graph,
                                                   max_triplets=args.maxTriplets, seed=args.seed + 2)
        summaries.append(summary)
        write_table(triplet_table, args.outputPrefix + '_triplets.tsv.gz')

    if not args.skipCommunities:
        members = defaultdict(set)
        for node, community in graph.nodes(data='community'):
            if community is not None and int(community) >= 0:
                members[int(community)].add(node)
        communities = [members[k] for k in sorted(members)]
        baits = {n for n, b in graph.nodes(data='bait') if b}
        community_table, summary = validate_communities(communities, node_bins, data,
                                                        min_members=args.communityMinMembers, seed=args.seed + 3,
                                                        bait_nodes=baits)
        summaries.append(summary)
        write_table(community_table, args.outputPrefix + '_communities.tsv')

    summaries = [s for s in summaries if len(s)]
    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    summary.insert(0, 'assay', data.kind)
    write_table(summary, args.outputPrefix + '_summary.tsv')


if __name__ == '__main__':
    main()
