"""egaCreateGraph: interaction graph, gene annotation and communities."""

import logging

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .differential import REPLICATE_MODES, consensus_scores
from .graph import (WEIGHT_CHOICES, build_graph, community_table, detect_communities, edge_table,
                    node_table, write_graph)
from .io import (INTERACTION_COLUMNS, INTERACTION_FORMATS, read_baits, read_interactions, read_tss_bed,
                 read_tss_from_gtf)

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Build the promoter capture Hi-C interaction graph, annotate baited promoters with genes '
                         'and partition the graph into communities.')
    parser.add_argument('--interactions', '-i', nargs='+', required=True,
                        help='Called interactions: CHiCAGO ibed or WashU text, BEDPE (bait first, score), the '
                             'Mifsud et al. table, or cool matrices read out at the baited bins')
    parser.add_argument('--format', '-f', choices=INTERACTION_FORMATS, default='ibed')
    parser.add_argument('--baits',
                        help="Baited fragments of the capture design (BED or CHiCAGO .baitmap); required for "
                             "--format cool, for other formats it sets the bait status and orientation of each pair")
    parser.add_argument('--coolField', default='count',
                        help="Pixel column of the cool matrices used as interaction score (e.g. count or a score)")
    parser.add_argument('--maxFdr', type=float, default=None,
                        help='Keep pairs with an FDR at most this value (inputs with an FDR column, e.g. HOMER)')
    parser.add_argument('--scoreThreshold', type=float, default=None,
                        help='Keep pairs with score >= threshold (e.g. CHiCAGO 5, or a minimum pixel value)')
    parser.add_argument('--minCount', type=int, default=0, help='Minimum read count of an edge (if available)')
    parser.add_argument('--replicates', choices=REPLICATE_MODES, default=None,
                        help="Treat the input files as replicates: 'all' keeps pairs present in every replicate and "
                             "reaching --scoreThreshold in at least one (others may be below it), 'intersect' pairs "
                             "reaching it in every replicate, 'union' in any; default: pool all files")
    parser.add_argument('--maxDistance', type=int, default=None, help='Maximum bait-other distance of an edge (bp)')
    parser.add_argument('--tss', help='Gene annotation: GTF(.gz) or BED of TSS/promoters')
    parser.add_argument('--geneTypes', nargs='*', default=['protein_coding', 'lincRNA', 'antisense'],
                        help='Gene types taken from a GTF')
    parser.add_argument('--weight', choices=WEIGHT_CHOICES, default='score', help='Edge weight')
    parser.add_argument('--communityMethod', choices=('leiden', 'louvain'), default='leiden')
    parser.add_argument('--resolution', type=float, default=1.0, help='Community resolution parameter')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>.graphml, _nodes.tsv, _edges.tsv, _communities.tsv')
    return parser.parse_args(args)


def load_interactions(paths, file_format, score_threshold=None, min_count=0, max_distance=None, baits=None,
                      cool_field='count', max_fdr=None):
    """Called interactions of all input files; pairs failing the thresholds are marked not significant."""
    table = read_interactions(paths, file_format, baits=baits, cool_field=cool_field, min_value=score_threshold)
    table['significant'] = True
    if score_threshold is not None:
        table['significant'] = table['significant'] & (table['score'] >= score_threshold)
    if max_fdr is not None and 'fdr' in table:
        table['significant'] = table['significant'] & (table['fdr'] <= max_fdr)
    if min_count:
        table['significant'] = table['significant'] & (table['count'].fillna(0) >= min_count)
    if max_distance is not None:
        table['significant'] = table['significant'] & (table['distance'] <= max_distance)
    missing = [c for c in INTERACTION_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError('interaction table misses columns {}'.format(missing))
    return table


def load_tss(path, gene_types):
    if path is None:
        return None
    lower = path.lower()
    if lower.endswith(('.gtf', '.gtf.gz', '.gff', '.gff3', '.gff.gz', '.gff3.gz')):
        return read_tss_from_gtf(path, gene_types)
    return read_tss_bed(path)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    baits = read_baits(args.baits) if args.baits else None
    if args.replicates and len(args.interactions) > 1:
        replicates = [read_interactions([path], args.format, baits=baits, cool_field=args.coolField)
                      for path in args.interactions]
        if args.maxFdr is not None:
            replicates = [t[t['fdr'] <= args.maxFdr] if 'fdr' in t else t for t in replicates]
        table = consensus_scores(replicates, args.replicates, args.scoreThreshold)
        if args.maxDistance is not None:
            table = table[table['distance'] <= args.maxDistance]
    else:
        table = load_interactions(args.interactions, args.format, args.scoreThreshold, args.minCount,
                                  args.maxDistance, baits=baits, cool_field=args.coolField, max_fdr=args.maxFdr)
    tss = load_tss(args.tss, args.geneTypes)
    graph = build_graph(table, tss=tss, weight=args.weight)
    communities = detect_communities(graph, method=args.communityMethod, resolution=args.resolution,
                                     seed=args.seed)
    prefix = args.outputPrefix
    write_graph(graph, ensure_parent(prefix + '.graphml'))
    write_table(node_table(graph), prefix + '_nodes.tsv')
    write_table(edge_table(graph), prefix + '_edges.tsv')
    write_table(community_table(graph, communities), prefix + '_communities.tsv')


if __name__ == '__main__':
    main()
