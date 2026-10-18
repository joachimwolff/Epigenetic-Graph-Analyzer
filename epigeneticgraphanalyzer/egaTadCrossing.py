"""egaTadCrossing: TAD / domain boundaries crossed by the gained, lost and shared contacts of a differential run."""

import logging
import sys

from . import tads as td
from ._cli import base_parser, setup_logging

log = logging.getLogger(__name__)


def add_tad_arguments(parser, required=False):
    """--tads and its options, shared by egaTadCrossing, egaDifferential, egaExplain, egaReport and egaExplorer."""
    parser.add_argument('--tads', required=required, default=None,
                        help='TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, '
                             '_tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools '
                             'insulation table, BED of boundaries or domains (optional score column), or egaBoundaries '
                             'output')
    parser.add_argument('--tadFormat', choices=td.FORMATS, default='auto', help='Format of --tads')
    parser.add_argument('--tadWindow', type=int, default=None,
                        help='cooltools: insulation window (default: the largest in the table)')
    parser.add_argument('--tadChain', default=None,
                        help='Chain file from the assembly of --tads to that of the graph (UCSC chromosome names)')
    parser.add_argument('--tadRegions', default=None,
                        help='BED of named regions (chrom start end name) used to name the partner regions')
    parser.add_argument('--tadMinStrength', type=float, default=0.0,
                        help='Only boundaries at or above this strength percentile (0-1) count as crossed')
    return parser


def tad_kwargs(args):
    return dict(fmt=args.tadFormat, chain=args.tadChain, window=args.tadWindow, regions=args.tadRegions,
                min_strength_pct=args.tadMinStrength)


def parse_arguments(args=None):
    parser = base_parser('Per edge of a differential run (egaDifferential): the TAD boundaries it crosses, the '
                         'strongest of them and same- vs cross-domain, per condition and by status (gained / lost / '
                         'shared); per promoter the crossing contacts, the boundary they cross, their partner region '
                         'and a sentence; genome-wide ranking of promoters (Fisher: gained / lost among crossing vs '
                         'non-crossing edges, BH) and of boundaries (gained / lost crossing edges vs the expectation '
                         'from edges of the same length).')
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential')
    add_tad_arguments(parser, required=True)
    parser.add_argument('--names', nargs=2, default=None, help='Condition names (default: from the run)')
    parser.add_argument('--gene', nargs='*', default=[], help='Print the sentence of these genes')
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>_edges.tsv.gz, _summary.tsv, _promoters.tsv, _boundaries.tsv, '
                             '_genes.tsv')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    analysis = td.Analysis(args.prefix, args.tads, names=args.names, **tad_kwargs(args))
    promoters, bounds = analysis.write(args.outputPrefix)
    log.info('%d promoters tested, %d with gained q < 0.05; %d boundaries, %d with gained q < 0.05',
             promoters['gained_p'].notna().sum(), (promoters['gained_q'] < 0.05).sum(), len(bounds),
             (bounds['gained_q'] < 0.05).sum())
    for gene in args.gene:
        for node in analysis.gene_nodes(gene):
            row = analysis.crossing.gene_row(node, gene)
            sys.stdout.write('{}\t{}\n'.format(node, row['sentence']))
    return analysis


if __name__ == '__main__':
    main()
