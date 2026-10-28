"""egaBuildMultiway: convert SPRITE clusters or GAM segregation tables into a binned incidence matrix."""

import logging

from ._cli import base_parser, ensure_parent, setup_logging
from .genome import IntervalIndex, read_chrom_sizes
from .multiway import read_gam_segregation, read_sprite_clusters

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Bin SPRITE clusters or GAM segregation tables into a sparse bin x group incidence '
                         'matrix (.npz) used by egaValidateMultiway.')
    parser.add_argument('--type', choices=('sprite', 'gam'), required=True)
    parser.add_argument('--input', nargs='+', required=True,
                        help='SPRITE cluster file, or one GAM segregation table per replicate')
    parser.add_argument('--chromSizes', required=True, help='Chromosome sizes of the data assembly')
    parser.add_argument('--binSize', type=int, default=25000, help='SPRITE bin size (GAM uses its windows)')
    parser.add_argument('--fragments', help='SPRITE: BED of intervals (e.g. restriction fragments of the capture '
                                            'Hi-C design); reads are assigned to the interval containing them '
                                            'instead of fixed bins')
    parser.add_argument('--minClusterSize', type=int, default=2, help='SPRITE: minimum reads per cluster')
    parser.add_argument('--maxClusterSize', type=int, default=1000, help='SPRITE: maximum reads per cluster')
    parser.add_argument('--output', '-o', required=True, help='Output .npz')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    chrom_sizes = read_chrom_sizes(args.chromSizes)
    if args.type == 'sprite':
        if len(args.input) != 1:
            raise ValueError('SPRITE expects one cluster file')
        index = IntervalIndex.from_bed(args.fragments) if args.fragments else None
        data = read_sprite_clusters(args.input[0], chrom_sizes, args.binSize,
                                    min_size=args.minClusterSize, max_size=args.maxClusterSize, index=index)
    else:
        data = read_gam_segregation(args.input, chrom_sizes)
    data.save(ensure_parent(args.output))


if __name__ == '__main__':
    main()
