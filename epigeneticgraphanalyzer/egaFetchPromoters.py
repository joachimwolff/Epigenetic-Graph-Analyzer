import argparse
import os

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Create promoter BED from a gene annotation table')
    parser.add_argument('--annotation', '-a', required=True, help='Input annotation TSV/CSV with gene coordinates')
    parser.add_argument('--output', '-o', required=True, help='Output BED file')
    parser.add_argument('--upstream', type=int, default=2000, help='Upstream window from TSS (default: 2000)')
    parser.add_argument('--downstream', type=int, default=500, help='Downstream window from TSS (default: 500)')
    parser.add_argument('--chromCol', default='Chromosome', help='Chromosome column name')
    parser.add_argument('--startCol', default='Start', help='Gene start column name')
    parser.add_argument('--endCol', default='End', help='Gene end column name')
    parser.add_argument('--strandCol', default='Strand', help='Strand column name')
    parser.add_argument('--nameCol', default='Name', help='Gene name column name')
    return parser.parse_args()


def main():
    args = parse_arguments()
    df = pd.read_csv(args.annotation, sep=None, engine='python')

    required = [args.chromCol, args.startCol, args.endCol, args.strandCol, args.nameCol]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'Missing required columns: {missing}')

    promoters = []
    for _, row in df.iterrows():
        chrom = str(row[args.chromCol])
        start = int(row[args.startCol])
        end = int(row[args.endCol])
        strand = str(row[args.strandCol])
        name = str(row[args.nameCol])

        tss = start if strand != '-' else end
        prom_start = max(0, tss - args.upstream)
        prom_end = max(prom_start + 1, tss + args.downstream)
        promoters.append((chrom, prom_start, prom_end, name, 0, strand))

    out_df = pd.DataFrame(promoters, columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    out_df.to_csv(args.output, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
