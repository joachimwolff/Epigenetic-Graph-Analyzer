import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description='Reconstruct restriction cut-site intervals from chromosome sizes')
    parser.add_argument('--chromSizes', '-c', required=True, help='Chromosome sizes file (chrom<TAB>size)')
    parser.add_argument('--resolution', '-r', required=True, type=int, help='Interval size')
    parser.add_argument('--output', '-o', required=True, help='Output BED path')
    return parser.parse_args()


def main():
    args = parse_arguments()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(args.chromSizes, 'r', encoding='utf-8') as in_handle, open(args.output, 'w', encoding='utf-8') as out_handle:
        for line in in_handle:
            line = line.strip()
            if not line:
                continue
            chrom, size_s = line.split('\t')[:2]
            size = int(size_s)
            for start in range(0, size, args.resolution):
                end = min(size, start + args.resolution)
                out_handle.write(f'{chrom}\t{start}\t{end}\n')


if __name__ == '__main__':
    main()
