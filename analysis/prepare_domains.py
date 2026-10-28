"""Convert 4DN compartment and boundary files into the BED inputs of the domain-matched controls.

compartments: the eigenvector bigWig (positive = A, as phased by the 4DN pipeline) is written as a bedGraph of its
              intervals (chrom, start, end, eigenvector).
domains:      the insulation boundaries (BED) define domains between consecutive boundaries of a chromosome
              (from the end of one boundary to the start of the next).

Usage: python prepare_domains.py <compartments.bw> <boundaries.bed.gz> <out compartments.bedGraph> <out domains.bed>
"""

import gzip
import sys

import numpy as np
import pyBigWig


def main():
    bigwig, boundaries, out_compartments, out_domains = sys.argv[1:5]
    handle = pyBigWig.open(bigwig)
    with open(out_compartments, 'w') as out:
        for chrom, size in handle.chroms().items():
            for start, end, value in handle.intervals(chrom) or []:
                if np.isfinite(value) and value != 0:
                    out.write('{}\t{}\t{}\t{:.6g}\n'.format(chrom, start, end, value))
    by_chrom = {}
    with gzip.open(boundaries, 'rt') as handle_b:
        for line in handle_b:
            fields = line.split('\t')
            by_chrom.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    n = 0
    with open(out_domains, 'w') as out:
        for chrom, items in by_chrom.items():
            items.sort()
            for (_, end), (start, _) in zip(items[:-1], items[1:]):
                if start > end:
                    out.write('{}\t{}\t{}\n'.format(chrom, end, start))
                    n += 1
    print('{} domains'.format(n))


if __name__ == '__main__':
    main()
