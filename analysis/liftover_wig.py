"""Lift a variableStep WIG (e.g. GEO DNase-seq tracks) to another assembly and write a sorted bedGraph.

Every step is lifted by its start with the UCSC chain and written immediately (constant memory); the lifted steps
are sorted with GNU sort, and steps that do not map, map to another chromosome, run past the chromosome end or
overlap the previous step after sorting are dropped.

Usage: python liftover_wig.py <in.wig.gz> <chain.gz> <chrom.sizes of the target> <out.bedGraph> [sort memory, 4G]
"""

import gzip
import os
import subprocess
import sys

from pyliftover import LiftOver


def main():
    wig, chain, sizes_path, output = sys.argv[1:5]
    memory = sys.argv[5] if len(sys.argv) > 5 else '4G'
    lift = LiftOver(chain)
    sizes = {}
    with open(sizes_path) as handle:
        for line in handle:
            name, size = line.split()[:2]
            sizes[name] = int(size)
    unsorted = output + '.unsorted'
    written = 0
    chrom, span = None, 1
    with gzip.open(wig, 'rt') as handle, open(unsorted, 'w') as out:
        for line in handle:
            if line.startswith('track'):
                continue
            if line.startswith('variableStep'):
                fields = dict(item.split('=') for item in line.split()[1:])
                chrom, span = fields['chrom'], int(fields.get('span', 1))
                continue
            position, value = line.split()
            hit = lift.convert_coordinate(chrom, int(position) - 1)
            if not hit or hit[0][0] != chrom:
                continue
            start = hit[0][1]
            if start + span > sizes.get(chrom, 0):
                continue
            out.write('{}\t{}\t{}\t{}\n'.format(chrom, start, start + span, value))
            written += 1
    # sort, then drop steps overlapping the previous step of the same chromosome
    command = ("LC_ALL=C sort -k1,1 -k2,2n -S {} -T {} {} | "
               "awk -F'\\t' -v OFS='\\t' '$1 != c {{c = $1; e = -1}} $2 >= e {{print; e = $3}}' > {}").format(
        memory, os.path.dirname(os.path.abspath(output)), unsorted, output)
    subprocess.run(command, shell=True, check=True)
    os.remove(unsorted)
    print('{} steps lifted'.format(written))


if __name__ == '__main__':
    main()
