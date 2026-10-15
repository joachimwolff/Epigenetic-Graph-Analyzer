"""Lift a variableStep WIG (e.g. GEO DNase-seq tracks) to another assembly and write a sorted bedGraph.

Every step is lifted by its start with the UCSC chain; steps that do not map, map to another chromosome or would
overlap the previous lifted step of the same chromosome after sorting are dropped.

Usage: python liftover_wig.py <in.wig.gz> <chain.gz> <chrom.sizes of the target> <out.bedGraph>
"""

import gzip
import sys

import numpy as np
import pandas as pd
from pyliftover import LiftOver


def main():
    wig, chain, sizes_path, output = sys.argv[1:5]
    lift = LiftOver(chain)
    sizes = dict(line.split()[:2] for line in open(sizes_path))
    chroms, starts, ends, values = [], [], [], []
    chrom, span = None, 1
    with gzip.open(wig, 'rt') as handle:
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
            if start + span > int(sizes.get(chrom, 0)):
                continue
            chroms.append(chrom)
            starts.append(start)
            ends.append(start + span)
            values.append(float(value))
    table = pd.DataFrame({'chrom': chroms, 'start': starts, 'end': ends, 'value': values})
    table = table.sort_values(['chrom', 'start'])
    previous_end = table.groupby('chrom')['end'].shift(1).fillna(-1).to_numpy()
    table = table[table['start'].to_numpy() >= previous_end]
    table.to_csv(output, sep='\t', header=False, index=False)
    print('{} steps lifted'.format(len(table)))


if __name__ == '__main__':
    main()
