"""Divide a coverage bigWig by its genome-wide mean (sumData / genome length of the bigWig header), so that tracks of
different sequencing depth share one scale (1 = genome average).

Usage: python scale_bigwig.py <in.bw> <out.bw>

The ChIP-seq bigWigs of Rubin et al. 2017 (GSE84657) are bedtools coverage without depth normalisation; their genome
means differ up to fourfold between donors and days, which the mark changes of egaDifferential (log2 of the node
means, pseudo count 1) would read as chromatin changes. Intervals are copied unchanged, values scaled.
"""

import sys

import numpy as np
import pyBigWig


def main():
    source, target = sys.argv[1], sys.argv[2]
    handle = pyBigWig.open(source)
    chroms = handle.chroms()
    header = handle.header()
    factor = sum(chroms.values()) / header['sumData']
    out = pyBigWig.open(target, 'w')
    out.addHeader(list(chroms.items()), maxZooms=10)
    for chrom in chroms:
        intervals = handle.intervals(chrom)
        if not intervals:
            continue
        array = np.array(intervals, dtype=np.float64)
        starts = array[:, 0].astype(np.int64)
        ends = array[:, 1].astype(np.int64)
        values = array[:, 2] * factor
        step = 1_000_000
        for i in range(0, len(starts), step):
            out.addEntries([chrom] * len(starts[i:i + step]), starts[i:i + step].tolist(),
                           ends=ends[i:i + step].tolist(), values=values[i:i + step].tolist())
    out.close()
    handle.close()
    print('{}\tfactor {:.4f}'.format(source, factor))


if __name__ == '__main__':
    main()
