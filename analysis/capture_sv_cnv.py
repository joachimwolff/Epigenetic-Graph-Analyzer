"""Structural variants and copy number from the read pairs of a capture Hi-C library (HiCUP BAM).

Translocations: inter-chromosomal read pairs are counted per pair of 1 Mb bins; each bin pair is ranked by observed
over expected pairs (product of the inter-chromosomal coverages of the two bins over the total), so that a
rearrangement joining two loci appears as a bin pair with many more pairs than its coverage predicts.
Copy number: for inter-chromosomal pairs with one end on a baited fragment, the other end samples the genome
largely independently of the capture; its coverage per chromosome and per 10 Mb, relative to the autosomal median,
estimates the relative copy number (one copy of X and none of Y in a male line with loss of Y).

Usage: python capture_sv_cnv.py <hicup.bam> <baitmap> <output prefix> [samtools] [max read pairs]
"""

import subprocess
import sys
from collections import Counter

import numpy as np
import pandas as pd

BIN = 1000000
CNV_BIN = 10000000
CHROMS = ['chr{}'.format(c) for c in list(range(1, 23)) + ['X', 'Y']]


def bait_index(path):
    baits = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    index = {}
    for chrom, group in baits.groupby('chrom'):
        group = group.sort_values('start')
        index[chrom] = (group['start'].to_numpy(), group['end'].to_numpy())
    return index


def is_bait(index, chrom, position):
    if chrom not in index:
        return False
    starts, ends = index[chrom]
    i = np.searchsorted(starts, position, side='right') - 1
    return i >= 0 and position <= ends[i]


def main():
    bam, baitmap, prefix = sys.argv[1:4]
    samtools = sys.argv[4] if len(sys.argv) > 4 else 'samtools'
    limit = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    baits = bait_index(baitmap)
    trans_pairs = Counter()
    trans_cov = Counter()
    other_end = Counter()
    n = 0
    process = subprocess.Popen([samtools, 'view', '-F', '0x900', bam], stdout=subprocess.PIPE, text=True,
                               bufsize=1 << 20)
    previous = None
    for line in process.stdout:
        fields = line.split('\t', 4)
        record = (fields[0], fields[2], int(fields[3]))
        if previous is None or previous[0] != record[0]:
            previous = record
            continue
        (_, c1, p1), (_, c2, p2) = previous, record
        previous = None
        if c1 not in CHROMS or c2 not in CHROMS or c1 == c2:
            continue
        b1, b2 = (c1, p1 // BIN), (c2, p2 // BIN)
        key = (b1, b2) if (CHROMS.index(c1), p1) < (CHROMS.index(c2), p2) else (b2, b1)
        trans_pairs[key] += 1
        trans_cov[b1] += 1
        trans_cov[b2] += 1
        bait1, bait2 = is_bait(baits, c1, p1), is_bait(baits, c2, p2)
        if bait1 and not bait2:
            other_end[(c2, p2 // CNV_BIN)] += 1
        elif bait2 and not bait1:
            other_end[(c1, p1 // CNV_BIN)] += 1
        n += 1
        if limit and n >= limit:
            process.terminate()
            break
    total = sum(trans_pairs.values())
    rows = []
    for ((c1, x1), (c2, x2)), observed in trans_pairs.items():
        if observed < 5:
            continue
        expected = trans_cov[(c1, x1)] * trans_cov[(c2, x2)] / (2.0 * total)
        rows.append((c1, x1 * BIN, c2, x2 * BIN, observed, expected, observed / expected))
    table = pd.DataFrame(rows, columns=['chrom1', 'start1', 'chrom2', 'start2', 'pairs', 'expected', 'obs_exp'])
    table = table.sort_values('obs_exp', ascending=False).reset_index(drop=True)
    table.to_csv(prefix + '_trans_hotspots.tsv', sep='\t', index=False)
    cnv = pd.DataFrame([(c, b * CNV_BIN, v) for (c, b), v in other_end.items()], columns=['chrom', 'start', 'reads'])
    autosomal = cnv[~cnv['chrom'].isin(['chrX', 'chrY'])]['reads'].median()
    cnv['relative'] = cnv['reads'] / autosomal
    cnv = cnv.sort_values(['chrom', 'start'], key=lambda c: c.map(CHROMS.index) if c.name == 'chrom' else c)
    cnv.to_csv(prefix + '_copy_number_10Mb.tsv', sep='\t', index=False)
    per_chrom = cnv.groupby('chrom')['reads'].sum() / cnv.groupby('chrom')['reads'].count()
    per_chrom = (per_chrom / per_chrom.drop(['chrX', 'chrY'], errors='ignore').median()).reindex(CHROMS)
    per_chrom.rename('relative_coverage').to_csv(prefix + '_copy_number_chromosomes.tsv', sep='\t')
    print('{} inter-chromosomal pairs'.format(total))
    print(table.head(15).to_string())
    print(per_chrom.round(2).to_string())


if __name__ == '__main__':
    main()
