"""Promoter-anchored loops of a Hi-C loop caller in the ibed format of the virtual capture calls.

Loops of Mustache (BIN1_CHR BIN1_START ... FDR) or Chromosight (chrom1 start1 ... qvalue) are kept if at least one
anchor lies in a baited bin of the virtual capture design (.baitmap, bins matched by their midpoint); that anchor
becomes the bait, and loops with two baited anchors are written in both orientations. The score column is
-log10 of the loop q-value (FDR), capped at 50, so that the benchmark can apply the same threshold to all loops.

Usage: python loops_to_ibed.py <loops.tsv> <mustache|chromosight> <design.baitmap> <bin size> <output.ibed>
"""

import sys

import numpy as np
import pandas as pd


def main():
    loops_path, tool, baitmap, bin_size, output = sys.argv[1:6]
    bin_size = int(bin_size)
    loops = pd.read_csv(loops_path, sep='\t')
    if tool == 'mustache':
        loops = loops.rename(columns={'BIN1_CHR': 'chrom1', 'BIN1_START': 'start1', 'BIN1_END': 'end1',
                                      'BIN2_CHROMOSOME': 'chrom2', 'BIN2_START': 'start2', 'BIN2_END': 'end2',
                                      'FDR': 'qvalue'})
    loops = loops[loops['chrom1'] == loops['chrom2']]
    baits = pd.read_csv(baitmap, sep='\t', header=None, usecols=[0, 1, 2, 3, 4],
                        names=['chrom', 'start', 'end', 'id', 'name'])
    key = {(c, (s + e) // 2 // bin_size): (s, e, str(n)) for c, s, e, n in
           zip(baits['chrom'], baits['start'], baits['end'], baits['name'])}
    score = np.minimum(-np.log10(np.clip(loops['qvalue'].astype(float).to_numpy(), 1e-50, 1)), 50)
    rows = []
    for (chrom, s1, e1, s2, e2), value in zip(loops[['chrom1', 'start1', 'end1', 'start2', 'end2']].itertuples(
            index=False, name=None), score):
        ends = [(s1, e1), (s2, e2)]
        for k in (0, 1):
            bait = key.get((chrom, (ends[k][0] + ends[k][1]) // 2 // bin_size))
            if bait is None:
                continue
            other = ends[1 - k]
            rows.append((chrom, bait[0], bait[1], bait[2], chrom, other[0] + 1, other[1], '.', 0, round(value, 3)))
    table = pd.DataFrame(rows, columns=['bait_chr', 'bait_start', 'bait_end', 'bait_name', 'otherEnd_chr',
                                        'otherEnd_start', 'otherEnd_end', 'otherEnd_name', 'N_reads', 'score'])
    table.to_csv(output, sep='\t', index=False)
    print('{} of {} loops with a baited anchor, {} rows'.format(
        int(len(set((r[0], r[1], r[5]) for r in rows))), len(loops), len(table)))


if __name__ == '__main__':
    main()
