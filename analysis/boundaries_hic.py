"""Hi-C insulation reference tracks for the boundary analysis (27_boundaries.sh); needs cooler (HiCExplorer env).

Diamond insulation (epigeneticgraphanalyzer.boundaries.diamond_insulation; square of window x window pixels
between the bins upstream and downstream of each bin, log2 over the chromosome mean, NaN with < 50 % valid pixels):
  4DN mcool files (balanced, GRCh38) of GM12878, K562 (in situ Hi-C) and H1 (Micro-C) at 10, 25 and 50 kb with a
  window of 200 kb (20 / 8 bins) at 10 / 25 kb and 250 kb (5 bins) at 50 kb.
Example regions (--regions chrom:start-end ...): balanced matrices are written as .npz for the figure.

Usage: python boundaries_hic.py <mcool dir> <out dir> [threads] [--regions cell:chrom:start-end ...]
Output: <out>/hic_insulation_<cell>_<res>kb.tsv.gz (chrom, start, end, insulation), <out>/hic_region_*.npz
"""

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer.boundaries import band_from_pixels, diamond_insulation  # noqa: E402

MCOOL = {'GM12878': '4DNFIIV4M7TF', 'K562': '4DNFI18UHVRO', 'H1': '4DNFI9GMP2J8'}
CHROMS = ['chr{}'.format(c) for c in list(range(1, 23)) + ['X']]
WINDOW_BP = {10000: 200000, 25000: 200000, 50000: 250000}


def mcool_insulation(args):
    import cooler
    path, cell, res, out = args
    target = os.path.join(out, 'hic_insulation_{}_{}kb.tsv.gz'.format(cell, res // 1000))
    if os.path.exists(target):
        return target
    clr = cooler.Cooler('{}::/resolutions/{}'.format(path, res))
    w = WINDOW_BP[res] // res
    rows = []
    for chrom in CHROMS:
        if chrom not in clr.chromnames:
            continue
        bins = clr.bins().fetch(chrom)
        n = len(bins)
        valid = np.isfinite(bins['weight'].values)
        pix = clr.matrix(balance=True, sparse=True, as_pixels=True).fetch(chrom)
        lo = clr.offset(chrom)
        i = pix.bin1_id.values - lo
        j = pix.bin2_id.values - lo
        v = np.nan_to_num(pix.balanced.values)
        band = band_from_pixels(n, i, j, v, 2 * w, valid=valid)
        ins = diamond_insulation(band, w)
        rows.append(pd.DataFrame({'chrom': chrom, 'start': bins.start.values, 'end': bins.end.values,
                                  'insulation': ins}))
        print(cell, res, chrom, n, int(np.isfinite(ins).sum()), flush=True)
    pd.concat(rows).to_csv(target, sep='\t', index=False, float_format='%.5g')
    return target


def region_matrix(path, cell, res, region, out):
    import cooler
    chrom, span = region.split(':')
    start, end = (int(x) for x in span.split('-'))
    clr = cooler.Cooler('{}::/resolutions/{}'.format(path, res))
    m = clr.matrix(balance=True).fetch((chrom, start, end))
    target = os.path.join(out, 'hic_region_{}_{}_{}_{}_{}kb.npz'.format(cell, chrom, start, end, res // 1000))
    np.savez_compressed(target, matrix=m, chrom=chrom, start=start, end=end, res=res)
    return target


def main():
    argv = sys.argv[1:]
    regions = []
    if '--regions' in argv:
        k = argv.index('--regions')
        regions = argv[k + 1:]
        argv = argv[:k]
    mcool_dir, out = argv[:2]
    threads = int(argv[2]) if len(argv) > 2 else 8
    os.makedirs(out, exist_ok=True)
    if regions:
        for spec in regions:
            cell, res, region = spec.split(':', 2)
            print(region_matrix(os.path.join(mcool_dir, MCOOL[cell] + '.mcool'), cell, int(res), region, out))
        return
    jobs = [(os.path.join(mcool_dir, acc + '.mcool'), cell, res, out) for cell, acc in MCOOL.items()
            for res in (10000, 25000, 50000)]
    with Pool(min(threads, len(jobs))) as pool:
        print(pool.map(mcool_insulation, jobs))


if __name__ == '__main__':
    main()
