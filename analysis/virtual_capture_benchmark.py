"""Agreement of virtual capture (CHiCAGO on bins of a Hi-C/Micro-C matrix) with real promoter capture Hi-C calls.

Both call sets are mapped to the virtual-capture bins (real calls lifted to the bin assembly when a chain file is
given; each end is assigned to the bin holding its midpoint). Only pairs whose bait bin is baited in both designs
are compared (a real bait counts when its midpoint falls into a virtual bait bin). Recall is the fraction of real
pairs found among the virtual calls, precision the fraction of virtual calls found among the real pairs, exactly
and with a tolerance of one bin at the other end, per distance stratum. Cis pairs within 2 Mb only.

Usage: python virtual_capture_benchmark.py <real calls> <real format> <real baits or -> <chain or -> <virtual calls>
       <virtual.baitmap> <bin size> <output prefix> [real score threshold=5]
       [<virtual format> <virtual baits or -> <virtual chain or -> <virtual score threshold>]
The optional arguments read the second call set like the first (default: CHiCAGO ibed on the bins, no liftover);
with two real call sets the result is the agreement between two capture experiments on the same bins.
"""

import sys

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.genome import LiftOver
from epigeneticgraphanalyzer.io import read_baits, read_ibed, read_interactions

STRATA = (0, 50000, 100000, 250000, 500000, 1000000, 2000001)


def lift(chroms, positions, liftover):
    out_chrom, out_pos = [], []
    for c, p in zip(chroms, positions):
        hit = liftover.position(c, int(p)) if liftover else (c, int(p))
        if hit is None or hit[0] != c:
            out_chrom.append(None)
            out_pos.append(-1)
        else:
            out_chrom.append(hit[0])
            out_pos.append(hit[1])
    return np.array(out_chrom, dtype=object), np.array(out_pos)


def keys(chrom, bait_bin, other_bin):
    return pd.Series(chrom).astype(str) + ':' + pd.Series(bait_bin).astype(str) + ':' + pd.Series(other_bin).astype(str)


def binned(path, file_format, bait_path, chain, threshold, bin_size):
    """Calls as (chrom, bait bin, other bin), lifted to the bin assembly when a chain file is given."""
    liftover = LiftOver(chain) if chain != '-' else None
    baits = read_baits(bait_path) if bait_path != '-' else None
    calls = read_ibed(path) if file_format == 'ibed' and baits is None else \
        read_interactions([path], file_format, baits=baits)
    calls = calls[(calls['score'] >= threshold) & (calls['distance'] <= 2000000)]
    chrom, bait_pos = lift(calls['chrom'], (calls['bait_start'] + calls['bait_end']) // 2, liftover)
    chrom2, other_pos = lift(calls['chrom'], (calls['other_start'] + calls['other_end']) // 2, liftover)
    ok = (chrom != None) & (chrom2 != None) & (chrom == chrom2)  # noqa: E711
    return pd.DataFrame({'chrom': chrom[ok], 'bait_bin': bait_pos[ok] // bin_size,
                         'other_bin': other_pos[ok] // bin_size})


def main():
    real_path, real_format, real_baits, chain, vc_path, vc_baitmap, bin_size, prefix = sys.argv[1:9]
    threshold = float(sys.argv[9]) if len(sys.argv) > 9 else 5.0
    vc_format, vc_bait_path, vc_chain = (sys.argv[10:13] if len(sys.argv) > 12 else ('ibed', '-', '-'))
    vc_threshold = float(sys.argv[13]) if len(sys.argv) > 13 else 5.0
    bin_size = int(bin_size)
    real = binned(real_path, real_format, real_baits, chain, threshold, bin_size)
    vc = binned(vc_path, vc_format, vc_bait_path, vc_chain, vc_threshold, bin_size)
    vc_baits = read_baits(vc_baitmap)
    vc_bait_bins = set(vc_baits['chrom'] + ':' + (vc_baits['start'] // bin_size).astype(str))
    real_bait_bins = set(real['chrom'] + ':' + real['bait_bin'].astype(str))
    shared = vc_bait_bins & real_bait_bins
    for table in (real, vc):
        table['bait_key'] = table['chrom'] + ':' + table['bait_bin'].astype(str)
        table['distance'] = (table['other_bin'] - table['bait_bin']).abs() * bin_size
    real = real[real['bait_key'].isin(shared) & (real['bait_bin'] != real['other_bin'])].drop_duplicates(
        ['chrom', 'bait_bin', 'other_bin'])
    vc = vc[vc['bait_key'].isin(shared)].drop_duplicates(['chrom', 'bait_bin', 'other_bin'])
    real_keys = set(keys(real['chrom'], real['bait_bin'], real['other_bin']))
    vc_keys = set(keys(vc['chrom'], vc['bait_bin'], vc['other_bin']))

    def found(table, reference, tolerance):
        hit = np.zeros(len(table), dtype=bool)
        for shift in range(-tolerance, tolerance + 1):
            hit |= keys(table['chrom'], table['bait_bin'], table['other_bin'] + shift).isin(reference).to_numpy()
        return hit

    rows = []
    stratum_real = np.searchsorted(STRATA, real['distance'].to_numpy(), side='right') - 1
    stratum_vc = np.searchsorted(STRATA, vc['distance'].to_numpy(), side='right') - 1
    for tolerance in (0, 1):
        recall_hit = found(real, vc_keys, tolerance)
        precision_hit = found(vc, real_keys, tolerance)
        for s in range(-1, len(STRATA) - 1):
            mr = np.ones(len(real), dtype=bool) if s < 0 else stratum_real == s
            mv = np.ones(len(vc), dtype=bool) if s < 0 else stratum_vc == s
            label = 'all' if s < 0 else '{}-{}kb'.format(STRATA[s] // 1000, STRATA[s + 1] // 1000)
            rows.append(dict(tolerance_bins=tolerance, stratum=label, real_pairs=int(mr.sum()),
                             virtual_pairs=int(mv.sum()),
                             recall=float(recall_hit[mr].mean()) if mr.any() else np.nan,
                             precision=float(precision_hit[mv].mean()) if mv.any() else np.nan))
    result = pd.DataFrame(rows)
    result.insert(0, 'shared_bait_bins', len(shared))
    result.to_csv(prefix + '_benchmark.tsv', sep='\t', index=False)
    print(result.to_string(index=False))


if __name__ == '__main__':
    main()
