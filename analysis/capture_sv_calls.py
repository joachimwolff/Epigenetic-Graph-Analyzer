"""Structural variant calls and copy number from the read pairs of capture Hi-C libraries (HiCUP BAM), per condition,
in the formats read by egaDifferential (--structuralVariants, --svContacts, --copyNumber).

Read pairs (first mate of every pair, primary alignments) are kept when their ends lie on different chromosomes or
at least --minIntra bp apart on one chromosome. Per condition (the libraries of a condition are pooled):

  joins          inter-chromosomal 1 Mb bin pairs are ranked by observed over expected pairs (expected: product of the
                 inter-chromosomal coverages of the two bins over twice the total, as capture_sv_cnv.py); intra-
                 chromosomal bin pairs at least --minIntra apart by observed over the mean of all bin pairs at the
                 same distance, scaled by the coverages of the two bins. A bin pair with at least --minPairs pairs and
                 an observed / expected of at least --minObsExp in a condition is a candidate; candidates of all
                 conditions within 2 Mb of each other on both sides form one join. Its breakend on each side is the
                 100 kb bin (within 2 Mb of the strongest 1 Mb bin pair) with the largest share of its
                 inter-chromosomal coverage in pairs with the 1 Mb around the other breakend (so that the enrichment
                 of baited bins cancels; conditions in which the join is present; one refinement). The join is present in a condition when its
                 strongest 1 Mb bin pair passes both thresholds there.
  contacts       for every baited fragment within --window of a breakend: the pairs to the 100 kb bins within
                 --window of the partner breakend, per condition (the promoter's own contacts across the junction).
  copy number    other ends of inter-chromosomal pairs with exactly one end on a baited fragment sample the genome
                 largely independently of the capture; their coverage per 1 Mb (10 Mb, chromosome) relative to the
                 autosomal median approximates the relative copy number (noisy: mappability and capture bias are not
                 corrected).

Outputs (<prefix>):
  _sv.bedpe                    one row per join: chrom1 start1 end1 chrom2 start2 end2 (the 100 kb breakpoint bins),
                               name (t(8;21) for inter-chromosomal joins, 'chr8 join 12.3-45.6 Mb' for intra-
                               chromosomal ones), score (largest observed / expected), strand1 strand2 ('.'), type
                               (inter / intra), conditions (comma-separated conditions where present); per condition
                               c: pairs_c (pairs in the 1 Mb bin pairs of the join), expected_c, obs_exp_c,
                               peak_pairs_c, peak_obs_exp_c (strongest 1 Mb bin pair), per_million_c (pairs per million
                               inter-chromosomal pairs of c), inter_total_c (inter-chromosomal pairs of c),
                               present_c (0/1)
  _sv_contacts.tsv             join, bait_chrom, bait_start, bait_end, partner_chrom, partner_start, partner_end,
                               pairs_c, per_million_c and inter_total_c per condition
  _<condition>_copy_number.bedGraph  relative copy number per 1 Mb bin (bins with at least 20 reads)
  _copy_number_10Mb.tsv, _copy_number_chromosomes.tsv  the same per 10 Mb and per chromosome, one column per condition
  _libraries.tsv               read pairs kept per library

Usage: python capture_sv_calls.py --bam siMM=a.bam,b.bam --bam siRE=c.bam --baits design.baitmap
           --chromSizes hg38.chrom.sizes -o prefix [--threads 12] [--samtools samtools] [--maxPairs N]
"""

import argparse
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

CHROMS = ['chr{}'.format(c) for c in list(range(1, 23)) + ['X', 'Y']]
BIN = 1000000
FINE = 100000
CNV_BIN = 10000000


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--bam', action='append', required=True,
                        help='condition=library.bam[,library2.bam], once per condition')
    parser.add_argument('--baits', required=True, help='Baited fragments (BED or CHiCAGO .baitmap)')
    parser.add_argument('--chromSizes', required=True)
    parser.add_argument('--outputPrefix', '-o', required=True)
    parser.add_argument('--minIntra', type=int, default=5000000,
                        help='Intra-chromosomal pairs at least this far apart are kept for long-range joins')
    parser.add_argument('--minPairs', type=int, default=50)
    parser.add_argument('--minObsExp', type=float, default=20.0)
    parser.add_argument('--window', type=int, default=2000000,
                        help='Baits within this distance of a breakend, partner bins within it of the partner breakend')
    parser.add_argument('--samtools', default='samtools')
    parser.add_argument('--threads', type=int, default=8, help='Threads in total (samtools decompression)')
    parser.add_argument('--maxPairs', type=int, default=0, help='Read at most this many kept pairs per library (test)')
    return parser.parse_args(args)


def read_baits(path):
    baits = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                        comment='#')
    baits = baits[baits['chrom'].isin(CHROMS)].sort_values(['chrom', 'start']).reset_index(drop=True)
    return baits


def bait_lookup(baits, chrom_index, positions):
    """Index of the bait fragment holding every (chrom, position), -1 elsewhere (vectorised per chromosome)."""
    out = np.full(len(positions), -1, dtype=np.int64)
    for c, group in baits.groupby('chrom', sort=False):
        if c not in CHROMS:
            continue
        mask = chrom_index == CHROMS.index(c)
        if not mask.any():
            continue
        starts = group['start'].to_numpy()
        ends = group['end'].to_numpy()
        ids = group.index.to_numpy()
        i = np.searchsorted(starts, positions[mask], side='right') - 1
        ok = (i >= 0) & (positions[mask] <= ends[np.clip(i, 0, None)])
        found = np.full(mask.sum(), -1, dtype=np.int64)
        found[ok] = ids[i[ok]]
        out[mask] = found
    return out


def read_library(task):
    """Kept pairs of one BAM: chrom index and position of both ends, ordered so that end 1 precedes end 2."""
    bam, samtools, threads, min_intra, limit = task
    expression = 'rnext != "=" || pnext - pos >= {0} || pos - pnext >= {0}'.format(min_intra)
    command = [samtools, 'view', '-@', str(max(1, threads - 1)), '-f', '0x40', '-F', '0x900', '-e', expression, bam]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1 << 22)
    code = {c: i for i, c in enumerate(CHROMS)}
    parts = []
    n = 0
    reader = pd.read_csv(process.stdout, sep='\t', header=None, usecols=[2, 3, 6, 7], names=['c1', 'p1', 'c2', 'p2'],
                         dtype={'c1': str, 'p1': np.int64, 'c2': str, 'p2': np.int64}, chunksize=2000000,
                         quoting=3, engine='c')
    for chunk in reader:
        c2 = chunk['c2'].where(chunk['c2'] != '=', chunk['c1'])
        a = chunk['c1'].map(code)
        b = c2.map(code)
        keep = a.notna() & b.notna()
        a = a[keep].to_numpy(dtype=np.int64)
        b = b[keep].to_numpy(dtype=np.int64)
        p = chunk.loc[keep, 'p1'].to_numpy()
        q = chunk.loc[keep, 'p2'].to_numpy()
        swap = (a > b) | ((a == b) & (p > q))
        parts.append(np.stack([np.where(swap, b, a), np.where(swap, q, p), np.where(swap, a, b),
                               np.where(swap, p, q)], axis=1))
        n += len(a)
        if limit and n >= limit:
            process.terminate()
            break
    process.wait()
    pairs = np.concatenate(parts) if parts else np.zeros((0, 4), dtype=np.int64)
    if limit:
        pairs = pairs[:limit]
    return bam, pairs


class Genome:
    def __init__(self, sizes, resolution):
        self.resolution = resolution
        self.sizes = np.array([sizes[c] for c in CHROMS], dtype=np.int64)
        self.n = (self.sizes + resolution - 1) // resolution
        self.offset = np.concatenate([[0], np.cumsum(self.n)[:-1]])
        self.total = int(self.n.sum())
        self.chrom_of = np.repeat(np.arange(len(CHROMS)), self.n)
        self.start_of = np.concatenate([np.arange(k) * resolution for k in self.n])

    def index(self, chrom, position):
        return self.offset[chrom] + np.minimum(position // self.resolution, self.n[chrom] - 1)


def bin_counts(i, j, size):
    keys = i * size + j
    unique, counts = np.unique(keys, return_counts=True)
    return unique // size, unique % size, counts


def condition_tables(pairs, coarse, min_intra):
    """1 Mb bin-pair counts with expected values: trans (coverage product) and intra >= min_intra (distance mean x
    coverage factors)."""
    i = coarse.index(pairs[:, 0], pairs[:, 1])
    j = coarse.index(pairs[:, 2], pairs[:, 3])
    trans = pairs[:, 0] != pairs[:, 2]
    coverage = np.bincount(i[trans], minlength=coarse.total) + np.bincount(j[trans], minlength=coarse.total)
    total = int(trans.sum())
    bi, bj, count = bin_counts(i[trans], j[trans], coarse.total)
    expected = coverage[bi] * coverage[bj] / (2.0 * max(total, 1))
    table_t = pd.DataFrame({'i': bi, 'j': bj, 'pairs': count, 'expected': expected, 'type': 'inter'})
    intra = ~trans
    ii, jj = i[intra], j[intra]
    min_d = int(np.ceil(min_intra / coarse.resolution))
    ok = (jj - ii) >= min_d
    bi, bj, count = bin_counts(ii[ok], jj[ok], coarse.total)
    distance = bj - bi
    # possible bin pairs at every distance over the genome
    max_d = int(coarse.n.max())
    possible = np.zeros(max_d + 1)
    for n in coarse.n:
        d = np.arange(min(n, max_d + 1))
        possible[d] += n - d
    observed_d = np.bincount(distance, weights=count, minlength=max_d + 1)
    mean_d = np.divide(observed_d, possible, out=np.zeros_like(observed_d), where=possible > 0)
    covered = coverage > 0
    factor = np.where(covered, coverage / coverage[covered].mean(), 0.0)
    expected = mean_d[distance] * factor[bi] * factor[bj]
    table_i = pd.DataFrame({'i': bi, 'j': bj, 'pairs': count, 'expected': expected, 'type': 'intra'})
    table = pd.concat([table_t, table_i], ignore_index=True)
    table['obs_exp'] = table['pairs'] / table['expected'].where(table['expected'] > 0)
    return table.set_index(['i', 'j']), total


def cluster(keys, gap=2):
    """Connected groups of bin pairs (i, j) within `gap` bins on both sides."""
    keys = sorted(keys)
    group = {}
    groups = []
    for k in keys:
        neighbours = {group[o] for o in group if abs(o[0] - k[0]) <= gap and abs(o[1] - k[1]) <= gap}
        if not neighbours:
            group[k] = len(groups)
            groups.append({k})
            continue
        target = min(neighbours)
        for other in sorted(neighbours - {target}):
            for member in groups[other]:
                group[member] = target
            groups[target] |= groups[other]
            groups[other] = set()
        group[k] = target
        groups[target].add(k)
    return [g for g in groups if g]


def chrom_label(c):
    return CHROMS[c][3:]


def main(args=None):
    args = parse_arguments(args)
    sizes = pd.read_csv(args.chromSizes, sep='\t', header=None, names=['chrom', 'size']).set_index('chrom')['size']
    coarse, fine = Genome(sizes, BIN), Genome(sizes, FINE)
    conditions = {}
    for item in args.bam:
        name, paths = item.split('=', 1)
        conditions[name] = paths.split(',')
    libraries = [p for paths in conditions.values() for p in paths]
    per = max(1, args.threads // max(1, len(libraries)))
    tasks = [(p, args.samtools, per, args.minIntra, args.maxPairs) for p in libraries]
    with Pool(min(len(tasks), args.threads)) as pool:
        read = dict(pool.map(read_library, tasks))
    baits = read_baits(args.baits)
    rows = []
    pooled = {}
    for name, paths in conditions.items():
        parts = []
        for p in paths:
            rows.append({'condition': name, 'library': p, 'pairs_kept': len(read[p]),
                         'inter_chromosomal': int((read[p][:, 0] != read[p][:, 2]).sum())})
            parts.append(read[p])
        pooled[name] = np.concatenate(parts)
    prefix = args.outputPrefix
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    pd.DataFrame(rows).to_csv(prefix + '_libraries.tsv', sep='\t', index=False)
    tables, totals = {}, {}
    for name, pairs in pooled.items():
        tables[name], totals[name] = condition_tables(pairs, coarse, args.minIntra)
        print('{}: {:,} pairs kept, {:,} inter-chromosomal'.format(name, len(pairs), totals[name]), flush=True)
    # candidates of every condition, clustered per chromosome pair
    candidates = set()
    for name, table in tables.items():
        hit = table[(table['pairs'] >= args.minPairs) & (table['obs_exp'] >= args.minObsExp)]
        candidates |= set(hit.index)
    by_pair = {}
    for i, j in candidates:
        by_pair.setdefault((coarse.chrom_of[i], coarse.chrom_of[j]), []).append((i, j))
    # fine bin pairs per condition for the breakpoint
    fine_index = {name: (fine.index(p[:, 0], p[:, 1]), fine.index(p[:, 2], p[:, 3])) for name, p in pooled.items()}
    coarse_index = {name: coarse.index(p[:, 0], p[:, 1]) * coarse.total + coarse.index(p[:, 2], p[:, 3])
                    for name, p in pooled.items()}
    # inter-chromosomal coverage of every 100 kb bin per condition
    fine_coverage = {}
    for name, p in pooled.items():
        trans = p[:, 0] != p[:, 2]
        fi, fj = fine_index[name]
        fine_coverage[name] = np.bincount(fi[trans], minlength=fine.total) + \
            np.bincount(fj[trans], minlength=fine.total)
    # pairs grouped by chromosome pair (ends are ordered, chromosome of end 1 <= chromosome of end 2)
    groups = {}
    for name, pairs in pooled.items():
        key = pairs[:, 0] * len(CHROMS) + pairs[:, 2]
        order = np.argsort(key, kind='stable')
        bounds = np.searchsorted(key[order], np.arange(len(CHROMS) ** 2 + 1))
        groups[name] = (order, bounds)
    joins = []
    for (c1, c2), keys in sorted(by_pair.items()):
        for group in cluster(keys):
            group = sorted(group)
            record = {}
            peak_score = {}
            for name, table in tables.items():
                part = table.reindex(group)
                pairs = float(part['pairs'].fillna(0).sum())
                expected = float(part['expected'].fillna(0).sum())
                best = part['obs_exp'].fillna(0).to_numpy()
                k = int(np.argmax(part['pairs'].fillna(0).to_numpy()))
                peak = part.iloc[k]
                record[name] = dict(pairs=pairs, expected=expected,
                                    obs_exp=pairs / expected if expected > 0 else np.nan,
                                    peak_pairs=float(np.nan_to_num(peak['pairs'])),
                                    peak_obs_exp=float(np.nan_to_num(peak['obs_exp'])),
                                    per_million=pairs / max(totals[name], 1) * 1e6, inter_total=totals[name],
                                    present=int(bool(((part['pairs'].fillna(0) >= args.minPairs) &
                                                      (part['obs_exp'].fillna(0) >= args.minObsExp)).any())))
                peak_score[name] = float(best.max()) if len(best) else 0.0
            # breakpoint: on each side the 100 kb bin with the largest share of its inter-chromosomal coverage in
            # pairs with the other side (the capture enrichment of baited bins cancels; bins with at least 5 such
            # pairs), in the conditions where the join is present (all when present in none). The other side is
            # first the window of 1 Mb around the strongest 1 Mb bin pair of the join (largest observed / expected
            # among those with at least --minPairs pairs), then 1 Mb around the breakend found there (one refinement)
            used = [n for n in pooled if record[n]['present']] or list(pooled)
            best, best_key = -1.0, group[0]
            for name in used:
                part = tables[name].reindex(group)
                ok = part['pairs'].fillna(0) >= args.minPairs
                if ok.any():
                    k = int(np.argmax(np.where(ok, part['obs_exp'].fillna(0), -1)))
                    if part['obs_exp'].iloc[k] > best:
                        best, best_key = float(part['obs_exp'].iloc[k]), group[k]
            positions = {}
            for name in used:
                pairs = pooled[name]
                order, bounds = groups[name]
                k = c1 * len(CHROMS) + c2
                sub = pairs[order[bounds[k]:bounds[k + 1]]]
                positions[name] = np.ones(len(sub), dtype=bool), sub[:, 1], sub[:, 3]
            coverage = sum(fine_coverage[n] for n in used)
            centre = [coarse.start_of[best_key[0]] + BIN // 2, coarse.start_of[best_key[1]] + BIN // 2]
            for _ in range(2):
                ends = []
                for side in (0, 1):
                    counts = np.zeros(fine.total)
                    for name in used:
                        on_pair, pos1, pos2 = positions[name]
                        own, other = (pos1, pos2) if side == 0 else (pos2, pos1)
                        chrom = c1 if side == 0 else c2
                        sel = on_pair & (np.abs(other - centre[1 - side]) <= BIN) & \
                            (np.abs(own - centre[side]) <= 2 * BIN)
                        counts += np.bincount(fine.index(np.full(int(sel.sum()), chrom), own[sel]),
                                              minlength=fine.total)
                    ratio = np.where((counts >= 5) & (coverage > 0), counts / np.maximum(coverage, 1), -1.0)
                    ends.append(int(np.argmax(ratio)) if ratio.max() > 0 else int(np.argmax(counts)))
                centre = [int(fine.start_of[ends[0]]) + FINE // 2, int(fine.start_of[ends[1]]) + FINE // 2]
            bx, by = ends
            start1, start2 = int(fine.start_of[bx]), int(fine.start_of[by])
            kind = 'inter' if c1 != c2 else 'intra'
            if kind == 'inter':
                name = 't({};{})'.format(chrom_label(c1), chrom_label(c2))
            else:
                name = '{} join {:.1f}-{:.1f} Mb'.format(CHROMS[c1], start1 / 1e6, start2 / 1e6)
            row = dict(chrom1=CHROMS[c1], start1=start1, end1=start1 + FINE, chrom2=CHROMS[c2], start2=start2,
                       end2=start2 + FINE, name=name, score=round(max(peak_score.values()), 1), strand1='.',
                       strand2='.', type=kind, conditions=','.join(n for n in record if record[n]['present']))
            for cond, values in record.items():
                for key, value in values.items():
                    row['{}_{}'.format(key, cond)] = value
            joins.append(row)
    calls = pd.DataFrame(joins)
    if len(calls):
        calls = calls.sort_values('score', ascending=False).reset_index(drop=True)
        # several joins between the same chromosome pair: numbered
        counts = calls.groupby('name').cumcount()
        calls['name'] = [n if k == 0 else '{} #{}'.format(n, k + 1) for n, k in zip(calls['name'], counts)]
    columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2', 'type',
               'conditions']
    with open(prefix + '_sv.bedpe', 'w') as handle:
        handle.write('#' + '\t'.join(list(calls.columns) if len(calls) else columns) + '\n')
        if len(calls):
            calls.to_csv(handle, sep='\t', index=False, header=False, float_format='%.6g')
    print(calls.head(30).to_string() if len(calls) else 'no joins', flush=True)
    # contacts of the baits near every breakend with the partner side
    contact_rows = []
    window = args.window
    for _, call in calls.iterrows():
        sides = [(call['chrom1'], call['start1'] + FINE // 2, call['chrom2'], call['start2'] + FINE // 2),
                 (call['chrom2'], call['start2'] + FINE // 2, call['chrom1'], call['start1'] + FINE // 2)]
        found = {}
        a_index, b_index = sorted((CHROMS.index(call['chrom1']), CHROMS.index(call['chrom2'])))
        for name, all_pairs in pooled.items():
            order, bounds = groups[name]
            k = a_index * len(CHROMS) + b_index
            pairs = all_pairs[order[bounds[k]:bounds[k + 1]]]
            for end in (0, 1):
                ca, pa = pairs[:, 2 * end], pairs[:, 2 * end + 1]
                cb, pb = pairs[:, 2 - 2 * end], pairs[:, 3 - 2 * end]
                for bait_chrom, bait_pos, partner_chrom, partner_pos in sides:
                    sel = (ca == CHROMS.index(bait_chrom)) & (np.abs(pa - bait_pos) <= window) & \
                          (cb == CHROMS.index(partner_chrom)) & (np.abs(pb - partner_pos) <= window)
                    if call['type'] == 'intra':
                        # both breakends on one chromosome: the two windows must be different
                        sel &= np.abs(pa - pb) >= args.minIntra
                    if not sel.any():
                        continue
                    ids = bait_lookup(baits, ca[sel], pa[sel])
                    partner_bin = pb[sel] // FINE * FINE
                    keep = ids >= 0
                    for bait, start in zip(ids[keep], partner_bin[keep]):
                        key = (int(bait), partner_chrom, int(start))
                        found.setdefault(key, {n: 0 for n in pooled})[name] += 1
        for (bait, partner_chrom, start), counts in found.items():
            row = dict(join=call['name'], bait_chrom=baits.loc[bait, 'chrom'], bait_start=int(baits.loc[bait, 'start']),
                       bait_end=int(baits.loc[bait, 'end']), partner_chrom=partner_chrom, partner_start=start,
                       partner_end=start + FINE)
            for name in pooled:
                row['pairs_' + name] = counts[name]
                row['per_million_' + name] = counts[name] / max(totals[name], 1) * 1e6
                row['inter_total_' + name] = totals[name]
            contact_rows.append(row)
    contacts = pd.DataFrame(contact_rows)
    if len(contacts):
        contacts = contacts.sort_values(['join', 'bait_chrom', 'bait_start', 'partner_chrom', 'partner_start'])
    contacts.to_csv(prefix + '_sv_contacts.tsv', sep='\t', index=False, float_format='%.6g')
    # copy number from the other ends of captured inter-chromosomal pairs
    per_10mb, per_chrom = {}, {}
    for name, pairs in pooled.items():
        trans = pairs[:, 0] != pairs[:, 2]
        t = pairs[trans]
        b1 = bait_lookup(baits, t[:, 0], t[:, 1]) >= 0
        b2 = bait_lookup(baits, t[:, 2], t[:, 3]) >= 0
        other_c = np.concatenate([t[b1 & ~b2, 2], t[b2 & ~b1, 0]])
        other_p = np.concatenate([t[b1 & ~b2, 3], t[b2 & ~b1, 1]])
        counts = np.bincount(coarse.index(other_c, other_p), minlength=coarse.total).astype(float)
        autosomal = coarse.chrom_of < 22
        median = np.median(counts[autosomal & (counts >= 20)])
        relative = counts / median
        keep = counts >= 20
        bedgraph = pd.DataFrame({'chrom': [CHROMS[c] for c in coarse.chrom_of[keep]],
                                 'start': coarse.start_of[keep],
                                 'end': np.minimum(coarse.start_of[keep] + BIN, coarse.sizes[coarse.chrom_of[keep]]),
                                 'relative': relative[keep]})
        bedgraph.to_csv('{}_{}_copy_number.bedGraph'.format(prefix, name), sep='\t', index=False, header=False,
                        float_format='%.4f')
        frame = pd.DataFrame({'chrom': coarse.chrom_of, 'start': coarse.start_of // CNV_BIN * CNV_BIN,
                              'reads': counts})
        ten = frame.groupby(['chrom', 'start'])['reads'].sum()
        ten = ten / ten[ten.index.get_level_values(0) < 22].median()
        per_10mb[name] = ten
        chrom_mean = frame[frame['reads'] > 0].groupby('chrom')['reads'].mean()
        per_chrom[name] = chrom_mean / chrom_mean[chrom_mean.index < 22].median()
    ten = pd.DataFrame(per_10mb).reset_index()
    ten['chrom'] = ten['chrom'].map(lambda c: CHROMS[c])
    ten.to_csv(prefix + '_copy_number_10Mb.tsv', sep='\t', index=False, float_format='%.4f')
    chrom = pd.DataFrame(per_chrom).reindex(range(len(CHROMS)))
    chrom.index = CHROMS
    chrom.index.name = 'chrom'
    chrom.to_csv(prefix + '_copy_number_chromosomes.tsv', sep='\t', float_format='%.4f')
    print(chrom.round(2).to_string(), flush=True)


if __name__ == '__main__':
    sys.exit(main())
