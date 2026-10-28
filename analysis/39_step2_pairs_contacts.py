"""Contacts of the explaining step-2 pairs (promoter - element) in pairwise contact data of one condition.

Promoter capture Hi-C never measures the promoter and the element of a step-2 path as a pair. This script takes the
explaining and background triplets of 36_gm12878_multiway_paths.py (genes without an explanation at the promoter or
step 1, explained at step 2) and counts, for every pair promoter - element, the contacts of the pair in contact maps
(cooler, unbalanced counts summed over the bins of the two fragments) and in pair lists with counts (ChIA-PET PET
clusters: chrom1 start1 end1 chrom2 start2 end2 count; the counts of the pairs whose anchors overlap the two
fragments). Controls: the same pair shifted along the chromosome by --shifts random offsets (both fragments by one
offset, a whole number of bins of the contact map, larger than the distance of the pair). Empirical
P = (1 + controls >= observed) / (1 + controls).

Outputs: <out>_pairs.tsv (one row per pair and data set) and <out>_summary.tsv (per data set, set and presence of the
path in the condition: pairs, fraction with P < 0.05, median log2 (observed + 1) / (control mean + 1)) and
<out>_per_gene.tsv (per explained gene: explaining pairs and pairs with P < 0.05 per data set).

Coordinates of a contact map on another assembly are lifted with --coolChain LABEL=CHAIN (start and end of each
fragment; pairs whose fragments do not lift to one chromosome in order are skipped).

Usage: python 39_step2_pairs_contacts.py <multiway_paths prefix> <graph.graphml> <out prefix>
           --cool LABEL=FILE.cool ... --pairs LABEL=FILE.txt.gz ... [--coolChain LABEL=CHAIN] [--shifts 50]
           [--seed 0]
"""

import argparse
import gzip

import cooler
import networkx as nx
import numpy as np
import pandas as pd


def norm_chrom(c):
    return c[3:] if c.startswith('chr') else c


def span(graph, node):
    return int(graph.nodes[node]['start']), int(graph.nodes[node]['end'])


def lift_pairs(pairs, chain):
    """The pairs with start and end of both fragments lifted (pairs that do not lift in order dropped)."""
    from pyliftover import LiftOver
    lift = LiftOver(chain)

    def one(chrom, start, end):
        a, b = lift.convert_coordinate(chrom, int(start)), lift.convert_coordinate(chrom, int(end))
        if not a or not b or a[0][0] != b[0][0] or b[0][1] <= a[0][1]:
            return None
        return a[0][0], int(a[0][1]), int(b[0][1])

    rows = {}
    for i, r in pairs.iterrows():
        p, t = one(r.chrom, r.p_start, r.p_end), one(r.chrom, r.t_start, r.t_end)
        if p and t and p[0] == t[0]:
            rows[i] = (p[0], p[1], p[2], t[1], t[2])
    lifted = pairs.loc[list(rows)].copy()
    lifted[['chrom', 'p_start', 'p_end', 't_start', 't_end']] = pd.DataFrame.from_dict(rows, orient='index')
    return lifted


def cool_counts(path, pairs, shifts, rng):
    """observed, control mean, P for each pair (dict index -> tuple) from a cooler file."""
    c = cooler.Cooler(path)
    binsize = c.binsize
    names = {norm_chrom(n): n for n in c.chromnames}
    out = {}
    for chrom, sub in pairs.groupby('chrom'):
        name = names.get(norm_chrom(chrom))
        if name is None:
            continue
        matrix = c.matrix(balance=False, sparse=True).fetch(name).tocsr()
        n = matrix.shape[0]

        def value(a0, a1, b0, b1):
            return float(matrix[a0:a1 + 1, b0:b1 + 1].sum())

        for i, r in sub.iterrows():
            a0, a1 = r.p_start // binsize, (r.p_end - 1) // binsize
            b0, b1 = r.t_start // binsize, (r.t_end - 1) // binsize
            if max(a1, b1) >= n:
                continue
            observed = value(a0, a1, b0, b1)
            distance = abs((a0 + a1) - (b0 + b1)) // 2 + max(a1 - a0, b1 - b0) + 1
            low, high = -min(a0, b0), n - 1 - max(a1, b1)
            controls, tries = [], 0
            while len(controls) < shifts and tries < shifts * 20 and high - low > 2 * distance:
                tries += 1
                k = int(rng.integers(low, high + 1))
                if abs(k) <= distance:
                    continue
                controls.append(value(a0 + k, a1 + k, b0 + k, b1 + k))
            if controls:
                controls = np.array(controls)
                out[i] = (observed, float(controls.mean()),
                          (1 + int(np.sum(controls >= observed))) / (1 + len(controls)))
    return out


def pair_list_counts(path, pairs, shifts, rng, step=5000):
    """observed, control mean, P for each pair from a list of pairs with counts (both anchors overlapping)."""
    table = pd.read_csv(path, sep='\t', header=None, usecols=range(7),
                        names=['c1', 's1', 'e1', 'c2', 's2', 'e2', 'count'])
    table = table[table.c1 == table.c2]
    table['chrom'] = table.c1.astype(str).map(norm_chrom)
    out = {}
    for chrom, sub in pairs.groupby('chrom'):
        t = table[table.chrom == norm_chrom(chrom)]
        if not len(t):
            continue
        s1, e1, s2, e2, cnt = (t[x].to_numpy() for x in ('s1', 'e1', 's2', 'e2', 'count'))
        size = int(max(e1.max(), e2.max()))

        def value(ps, pe, ts, te):
            hit = ((s1 < pe) & (e1 > ps) & (s2 < te) & (e2 > ts)) | ((s1 < te) & (e1 > ts) & (s2 < pe) & (e2 > ps))
            return float(cnt[hit].sum())

        for i, r in sub.iterrows():
            observed = value(r.p_start, r.p_end, r.t_start, r.t_end)
            distance = abs((r.p_start + r.p_end) - (r.t_start + r.t_end)) // 2 + max(r.p_end - r.p_start,
                                                                                     r.t_end - r.t_start)
            low = -min(r.p_start, r.t_start) // step
            high = (size - max(r.p_end, r.t_end)) // step
            controls, tries = [], 0
            while len(controls) < shifts and tries < shifts * 20 and high - low > 2:
                tries += 1
                k = int(rng.integers(low, high + 1)) * step
                if abs(k) <= distance:
                    continue
                controls.append(value(r.p_start + k, r.p_end + k, r.t_start + k, r.t_end + k))
            if controls:
                controls = np.array(controls)
                out[i] = (observed, float(controls.mean()),
                          (1 + int(np.sum(controls >= observed))) / (1 + len(controls)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('paths')
    parser.add_argument('graph')
    parser.add_argument('output')
    parser.add_argument('--cool', nargs='*', default=[])
    parser.add_argument('--pairs', nargs='*', default=[])
    parser.add_argument('--coolChain', nargs='*', default=[], help='LABEL=CHAIN: lift the pairs for this map')
    parser.add_argument('--shifts', type=int, default=50)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    graph = nx.read_graphml(args.graph)
    triplets = pd.read_csv(args.paths + '_triplets.tsv', sep='\t')
    pairs = triplets.drop_duplicates(['set', 'promoter', 'target'])[['set', 'genes', 'promoter', 'target', 'present']]
    pairs = pairs.reset_index(drop=True)
    pairs['chrom'] = pairs.promoter.map(lambda n: graph.nodes[n]['chrom'])
    pairs[['p_start', 'p_end']] = pd.DataFrame([span(graph, n) for n in pairs.promoter], index=pairs.index)
    pairs[['t_start', 't_end']] = pd.DataFrame([span(graph, n) for n in pairs.target], index=pairs.index)
    pairs['distance_kb'] = (((pairs.t_start + pairs.t_end) - (pairs.p_start + pairs.p_end)).abs() / 2000).round(1)
    print('{} explaining and {} background pairs'.format(int((pairs.set == 'explaining').sum()),
                                                          int((pairs.set == 'background').sum())), flush=True)
    rows = []
    chains = dict(item.split('=', 1) for item in args.coolChain)
    sources = [(item.split('=', 1), cool_counts) for item in args.cool] + \
        [(item.split('=', 1), pair_list_counts) for item in args.pairs]
    for (label, path), function in sources:
        query = lift_pairs(pairs, chains[label]) if label in chains else pairs
        counts = function(path, query, args.shifts, np.random.default_rng(args.seed))
        table = pairs.copy()
        table['data_set'] = label
        table[['observed', 'control_mean', 'p_empirical']] = pd.DataFrame(
            [counts.get(i, (np.nan, np.nan, np.nan)) for i in table.index], index=table.index)
        table['log2_ratio'] = np.log2((table.observed + 1) / (table.control_mean + 1))
        rows.append(table)
        print(label, 'done', flush=True)
    table = pd.concat(rows, ignore_index=True)
    table.to_csv(args.output + '_pairs.tsv', sep='\t', index=False)
    tested = table.dropna(subset=['p_empirical'])
    summary = tested.groupby(['data_set', 'set', 'present']).agg(
        pairs=('p_empirical', 'size'), p_below_005=('p_empirical', lambda x: float((x < 0.05).mean())),
        median_log2_ratio=('log2_ratio', 'median'), median_distance_kb=('distance_kb', 'median')).reset_index()
    summary.to_csv(args.output + '_summary.tsv', sep='\t', index=False)
    per_gene = tested[tested.set == 'explaining'].assign(significant=lambda x: x.p_empirical < 0.05).groupby(
        ['genes', 'promoter', 'present', 'data_set']).significant.agg(['size', 'sum']).unstack()
    per_gene.to_csv(args.output + '_per_gene.tsv', sep='\t')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
