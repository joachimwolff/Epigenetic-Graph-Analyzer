"""Co-occurrence of explaining paths in multi-way data of one condition, by the presence of the paths in it.

For every gene explained at step 2 (egaDifferential, explains_declared_mean) and not explained at the promoter (marks,
contacts) or at step 1 (marks), the explaining paths are the triplets
promoter - mediator - target with a target two steps away whose marks change at least two-fold (log2 >= --minChange)
in the direction of the expression change (active marks with the sign of the change, repressive marks against it).
A triplet is present in the condition of the multi-way data (--condition B or A) when both of its interactions are
called there (status shared or gained in B; shared or lost in A), and absent otherwise. For every data set (multi-way
groups on restriction fragments, MultiwayData npz) the script counts the groups holding all three fragments and the
same triplet shifted by --shifts random offsets along the chromosome (all three nodes by one whole number of
restriction fragments, so that every shifted node covers as many fragments as the node; offsets larger than the span
of the triplet in fragments; every shifted node covered by a group). A random sample of all step-2
triplets of the graph (--background) gives the same counts for paths that explain nothing.

Outputs: <out>_triplets.tsv (one row per triplet and data set: observed, control mean, empirical P, log2 ratio),
<out>_per_gene.tsv (per explained gene and presence: whether at least one explaining triplet has P < 0.05 in each data
set, and in how many data sets) and
<out>_summary.tsv (per data set, set (explaining or background) and presence: triplets, median log2 ratio, fraction
with observed above every control, fraction with P < 0.05) with Mann-Whitney tests of present against absent.

Usage: python 36_gm12878_multiway_paths.py <cmp prefix> <out prefix> --data LABEL=FILE.npz ... [--condition B]
           [--shifts 50] [--background 20000] [--minChange 1] [--seed 0]
"""

import argparse
import multiprocessing
import os
import sys
import zlib

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_gene_multiway import ChromGroups, intersect  # noqa: E402

from epigeneticgraphanalyzer.multiway import MultiwayData  # noqa: E402


def change(attrs, mark, a, b):
    return float(np.log2((float(attrs.get('{}_{}'.format(mark, b), 0) or 0) + 1) /
                         (float(attrs.get('{}_{}'.format(mark, a), 0) or 0) + 1)))


def explaining_triplets(graph, levels, roles, a, b, min_change):
    """(gene, promoter, mediator, target) of the genes explained at step 2."""
    rows = []
    explains = levels['explains_declared_mean'].astype(str) == 'True'
    earlier = set(levels[explains & levels.level.isin(['promoter_marks', 'promoter_contacts', 'step1_marks'])].node)
    explained = levels[(levels.level == 'step2') & explains & ~levels.node.isin(earlier)]
    for row in explained.itertuples():
        source, sign = row.node, np.sign(row.rna_log2FC)
        if source not in graph:
            continue
        distance = nx.single_source_shortest_path_length(graph, source, cutoff=2)
        for target, k in distance.items():
            if k != 2:
                continue
            attrs = graph.nodes[target]
            if not any(sign * role * change(attrs, mark, a, b) >= min_change for mark, role in roles.items()):
                continue
            for path in nx.all_shortest_paths(graph, source, target):
                rows.append((row.genes, source, path[1], target))
    return pd.DataFrame(rows, columns=['genes', 'promoter', 'mediator', 'target']).drop_duplicates()


def background_triplets(graph, n, seed):
    rng = np.random.default_rng(seed)
    nodes = [v for v in graph.nodes if int(graph.nodes[v].get('bait', 0) or 0) == 1]
    rows = []
    for source in rng.permutation(nodes):
        neighbours = list(graph.neighbors(source))
        if not neighbours:
            continue
        mediator = neighbours[rng.integers(len(neighbours))]
        targets = [t for t in graph.neighbors(mediator) if t != source and not graph.has_edge(source, t)]
        if targets:
            rows.append(('', source, mediator, targets[rng.integers(len(targets))]))
        if len(rows) >= n:
            break
    return pd.DataFrame(rows, columns=['genes', 'promoter', 'mediator', 'target'])


def presence(graph, triplets, condition):
    keep = ('shared', 'gained') if condition == 'B' else ('shared', 'lost')
    return [all(graph.edges[u, v].get('status', 'shared') in keep for u, v in ((p, m), (m, t)))
            for p, m, t in zip(triplets.promoter, triplets.mediator, triplets.target)]


# shared with the worker processes (fork); each worker builds the groups of a chromosome once
DATA = None
GROUPS = {}


def count_chunk(task):
    chrom, chunk, rows, shifts, seed = task
    rng = np.random.default_rng([seed, zlib.crc32(chrom.encode()), chunk])
    out = []
    if chrom not in DATA.bins.offsets:
        return [(i, np.nan, np.nan, np.nan, 0) for i, _ in rows]
    if chrom not in GROUPS:
        GROUPS.clear()
        GROUPS[chrom] = ChromGroups(DATA, chrom)
    groups = GROUPS[chrom]
    for i, spans in rows:
        # the fragments of each node; controls shift all three by one whole number of fragments, so that every
        # shifted node covers as many fragments as the node itself
        units = [groups.overlapping(s + 1, e - 1) for s, e in spans]
        if any(not u for u in units):
            out.append((i, np.nan, np.nan, np.nan, 0))
            continue

        def holding(offset):
            sets = [groups.groups([f + offset for f in u]) for u in units]
            return None if any(len(x) == 0 for x in sets) else len(intersect(*sets))

        observed = holding(0)
        if observed is None:
            out.append((i, np.nan, np.nan, np.nan, 0))
            continue
        first, last = min(min(u) for u in units), max(max(u) for u in units)
        span = last - first
        low, high = -first, groups.n - 1 - last
        controls, tries = [], 0
        while len(controls) < shifts and tries < shifts * 20 and high - low > 2 * span:
            tries += 1
            offset = int(rng.integers(low, high + 1))
            if abs(offset) <= span:
                continue
            value = holding(offset)
            if value is not None:
                controls.append(value)
        controls = np.array(controls)
        if not len(controls):
            out.append((i, observed, np.nan, np.nan, 0))
            continue
        p = (1 + int(np.sum(controls >= observed))) / (1 + len(controls))
        out.append((i, observed, float(controls.mean()), p, len(controls)))
    return out


def count_triplets(graph, triplets, data, shifts, seed, threads, chunk_size=100):
    """Observed and shifted counts per triplet, in chunks of triplets of one chromosome per task."""
    global DATA
    DATA = data
    chroms = triplets.promoter.map(lambda n: graph.nodes[n]['chrom'])
    tasks = []
    for chrom, sub in triplets.groupby(chroms):
        rows = [(i, [(int(graph.nodes[n]['start']), int(graph.nodes[n]['end'])) for n in (p, m, t)])
                for i, p, m, t in zip(sub.index, sub.promoter, sub.mediator, sub.target)]
        for chunk, first in enumerate(range(0, len(rows), chunk_size)):
            tasks.append((chrom, chunk, rows[first:first + chunk_size], shifts, seed))
    # chunks of one chromosome in sequence, so that a worker mostly reuses the groups it built
    with multiprocessing.get_context('fork').Pool(threads) as pool:
        out = [row for rows in pool.imap(count_chunk, tasks, chunksize=1) for row in rows]
    table = pd.DataFrame(out, columns=['index', 'observed', 'control_mean', 'p_empirical', 'n_controls'])
    return table.set_index('index').reindex(triplets.index)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('prefix')
    parser.add_argument('output')
    parser.add_argument('--data', nargs='+', required=True, help='LABEL=FILE.npz')
    parser.add_argument('--conditions', nargs=2, default=['naiveB', 'GM12878'])
    parser.add_argument('--condition', choices=['A', 'B'], default='B', help='Condition of the multi-way data')
    parser.add_argument('--roles', nargs='+', default=['H3K27ac=1', 'H3K4me3=1', 'H3K4me1=1', 'H3K27me3=-1'])
    parser.add_argument('--shifts', type=int, default=50)
    parser.add_argument('--background', type=int, default=20000)
    parser.add_argument('--minChange', type=float, default=1.0)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threads', type=int, default=12)
    args = parser.parse_args()
    a, b = args.conditions
    roles = {k: float(v) for k, v in (r.split('=') for r in args.roles)}
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    explaining = explaining_triplets(graph, levels, roles, a, b, args.minChange).assign(set='explaining')
    background = background_triplets(graph, args.background, args.seed).assign(set='background')
    triplets = pd.concat([explaining, background], ignore_index=True)
    triplets['present'] = presence(graph, triplets, args.condition)
    print('{} explaining triplets of {} genes ({} present in {}), {} background triplets'.format(
        len(explaining), explaining.promoter.nunique(), int(triplets[triplets.set == 'explaining'].present.sum()),
        b if args.condition == 'B' else a, len(background)), flush=True)
    rows, summary = [], []
    for item in args.data:
        label, path = item.split('=', 1)
        data = MultiwayData.load(path)
        counts = count_triplets(graph, triplets, data, args.shifts, args.seed, args.threads)
        del data
        table = pd.concat([triplets, counts], axis=1).assign(data_set=label)
        table['log2_ratio'] = np.log2((table.observed + 1) / (table.control_mean + 1))
        rows.append(table)
        tested = table.dropna(subset=['control_mean'])
        for (name, present), sub in tested.groupby(['set', 'present']):
            summary.append(dict(data_set=label, set=name, present=present, triplets=len(sub),
                                median_log2_ratio=sub.log2_ratio.median(), mean_log2_ratio=sub.log2_ratio.mean(),
                                above_all_controls=float((sub.p_empirical <= 1 / (1 + sub.n_controls)).mean()),
                                p_below_005=float((sub.p_empirical < 0.05).mean())))
        for name in ('explaining', 'background'):
            sub = tested[tested.set == name]
            x, y = sub[sub.present].log2_ratio, sub[~sub.present].log2_ratio
            if len(x) and len(y):
                test = mannwhitneyu(x, y, alternative='greater')
                summary.append(dict(data_set=label, set=name, present='present > absent', triplets=len(sub),
                                    median_log2_ratio=x.median() - y.median(), p_mannwhitney=test.pvalue))
        print(label, 'done', flush=True)
    table = pd.concat(rows)
    table.to_csv(args.output + '_triplets.tsv', sep='\t', index=False)
    explaining = table[table.set == 'explaining'].assign(significant=lambda x: x.p_empirical < 0.05)
    per_gene = explaining.groupby(['genes', 'promoter', 'present', 'data_set']).significant.any().unstack()
    per_gene['data_sets_p_lt_005'] = per_gene.sum(axis=1)
    per_gene.reset_index().to_csv(args.output + '_per_gene.tsv', sep='\t', index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(args.output + '_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
