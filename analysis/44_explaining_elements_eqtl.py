"""eQTL of the gene in its explaining elements beyond step 1, against fragments at the same distance with the same
chromatin change but without a graph connection (window), and against graph neighbours at the same step without a
concordant change.

Genes: promoter fragments without an explanation at the promoter (marks, contacts, state) or at step 1 (marks, state)
and with a candidate explanation at step 2 or, failing that, step 3 (explains_declared_mean); optionally only those
with a per-gene P below --maxP beyond the mirrored loci. Explaining elements: the nodes at the explaining step whose
explaining marks change concordantly (as in 41_replication.py). For every explaining element at distance d from the
promoter, the controls are the HindIII fragments at a distance between d / --distanceFactor and d * --distanceFactor
(either side) that are
  window   not within three steps of the promoter in the union graph, with a concordant change of one of the same
           explaining marks (sign x expression sign x log2 change >= --minChange, the same signal as the graph nodes:
           exact mean of the bigWig over the fragment);
  distance not within three steps, any chromatin;
  graph    nodes at the same step without a concordant change.
An element (or a control fragment) holds an eQTL when a significant variant of one of the genes of the promoter
fragment (gene_ids) lies inside it. Only elements within --maxDistance of the promoter are evaluated (the cis window of
the eQTL mapping). Reported: the share of explaining elements with an eQTL and, per control, the expected share (mean
over elements of the share among its controls) and the difference to the explaining elements (over the elements with
controls), with a bootstrap over genes. Per gene (<out>_genes.tsv), the explaining elements within --maxDistance are
set against all window fragments within --maxDistance (not within three steps, concordant change of an explaining
mark): their numbers and how many of each hold a link to the gene; the summary adds the pooled shares and the median
number of window fragments per explaining element.

Usage: python 44_explaining_elements_eqtl.py <run prefix> <rmap> <eQTL significant tsv(.gz)> <out prefix>
           --tracksA mark=bw ... --tracksB mark=bw ... [--maxP 1] [--minChange 1] [--distanceFactor 1.25]
           [--maxDistance 1000000] [--bootstrap 2000] [--seed 0]
"""

import argparse
import bisect
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module  # noqa: E402

replication = import_module('41_replication')

EARLY = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state']


def tracks(values):
    return dict(v.split('=', 1) for v in values)


class Signal:
    """Exact bigWig mean over a fragment, per mark and condition, cached."""

    def __init__(self, tracks_a, tracks_b):
        self.handles = {(m, 0): pyBigWig.open(p) for m, p in tracks_a.items()}
        self.handles.update({(m, 1): pyBigWig.open(p) for m, p in tracks_b.items()})
        self.cache = {}

    def mean(self, mark, condition, chrom, start, end):
        key = (mark, condition, chrom, start)
        if key not in self.cache:
            try:
                value = self.handles[(mark, condition)].stats(chrom, start, end, type='mean', exact=True)[0]
            except RuntimeError:
                value = None
            self.cache[key] = float(value or 0.0)
        return self.cache[key]

    def change(self, mark, chrom, start, end):
        return float(np.log2((self.mean(mark, 1, chrom, start, end) + 1) / (self.mean(mark, 0, chrom, start, end) + 1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('run', 'rmap', 'eqtl', 'output'):
        parser.add_argument(name)
    parser.add_argument('--tracksA', nargs='+', required=True)
    parser.add_argument('--tracksB', nargs='+', required=True)
    parser.add_argument('--maxP', type=float, default=1.0)
    parser.add_argument('--minChange', type=float, default=1.0)
    parser.add_argument('--distanceFactor', type=float, default=1.25)
    parser.add_argument('--maxDistance', type=int, default=1000000)
    parser.add_argument('--bootstrap', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    levels = pd.read_csv(args.run + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    levels['explains'] = levels.explains_declared_mean.astype(str) == 'True'
    graph = nx.read_graphml(args.run + '_union.graphml')
    conditions = replication.conditions_of(graph)
    early = set(levels[levels.explains & levels.level.isin(EARLY)].node)
    later = levels[~levels.node.isin(early) & levels.level.isin(['step2', 'step3']) & levels.explains &
                   (levels.p_declared < args.maxP)].sort_values('level').drop_duplicates('node')
    rmap = pd.read_csv(args.rmap, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
                       dtype={'chrom': str})
    rmap['chrom'] = 'chr' + rmap.chrom.astype(str).str.replace('^chr', '', regex=True)
    rmap['mid'] = (rmap.start + rmap.end) // 2
    fragments = {c: sub.sort_values('mid').reset_index(drop=True) for c, sub in rmap.groupby('chrom')}
    eqtl = pd.read_csv(args.eqtl, sep='\t')
    eqtl['gene_id'] = eqtl.gene_id.str.split('.').str[0]
    variants = {g: sub for g, sub in eqtl.groupby('gene_id')}
    signal = Signal(tracks(args.tracksA), tracks(args.tracksB))

    def holds(gene_ids, chrom, start, end):
        for g in gene_ids:
            sub = variants.get(g)
            if sub is not None and ((sub.chrom == chrom) & (sub.pos > start) & (sub.pos <= end)).any():
                return True
        return False

    rows, gene_rows = [], []
    for row in later.itertuples():
        node = row.node
        attrs = graph.nodes[node]
        gene_ids = [g.split('.')[0] for g in str(attrs.get('gene_ids', '')).split(',') if g]
        gene_ids = [g for g in gene_ids if g in variants]
        if not gene_ids:
            continue
        chrom = attrs['chrom']
        centre = (int(attrs['start']) + int(attrs['end'])) // 2
        sign = np.sign(float(row.rna_log2FC))
        step = int(row.level.replace('step', ''))
        signs = dict(x.split(':') for x in str(row.signs_declared).split(';') if ':' in x)
        marks = [m for m in str(row.concordant_marks_declared).split(',') if m in signs]
        elements = replication.explaining_elements(graph, levels, node, row.level, sign, conditions,
                                                   args.minChange)
        distance = nx.single_source_shortest_path_length(graph, node, cutoff=3)
        near = {(graph.nodes[n]['chrom'], int(graph.nodes[n]['start'])) for n in distance}
        same_step = [n for n, d in distance.items() if d == step and n not in elements]
        table = fragments[chrom]

        def concordant(f):
            return any(float(signs[m]) * sign * signal.change(m, chrom, int(f.start), int(f.end)) >= args.minChange
                       for m in marks)
        kept = [n for n in elements if graph.nodes[n]['chrom'] == chrom and
                abs((int(graph.nodes[n]['start']) + int(graph.nodes[n]['end'])) // 2 - centre) <= args.maxDistance]
        if kept:
            i = bisect.bisect_left(table.mid.values, centre - args.maxDistance)
            j = bisect.bisect_right(table.mid.values, centre + args.maxDistance)
            around = table.iloc[i:j]
            around = [f for f in around.itertuples() if (chrom, int(f.start)) not in near and concordant(f)]
            gene_rows.append(dict(
                gene_node=node, genes=row.genes, step=step, p=row.p_declared, elements=len(kept),
                elements_linked=sum(holds(gene_ids, chrom, int(graph.nodes[n]['start']), int(graph.nodes[n]['end']))
                                    for n in kept),
                window=len(around), window_linked=sum(holds(gene_ids, chrom, int(f.start), int(f.end))
                                                      for f in around)))
        for element in elements:
            e = graph.nodes[element]
            d = abs((int(e['start']) + int(e['end'])) // 2 - centre)
            if e['chrom'] != chrom or d > args.maxDistance:
                continue
            low, high = d / args.distanceFactor, d * args.distanceFactor
            candidates = []
            for lo, hi in ((centre + low, centre + high), (centre - high, centre - low)):
                i, j = bisect.bisect_left(table.mid.values, lo), bisect.bisect_right(table.mid.values, hi)
                candidates.append(table.iloc[i:j])
            candidates = pd.concat(candidates)
            candidates = candidates[[(chrom, int(s)) not in near for s in candidates.start]]
            window = [f for f in candidates.itertuples() if concordant(f)]
            graph_controls = [graph.nodes[n] for n in same_step
                              if low <= abs((int(graph.nodes[n]['start']) + int(graph.nodes[n]['end'])) // 2 -
                                            centre) <= high]

            def share(items):
                if not len(items):
                    return np.nan
                return float(np.mean([holds(gene_ids, chrom, int(f[0]), int(f[1])) for f in items]))
            rows.append(dict(
                gene_node=node, genes=row.genes, step=step, p=row.p_declared, element=element, distance=d,
                eqtl=holds(gene_ids, chrom, int(e['start']), int(e['end'])),
                window_n=len(window), window=share([(f.start, f.end) for f in window]),
                distance_n=len(candidates), distance_share=share(list(zip(candidates.start, candidates.end))),
                graph_n=len(graph_controls),
                graph=share([(g['start'], g['end']) for g in graph_controls])))
    table = pd.DataFrame(rows)
    table.to_csv(args.output + '_elements.tsv', sep='\t', index=False)
    rng = np.random.default_rng(args.seed)
    genes = table.gene_node.unique()
    by_gene = {g: sub for g, sub in table.groupby('gene_node')}
    columns = ['eqtl', 'window', 'distance_share', 'graph']

    def statistics(sub):
        out = {c: sub[c].astype(float).mean() for c in columns}
        for c in columns[1:]:
            paired = sub[sub[c].notna()]
            out['eqtl_minus_' + c] = paired.eqtl.astype(float).mean() - paired[c].mean()
        return out
    point = statistics(table)
    boot = {c: [] for c in point}
    for _ in range(args.bootstrap):
        sample = statistics(pd.concat([by_gene[g] for g in rng.choice(genes, size=len(genes))]))
        for c in point:
            boot[c].append(sample[c])
    summary = pd.DataFrame([dict(set=c, share=point[c], low=np.nanpercentile(boot[c], 2.5),
                                 high=np.nanpercentile(boot[c], 97.5)) for c in point])
    per_gene = pd.DataFrame(gene_rows)
    per_gene.to_csv(args.output + '_genes.tsv', sep='\t', index=False)
    extra = [dict(set='pooled_elements_linked', share=per_gene.elements_linked.sum() / per_gene.elements.sum()),
             dict(set='pooled_window_linked', share=per_gene.window_linked.sum() / max(1, per_gene.window.sum())),
             dict(set='median_elements_per_gene', share=per_gene.elements.median()),
             dict(set='median_window_per_gene', share=per_gene.window.median()),
             dict(set='genes_with_linked_element', share=(per_gene.elements_linked > 0).mean()),
             dict(set='genes_with_linked_window', share=(per_gene.window_linked > 0).mean())]
    summary = pd.concat([summary, pd.DataFrame(extra)], ignore_index=True)
    summary.insert(0, 'elements', len(table))
    summary.insert(0, 'genes', len(genes))
    summary.to_csv(args.output + '_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
