"""Genome-browser panels of one gene report, one PDF per explanation level, all over the same region.

<prefix>_level_step1.pdf     the promoter and its interactions with its direct contacts (step 1); the promoter shaded when
                             it explains the gene, the direct contacts with a concordant change when step 1 explains it
<prefix>_level_step<k>.pdf   for k = 2 .. --maxStep: the interactions of the shortest paths to the nodes k steps away with
                             a concordant change, coloured by the step they reach; only these nodes shaded

Every panel shows the genes (the gene of the report in orange, its body lightly and the promoter fragment of the report
strongly shaded), the interactions of the level in both
conditions and one log2 track (condition B over A, centred at its median in the region; red, higher in B; blue, higher
in A) per data type. Nodes are shaded only at a level that explains the gene in the egaDifferential levels table
(explains_declared_mean of promoter_marks, step1_marks, step<k>), and only for the marks that explain it there
(concordant_marks_declared): a node is shaded when, for one of these marks, sign * sign(expression change) *
log2 change >= --minChange (signs from signs_declared).

Usage: python figure_gene_level_tracks.py <cmp prefix> <promoter node> <genes.gtf.gz> <output prefix>
           --tracks NAME=A.bw,B.bw ... [--rna A.bw,B.bw] [--conditions A B] [--gene ASUN] [--maxStep 3]
           [--roles H3K27ac=1 H3K4me1=1 H3K4me3=1 H3K27me3=-1] [--minChange 1]
"""

import argparse
import os
import subprocess
import sys

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_gene_tracks import (HIGHLIGHT_COLOUR, STEP_COLOURS, difference_bedgraph, log2_change,  # noqa: E402
                                write_genes, write_links)


def concordant(attrs, roles, a, b, sign, min_change):
    return any(role * sign * log2_change(attrs, mark, a, b) >= min_change for mark, role in roles.items())


def write_bed(path, chrom, intervals):
    with open(path, 'w') as out:
        for s, e in intervals:
            out.write('{}\t{}\t{}\n'.format(chrom, s, e))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('prefix', 'promoter', 'gtf', 'output'):
        parser.add_argument(name)
    parser.add_argument('--tracks', nargs='+', required=True, help='NAME=A.bw,B.bw')
    parser.add_argument('--rna', help='RNA-seq coverage A.bw,B.bw')
    parser.add_argument('--conditions', nargs=2, default=['hESC', 'hNEC'])
    parser.add_argument('--gene', default=None)
    parser.add_argument('--geneTypes', nargs='+', default=['protein_coding'])
    parser.add_argument('--maxStep', type=int, default=3)
    parser.add_argument('--minChange', type=float, default=1.0)
    parser.add_argument('--padding', type=int, default=60000)
    parser.add_argument('--bins', type=int, default=1000)
    parser.add_argument('--pseudocount', nargs='+', default=['RNA-seq=0.1'])
    parser.add_argument('--pgt', default=os.path.expanduser('~/miniconda3/envs/hicexplorer/bin/pyGenomeTracks'))
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    sign = np.sign(float(levels.loc[levels.node == args.promoter, 'rna_log2FC'].iloc[0]))
    rows = levels[levels.node == args.promoter].set_index('level')

    def explaining_roles(level):
        """{mark: sign} of the marks that explain the gene at this level, or {} when the level does not explain it."""
        if level not in rows.index or str(rows.loc[level, 'explains_declared_mean']) != 'True':
            return {}
        signs = dict(x.split(':') for x in str(rows.loc[level, 'signs_declared']).split(';') if ':' in x)
        marks = [m for m in str(rows.loc[level, 'concordant_marks_declared']).split(',') if m in signs]
        return {m: float(signs[m]) for m in marks}

    gene = args.gene or graph.nodes[args.promoter]['genes'].split(',')[0]
    source = args.promoter
    chrom = graph.nodes[source]['chrom']
    span = lambda n: (int(graph.nodes[n]['start']), int(graph.nodes[n]['end']))  # noqa: E731

    # per level: the interactions (edge sets by the step they reach) and the explaining nodes
    distance = nx.single_source_shortest_path_length(graph, source, cutoff=args.maxStep)
    direct = [n for n, k in distance.items() if k == 1]
    roles = explaining_roles('promoter_marks')
    promoter = [source] if roles and concordant(graph.nodes[source], roles, a, b, sign, args.minChange) else []
    roles = explaining_roles('step1_marks')
    contacts = [n for n in direct if roles and concordant(graph.nodes[n], roles, a, b, sign, args.minChange)]
    panels = [('step1', {1: {frozenset((source, n)) for n in direct}}, promoter + contacts, 1)]
    for k in range(2, args.maxStep + 1):
        roles = explaining_roles('step{}'.format(k))
        # the paths to the concordant nodes when the level explains, otherwise to all nodes of the step, unshaded
        reached = [n for n, d in distance.items() if d == k]
        targets = [n for n in reached if roles and concordant(graph.nodes[n], roles, a, b, sign, args.minChange)]
        drawn = targets or reached
        edges = {}
        for t in drawn:
            for path in nx.all_shortest_paths(graph, source, t):
                for i, e in enumerate(zip(path[:-1], path[1:])):
                    edges.setdefault(i + 1, set()).add(frozenset(e))
        panels.append(('step{}'.format(k), edges, targets, k))
    # the region: the promoter, its direct contacts and the levels that explain the gene; a level that explains
    # nothing shows its interactions within this region only
    nodes = {n for label, edges, highlight, _ in panels if label == 'step1' or highlight
             for es in edges.values() for e in es for n in e}
    start = min(span(n)[0] for n in nodes) - args.padding
    end = max(span(n)[1] for n in nodes) + args.padding
    out_dir = os.path.dirname(os.path.abspath(args.output))
    stem = os.path.basename(args.output)
    body = write_genes(args.gtf, chrom, start, end, os.path.join(out_dir, stem + '_genes.bed'), args.geneTypes, gene)
    write_bed(os.path.join(out_dir, stem + '_gene.bed'), chrom, [body or span(source)])
    write_bed(os.path.join(out_dir, stem + '_promoter.bed'), chrom, [span(source)])
    pseudo = {k: float(v) for k, v in (x.split('=') for x in args.pseudocount)}
    tracks = []
    for spec in (['RNA-seq=' + args.rna] if args.rna else []) + args.tracks:
        name, files = spec.split('=')
        fa, fb = files.split(',')
        path, top = difference_bedgraph(fa, fb, chrom, start, end, args.bins, pseudo.get(name, 1.0),
                                        os.path.join(out_dir, '{}_{}_log2.bedgraph'.format(stem, name)))
        tracks += ['[{} log2]'.format(name), 'file = {}'.format(path), 'title = {}'.format(name), 'height = 0.7',
                   'color = #b2182b', 'negative_color = #2166ac', 'min_value = {:.2f}'.format(-top),
                   'max_value = {:.2f}'.format(top), 'fontsize = 8', '', '[spacer]', 'height = 0.05', '']
    for label, edges, highlight, k in panels:
        ini = ['[x-axis]', 'where = top', 'fontsize = 9', '',
               '[genes]', 'file = {}_genes.bed'.format(stem), 'title = genes', 'height = 1.1', 'fontsize = 8',
               'style = flybase', 'labels = true', 'color = bed_rgb', 'color_utr = bed_rgb', '',
               '[spacer]', 'height = 0.2', '']
        for cond, key in ((a, 'a'), (b, 'b')):
            first = True
            for step, es in sorted(edges.items()):
                path = os.path.join(out_dir, '{}_{}_{}_reach{}.links'.format(stem, label, cond, step))
                write_links(graph, {n for e in es for n in e}, key, start, end, path, es)
                head = ['title = pcHi-C {}'.format(cond), 'height = 1.1', 'fontsize = 8'] if first else \
                    ['overlay_previous = share-y']
                ini += ['[links {} {}]'.format(cond, step), 'file = {}'.format(os.path.basename(path))] + head + \
                    ['color = {}'.format(STEP_COLOURS[step]), 'line_width = 1.1', 'links_type = arcs',
                     'orientation = inverted' if cond == a else '', 'compact_arcs_level = 2', '']
                first = False
        ini += tracks
        # the gene body lightly, the promoter fragment of the report strongly
        ini += ['[gene]', 'file = {}_gene.bed'.format(stem), 'type = vhighlight',
                'color = {}'.format(HIGHLIGHT_COLOUR), 'alpha = 0.12', '',
                '[promoter]', 'file = {}_promoter.bed'.format(stem), 'type = vhighlight',
                'color = {}'.format(HIGHLIGHT_COLOUR), 'alpha = 0.45', '']
        if highlight:
            bed = os.path.join(out_dir, '{}_{}_explaining.bed'.format(stem, label))
            write_bed(bed, chrom, sorted(span(n) for n in highlight))
            ini += ['[explaining]', 'file = {}'.format(os.path.basename(bed)), 'type = vhighlight',
                    'color = {}'.format(STEP_COLOURS[k]), 'alpha = 0.35', '']
        ini_path = os.path.join(out_dir, '{}_{}.ini'.format(stem, label))
        open(ini_path, 'w').write('\n'.join(ini))
        subprocess.run([args.pgt, '--tracks', ini_path, '--region', '{}:{}-{}'.format(chrom, start, end),
                        '--outFileName', '{}_level_{}.pdf'.format(args.output, label), '--width', '18',
                        '--fontSize', '9'], check=True, cwd=out_dir, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        print('{} {}: {} interactions, {} explaining nodes'.format(gene, label, sum(len(v) for v in edges.values()),
                                                                     len(highlight)))


if __name__ == '__main__':
    main()
