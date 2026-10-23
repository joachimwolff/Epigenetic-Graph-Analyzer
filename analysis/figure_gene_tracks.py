"""Genome-browser panel of one gene report with pyGenomeTracks, and its expression, two PDFs.

<prefix>_tracks.pdf      genes (the gene of the report in colour, its promoter fragment shaded across all tracks),
                         promoter capture Hi-C arcs of both conditions (interactions of the promoter fragment and of
                         the fragments on the paths of the report), RNA-seq coverage (--rna) and chromatin tracks of
                         both conditions, and the fragments of each step of the report
<prefix>_expression.pdf  counts per million of the gene in the replicates of both conditions (--counts)

The fragments of each step are those of figure_gene_report.py: the promoter, the nodes on the shortest paths to the
nodes at --step whose mark changes at least two-fold in the direction of the expression change, and these nodes.

The arcs are the edges of the union graph of the comparison at the fragments of the report; the edges of the
explaining paths (promoter to the concordant nodes at --step) are drawn in colour over the others. With --zoom, a second
PDF <prefix>_zoom.pdf shows the genes, RNA-seq, chromatin tracks and fragments of the report in a smaller region.

Usage: python figure_gene_tracks.py <cmp prefix> <promoter node> <mark> <genes.gtf.gz> <output prefix>
           --tracks NAME=A.bw,B.bw ... [--rna A.bw,B.bw] [--counts counts.tsv.gz] [--conditions hESC hNEC]
           [--step 3] [--gene PRRX2] [--geneTypes protein_coding lncRNA] [--pgt pyGenomeTracks] [--padding 60000]
"""

import argparse
import gzip
import os
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig

STEP_COLOURS = {0: '#222222', 1: '#7f7f7a', 2: '#b07f2a', 3: '#b2182b'}
GENE_COLOUR, HIGHLIGHT_COLOUR, PATH_COLOUR = '#7f7f7a', '#d95f02', '#d95f02'


def log2_change(data, mark, a, b):
    return float(np.log2((data.get('{}_{}'.format(mark, b), 0.0) + 1) / (data.get('{}_{}'.format(mark, a), 0.0) + 1)))


def report_nodes(graph, source, mark, a, b, step, sign):
    distance = nx.single_source_shortest_path_length(graph, source, cutoff=step)
    targets = [n for n, k in distance.items() if k == step and sign * log2_change(graph.nodes[n], mark, a, b) >= 1]
    by_step, edges = {0: {source}, step: set(targets)}, set()
    for t in targets:
        for path in nx.all_shortest_paths(graph, source, t):
            for k in range(1, step):
                by_step.setdefault(k, set()).add(path[k])
            edges.update(frozenset(e) for e in zip(path[:-1], path[1:]))
    return by_step, edges


def write_links(graph, nodes, condition, start, end, path, edges=None):
    """Edges of the union graph with at least one end on a fragment of the report (or only the given edges) present
    in one condition (shared and lost in A, shared and gained in B) with both ends in the region, as pyGenomeTracks
    links with the interaction strength."""
    keep = {'a': ('shared', 'lost'), 'b': ('shared', 'gained')}[condition]
    rows = []
    for u, v, data in graph.edges(nodes, data=True):
        if data.get('status', 'shared') not in keep or (edges is not None and frozenset((u, v)) not in edges):
            continue
        (s1, e1), (s2, e2) = sorted((int(graph.nodes[x]['start']), int(graph.nodes[x]['end'])) for x in (u, v))
        if s1 < start or e2 > end or (s1, s2) in {(r[1], r[4]) for r in rows}:
            continue
        rows.append((graph.nodes[u]['chrom'], s1, e1, graph.nodes[v]['chrom'], s2, e2,
                     float(data.get('r_' + condition, 1.0))))
    with open(path, 'w') as out:
        for r in rows:
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n'.format(*r))
    return len(rows)


def write_genes(gtf, chrom, start, end, path, types, highlight):
    """Genes of the given types overlapping the region as BED9, the highlighted gene in its own colour."""
    with gzip.open(gtf, 'rt') as handle, open(path, 'w') as out:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] != 'gene' or f[0] != chrom or int(f[4]) < start or int(f[3]) > end:
                continue
            attrs = dict(x.strip().split(' ', 1) for x in f[8].strip().rstrip(';').split(';') if x.strip())
            name = attrs['gene_name'].strip('"')
            if attrs.get('gene_type', '').strip('"') not in types and name != highlight:
                continue
            colour = (HIGHLIGHT_COLOUR if name == highlight else GENE_COLOUR).lstrip('#')
            rgb = ','.join(str(int(colour[i:i + 2], 16)) for i in (0, 2, 4))
            out.write('{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\n'.format(chrom, int(f[3]) - 1, f[4], name, f[6],
                                                                   int(f[3]) - 1, int(f[3]) - 1, rgb))


def panel_expression(counts, gene, conditions, path, gene_column, prefix, labels):
    table = pd.read_csv(counts, sep=',' if counts.endswith(('.csv', '.csv.gz')) else '\t')
    # GEO tables may repeat a column; pandas renames the copies to <name>.1
    samples = [c for c in table.columns if c.startswith(prefix) and '.' not in c.split('_')[-1]]
    cpm = table[samples] / table[samples].sum() * 1e6
    row = cpm[table[gene_column] == gene].sum()
    groups = {c: [s for s in samples if '_{}_'.format(c) in s] for c in conditions}
    fig, ax = plt.subplots(figsize=(1.5, 2.0))
    for i, (cond, cols) in enumerate(groups.items()):
        values = [row[c] for c in cols]
        ax.bar(i, np.mean(values), color=['#9a9a96', '#2a6fb0'][i], width=0.6)
        ax.scatter([i] * len(values), values, color='black', s=10, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('{} counts per million'.format(gene), fontsize=8)
    ax.tick_params(axis='y', labelsize=7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return {c: [round(row[s], 1) for s in cols] for c, cols in groups.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('prefix', 'promoter', 'mark', 'gtf', 'output'):
        parser.add_argument(name)
    parser.add_argument('--tracks', nargs='+', required=True, help='NAME=A.bw,B.bw')
    parser.add_argument('--rna', help='RNA-seq coverage A.bw,B.bw')
    parser.add_argument('--counts', help='Count table for the expression panel')
    parser.add_argument('--countsGeneColumn', default='Probe')
    parser.add_argument('--countsPrefix', default='RNA_Seq', help='Prefix of the sample columns of the count table')
    parser.add_argument('--geneTypes', nargs='+', default=['protein_coding', 'lincRNA', 'lncRNA'])
    parser.add_argument('--conditions', nargs=2, default=['hESC', 'hNEC'])
    parser.add_argument('--expressionConditions', nargs=2, default=['hESC', 'hNEC'],
                        help='Condition names as they appear in the column names of the count table')
    parser.add_argument('--step', type=int, default=3)
    parser.add_argument('--gene', default=None)
    parser.add_argument('--pgt', default=os.path.expanduser('~/miniconda3/envs/hicexplorer/bin/pyGenomeTracks'))
    parser.add_argument('--padding', type=int, default=60000)
    parser.add_argument('--zoom', help='start-end of a second, smaller region on the same chromosome')
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    sign = np.sign(float(levels.loc[levels.node == args.promoter, 'rna_log2FC'].iloc[0]))
    gene = args.gene or graph.nodes[args.promoter]['genes'].split(',')[0]
    by_step, path_edges = report_nodes(graph, args.promoter, args.mark, a, b, args.step, sign)
    nodes = [n for v in by_step.values() for n in v]
    chrom = graph.nodes[args.promoter]['chrom']
    start = min(int(graph.nodes[n]['start']) for n in nodes) - args.padding
    end = max(int(graph.nodes[n]['end']) for n in nodes) + args.padding
    out_dir = os.path.dirname(os.path.abspath(args.output))
    stem = os.path.basename(args.output)
    counts = {}
    for cond, key in ((a, 'a'), (b, 'b')):
        counts[cond] = (write_links(graph, nodes, key, start, end, os.path.join(out_dir, '{}_{}.links'.format(stem, cond))),
                        write_links(graph, nodes, key, start, end,
                                    os.path.join(out_dir, '{}_{}_paths.links'.format(stem, cond)), path_edges))
    write_genes(args.gtf, chrom, start, end, os.path.join(out_dir, stem + '_genes.bed'), args.geneTypes, gene)
    with open(os.path.join(out_dir, stem + '_promoter.bed'), 'w') as out:
        out.write('{}\t{}\t{}\n'.format(chrom, graph.nodes[args.promoter]['start'], graph.nodes[args.promoter]['end']))
    with open(os.path.join(out_dir, stem + '_steps.bed'), 'w') as out:
        for k in sorted(by_step):
            for n in sorted(by_step[k], key=lambda n: int(graph.nodes[n]['start'])):
                c = STEP_COLOURS[k].lstrip('#')
                rgb = ','.join(str(int(c[i:i + 2], 16)) for i in (0, 2, 4))
                label = 'promoter' if k == 0 else 'step {}'.format(k)
                out.write('{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t{}\n'.format(chrom, graph.nodes[n]['start'],
                                                                     graph.nodes[n]['end'], label,
                                                                     graph.nodes[n]['start'], graph.nodes[n]['end'],
                                                                     rgb))
    ini = ['[x-axis]', 'where = top', 'fontsize = 9', '',
           '[genes]', 'file = {}_genes.bed'.format(stem), 'title = genes', 'height = 1.6', 'fontsize = 8',
           'style = flybase', 'labels = true', 'color = bed_rgb', 'color_utr = bed_rgb', '',
           '[spacer]', 'height = 0.2', '']
    for cond, colour in ((a, '#c4c4c0'), (b, '#a9c6e3')):
        common = ['links_type = arcs', 'orientation = inverted' if cond == a else '', 'compact_arcs_level = 2']
        ini += ['[links {}]'.format(cond), 'file = {}_{}.links'.format(stem, cond),
                'title = pcHi-C {}'.format(cond), 'height = 1.8', 'color = {}'.format(colour),
                'line_width = 0.6', 'fontsize = 8'] + common + ['',
                '[links {} paths]'.format(cond), 'file = {}_{}_paths.links'.format(stem, cond),
                'overlay_previous = share-y', 'color = {}'.format(PATH_COLOUR), 'line_width = 1.0'] + common + ['']
    specs = (['RNA-seq=' + args.rna] if args.rna else []) + args.tracks
    for spec in specs:
        name, files = spec.split('=')
        fa, fb = files.split(',')
        # one scale for both conditions of a mark: the larger 99th percentile of 600 bins of the region
        top = 0.0
        for path in (fa, fb):
            bw = pyBigWig.open(path)
            values = np.nan_to_num(np.array(bw.stats(chrom, start, end, type='mean', nBins=600), dtype=float))
            bw.close()
            top = max(top, float(np.percentile(values, 99)))
        for cond, path, colour in ((a, fa, '#7f7f7a'), (b, fb, '#2a6fb0')):
            ini += ['[{} {}]'.format(name, cond), 'file = {}'.format(path), 'title = {} {}'.format(name, cond),
                    'height = 1.0', 'color = {}'.format(colour), 'summary_method = mean', 'number_of_bins = 600',
                    'min_value = 0', 'max_value = {:.1f}'.format(top), 'fontsize = 8', '']
        ini += ['[spacer]', 'height = 0.15', '']
    ini += ['[steps]', 'file = {}_steps.bed'.format(stem), 'title = report', 'height = 1.0', 'color = bed_rgb',
            'labels = false', 'display = collapsed', 'fontsize = 8', '',
            '[promoter]', 'file = {}_promoter.bed'.format(stem), 'type = vhighlight',
            'color = {}'.format(HIGHLIGHT_COLOUR), 'alpha = 0.3', '']
    ini_path = os.path.join(out_dir, stem + '_tracks.ini')
    open(ini_path, 'w').write('\n'.join(ini))
    subprocess.run([args.pgt, '--tracks', ini_path, '--region', '{}:{}-{}'.format(chrom, start, end),
                    '--outFileName', args.output + '_tracks.pdf', '--width', '18', '--fontSize', '9'],
                   check=True, cwd=out_dir)
    if args.zoom:
        # the same tracks without the arcs, which leave a small region
        lines, skip = [], False
        for line in ini:
            skip = line.startswith('[links') or (skip and not line.startswith('['))
            if not skip:
                lines.append('fontsize = 7' if line == 'fontsize = 9' else line)
        zoom_ini = os.path.join(out_dir, stem + '_zoom.ini')
        open(zoom_ini, 'w').write('\n'.join(lines))
        subprocess.run([args.pgt, '--tracks', zoom_ini, '--region', '{}:{}'.format(chrom, args.zoom),
                        '--outFileName', args.output + '_zoom.pdf', '--width', '10', '--fontSize', '9',
                        '--trackLabelFraction', '0.3'],
                       check=True, cwd=out_dir)
    expression = None
    if args.counts:
        expression = panel_expression(args.counts, gene, args.expressionConditions, args.output + '_expression.pdf',
                                      args.countsGeneColumn, args.countsPrefix, args.conditions)
    print('{}: region {}:{}-{}, interactions {}, steps {}, CPM {}'.format(
        gene, chrom, start, end, counts, {k: len(v) for k, v in by_step.items()}, expression))


if __name__ == '__main__':
    main()
