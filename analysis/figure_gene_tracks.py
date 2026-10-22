"""Genome-browser panel of one gene report with pyGenomeTracks, and its expression, two PDFs.

<prefix>_tracks.pdf      genes, promoter capture Hi-C arcs of both conditions (interactions of the promoter fragment
                         and of the fragments on the paths of the report), chromatin tracks of both conditions and
                         the fragments of each step of the report
<prefix>_expression.pdf  counts per million of the gene in the replicates of both conditions

The fragments of each step are those of figure_gene_report.py: the promoter, the nodes on the shortest paths to the
nodes at --step whose mark changes at least two-fold in the direction of the expression change, and these nodes.

Usage: python figure_gene_tracks.py <cmp prefix> <promoter node> <mark> <washU A> <washU B> <genes.gtf.gz>
           <counts.tsv.gz> <output prefix> --tracks NAME=A.bw,B.bw ... [--conditions hESC hNEC] [--step 3]
           [--score 12] [--gene PRRX2] [--pgt pyGenomeTracks] [--padding 40000]
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


def log2_change(data, mark, a, b):
    return float(np.log2((data.get('{}_{}'.format(mark, b), 0.0) + 1) / (data.get('{}_{}'.format(mark, a), 0.0) + 1)))


def report_nodes(graph, source, mark, a, b, step, sign):
    distance = nx.single_source_shortest_path_length(graph, source, cutoff=step)
    targets = [n for n, k in distance.items() if k == step and sign * log2_change(graph.nodes[n], mark, a, b) >= 1]
    by_step = {0: {source}, step: set(targets)}
    for t in targets:
        for path in nx.all_shortest_paths(graph, source, t):
            for k in range(1, step):
                by_step.setdefault(k, set()).add(path[k])
    return by_step


def write_links(washu, chrom, start, end, anchors, score, path):
    """Interactions within the region with at least one end on an anchor fragment, as pyGenomeTracks links."""
    rows = []
    with gzip.open(washu, 'rt') as handle:
        for line in handle:
            fields = line.split('\t')
            if fields[0].replace('chr', '') != chrom.replace('chr', ''):
                continue
            s1, e1 = int(fields[1]), int(fields[2])
            other, value = fields[3].split(',')
            c2, span = other.split(':')
            s2, e2 = (int(x) for x in span.split('-'))
            if float(value) < score or min(s1, s2) < start or max(e1, e2) > end:
                continue
            if not any(s <= (s1 + e1) // 2 <= e or s <= (s2 + e2) // 2 <= e for s, e in anchors):
                continue
            if (s2, e2, s1, e1) in {(r[1], r[2], r[4], r[5]) for r in rows}:
                continue
            rows.append((chrom, s1, e1, chrom, s2, e2, float(value)))
    with open(path, 'w') as out:
        for r in rows:
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n'.format(*r))
    return len(rows)


def write_genes(gtf, chrom, start, end, path):
    """Protein-coding and lncRNA genes overlapping the region as BED6."""
    with gzip.open(gtf, 'rt') as handle, open(path, 'w') as out:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            if f[2] != 'gene' or f[0] != chrom or int(f[4]) < start or int(f[3]) > end:
                continue
            attrs = dict(x.strip().split(' ', 1) for x in f[8].strip().rstrip(';').split(';') if x.strip())
            if attrs.get('gene_type', '').strip('"') not in ('protein_coding', 'lincRNA'):
                continue
            out.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(chrom, int(f[3]) - 1, f[4], attrs['gene_name'].strip('"'),
                                                       f[6]))


def panel_expression(counts, gene, conditions, path):
    table = pd.read_csv(counts, sep='\t')
    samples = [c for c in table.columns if c.startswith('RNA_Seq')]
    cpm = table[samples] / table[samples].sum() * 1e6
    row = cpm[table['Probe'] == gene].iloc[0]
    groups = {c: [s for s in samples if '_{}_'.format(c) in s] for c in conditions}
    fig, ax = plt.subplots(figsize=(1.5, 2.0))
    for i, (cond, cols) in enumerate(groups.items()):
        values = [row[c] for c in cols]
        ax.bar(i, np.mean(values), color=['#9a9a96', '#2a6fb0'][i], width=0.6)
        ax.scatter([i] * len(values), values, color='black', s=10, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(list(groups), fontsize=8)
    ax.set_ylabel('{} counts per million'.format(gene), fontsize=8)
    ax.tick_params(axis='y', labelsize=7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return {c: [round(row[s], 1) for s in cols] for c, cols in groups.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('prefix', 'promoter', 'mark', 'washu_a', 'washu_b', 'gtf', 'counts', 'output'):
        parser.add_argument(name)
    parser.add_argument('--tracks', nargs='+', required=True, help='NAME=A.bw,B.bw')
    parser.add_argument('--conditions', nargs=2, default=['hESC', 'hNEC'])
    parser.add_argument('--expressionConditions', nargs=2, default=['hESC', 'hNEC'],
                        help='Condition names as they appear in the column names of the count table')
    parser.add_argument('--step', type=int, default=3)
    parser.add_argument('--score', type=float, default=12)
    parser.add_argument('--gene', default=None)
    parser.add_argument('--pgt', default=os.path.expanduser('~/miniconda3/envs/hicexplorer/bin/pyGenomeTracks'))
    parser.add_argument('--padding', type=int, default=60000)
    args = parser.parse_args()
    a, b = args.conditions
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    sign = np.sign(float(levels.loc[levels.node == args.promoter, 'rna_log2FC'].iloc[0]))
    gene = args.gene or graph.nodes[args.promoter]['genes'].split(',')[0]
    by_step = report_nodes(graph, args.promoter, args.mark, a, b, args.step, sign)
    nodes = [n for v in by_step.values() for n in v]
    chrom = graph.nodes[args.promoter]['chrom']
    start = min(int(graph.nodes[n]['start']) for n in nodes) - args.padding
    end = max(int(graph.nodes[n]['end']) for n in nodes) + args.padding
    out_dir = os.path.dirname(os.path.abspath(args.output))
    stem = os.path.basename(args.output)
    anchors = [(int(graph.nodes[n]['start']), int(graph.nodes[n]['end'])) for n in nodes]
    counts = {}
    for cond, washu in ((a, args.washu_a), (b, args.washu_b)):
        counts[cond] = write_links(washu, chrom, start, end, anchors, args.score,
                                   os.path.join(out_dir, '{}_{}.links'.format(stem, cond)))
    write_genes(args.gtf, chrom, start, end, os.path.join(out_dir, stem + '_genes.bed'))
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
           'style = flybase', 'labels = true', '', '[spacer]', 'height = 0.2', '']
    for cond, colour in ((a, '#7f7f7a'), (b, '#2a6fb0')):
        ini += ['[links {}]'.format(cond), 'file = {}_{}.links'.format(stem, cond),
                'title = pcHi-C {}'.format(cond), 'height = 1.8', 'color = {}'.format(colour),
                'line_width = 0.8', 'links_type = arcs', 'orientation = inverted' if cond == a else '',
                'compact_arcs_level = 2', 'fontsize = 8', '']
    for spec in args.tracks:
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
            'labels = false', 'display = collapsed', 'fontsize = 8', '']
    ini_path = os.path.join(out_dir, stem + '_tracks.ini')
    open(ini_path, 'w').write('\n'.join(ini))
    subprocess.run([args.pgt, '--tracks', ini_path, '--region', '{}:{}-{}'.format(chrom, start, end),
                    '--outFileName', args.output + '_tracks.pdf', '--width', '18', '--fontSize', '9'],
                   check=True, cwd=out_dir)
    expression = panel_expression(args.counts, gene, args.expressionConditions, args.output + '_expression.pdf')
    print('{}: region {}:{}-{}, interactions {}, steps {}, CPM {}'.format(
        gene, chrom, start, end, counts, {k: len(v) for k, v in by_step.items()}, expression))


if __name__ == '__main__':
    main()
