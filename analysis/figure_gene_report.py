"""Report of one gene explained at graph neighbours, three PDFs (one per panel).

<prefix>_levels.pdf  change of one mark at each level (promoter, step 1, 2, 3) next to the mirrored loci, with the
                     two-fold line and the per-gene P-value of every step (egaDifferential <cmp>_explanations_all_levels)
<prefix>_paths.pdf   paths from the promoter to the nodes at the explaining step whose mark changes at least two-fold
                     in the direction of the expression change, one column per step; fragments of a column within
                     --merge bp are drawn as one node, coloured by the mean mark change
<prefix>_track.pdf   mark signal of both conditions along the region of these paths, with the fragments marked by step

Usage: python figure_gene_report.py <cmp prefix> <promoter node> <mark> <bigwig A> <bigwig B> <output prefix>
           [--conditions hESC hNEC] [--step 3] [--merge 50000] [--name PRRX2]
<cmp prefix> is the output prefix of egaDifferential (<prefix>_union.graphml, <prefix>_explanations_all_levels.tsv).
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

INK, MUTED, GRAPH, MIRROR = '#222222', '#888888', '#2a6fb0', '#c8c8c4'
STEP_COLOURS = {0: '#222222', 1: '#7f7f7a', 2: '#b07f2a', 3: '#b2182b'}


def parse_changes(text, mark):
    if not isinstance(text, str):
        return np.nan
    for item in text.split(';'):
        name, value = item.split(':')
        if name == mark:
            return float(value)
    return np.nan


def log2_change(data, mark, a, b):
    return float(np.log2((data.get('{}_{}'.format(mark, b), 0.0) + 1) / (data.get('{}_{}'.format(mark, a), 0.0) + 1)))


def middle(data):
    return (int(data['start']) + int(data['end'])) // 2


def panel_levels(rows, mark, path):
    levels = [('promoter_marks', 'promoter'), ('step1_marks', 'step 1'), ('step2', 'step 2'), ('step3', 'step 3')]
    fig, ax = plt.subplots(figsize=(3.0, 2.4))
    for y, (level, label) in enumerate(levels):
        row = rows[rows.level == level]
        if row.empty:
            continue
        row = row.iloc[0]
        graph = parse_changes(row['mark_changes'], mark)
        mirror = parse_changes(row.get('mirror_changes'), mark)
        ax.barh(y - 0.18, graph, height=0.34, color=GRAPH, lw=0)
        if not np.isnan(mirror):
            ax.barh(y + 0.18, mirror, height=0.34, color=MIRROR, lw=0)
        p = row.get('p_declared')
        if isinstance(p, float) and not np.isnan(p):
            text = 'P = {:.2g}'.format(p)
            q = row.get('q_declared')
            if isinstance(q, float) and not np.isnan(q):
                text += ', q = {:.2g}'.format(q)
            ax.text(max(graph, mirror if not np.isnan(mirror) else 0) + 0.08, y, text, va='center', fontsize=8,
                    color=INK)
    ax.axvline(1, color=MUTED, lw=0.8, ls='--')
    ax.axvline(0, color=INK, lw=0.6)
    ax.set_yticks(range(len(levels)))
    ax.set_yticklabels([label for _, label in levels], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('{} change, log$_2$'.format(mark), fontsize=9)
    ax.tick_params(axis='x', labelsize=8)
    ax.set_xlim(min(-0.3, ax.get_xlim()[0]), 3.6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=GRAPH), plt.Rectangle((0, 0), 1, 1, color=MIRROR)],
              ['graph neighbours', 'mirrored loci'], fontsize=8, frameon=False, loc='lower center',
              bbox_to_anchor=(0.5, 1.0), ncol=2)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def group(nodes, graph, merge):
    groups, current = [], []
    for n in sorted(nodes, key=lambda n: middle(graph.nodes[n])):
        if current and middle(graph.nodes[n]) - middle(graph.nodes[current[-1]]) > merge:
            groups.append(current)
            current = []
        current.append(n)
    if current:
        groups.append(current)
    return groups


def panel_paths(graph, source, targets, step, mark, a, b, merge, name, path):
    columns = {0: [[source]]}
    on_paths = {k: set() for k in range(1, step)}
    for t in targets:
        for p in nx.all_shortest_paths(graph, source, t):
            for k in range(1, step):
                on_paths[k].add(p[k])
    for k in range(1, step):
        columns[k] = group(on_paths[k], graph, merge)
    columns[step] = group(targets, graph, merge)
    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
    cmap = plt.get_cmap('RdBu_r')
    pos, members = {}, {}
    for k, groups in columns.items():
        ys = np.linspace(1, 0, len(groups)) if len(groups) > 1 else [0.5]
        for i, (grp, y) in enumerate(zip(groups, ys)):
            pos[(k, i)] = (k, y)
            members[(k, i)] = grp
    node_key = {n: key for key, grp in members.items() for n in grp}
    fig, ax = plt.subplots(figsize=(1.7 * (step + 1), 3.2))
    drawn = set()
    for t in targets:
        for p in nx.all_shortest_paths(graph, source, t):
            for u, v in zip(p[:-1], p[1:]):
                ku, kv = node_key.get(u), node_key.get(v)
                if ku is None or kv is None or (ku, kv) in drawn:
                    continue
                drawn.add((ku, kv))
                status = graph.edges[u, v].get('status', 'shared')
                colour = {'gained': '#b2182b', 'lost': '#2166ac'}.get(status, '#9a9a96')
                style = '--' if status in ('gained', 'lost') else '-'
                (x0, y0), (x1, y1) = pos[ku], pos[kv]
                ax.plot([x0, x1], [y0, y1], color=colour, ls=style, lw=1.0, alpha=0.8, zorder=1)
    for key, (x, y) in pos.items():
        grp = members[key]
        change = float(np.mean([log2_change(graph.nodes[n], mark, a, b) for n in grp]))
        bait = any(int(graph.nodes[n].get('bait', 0) or 0) == 1 for n in grp)
        ax.scatter([x], [y], s=60 + 20 * len(grp) if key[0] else 160, marker='s' if bait else 'o',
                   c=[cmap(norm(change))], edgecolors='black', linewidths=0.6, zorder=3)
        if key[0] == 0:
            ax.text(x, y - 0.09, name, ha='center', va='top', fontsize=9, fontweight='bold')
        else:
            lo = min(int(graph.nodes[n]['start']) for n in grp) / 1e6
            hi = max(int(graph.nodes[n]['end']) for n in grp) / 1e6
            span = '{:.2f}-{:.2f} Mb'.format(lo, hi) if hi - lo > 0.005 else '{:.2f} Mb'.format(lo)
            ax.text(x, y - 0.07, '{}\n{} fragment{}'.format(span, len(grp), 's' if len(grp) > 1 else ''),
                    ha='center', va='top', fontsize=6.5, color='#333333')
    for k in range(step + 1):
        ax.text(k, 1.16, 'promoter' if k == 0 else 'step {}'.format(k), ha='center', va='bottom', fontsize=9)
    ax.set_xlim(-0.5, step + 0.5)
    ax.set_ylim(-0.25, 1.25)
    ax.axis('off')
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', fraction=0.04,
                       pad=0.04, shrink=0.5)
    bar.set_label('{} change {}/{}, log$_2$'.format(mark, b, a), fontsize=8)
    bar.ax.tick_params(labelsize=7)
    handles = [Line2D([0], [0], color='#9a9a96', label='shared'),
               Line2D([0], [0], color='#b2182b', ls='--', label='gained in ' + b),
               Line2D([0], [0], color='#2166ac', ls='--', label='lost in ' + b)]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.08), ncol=3, fontsize=7.5, frameon=False)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return columns


def panel_track(graph, columns, bigwig_a, bigwig_b, a, b, mark, path, bin_size=2000):
    nodes = [n for groups in columns.values() for grp in groups for n in grp]
    chrom = graph.nodes[nodes[0]]['chrom']
    start = min(int(graph.nodes[n]['start']) for n in nodes) - 30000
    end = max(int(graph.nodes[n]['end']) for n in nodes) + 30000
    edges = np.arange(start, end + bin_size, bin_size)
    fig, (ax, lane) = plt.subplots(2, 1, figsize=(6.4, 2.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    for bw_path, label, colour in ((bigwig_a, a, '#9a9a96'), (bigwig_b, b, GRAPH)):
        bw = pyBigWig.open(bw_path)
        values = [bw.stats(chrom, int(s), int(min(s + bin_size, bw.chroms()[chrom])), type='mean')[0] or 0.0
                  for s in edges[:-1]]
        bw.close()
        ax.plot((edges[:-1] + bin_size / 2) / 1e6, values, color=colour, lw=0.9, label=label)
    ax.set_ylabel('{} fold\nenrichment'.format(mark), fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7.5, frameon=False, loc='upper left')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for k, groups in columns.items():
        for grp in groups:
            for n in grp:
                s, e = int(graph.nodes[n]['start']) / 1e6, int(graph.nodes[n]['end']) / 1e6
                lane.add_patch(plt.Rectangle((s, 3 - k - 0.35), max(e - s, 0.002), 0.7, color=STEP_COLOURS[k]))
                ax.axvspan(s, e, color=STEP_COLOURS[k], alpha=0.12, lw=0)
    lane.set_ylim(-0.6, 3.6)
    lane.set_yticks([3, 2, 1, 0])
    lane.set_yticklabels(['promoter', 'step 1', 'step 2', 'step 3'], fontsize=7)
    lane.set_xlabel('{} (Mb)'.format(chrom), fontsize=8)
    lane.tick_params(axis='x', labelsize=7)
    for side in ('top', 'right'):
        lane.spines[side].set_visible(False)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('prefix')
    parser.add_argument('promoter')
    parser.add_argument('mark')
    parser.add_argument('bigwig_a')
    parser.add_argument('bigwig_b')
    parser.add_argument('output')
    parser.add_argument('--conditions', nargs=2, default=['hESC', 'hNEC'])
    parser.add_argument('--step', type=int, default=3)
    parser.add_argument('--merge', type=int, default=50000)
    parser.add_argument('--name', default=None)
    args = parser.parse_args()
    a, b = args.conditions
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    rows = levels[levels.node == args.promoter]
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    name = args.name or (graph.nodes[args.promoter].get('genes') or args.promoter).split(',')[0]
    sign = np.sign(float(rows['rna_log2FC'].iloc[0]))
    distance = nx.single_source_shortest_path_length(graph, args.promoter, cutoff=args.step)
    targets = [n for n, k in distance.items() if k == args.step
               and sign * log2_change(graph.nodes[n], args.mark, a, b) >= 1]
    panel_levels(rows, args.mark, args.output + '_levels.pdf')
    columns = panel_paths(graph, args.promoter, targets, args.step, args.mark, a, b, args.merge, name,
                          args.output + '_paths.pdf')
    panel_track(graph, columns, args.bigwig_a, args.bigwig_b, a, b, args.mark, args.output + '_track.pdf')
    print('{}: {} nodes at step {} with a concordant {} change; columns {}'.format(
        name, len(targets), args.step, args.mark, {k: [len(g) for g in v] for k, v in columns.items()}))


if __name__ == '__main__':
    main()
