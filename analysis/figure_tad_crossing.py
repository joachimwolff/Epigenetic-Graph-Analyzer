"""Locus figure of the TAD-crossing analysis (28_tad_crossing.sh): Hoxa2 and the boundary at chr6:51.175 Mb (mm10).

usage: figure_tad_crossing.py <tad_crossing results dir> <kessler2023_hic_ctcf dir> <egaDifferential prefix>
       [region chr6:50400000-52700000] [gene Hoxa2]

One PDF per panel in <results>/figures/ (no panel letters), all on the same genomic axis:
  hic_Md.pdf, hic_PA2.pdf      balanced Hi-C (25 kb, Kessler et al. 2023) as triangles, log10, one colour scale
  insulation.pdf               hicFindTADs TAD-separation score of Md and PA2 with the called boundaries
  ctcf.pdf                     CTCF ChIP-seq of Md and PA2 with the cohesin-assumed CTCF summits of egaBoundaries
                               and their motif orientation
  egaboundaries.pdf            capture insulation of each condition graph and the egaBoundaries rule calls
  arcs_Md.pdf, arcs_PA2.pdf    EGA arcs of the gene's promoter in each condition: squares baited promoters, circles
                               other ends, width = strength, gained red / lost blue dashed, shared grey; contacts
                               crossing a boundary (at or above the strength percentile of the run) drawn in full,
                               the others faded
HIRE1 / HIRE2 (Kessler et al. 2023) are shaded on every panel.
"""

import os
import sys

import h5py
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyBigWig  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import tads as td  # noqa: E402
from epigeneticgraphanalyzer.egaPlotDifferential import (INK, ONLY_DASH, STATUS_COLORS, edge_width,  # noqa: E402
                                                          width_handles)
from epigeneticgraphanalyzer.genome import parse_node_id  # noqa: E402

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'ps.fonttype': 42})
COND_COLOUR = {'Md': '#2a78d6', 'PA2': '#eb6834'}
HIRE = [('HIRE2', 50789172, 50828639), ('HIRE1', 50913170, 51087888)]
HIRE_COLOUR = '#e6e3d6'
BOUNDARY = '#5f5e58'
WIDTH = 6.4


def mb_axis(ax, start, end, label=True):
    ax.set_xlim(start, end)
    ticks = np.arange(np.ceil(start / 2e5) * 2e5, end + 1, 2e5)
    ax.set_xticks(ticks)
    ax.set_xticklabels(['{:.1f}'.format(t / 1e6) for t in ticks] if label else [])
    if label:
        ax.set_xlabel('chr6 (Mb, mm10)')


def shade_hire(ax, label=False, y=None):
    for name, s, e in HIRE:
        ax.axvspan(s, e, color=HIRE_COLOUR, zorder=0, lw=0)
        if label:
            ax.text((s + e) / 2, 0.97 if y is None else y, name, ha='center', va='top', fontsize=6.5, color=INK,
                    transform=ax.get_xaxis_transform() if y is None else ax.transData, zorder=7,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.5, alpha=0.8))


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print(path)


# ------------------------------------------------------------------------------------------------------------------
# Hi-C

def cool_region(path, chrom, start, end):
    """Balanced matrix of a region from a cooler (h5py only)."""
    with h5py.File(path, 'r') as f:
        names = [x.decode() if isinstance(x, bytes) else x for x in f['chroms/name'][:]]
        k = names.index(chrom)
        offset = f['indexes/chrom_offset'][k]
        res = int(f.attrs['bin-size'])
        lo = offset + start // res
        hi = offset + (end + res - 1) // res
        weight = f['bins/weight'][lo:hi]
        b1 = f['indexes/bin1_offset'][lo:hi + 1]
        p1 = f['pixels/bin1_id'][b1[0]:b1[-1]]
        p2 = f['pixels/bin2_id'][b1[0]:b1[-1]]
        v = f['pixels/count'][b1[0]:b1[-1]].astype(float)
    keep = (p2 >= lo) & (p2 < hi)
    n = hi - lo
    m = np.full((n, n), np.nan)
    i, j = p1[keep] - lo, p2[keep] - lo
    val = v[keep] * weight[i] * weight[j]
    m[i, j] = val
    m[j, i] = val
    return m, (lo - offset) * res, res


def hic_panel(matrix, first, res, start, end, vmin, vmax, title, depth=1.5e6):
    n = matrix.shape[0]
    edges = first + np.arange(n + 1) * res
    X, Y = np.meshgrid(edges, edges)
    x = (X + Y) / 2.0
    y = (Y - X) / 2.0
    z = np.log10(np.where(matrix > 0, matrix, np.nan))
    z = np.ma.masked_where(~np.isfinite(z) | (np.arange(n)[None, :] < np.arange(n)[:, None]), z)
    fig, ax = plt.subplots(figsize=(WIDTH, 2.3))
    mesh = ax.pcolormesh(x, y, z.T, cmap='Reds', vmin=vmin, vmax=vmax, rasterized=True, shading='flat')
    ax.set_ylim(0, depth / 2)
    ax.set_yticks([0, depth / 4, depth / 2])
    ax.set_yticklabels(['0', '{:.2g}'.format(depth / 2e6), '{:.2g}'.format(depth / 1e6)])
    ax.set_ylabel('distance (Mb)')
    mb_axis(ax, start, end)
    ax.set_title(title, loc='left', fontsize=8)
    for name, s, e in HIRE:
        ax.plot([s, e], [-0.035 * depth / 2] * 2, color='#6b6a64', lw=3, clip_on=False, solid_capstyle='butt')
        ax.text(e + 0.004 * (end - start), -0.035 * depth / 2, name, ha='left', va='center', fontsize=6,
                clip_on=False)
    ax.tick_params(axis='x', pad=11)
    bar = fig.colorbar(mesh, ax=ax, fraction=0.025, pad=0.01)
    bar.set_label('log10 balanced contacts', fontsize=7)
    bar.ax.tick_params(labelsize=6)
    return fig


# ------------------------------------------------------------------------------------------------------------------
# tracks

def bedgraph_region(path, chrom, start, end):
    t = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'value'],
                    dtype={'chrom': str})
    return t[(t.chrom == chrom) & (t.end > start) & (t.start < end)]


def insulation_panel(results, start, end):
    fig, ax = plt.subplots(figsize=(WIDTH, 1.6))
    shade_hire(ax)
    marks, total = [], 0
    for cond in ('Md', 'PA2'):
        prefix = os.path.join(results, 'hicfindtads', '{}_25kb'.format(cond))
        track = bedgraph_region(prefix + '_score.bedgraph', '6', start, end)
        ax.plot((track.start + track.end) / 2, track.value, color=COND_COLOUR[cond], lw=1.0, label=cond)
        tads = td.read_hicexplorer(prefix)
        b = tads.region('6', start, end)
        ax.scatter(b.pos, b.score, marker='v', s=34, color=COND_COLOUR[cond], zorder=3, edgecolor='none')
        marks.append(b)
        total = len(tads)
    strongest = marks[1].sort_values('strength', ascending=False).iloc[0]
    ax.annotate('{:.3f} Mb, PA2 strength rank {} of {}'.format(strongest.pos / 1e6, int(strongest.strength_rank),
                                                               total),
                (strongest.pos, strongest.score), xytext=(8, -2), textcoords='offset points', fontsize=6.5,
                va='center')
    ax.axhline(0, color='#c9c8c2', lw=0.5, zorder=0)
    ax.set_ylabel('TAD-separation\nscore')
    mb_axis(ax, start, end)
    ax.legend(handles=[plt.Line2D([0], [0], color=COND_COLOUR[c], lw=1.2, label=c) for c in ('Md', 'PA2')] +
              [plt.Line2D([0], [0], marker='v', ls='', color='#6b6a64', label='hicFindTADs boundary')],
              fontsize=6.5, frameon=False, ncol=3, loc='lower left', bbox_to_anchor=(0, 1.0))
    return fig


def ctcf_panel(raw, results, start, end):
    files = {'Md': 'GSM6505197_E10_Md_CTCF_normalized_by_autosomal_lib_1e6.scaled_0100.bw',
             'PA2': 'GSM6505196_E10_PA2_CTCF_normalized_by_autosomal_lib_1e6.scaled_0100.bw'}
    fig, axes = plt.subplots(2, 1, figsize=(WIDTH, 1.9), sharex=True)
    for ax, cond in zip(axes, ('Md', 'PA2')):
        shade_hire(ax)
        bw = pyBigWig.open(os.path.join(raw, files[cond]))
        n = (end - start) // 2000
        v = np.nan_to_num(np.array(bw.stats('chr6', start, end, type='max', nBins=n), dtype=float))
        x = start + (np.arange(n) + 0.5) * 2000
        ax.fill_between(x, 0, v, color=COND_COLOUR[cond], lw=0, step='mid')
        sites = os.path.join(results, 'egaboundaries', '{}_ega_ctcf_sites.bed'.format(cond))
        if os.path.exists(sites):
            s = pd.read_csv(sites, sep='\t', header=None, names=['chrom', 's', 'e', 'name', 'ctcf', 'strand'])
            s = s[(s.chrom == 'chr6') & (s.s > start) & (s.e < end)]
            top = v.max() * 1.05 if v.max() > 0 else 1
            for r in s.itertuples():
                marker = '>' if r.strand == '+' else '<' if r.strand == '-' else 'o'
                ax.plot([(r.s + r.e) / 2], [top], marker=marker, ms=3.5, color=INK, ls='', clip_on=False)
        ax.set_ylabel('CTCF\n{}'.format(cond))
        ax.set_ylim(0, None)
        mb_axis(ax, start, end)
    axes[0].legend(handles=[plt.Line2D([0], [0], marker='>', ls='', color=INK, ms=4, label='forward motif'),
                            plt.Line2D([0], [0], marker='<', ls='', color=INK, ms=4, label='reverse motif'),
                            plt.Line2D([0], [0], marker='o', ls='', color=INK, ms=3, label='both / none')],
                   fontsize=6.5, frameon=False, ncol=3, loc='lower left', bbox_to_anchor=(0, 1.05))
    return fig


def egaboundaries_panel(results, start, end):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(WIDTH, 2.4), sharex=True, gridspec_kw={'height_ratios': [1.6, 1]})
    rules = (('call_ctcf_pattern', 'CTCF pattern'), ('call_ctcf', 'CTCF / cohesin'), ('call_dip', 'capture dip'),
             ('call_dip_and_ctcf_pattern', 'dip AND pattern (default)'))
    shade_hire(ax)
    shade_hire(rx)
    for k, cond in enumerate(('Md', 'PA2')):
        bins = pd.read_csv(os.path.join(results, 'egaboundaries', '{}_ega_bins.tsv'.format(cond)), sep='\t')
        bins = bins[(bins.chrom == 'chr6') & (bins.end > start) & (bins.start < end)]
        covered = bins.covered.astype(str).str.lower() == 'true'
        y = np.where(covered, bins.capture_insulation.astype(float), np.nan)
        ax.plot((bins.start + bins.end) / 2, y, color=COND_COLOUR[cond], lw=1.0, marker='.', ms=2,
                label='capture insulation {}'.format(cond))
        for j, (col, label) in enumerate(rules):
            hits = bins[bins[col].astype(str).str.lower() == 'true']
            row = len(rules) - 1 - j + (0.15 if k == 0 else -0.15)
            rx.scatter((hits.start + hits.end) / 2, [row] * len(hits), marker='|', s=30, color=COND_COLOUR[cond])
    rx.set_yticks(range(len(rules)))
    rx.set_yticklabels([label for _, label in rules][::-1], fontsize=6.5)
    rx.set_ylim(-0.6, len(rules) - 0.4)
    ax.axhline(0, color='#c9c8c2', lw=0.5, zorder=0)
    ax.set_ylabel('capture insulation\n(log2 obs / exp)')
    mb_axis(rx, start, end)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc='lower left', bbox_to_anchor=(0, 1.0))
    return fig


# ------------------------------------------------------------------------------------------------------------------
# arcs

def arcs_panel(analysis, gene_nodes, cond, start, end, min_pct):
    a, b = analysis.names
    ann = analysis.annotated
    sel = (ann.node_a.isin(gene_nodes) | ann.node_b.isin(gene_nodes)) & ann['in_{}'.format(cond)]
    part = ann[sel]
    graph = analysis.graph
    fig, ax = plt.subplots(figsize=(WIDTH, 1.9))
    shade_hire(ax, label=True)
    longest = max(float((part.b_mid - part.a_mid).max()), 1.0) if len(part) else 1.0
    counts = {'cross': 0, 'within': 0}
    for status in ('shared', 'lost', 'gained'):
        for r in part[part.status == status].itertuples():
            data = graph.edges[r.node_a, r.node_b] if graph.has_edge(r.node_a, r.node_b) else {}
            strength = float(data.get('r_a' if cond == a else 'r_b', 0) or 0)
            x0, x1 = r.a_mid, r.b_mid
            counts['cross' if r.crosses else 'within'] += 1
            ax.add_patch(Arc(((x0 + x1) / 2, 0), x1 - x0, 2 * np.sqrt((x1 - x0) / longest), theta1=0, theta2=180,
                             color=STATUS_COLORS[status], lw=edge_width(strength),
                             ls='-' if status == 'shared' else ONLY_DASH, alpha=0.95 if r.crosses else 0.3,
                             zorder=3 if r.crosses else 2))
    nodes = set(part.node_a) | set(part.node_b)
    for n in nodes:
        c, s, e = parse_node_id(n)
        bait = bool(graph.nodes[n].get('bait')) if n in graph else False
        if n in gene_nodes:
            ax.plot([(s + e) / 2], [0], marker='s', ms=6, color=INK, zorder=6)
        else:
            ax.plot([(s + e) / 2], [0], marker='s' if bait else 'o', ms=3.2 if bait else 2.8,
                    mfc=INK if bait else 'white', mec=INK, mew=0.5, zorder=5, ls='')
    b_all = analysis.tads.region('chr6', start, end)
    b_all = b_all[b_all.strength_pct.fillna(1) >= min_pct]
    for r in b_all.itertuples():
        ax.axvline(r.pos, color=BOUNDARY, lw=1.3 if r.strength_pct >= 0.9 else 0.5, ls=(0, (2, 2)), zorder=1)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_ylim(-0.08, 1.15)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    mb_axis(ax, start, end)
    ax.text(0.0, 1.02, '{}: {} contacts of the promoter, {} crossing a boundary'.format(
        cond, len(part), counts['cross']), transform=ax.transAxes, fontsize=7, va='bottom')
    status = ('shared', 'lost') if cond == a else ('shared', 'gained')
    handles = [plt.Line2D([0], [0], color=STATUS_COLORS[s], lw=1.6, ls='-' if s == 'shared' else ONLY_DASH,
                          label={'shared': 'shared', 'gained': 'gained (only in {})'.format(b),
                                 'lost': 'lost (only in {})'.format(a)}[s]) for s in status]
    handles += [plt.Line2D([0], [0], color='#9a9994', lw=1.6, alpha=0.3, label='within its domain (faded)')]
    handles += width_handles()
    handles += [plt.Line2D([0], [0], marker='s', ls='', color=INK, ms=4, label='baited promoter'),
                plt.Line2D([0], [0], marker='o', ls='', mfc='white', mec=INK, ms=3.5, label='other end'),
                plt.Line2D([0], [0], color=BOUNDARY, lw=1.3, ls=(0, (2, 2)), label='TAD boundary (thick: top decile)')]
    ax.legend(handles=handles, fontsize=6, frameon=False, ncol=4, loc='upper center',
              bbox_to_anchor=(0.5, -0.32), handlelength=2.4)
    return fig


def main(results, raw, prefix, region='chr6:50400000-52700000', gene='Hoxa2'):
    chrom, span = region.split(':')
    start, end = (int(x) for x in span.split('-'))
    out = os.path.join(results, 'figures')
    os.makedirs(out, exist_ok=True)
    mats = {}
    for cond in ('Md', 'PA2'):
        cool = os.path.join(os.path.dirname(os.path.abspath(results)), 'kessler_tad_crossing',
                            '{}_25kb_25000.cool'.format(cond))
        mats[cond] = cool_region(cool, chrom.replace('chr', ''), start, end)
    values = np.concatenate([np.log10(m[0][m[0] > 0]) for m in mats.values()])
    vmin, vmax = np.nanpercentile(values, [5, 99.5])
    for cond, (m, first, res) in mats.items():
        save(hic_panel(m, first, res, start, end, vmin, vmax, 'Hi-C {} (25 kb)'.format(cond)), out,
             'hic_{}.pdf'.format(cond))
    save(insulation_panel(results, start, end), out, 'insulation.pdf')
    save(ctcf_panel(raw, results, start, end), out, 'ctcf.pdf')
    if os.path.exists(os.path.join(results, 'egaboundaries', 'Md_ega_bins.tsv')):
        save(egaboundaries_panel(results, start, end), out, 'egaboundaries.pdf')
    min_pct = 0.5
    analysis = td.Analysis(prefix, os.path.join(results, 'hicfindtads', 'PA2_25kb'),
                           regions=os.path.join(results, 'kessler2023_HIRE_mm10.bed'), min_strength_pct=min_pct)
    nodes = analysis.gene_nodes(gene)
    for cond in analysis.names:
        save(arcs_panel(analysis, nodes, cond, start, end, min_pct), out, 'arcs_{}.pdf'.format(cond))


if __name__ == '__main__':
    main(*sys.argv[1:])
