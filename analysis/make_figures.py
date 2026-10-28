"""Figures and tables of the EGA manuscript from the results of analysis/run_all.sh.

Usage: python make_figures.py <results_dir> <output_dir>
"""

import glob
import json
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# categorical slots 1-4 of the validated reference palette, assigned in fixed order
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#d9d8d2'
TIERS = ['1e-2', '1e-3', '1e-4']

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.edgecolor': MUTED,
    'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.dpi': 300,
})


def summary(results, name):
    path = os.path.join(results, 'validation', name + '_summary.tsv')
    return pd.read_csv(path, sep='\t') if os.path.exists(path) else None


def pick(table, set_name, pair_class, metric, stratum='all'):
    if table is None:
        return None
    row = table[(table['set'] == set_name) & (table['pair_class'] == pair_class) & (table['metric'] == metric) &
                (table['stratum'] == stratum)]
    return row.iloc[0] if len(row) else None


def primary_metric(assay):
    return 'npmi'


def panel_label(ax, label):
    ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontweight='bold', fontsize=9, color=INK, va='bottom')


def category_value(table, set_name, pair_class, assay):
    if set_name == 'triplets':
        return pick(table, 'triplets', pair_class, 'c_uvw')
    if set_name == 'communities':
        rows = table[(table['set'] == 'communities') & (table['pair_class'] == pair_class)]
        return rows.iloc[0] if len(rows) else None
    return pick(table, set_name, pair_class, primary_metric(assay))


def figure_validation(results, output, runs):
    fig = plt.figure(figsize=(7.1, 5.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.55, wspace=0.3)
    axes = [fig.add_subplot(grid[0, :])]
    bottom = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.55, wspace=0.3, left=0.09, right=0.97,
                              top=0.92, bottom=0.09)
    axes += [fig.add_subplot(bottom[1, 0]), fig.add_subplot(bottom[1, 1])]
    ax = axes[0]
    y = np.arange(len(CATEGORIES))[::-1]
    for offset, (name, assay, label, color, marker) in zip((-0.17, 0.17), runs):
        table = summary(results, name)
        for yi, (set_name, pair_class, _) in zip(y, CATEGORIES):
            row = category_value(table, set_name, pair_class, assay) if table is not None else None
            if row is None:
                continue
            ax.scatter(row['fraction_above_control'], yi + offset, s=24, color=color, marker=marker, zorder=3,
                       edgecolor='white', lw=0.6)
            ax.text(0.985, yi + offset, '{:,}'.format(int(row['n'])), fontsize=5.5, color=MUTED, va='center',
                    ha='right', transform=ax.get_yaxis_transform())
        ax.scatter([], [], color=color, marker=marker, s=24, label=label)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([c[2] for c in CATEGORIES])
    ax.set_xlim(0.4, 1.0)
    ax.set_xlabel('P(tested > matched control)')
    ax.grid(axis='x', color=GRID, lw=0.5)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=6.5)
    panel_label(ax, 'a')

    strata = ['0-100kb', '100-250kb', '250-500kb', '500-1000kb', '1000-2000kb']
    for ax, (name, assay, label, color, marker), letter in zip(axes[1:], runs, 'bc'):
        table = summary(results, name)
        for pair_class, colour, lab in (('oe-oe', BLUE, 'indirect, other-other'),
                                        ('bait-bait', AQUA, 'indirect, promoter-promoter')):
            values = []
            for stratum in strata:
                row = pick(table, 'indirect', pair_class, primary_metric(assay), stratum) if table is not None else None
                values.append(np.nan if row is None or row['n'] < 30 else row['fraction_above_control'])
            ax.plot(range(len(strata)), values, color=colour, lw=2, marker='o', ms=4, label=lab)
        direct = []
        for stratum in strata:
            row = pick(table, 'direct', 'bait-oe', primary_metric(assay), stratum) if table is not None else None
            direct.append(np.nan if row is None or row['n'] < 30 else row['fraction_above_control'])
        ax.plot(range(len(strata)), direct, color=MUTED, lw=1.2, ls=(0, (3, 2)), marker='s', ms=3,
                label='direct, promoter-other')
        ax.axhline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
        ax.set_xticks(range(len(strata)))
        ax.set_xticklabels(['<0.1', '0.1-0.25', '0.25-0.5', '0.5-1', '1-2'])
        ax.set_xlabel('distance (Mb)')
        ax.set_ylabel('P(tested > matched control)')
        ax.set_ylim(0.4, 1.0)
        ax.grid(axis='y', color=GRID, lw=0.5)
        ax.set_title(label, loc='left')
        if letter == 'b':
            ax.legend(frameon=False, fontsize=6, loc='lower right')
        panel_label(ax, letter)
    fig.subplots_adjust(left=0.3, right=0.97, top=0.92, bottom=0.09)
    fig.savefig(os.path.join(output, 'figure2_validation.pdf'))
    plt.close(fig)


# ----------------------------------------------------------------------------------- figure 3
def figure_stringency(results, output, tiers_runs, crosscell, fdr):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.6), sharey=True)
    series = [('direct', 'bait-oe', 'direct, promoter-other', MUTED, 's'),
              ('indirect', 'oe-oe', 'indirect, other-other', BLUE, 'o'),
              ('triplets', 'all', 'two-hop triplets', YELLOW, 'D'),
              ('communities', 'null_locus', 'communities, shifted captured sets', AQUA, '^')]
    for ax, (sample, assay, title), letter in zip(axes, tiers_runs, 'ab'):
        for set_name, pair_class, label, color, marker in series:
            values = []
            for tier in TIERS:
                table = summary(results, '{}_p{}_{}'.format(sample, tier, assay))
                row = category_value(table, set_name, pair_class, assay) if table is not None else None
                values.append(np.nan if row is None else row['fraction_above_control'])
            ax.plot(range(len(TIERS)), values, color=color, marker=marker, lw=2, ms=5, label=label)
        cross_name, cross_assay = crosscell[sample]
        table = summary(results, cross_name)
        if table is not None:
            for set_name, pair_class, label, color, marker in series:
                row = category_value(table, set_name, pair_class, cross_assay)
                if row is not None:
                    ax.scatter(1.25, row['fraction_above_control'], marker=marker, s=26, facecolor='white',
                               edgecolor=color, lw=1.2, zorder=3)
            ax.text(1.25, 0.415, 'graph of the\nother cell type', ha='center', va='bottom', fontsize=5.5, color=MUTED)
        ax.axhline(0.5, color=MUTED, lw=0.8, ls=(0, (2, 2)))
        labels = []
        for tier in TIERS:
            estimate = fdr.get((sample, tier))
            labels.append('p < {}\nFDR {:.0%}'.format(tier.replace('1e-', '10$^{-') + '}$', estimate)
                          if estimate is not None else 'p < ' + tier)
        ax.set_xticks(range(len(TIERS)))
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.4, 2.4)
        ax.set_ylim(0.4, 0.9)
        ax.set_xlabel('edge threshold (estimated FDR of the edge set)')
        ax.grid(axis='y', color=GRID, lw=0.5)
        ax.set_title(title, loc='left')
        panel_label(ax, letter)
    axes[0].set_ylabel('P(tested > matched control)')
    axes[1].legend(frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(output, 'figure6_stringency.pdf'))
    plt.close(fig)


# ----------------------------------------------------------------------------------- figure 4
def figure_differential(results, output, contrasts, tier):
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.9), gridspec_kw={'width_ratios': [1, 1.2, 1]})
    names, promoter, direct, indirect, rows = [], [], [], [], []
    for sample, label in contrasts:
        path = os.path.join(results, 'differential', 'H1_vs_{}_p{}_permutation.json'.format(sample, tier))
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            result = json.load(handle)
        mirror = result['mirror_comparison']
        names.append(label)
        promoter.append(mirror['r2_promoter'])
        direct.append(mirror['gain_direct_graph'])
        indirect.append(mirror['gain_indirect_graph'])
        rows.append((label, mirror, result))
    x = np.arange(len(names))
    ax = axes[0]
    ax.bar(x, promoter, color=BLUE, width=0.6, label='promoter')
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('R$^2$ of RNA log$_2$FC')
    ax.set_title('Promoter marks', loc='left')
    ax.grid(axis='y', color=GRID, lw=0.5)
    panel_label(ax, 'a')

    ax = axes[1]
    width = 0.2
    for k, (label, mirror, _) in enumerate(rows):
        for j, (level, color) in enumerate((('direct', BLUE), ('indirect', ORANGE))):
            graph = mirror['gain_{}_graph'.format(level)]
            control = mirror['gain_{}_mirror'.format(level)]
            ax.bar(k + (j - 0.5) * 2 * width - width / 2, graph, width=width * 0.9, color=color)
            ax.bar(k + (j - 0.5) * 2 * width + width / 2, control, width=width * 0.9, facecolor='white',
                   edgecolor=color, hatch='////', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('gain in R$^2$')
    ax.set_title('Neighbours versus mirrored loci', loc='left')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.legend(handles=[mpatches.Patch(color=BLUE, label='direct, graph'),
                       mpatches.Patch(facecolor='white', edgecolor=BLUE, hatch='////', label='direct, mirrored'),
                       mpatches.Patch(color=ORANGE, label='indirect, graph'),
                       mpatches.Patch(facecolor='white', edgecolor=ORANGE, hatch='////', label='indirect, mirrored')],
              frameon=False, fontsize=5.8, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
    panel_label(ax, 'b')

    ax = axes[2]
    for j, (level, color, marker) in enumerate((('direct', BLUE, 'o'), ('indirect', ORANGE, 's'))):
        values = [m['{}_graph_minus_mirror'.format(level)] for _, m, _ in rows]
        low = [m['{}_ci_low'.format(level)] for _, m, _ in rows]
        high = [m['{}_ci_high'.format(level)] for _, m, _ in rows]
        xs = x + (j - 0.5) * 0.3
        ax.errorbar(xs, values, yerr=[np.subtract(values, low), np.subtract(high, values)], fmt=marker, color=color,
                    ms=4, lw=1.2, capsize=0, label=level)
    ax.axhline(0, color=MUTED, lw=0.8, ls=(0, (2, 2)))
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('graph $-$ mirrored gain in R$^2$')
    ax.set_title('Graph minus mirrored loci', loc='left')
    ax.legend(frameon=False, fontsize=6, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
    ax.grid(axis='y', color=GRID, lw=0.5)
    panel_label(ax, 'c')
    fig.tight_layout()
    fig.savefig(os.path.join(output, 'figure7_differential.pdf'))
    plt.close(fig)


# ----------------------------------------------------------------------------------- figure 5
def figure_example(results, output, graph_name, multiway_name, run_name, n_shifts=50, seed=0):
    """A multi-promoter community: interaction arcs and SPRITE clusters by number of community loci."""
    from epigeneticgraphanalyzer.graph import read_graph
    from epigeneticgraphanalyzer.multiway import MultiwayData
    from epigeneticgraphanalyzer.validation import map_nodes_to_bins

    table_path = os.path.join(results, 'validation', run_name + '_communities.tsv')
    if not os.path.exists(table_path):
        return None
    table = pd.read_csv(table_path, sep='\t')
    candidates = table[(table['n_baits'] >= 3) & (table['n_bins'].between(6, 25)) &
                       (table['null_locus_n'] >= 10)].sort_values('log2_enrichment_locus', ascending=False)
    if not len(candidates):
        return None
    chosen = candidates.iloc[0]
    graph = read_graph(os.path.join(results, 'graphs', graph_name + '.graphml'))
    data = MultiwayData.load(os.path.join(results, 'multiway', multiway_name + '.npz'))
    node_bins = map_nodes_to_bins(graph, data)
    members = [n for n, c in graph.nodes(data='community') if int(c) == int(chosen['community'])]
    bins = np.unique([node_bins[n] for n in members if n in node_bins])

    def hit_distribution(bin_set):
        groups, hits = np.unique(data.incidence[bin_set].indices, return_counts=True)
        return np.bincount(hits, minlength=len(bins) + 1)[:len(bins) + 1]

    observed = hit_distribution(bins)
    chrom = data.bins.chrom_names[data.bins.chrom_index_of_bins(bins[:1])[0]]
    first, n_chrom = data.bins.offsets[chrom], data.bins.n_chrom_bins(chrom)
    covered = data.coverage_raw > 0
    locus = np.zeros(len(covered), dtype=bool)
    locus[np.fromiter(node_bins.values(), dtype=np.int64)] = True
    rng = np.random.default_rng(seed)
    null = []
    while len(null) < n_shifts:
        shifted = (bins - first + int(rng.integers(1, n_chrom))) % n_chrom + first
        if covered[shifted].mean() >= 0.8 and locus[shifted].mean() >= 0.5:
            null.append(hit_distribution(shifted))
    null = np.array(null)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.4), gridspec_kw={'width_ratios': [1.5, 1]})
    ax = axes[0]
    sub = graph.subgraph(members)
    mid = {n: (graph.nodes[n]['start'] + graph.nodes[n]['end']) / 2e6 for n in members}
    span = max(mid.values()) - min(mid.values())
    for u, v in sub.edges:
        xs = np.linspace(mid[u], mid[v], 60)
        c, h = (mid[u] + mid[v]) / 2, abs(mid[u] - mid[v]) / 2
        ax.plot(xs, abs(mid[u] - mid[v]) / span * np.sqrt(np.clip(1 - ((xs - c) / max(h, 1e-9)) ** 2, 0, 1)),
                color=BLUE, lw=0.8)
    for n in members:
        bait = graph.nodes[n].get('bait')
        ax.scatter(mid[n], 0, marker='s' if bait else 'o', s=18 if bait else 10, zorder=3,
                   color=INK if bait else 'white', edgecolor=INK, lw=0.6)
    genes = sorted({g for n in members for g in str(graph.nodes[n].get('genes', '')).split(',') if g})
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('{} (Mb)'.format(chrom))
    ax.set_title('Community of {} promoters, {} other ends'.format(
        sum(1 for n in members if graph.nodes[n].get('bait')), sum(1 for n in members if not graph.nodes[n].get('bait'))),
        loc='left')
    ax.text(0.0, -0.32, ', '.join(genes[:12]) + (' ...' if len(genes) > 12 else ''), transform=ax.transAxes,
            fontsize=5.8, color=MUTED, style='italic')
    panel_label(ax, 'a')

    ax = axes[1]
    k = np.arange(2, len(bins) + 1)
    ax.bar(k - 0.2, observed[2:], width=0.4, color=BLUE, label='community')
    ax.bar(k + 0.2, null[:, 2:].mean(axis=0), width=0.4, facecolor='white', edgecolor=MUTED, hatch='////', lw=0.6,
           label='shifted captured sets')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1)
    ax.set_xlabel('community bins in one SPRITE cluster')
    ax.set_ylabel('clusters')
    ax.set_xlim(1.4, min(len(bins), 10) + 0.6)
    ax.legend(frameon=False, fontsize=6)
    ax.grid(axis='y', color=GRID, lw=0.5)
    panel_label(ax, 'b')
    fig.tight_layout()
    fig.savefig(os.path.join(output, 'figure3_example.pdf'))
    plt.close(fig)
    return {'community': int(chosen['community']), 'chrom': chrom, 'n_bins': int(len(bins)), 'genes': genes,
            'observed_ge3': int(observed[3:].sum()), 'null_ge3_mean': float(null[:, 3:].sum(axis=1).mean()),
            'start': min(graph.nodes[n]['start'] for n in members), 'end': max(graph.nodes[n]['end'] for n in members)}


# ----------------------------------------------------------------------------------- figure 6
def figure_kessler(results, output, gene='Hoxa3', window=1500000):
    """Md versus PA2: explanation levels of differential genes and the indirect explanation of one gene."""
    import networkx as nx

    from epigeneticgraphanalyzer.graph import read_graph

    folder = os.path.join(results, 'kessler')
    classes = pd.read_csv(os.path.join(folder, 'Md_vs_PA2_gene_classes.tsv'), sep='\t')
    report = pd.read_csv(os.path.join(folder, 'Md_vs_PA2_indirect_only_genes.tsv'), sep='\t')
    graph = read_graph(os.path.join(folder, 'Md_vs_PA2_union.graphml'))

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), gridspec_kw={'width_ratios': [0.7, 2]})
    ax = axes[0]
    order = ['promoter', 'direct', 'indirect', 'unexplained']
    counts = [int((classes['explanation'] == level).sum()) for level in order]
    colors = [MUTED, BLUE, ORANGE, GRID]
    ax.bar(range(4), counts, color=colors, width=0.65)
    beyond = int(report['beyond_mirror'].sum())
    ax.bar(2, beyond, color='white', edgecolor=ORANGE, hatch='////', width=0.65, lw=0.8)
    for k, value in enumerate(counts):
        ax.text(k, value + 3, str(value), ha='center', va='bottom', fontsize=6.5, color=INK)
    ax.set_xticks(range(4))
    ax.set_xticklabels(order, rotation=30, ha='right')
    ax.set_ylabel('differential genes')
    ax.grid(axis='y', color=GRID, lw=0.5)
    panel_label(ax, 'a')

    ax = axes[1]
    row = report[report['genes'].astype(str).str.split(',').apply(lambda names: gene in names)]
    source = row.iloc[0]['node']
    first = set(graph.neighbors(source))
    second = set()
    for n in first:
        second |= set(graph.neighbors(n))
    second -= first | {source}
    drivers = set()
    if len(row):
        for entry in str(row.iloc[0]['drivers']).split('; '):
            drivers.add(entry.split(':', 1)[1].split('[')[0])
    centre = (graph.nodes[source]['start'] + graph.nodes[source]['end']) / 2
    shown = {n for n in {source} | first | {d for d in drivers if d in graph}
             if abs((graph.nodes[n]['start'] + graph.nodes[n]['end']) / 2 - centre) <= window}
    mediators = {m for d in drivers if d in graph for m in set(graph.neighbors(d)) & first}
    mid = {n: (graph.nodes[n]['start'] + graph.nodes[n]['end']) / 2e6 for n in shown}
    span = max(mid.values()) - min(mid.values())
    status_color = {'gained': ORANGE, 'lost': BLUE, 'shared': MUTED}
    for u, v, data in graph.subgraph(shown).edges(data=True):
        x1, x2 = sorted((mid[u], mid[v]))
        xs = np.linspace(x1, x2, 60)
        c, h = (x1 + x2) / 2, max((x2 - x1) / 2, 1e-9)
        on_path = (u in drivers or v in drivers) and (u in mediators or v in mediators)
        sign = -1 if on_path else 1
        ax.plot(xs, sign * (x2 - x1) / span * np.sqrt(np.clip(1 - ((xs - c) / h) ** 2, 0, 1)),
                color=status_color.get(data.get('status'), MUTED), lw=1.2 if (source in (u, v) or on_path) else 0.6)
    change = {n: np.log2((graph.nodes[n].get('H3K27ac_PA2', 0) + 1) / (graph.nodes[n].get('H3K27ac_Md', 0) + 1))
              for n in shown}
    limit = max(1.0, max(abs(v) for v in change.values()))
    scatter = None
    for n in shown:
        bait = graph.nodes[n].get('bait')
        scatter = ax.scatter(mid[n], 0, marker='s' if bait else 'o', s=36 if bait else 14, c=[change[n]],
                             cmap='RdBu_r', vmin=-limit, vmax=limit, edgecolor=INK, lw=0.6, zorder=3)
    driver_names = sorted({g for n in drivers & shown for g in str(graph.nodes[n].get('genes', '')).split(',')
                           if g and not g.startswith(('Gm', 'LOC', 'Hoxaas', 'Hotair'))})
    ax.annotate('{} (gene)\n{} (explaining promoters)'.format(gene, ', '.join(driver_names)),
                xy=(mid[source], 0.02), xytext=(mid[source] - 0.12 * span, 0.75), fontsize=6, style='italic',
                ha='right', color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.5))
    fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.01, label='H3K27ac log$_2$(PA2/Md)')
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('{} (Mb)'.format(graph.nodes[source]['chrom']))
    ax.legend(handles=[plt.Line2D([], [], color=ORANGE, lw=1.2, label='gained in PA2'),
                       plt.Line2D([], [], color=BLUE, lw=1.2, label='lost in PA2'),
                       plt.Line2D([], [], color=MUTED, lw=1.2, label='shared')],
              frameon=False, fontsize=6, loc='upper left')
    ax.set_title('{} promoter (RNA log$_2$FC {:+.1f})'.format(gene, row.iloc[0]['rna_log2FC']), loc='left',
                 fontsize=7.5)
    panel_label(ax, 'b')
    fig.tight_layout()
    fig.savefig(os.path.join(output, 'figure8_md_pa2.pdf'))
    plt.close(fig)


# ----------------------------------------------------------------------------------- tables
def call_table(results, samples):
    rows = []
    for sample, label in samples:
        path = os.path.join(results, 'called', sample + '.summary.json')
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            summary_json = json.load(handle)
        row = {'sample': label, 'tested_pairs': summary_json['n_tested'], 'nb_size': summary_json['nb_size']}
        for entry in summary_json.get('fdr_curve', []):
            tier = '{:.0e}'.format(entry['pvalue_threshold']).replace('e-0', 'e-')
            row['calls_p' + tier] = entry['n_calls']
            row['fdr_p' + tier] = entry['estimated_fdr']
        for tier in TIERS:
            nodes = os.path.join(results, 'graphs', '{}_p{}_nodes.tsv'.format(sample, tier))
            communities = os.path.join(results, 'graphs', '{}_p{}_communities.tsv'.format(sample, tier))
            if os.path.exists(nodes):
                n = pd.read_csv(nodes, sep='\t')
                c = pd.read_csv(communities, sep='\t')
                row['nodes_p' + tier] = len(n)
                row['communities_p' + tier] = len(c)
                row['multi_promoter_communities_p' + tier] = int((c['n_baits'] > 1).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def validation_table(results):
    rows = []
    for path in sorted(glob.glob(os.path.join(results, 'validation', '*_summary.tsv'))):
        name = os.path.basename(path).replace('_summary.tsv', '')
        table = pd.read_csv(path, sep='\t')
        table.insert(0, 'run', name)
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def differential_table(results, contrasts):
    rows = []
    for sample, label in contrasts:
        for tier in ('1e-2', '1e-3'):
            path = os.path.join(results, 'differential', 'H1_vs_{}_p{}_permutation.json'.format(sample, tier))
            if not os.path.exists(path):
                continue
            with open(path) as handle:
                result = json.load(handle)
            mirror = result['mirror_comparison']
            row = {'contrast': 'H1-hESC vs ' + label, 'edge_threshold': '$P<10^{{-{}}}$'.format(tier.split('-')[1]), 'genes': mirror['n_genes'],
                   'differential_genes': result['n_differential_genes']}
            for key in ('r2_promoter', 'gain_direct_graph', 'gain_direct_mirror', 'direct_graph_minus_mirror',
                        'direct_ci_low', 'direct_ci_high', 'gain_indirect_graph', 'gain_indirect_mirror',
                        'indirect_graph_minus_mirror', 'indirect_ci_low', 'indirect_ci_high'):
                row[key] = mirror[key]
            row['permutation_p_indirect'] = result['empirical_p_delta_r2']
            row['fraction_indirect'] = result['observed_fraction_indirect']
            row['fraction_indirect_null'] = result['null_fraction_indirect_mean']
            row['permutation_p_fraction'] = result['empirical_p_fraction']
            for level in ('promoter', 'direct', 'indirect', 'unexplained'):
                row['class_' + level] = result['class_counts'].get(level, 0)
            models = pd.read_csv(path.replace('_permutation.json', '_nested_models.tsv'), sep='\t')
            row['f_pvalue_direct'] = models['f_pvalue'].iloc[1]
            row['f_pvalue_indirect'] = models['f_pvalue'].iloc[2]
            rows.append(row)
    return pd.DataFrame(rows)


def latex_number(value, digits=3):
    if isinstance(value, str):
        return value if '$' in value else value.replace('_', '\\_')
    if pd.isna(value):
        return '--'
    if isinstance(value, (int, np.integer)) or float(value).is_integer() and abs(value) >= 1:
        return '{:,}'.format(int(value))
    if abs(value) < 1e-3 and value != 0:
        mantissa, exponent = '{:.1e}'.format(value).split('e')
        return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))
    return '{:.{}f}'.format(value, digits)


def write_latex_table(df, path, columns, headers, caption, label, size='\\scriptsize', longtable=False):
    spec = '@{}' + 'l' * 2 + 'r' * (len(columns) - 2) + '@{}'
    header = [' & '.join(headers) + ' \\\\', '\\midrule']
    def cell(row, column):
        value = row[column]
        # p-values below the double precision range are reported as an upper bound
        if 'p' in column.split('_') and not isinstance(value, str) and value == 0:
            return '$<10^{-300}$'
        return latex_number(value)

    body = [' & '.join(cell(row, c) for c in columns) + ' \\\\' for _, row in df.iterrows()]
    if longtable:
        lines = ['{' + size, '\\begin{longtable}{' + spec + '}', '\\toprule'] + header + ['\\endhead'] + body + \
            ['\\bottomrule', '\\caption{' + caption + '}\\label{' + label + '}', '\\end{longtable}}']
    else:
        lines = ['\\begin{table}[!htb]', '\\centering', size, '\\begin{tabular}{' + spec + '}', '\\toprule'] + \
            header + body + ['\\bottomrule', '\\end{tabular}', '\\caption{' + caption + '}',
                             '\\label{' + label + '}', '\\end{table}']
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def fdr_lookup(results):
    lookup = {}
    for path in glob.glob(os.path.join(results, 'called', '*.summary.json')):
        sample = os.path.basename(path).replace('.summary.json', '')
        with open(path) as handle:
            for entry in json.load(handle).get('fdr_curve', []):
                tier = '{:.0e}'.format(entry['pvalue_threshold']).replace('e-0', 'e-')
                lookup[(sample, tier)] = entry['estimated_fdr']
    return lookup


def main():
    results, output = sys.argv[1], sys.argv[2]
    os.makedirs(output, exist_ok=True)
    figure_validation(results, output, [
        ('GM_mifsud_c20_sprite', 'sprite', 'GM12878, SPRITE', BLUE, 'o'),
        ('H1_p1e-3_gam', 'gam', 'H1-hESC, GAM', ORANGE, 's')])
    figure_stringency(results, output, [('GM', 'sprite', 'GM12878 SPRITE'), ('H1', 'gam', 'H1-hESC GAM')],
                      {'GM': ('H1_p1e-3_sprite_crosscell', 'sprite'), 'H1': ('GM_p1e-3_gam_crosscell', 'gam')},
                      fdr_lookup(results))
    contrasts = [('ME', 'ME'), ('TB', 'TB'), ('msc', 'MSC'), ('npc', 'NPC')]
    for tier in ('1e-2', '1e-3'):
        figure_differential(results, output, contrasts, tier)
        os.replace(os.path.join(output, 'figure7_differential.pdf'),
                   os.path.join(output, 'figure7_differential_p{}.pdf'.format(tier)))
    if os.path.exists(os.path.join(results, 'kessler', 'Md_vs_PA2_union.graphml')):
        figure_kessler(results, output)
    samples = [('H1', 'H1-hESC'), ('GM', 'GM12878'), ('ME', 'mesendoderm'), ('TB', 'trophoblast'),
               ('msc', 'mesenchymal stem cells'), ('npc', 'neural progenitors')]
    calls = call_table(results, samples)
    calls.to_csv(os.path.join(output, 'table_calls.tsv'), sep='\t', index=False)
    validation = validation_table(results)
    validation.to_csv(os.path.join(output, 'table_validation.tsv'), sep='\t', index=False)
    differential = differential_table(results, contrasts)
    differential.to_csv(os.path.join(output, 'table_differential.tsv'), sep='\t', index=False)

    long_calls = []
    for _, row in calls.iterrows():
        for tier in TIERS:
            long_calls.append({'sample': row['sample'], 'threshold': '$P<10^{{-{}}}$'.format(tier.split('-')[1]), 'calls': row['calls_p' + tier],
                               'fdr': row['fdr_p' + tier], 'nodes': row.get('nodes_p' + tier),
                               'communities': row.get('communities_p' + tier),
                               'multi': row.get('multi_promoter_communities_p' + tier)})
    write_latex_table(pd.DataFrame(long_calls), os.path.join(output, 'table_s1.tex'),
                      ['sample', 'threshold', 'calls', 'fdr', 'nodes', 'communities', 'multi'],
                      ['Sample', 'Edge threshold', 'Interactions', 'Estimated FDR', 'Nodes', 'Communities',
                       'Multi-promoter communities'],
                      'Interaction calls and graphs of the samples of Jung et al. (2019) at three nominal $P$-value '
                      'thresholds (at least six reads, 15\\,kb to 2\\,Mb). The estimated FDR is the expected number '
                      'of background pairs passing the thresholds divided by the number of calls.', 'tab:calls')
    keep = validation[(validation['stratum'] == 'all') & ~validation['metric'].isin(['npmi_delta', 'cooc'])].copy()
    keep['test'] = keep['set'] + ', ' + keep['pair_class']
    write_latex_table(keep, os.path.join(output, 'table_s2.tex'),
                      ['run', 'test', 'metric', 'n', 'test_mean', 'control_mean', 'fraction_above_control', 'wilcoxon_p'],
                      ['Run', 'Test', 'Metric', '$n$', 'Mean, test', 'Mean, control', 'P(test > control)',
                       '$P$ (Wilcoxon)'],
                      'Validation of all graphs against SPRITE (on restriction fragments) and GAM (30\\\\,kb windows). '
                      'Runs are named graph, edge threshold and multi-way assay; cross-cell runs test a graph against '
                      'the data of the other cell type. P(test > control): probability that a tested item exceeds '
                      'one of its matched controls, ties counting one half. '
                      'cooc\\_oe: $\\log_2$ ratio of co-occurrence to the distance expectation; npmi: normalised '
                      'pointwise mutual information; c\\_uvw: groups containing all three loci of a triplet; '
                      'groups\\_with\\_3\\_members: groups containing at least three bins of a community.',
                      'tab:validation', size='\\tiny', longtable=True)
    write_latex_table(differential, os.path.join(output, 'table_s3.tex'),
                      ['contrast', 'edge_threshold', 'genes', 'r2_promoter', 'gain_direct_graph', 'gain_direct_mirror',
                       'direct_ci_low', 'direct_ci_high', 'gain_indirect_graph', 'gain_indirect_mirror',
                       'indirect_ci_low', 'indirect_ci_high', 'permutation_p_indirect'],
                      ['Contrast', 'Edges', 'Genes', '$R^2$ promoter', 'Direct, graph', 'Direct, mirrored',
                       'CI low', 'CI high', 'Indirect, graph', 'Indirect, mirrored', 'CI low', 'CI high',
                       '$P$ permutation'],
                      'Differential mode. Gains in explained variance of RNA $\\log_2$ fold changes by the direct and '
                      'the indirect level for graph neighbours and mirrored loci; CI: 95\\% bootstrap interval of '
                      'graph minus mirrored gain; $P$ permutation: indirect level against permuted node signals.',
                      'tab:differential', size='\\tiny')
    write_latex_table(differential, os.path.join(output, 'table_s4.tex'),
                      ['contrast', 'edge_threshold', 'class_promoter', 'class_direct', 'class_indirect',
                       'class_unexplained', 'fraction_indirect', 'fraction_indirect_null', 'permutation_p_fraction'],
                      ['Contrast', 'Edges', 'Promoter', 'Direct', 'Indirect', 'Unexplained',
                       'Indirect (informative genes)', 'Indirect, permuted', '$P$ permutation'],
                      'Hierarchical explanation of genes with an at least two-fold expression change by a concordant '
                      'two-fold change of a chromatin mark at the promoter, at direct or at indirect neighbours. '
                      'The last three columns refer to genes with direct and indirect neighbours and compare the '
                      'fraction explained at the indirect level with permuted node signals.',
                      'tab:classes', size='\\scriptsize')


if __name__ == '__main__':
    main()
