"""egaPlotDifferential: where the differential graph changes, as arcs along the genome and as a graph.

Reads the union graph of egaDifferential (<prefix>_union.graphml). The upper panel draws the interactions of a region
or of the neighbourhood of a gene as arcs along the chromosome, gained interactions above the axis, lost
interactions below and shared interactions in grey, with line widths by edge strength. The lower panel draws the same
nodes as a promoter-centred graph: edges coloured by status, other ends coloured by the change of a chromatin mark,
promoters by the expression change of their genes, and candidates for regulation through indirect contacts
(egaDifferential <prefix>_indirect_only_genes.tsv) outlined. With --sideBySide the graph is drawn for condition A,
condition B and their union on one layout.
"""

import logging

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

from ._cli import base_parser, ensure_parent, setup_logging  # noqa: E402
from .egaPlotRegion import parse_region, region_subgraph  # noqa: E402
from .graph import read_graph  # noqa: E402
from .egaPlotRegion import ordered_subgraph  # noqa: E402

log = logging.getLogger(__name__)

STATUS_COLORS = {'gained': '#d95f02', 'lost': '#1b9e77', 'shared': '#b9b8b2'}
INK = '#1f1f1e'
IGNORED_GENE_PREFIXES = ('Gm', 'LOC', 'RP1', 'RP3', 'RP4', 'RP5', 'RP6', 'AC0', 'AL0', 'CTD', 'CTC', 'CTA')


def parse_arguments(args=None):
    parser = base_parser('Plot the differential interaction graph of egaDifferential for a region or the neighbourhood '
                         'of a gene: gained, lost and shared interactions as arcs and as a promoter-centred graph '
                         'with chromatin and expression changes.')
    parser.add_argument('--graph', '-g', required=True, help='<prefix>_union.graphml of egaDifferential')
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--region', help='chrom:start-end')
    selection.add_argument('--gene', help='Gene name; draws its neighbourhood up to --steps')
    parser.add_argument('--steps', type=int, default=2, help='Neighbourhood size around --gene (graph steps)')
    parser.add_argument('--names', nargs=2, default=None, help='Condition names A B (default: from the graph)')
    parser.add_argument('--mark', help='Chromatin mark whose log2 change colours the nodes (default: the first one)')
    parser.add_argument('--candidates', help='<prefix>_indirect_only_genes.tsv of egaDifferential; outlines the '
                                             'candidates and the nodes that explain them')
    parser.add_argument('--sideBySide', action='store_true', help='Also draw condition A and condition B')
    parser.add_argument('--maxLog2', type=float, default=2.0, help='Colour scale limit of log2 changes')
    parser.add_argument('--maxDistance', type=int, default=2000000,
                        help='Interactions longer than this (bp) are not drawn (default: 2 Mb, as in the graphs)')
    parser.add_argument('--title', default=None)
    parser.add_argument('--outFileName', '-o', required=True,
                        help='Figure file; the format follows the extension: pdf, svg, eps, ps (vector) or png, jpg, '
                             'tif (raster)')
    parser.add_argument('--dpi', type=int, default=300, help='Resolution of raster formats (default: 300)')
    return parser.parse_args(args)


def condition_names(graph):
    """Condition names in the order egaDifferential wrote them (mark_A before mark_B on every node)."""
    names = []
    for _, data in graph.nodes(data=True):
        for key in data:
            if '_' in key and key != 'gene_ids' and not key.startswith('rna'):
                name = key.rsplit('_', 1)[1]
                if name not in names:
                    names.append(name)
        if len(names) >= 2:
            break
    return names


def marks_of(graph, names):
    marks = set()
    for _, data in graph.nodes(data=True):
        for key in data:
            for name in names:
                if key.endswith('_' + name):
                    marks.add(key[:-len(name) - 1])
        if marks:
            break
    return sorted(marks)


def neighbourhood(graph, gene, steps):
    seeds = [n for n, d in graph.nodes(data=True) if gene in str(d.get('genes', '')).split(',')]
    if not seeds:
        raise ValueError('Gene {} not found in the graph'.format(gene))
    nodes = set()
    for seed in seeds:
        nodes |= set(nx.single_source_shortest_path_length(graph, seed, cutoff=steps))
    return ordered_subgraph(graph, nodes), seeds


def mark_change(data, mark, names, pseudo=1.0):
    a, b = data.get('{}_{}'.format(mark, names[0])), data.get('{}_{}'.format(mark, names[1]))
    if a is None or b is None:
        return np.nan
    return float(np.log2((float(b) + pseudo) / (float(a) + pseudo)))


def gene_label(data):
    names = [g for g in str(data.get('genes', '')).split(',') if g and not g.startswith(IGNORED_GENE_PREFIXES)]
    return names[0] if names else ''


def candidate_nodes(path, nodes):
    """Candidate promoters and the two-step nodes named as their drivers."""
    if not path:
        return set(), set()
    table = pd.read_csv(path, sep='\t')
    candidates = set(table['node']) & set(nodes)
    drivers = set()
    for text in table.loc[table['node'].isin(candidates), 'drivers'].dropna():
        for item in str(text).split('; '):
            parts = item.split(':')
            if len(parts) >= 3:
                drivers.add('{}:{}'.format(parts[1], parts[2].split('[')[0]))
    return candidates, drivers & set(nodes)


def spread_layout(sub, aspect=1.6, seed=0, iterations=400):
    """Layout with little overlap: ForceAtlas2 per connected part (node sizes respected, hubs dissuaded from the
    periphery), parts placed in genomic order on a grid, then remaining node overlaps pushed apart."""
    components = sorted(nx.connected_components(sub), key=lambda c: min(int(sub.nodes[n]['start']) for n in c))
    columns = int(np.ceil(np.sqrt(len(components) * aspect))) or 1
    largest = max(len(c) for c in components)
    positions = {}
    for k, nodes in enumerate(components):
        part = ordered_subgraph(sub, nodes)
        if len(nodes) <= 2:
            layout = {n: np.array([i - 0.5, 0.0]) for i, n in enumerate(nodes)}
        else:
            size = {n: 3.0 if part.nodes[n].get('bait') else 1.0 for n in part}
            layout = nx.forceatlas2_layout(part, max_iter=iterations, node_size=size, distributed_action=True,
                                           scaling_ratio=4.0, gravity=1.0, seed=seed)
        xy = np.array([layout[n] for n in nodes], dtype=float)
        xy -= xy.mean(axis=0)
        scale = np.abs(xy).max() or 1.0
        extent = 0.45 * np.sqrt(len(nodes) / largest) + 0.12
        row, column = divmod(k, columns)
        for n, p in zip(nodes, xy):
            positions[n] = p / scale * extent + np.array([column, -row])
    return separate(positions, sub)


def separate(positions, sub, min_distance=0.035, rounds=200):
    """Push nodes closer than min_distance apart (baits count twice as large)."""
    nodes = list(positions)
    xy = np.array([positions[n] for n in nodes], dtype=float)
    radius = np.array([2.0 if sub.nodes[n].get('bait') else 1.0 for n in nodes]) * min_distance / 2
    for _ in range(rounds):
        delta = xy[:, None, :] - xy[None, :, :]
        distance = np.sqrt((delta ** 2).sum(axis=2)) + 1e-9
        need = radius[:, None] + radius[None, :]
        overlap = np.clip(need - distance, 0, None)
        np.fill_diagonal(overlap, 0)
        if overlap.max() < 1e-4:
            break
        xy += (delta / distance[:, :, None] * overlap[:, :, None] / 2).sum(axis=1)
    return {n: xy[i] for i, n in enumerate(nodes)}


def place_labels(ax, labelled, positions, fontsize=6, obstacles=None, marker_points=9.0, **_):
    """Gene labels next to their nodes without overlapping each other or any node marker.

    Labels are placed one after another (top to bottom); each takes the first free position from a list of
    candidates around its node (above, below, right, left, the diagonals, then the same at twice the distance). A position is free if the label box overlaps neither a placed label nor a node marker (`obstacles`,
    default all of `positions`, squares of `marker_points` points). Labels placed in the second ring get a leader line;
    if no candidate is free, the label takes the least overlapping position next to its node."""
    if not labelled:
        return
    fig = ax.figure
    fig.canvas.draw()
    inverse = ax.transData.inverted()
    renderer = fig.canvas.get_renderer()
    origin = inverse.transform([[0, 0]])[0]
    pixels = marker_points * fig.dpi / 72.0
    half = np.abs(inverse.transform([[pixels, pixels]])[0] - origin) / 2
    blocks = np.array(list((obstacles if obstacles is not None else positions).values()), dtype=float)
    texts, sizes = [], []
    for n, label in labelled:
        t = ax.text(0, 0, label, fontsize=fontsize, style='italic', ha='center', va='center', color=INK, zorder=5,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=0.15))
        box = t.get_window_extent(renderer)
        (x0, y0), (x1, y1) = inverse.transform([[box.x0, box.y0], [box.x1, box.y1]])
        texts.append(t)
        sizes.append((abs(x1 - x0), abs(y1 - y0)))
    sizes = np.array(sizes)
    placed = []   # (cx, cy, w, h)

    def overlap(cx, cy, w, h):
        total = 0.0
        for px, py, pw, ph in placed:
            ox = (w + pw) / 2 - abs(cx - px)
            oy = (h + ph) / 2 - abs(cy - py)
            if ox > 0 and oy > 0:
                total += ox * oy
        ox = w / 2 + half[0] - np.abs(cx - blocks[:, 0])
        oy = h / 2 + half[1] - np.abs(cy - blocks[:, 1])
        hit = (ox > 0) & (oy > 0)
        return total + float((ox[hit] * oy[hit]).sum())

    order = sorted(range(len(labelled)), key=lambda i: -positions[labelled[i][0]][1])
    for i in order:
        x, y = positions[labelled[i][0]]
        w, h = sizes[i]
        candidates = []
        for ring in (1, 2):
            dx = (w / 2 + half[0]) * 1.05 + (ring - 1) * w * 0.6
            dy = (h / 2 + half[1]) * 1.1 + (ring - 1) * h
            candidates += [(ring, x, y + dy), (ring, x, y - dy), (ring, x + dx, y), (ring, x - dx, y),
                           (ring, x + dx, y + dy), (ring, x - dx, y + dy), (ring, x + dx, y - dy),
                           (ring, x - dx, y - dy)]
        # the nearest free position; without one, the least overlapping position of the first ring, so that a label
        # never moves far from its node
        scored = [(overlap(cx, cy, w, h), k, ring, cx, cy) for k, (ring, cx, cy) in enumerate(candidates)]
        free = [c for c in scored if c[0] == 0]
        _, _, ring, cx, cy = free[0] if free else min(c for c in scored if c[2] == 1)
        texts[i].set_position((cx, cy))
        placed.append((cx, cy, w, h))
        if ring > 1:
            ax.plot([x, cx], [y, cy], color='#8c8c88', lw=0.4, zorder=4)


def draw_arcs(ax, sub, names):
    chrom = {d['chrom'] for _, d in sub.nodes(data=True)}
    middle = {n: (int(d['start']) + int(d['end'])) / 2 for n, d in sub.nodes(data=True)}
    span = (min(middle.values()), max(middle.values()))
    longest = max((abs(middle[u] - middle[v]) for u, v in sub.edges), default=1.0)
    for status, sign, order in (('shared', 1, 1), ('gained', 1, 2), ('lost', -1, 2)):
        for u, v, data in sub.edges(data=True):
            if data.get('status') != status:
                continue
            x0, x1 = sorted((middle[u], middle[v]))
            strength = max(float(data.get('r_a', 0)), float(data.get('r_b', 0)))
            height = 2 * np.sqrt((x1 - x0) / longest)
            ax.add_patch(Arc(((x0 + x1) / 2, 0), x1 - x0, height, theta1=0 if sign > 0 else 180,
                             theta2=180 if sign > 0 else 360, color=STATUS_COLORS[status],
                             lw=0.4 + 0.5 * min(strength, 4), alpha=0.9, zorder=order))
    last = -np.inf
    for n in sorted(sub.nodes, key=lambda n: middle[n]):
        data = sub.nodes[n]
        if data.get('bait'):
            ax.plot([middle[n]], [0], marker='s', ms=3, color=INK, zorder=3)
            label = gene_label(data)
            # labels closer than 1.5% of the window to the previous one are left out
            if label and middle[n] - last > 0.015 * (span[1] - span[0] or 1.0):
                ax.text(middle[n], -1.12, label, rotation=90, ha='center', va='top', fontsize=5.5, style='italic',
                        color=INK)
                last = middle[n]
    ax.axhline(0, color=INK, lw=0.6)
    pad = 0.03 * (span[1] - span[0] or 1.0)
    ax.set_xlim(span[0] - pad, span[1] + pad)
    ax.set_ylim(-1.6, 1.15)
    ax.set_yticks([])
    for side in ('left', 'right', 'top'):
        ax.spines[side].set_visible(False)
    ax.set_xlabel('{} (Mb)'.format(', '.join(sorted(chrom))), fontsize=7)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: '{:.2f}'.format(x / 1e6)))
    ax.tick_params(labelsize=6)
    ax.text(0.0, 0.98, 'gained in {}'.format(names[1]), transform=ax.transAxes, fontsize=6.5,
            color=STATUS_COLORS['gained'], va='top')
    ax.text(0.0, 0.02, 'lost in {}'.format(names[1]), transform=ax.transAxes, fontsize=6.5,
            color=STATUS_COLORS['lost'], va='bottom')


def draw_graph(ax, sub, positions, names, mark, max_log2, edges_of='union', candidates=(), drivers=(), title=''):
    norm = mcolors.TwoSlopeNorm(vmin=-max_log2, vcenter=0, vmax=max_log2)
    cmap = plt.get_cmap('RdBu_r')
    labelled, seen = [], set()
    for u, v, data in sub.edges(data=True):
        status = data.get('status', 'shared')
        if edges_of == 'A' and status == 'gained' or edges_of == 'B' and status == 'lost':
            continue
        color = STATUS_COLORS[status] if edges_of == 'union' else '#7f7f7a'
        strength = float(data.get('r_a' if edges_of == 'A' else 'r_b' if edges_of == 'B' else 'r_b', 0) or 0)
        if edges_of == 'union':
            strength = max(float(data.get('r_a', 0)), float(data.get('r_b', 0)))
        ax.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]], color=color,
                lw=0.4 + 0.4 * min(strength, 4), alpha=0.85, zorder=1)
    for n, data in sub.nodes(data=True):
        if n not in positions:
            continue
        bait = bool(data.get('bait'))
        value = data.get('rna_log2FC') if bait else mark_change(data, mark, names) if mark else np.nan
        value = np.nan if value is None else float(value)
        face = '#ffffff' if np.isnan(value) else cmap(norm(np.clip(value, -max_log2, max_log2)))
        outline = INK
        width = 0.5
        if n in candidates:
            outline, width = '#000000', 1.8
        elif n in drivers:
            outline, width = '#6a3d9a', 1.6
        ax.scatter(*positions[n], marker='s' if bait else 'o', s=70 if bait else 18, color=face, edgecolor=outline,
                   lw=width, zorder=3 if bait else 2)
        label = gene_label(data) if bait else ''
        if label and label not in seen:
            seen.add(label)
            labelled.append((n, label))
    ax.set_axis_off()
    xy = np.array(list(positions.values()))
    ax.set_xlim(xy[:, 0].min() - 0.12, xy[:, 0].max() + 0.12)
    ax.set_ylim(xy[:, 1].min() - 0.08, xy[:, 1].max() + 0.15)
    place_labels(ax, labelled, positions)
    ax.set_title(title, loc='left', fontsize=7.5)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    names = args.names or condition_names(graph)[:2]
    marks = marks_of(graph, names)
    mark = args.mark or (marks[0] if marks else None)
    graph = graph.copy()
    middle = {n: (int(d['start']) + int(d['end'])) // 2 for n, d in graph.nodes(data=True)}
    graph.remove_edges_from([(u, v) for u, v in graph.edges if abs(middle[u] - middle[v]) > args.maxDistance])
    if args.region:
        sub = region_subgraph(graph, parse_region(args.region))
    else:
        sub, _ = neighbourhood(graph, args.gene, args.steps)
    sub = ordered_subgraph(sub, [n for n in sub.nodes if sub.degree[n] > 0])
    if sub.number_of_edges() == 0:
        raise ValueError('No interactions in the selection')
    status = pd.Series([d.get('status') for _, _, d in sub.edges(data=True)]).value_counts().to_dict()
    log.info('%d nodes, %d edges: %s', sub.number_of_nodes(), sub.number_of_edges(), status)
    candidates, drivers = candidate_nodes(args.candidates, sub.nodes)
    positions = spread_layout(sub, aspect=1.6)

    panels = ['A', 'B', 'union'] if args.sideBySide else ['union']
    fig = plt.figure(figsize=(7.2 if not args.sideBySide else 10.5, 7.2))
    grid = fig.add_gridspec(2, len(panels), height_ratios=[1, 1.9], hspace=0.25)
    draw_arcs(fig.add_subplot(grid[0, :]), sub, names)
    titles = {'A': names[0], 'B': names[1],
              'union': '{} vs {}: gained {}, lost {}, shared {}'.format(
                  names[0], names[1], status.get('gained', 0), status.get('lost', 0), status.get('shared', 0))}
    for k, panel in enumerate(panels):
        draw_graph(fig.add_subplot(grid[1, k]), sub, positions, names, mark, args.maxLog2, edges_of=panel,
                   candidates=candidates, drivers=drivers, title=titles[panel])
    norm = mcolors.TwoSlopeNorm(vmin=-args.maxLog2, vcenter=0, vmax=args.maxLog2)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), ax=fig.axes, orientation='horizontal',
                       fraction=0.025, pad=0.02, aspect=50)
    bar.set_label('log2 change {} / {}: expression (squares){}'.format(
        names[1], names[0], ', {} (circles)'.format(mark) if mark else ''), fontsize=7)
    bar.ax.tick_params(labelsize=6)
    handles = [plt.Line2D([0], [0], color=c, lw=2, label=s) for s, c in STATUS_COLORS.items()]
    handles += [plt.Line2D([0], [0], marker='s', ls='', color='#9a9994', markeredgecolor=INK, ms=6,
                           label='baited promoter'),
                plt.Line2D([0], [0], marker='o', ls='', color='#9a9994', markeredgecolor=INK, ms=4,
                           label='other end (OE)')]
    if candidates:
        handles.append(plt.Line2D([0], [0], marker='s', color='white', markeredgecolor='#000000', markeredgewidth=1.8,
                                  ms=7, lw=0, label='candidate (indirect only)'))
        handles.append(plt.Line2D([0], [0], marker='o', color='white', markeredgecolor='#6a3d9a', markeredgewidth=1.6,
                                  ms=6, lw=0, label='explaining node'))
    fig.legend(handles=handles, loc='upper right', fontsize=6.5, frameon=False, ncol=len(handles))
    if args.title:
        fig.suptitle(args.title, x=0.02, ha='left', fontsize=9)
    fig.savefig(ensure_parent(args.outFileName), bbox_inches='tight', dpi=args.dpi)
    plt.close(fig)
    log.info('wrote %s', args.outFileName)


if __name__ == '__main__':
    main()
