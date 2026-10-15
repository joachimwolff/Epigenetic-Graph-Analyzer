"""egaPlotLocus: a genomic region as capture Hi-C graph and as multi-way (SPRITE / GAM) contact map."""

import logging

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

from ._cli import base_parser, ensure_parent, setup_logging  # noqa: E402
from .genome import LiftOver, parse_node_id  # noqa: E402
from .graph import read_graph  # noqa: E402
from .indirect import indirect_contacts  # noqa: E402
from .multiway import MultiwayData, npmi_from_counts  # noqa: E402

log = logging.getLogger(__name__)

EDGE_COLOR = '#2a78d6'
COMMUNITY_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
OTHER_COLOR = '#b9b8b2'
INDIRECT_COLOR = '#eb6834'
INK = '#1f1f1e'


def parse_arguments(args=None):
    parser = base_parser('Plot a region of the interaction graph (arcs of direct and indirect contacts) above the '
                         'multi-way contact map of the same region with the graph contacts marked on it.')
    parser.add_argument('--graph', '-g', required=True)
    parser.add_argument('--multiway', '-m', required=True, help='.npz of egaBuildMultiway')
    parser.add_argument('--region', help='chrom:start-end in the assembly of the graph (map style)')
    parser.add_argument('--chain', help='Chain file from the graph assembly to the multi-way assembly')
    parser.add_argument('--metric', choices=('auto', 'cooc_oe', 'npmi'), default='auto',
                        help="Map values: log2 observed/expected co-occurrence or NPMI ('auto': SPRITE -> cooc_oe, "
                             'GAM -> npmi)')
    parser.add_argument('--highlightCommunity', type=int, default=None, help='Community id drawn in colour')
    parser.add_argument('--style', choices=('events', 'communities', 'map'), default='events',
                        help="'events': the SPRITE clusters / GAM tubes that contain several loci of the community, "
                             "'map': pairwise co-occurrence map of the region")
    parser.add_argument('--panels', choices=('all', 'evidence'), default='all',
                        help="communities: all panels or only the multi-way events and their purity")
    parser.add_argument('--minLoci', type=int, default=3, help="events: minimum community loci per cluster / tube")
    parser.add_argument('--maxRows', type=int, default=60, help='events: clusters / tubes drawn')
    parser.add_argument('--shifts', type=int, default=50, help='events: shifted control sets')
    parser.add_argument('--title', default=None)
    parser.add_argument('--output', '-o', required=True)
    return parser.parse_args(args)


def parse_region(region):
    chrom, span = region.split(':')
    start, end = span.replace(',', '').split('-')
    return chrom, int(start), int(end)


def region_matrix(data, first_bin, n_bins, metric):
    """Pairwise metric of all bins of a region (NaN on the diagonal and for undetected bins)."""
    index = np.arange(first_bin, first_bin + n_bins)
    a, b = np.triu_indices(n_bins, k=1)
    max_offset = n_bins
    data.build_band(max_offset)
    matrix = np.full((n_bins, n_bins), np.nan)
    if metric == 'npmi':
        values = npmi_from_counts(data.cooccurrence(index[a], index[b], weighted=False),
                                  data.coverage_raw[index[a]], data.coverage_raw[index[b]], data.n_groups)
    else:
        expected = data.expected_by_offset(max_offset, 'cooc')
        cooc = data.cooccurrence(index[a], index[b])
        pseudo = np.nanmin(expected[1:][expected[1:] > 0]) * 0.1
        values = np.log2((cooc + pseudo) / (expected[b - a] + pseudo))
    covered = data.coverage_raw[index] > 0
    values[~(covered[a] & covered[b])] = np.nan
    matrix[a, b] = values
    matrix[b, a] = values
    return matrix


def node_positions(graph, nodes, liftover):
    positions = {}
    for node in nodes:
        chrom, start, end = parse_node_id(node)
        mid = (start + end) // 2
        if liftover is not None:
            lifted = liftover.position(chrom, mid)
            if lifted is None:
                continue
            chrom, mid = lifted
        positions[node] = (chrom, mid)
    return positions


def plot_locus(graph, data, region, output, liftover=None, metric='auto', community=None, title=None):
    chrom, start, end = region
    original = region
    if liftover is not None:
        lifted_start = liftover.position(chrom, start)
        lifted_end = liftover.position(chrom, end)
        if lifted_start is None or lifted_end is None:
            raise ValueError('Region borders cannot be lifted')
        chrom, start, end = lifted_start[0], lifted_start[1], lifted_end[1]
    if metric == 'auto':
        metric = 'npmi' if data.kind == 'gam' else 'cooc_oe'
    if not data.bins.fixed:
        raise ValueError('--style map needs fixed-size bins')
    bin_size = data.bins.bin_size
    first_bin = data.bins.bin_of(chrom, start)
    n_bins = (end - start) // bin_size + 1
    matrix = region_matrix(data, first_bin, n_bins, metric)
    origin = (start // bin_size) * bin_size

    margin = (original[2] - original[1]) // 2
    in_region = [n for n, d in graph.nodes(data=True) if d.get('chrom') == original[0] and
                 original[1] - margin <= int(d.get('start', -1)) <= original[2] + margin]
    positions = node_positions(graph, in_region, liftover)
    positions = {n: p for n, p in positions.items() if p[0] == chrom and start <= p[1] < end}
    sub = graph.subgraph(positions)
    pairs = indirect_contacts(sub) if sub.number_of_edges() else None

    fig = plt.figure(figsize=(7.1, 5.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], width_ratios=[1, 0.025], hspace=0.08, wspace=0.03)
    ax_arcs = fig.add_subplot(grid[0, 0])
    ax_map = fig.add_subplot(grid[1, 0], sharex=ax_arcs)
    ax_bar = fig.add_subplot(grid[1, 1])

    def mb(position):
        return position / 1e6

    span = (end - start)
    for u, v in sub.edges:
        x1, x2 = sorted((mb(positions[u][1]), mb(positions[v][1])))
        xs = np.linspace(x1, x2, 60)
        height = (x2 - x1) / mb(span)
        c, h = (x1 + x2) / 2, max((x2 - x1) / 2, 1e-9)
        highlighted = community is not None and graph.nodes[u].get('community') == community and \
            graph.nodes[v].get('community') == community
        ax_arcs.plot(xs, height * np.sqrt(np.clip(1 - ((xs - c) / h) ** 2, 0, 1)),
                     color=EDGE_COLOR if (community is None or highlighted) else '#b9b8b2', lw=0.9)
    if pairs is not None:
        for row in pairs.itertuples():
            x1, x2 = sorted((mb(positions[row.node_a][1]), mb(positions[row.node_b][1])))
            xs = np.linspace(x1, x2, 60)
            height = (x2 - x1) / mb(span)
            c, h = (x1 + x2) / 2, max((x2 - x1) / 2, 1e-9)
            ax_arcs.plot(xs, -height * np.sqrt(np.clip(1 - ((xs - c) / h) ** 2, 0, 1)), color=INDIRECT_COLOR,
                         lw=0.6, ls=(0, (3, 2)), alpha=0.8)
    for node, (_, position) in positions.items():
        bait = graph.nodes[node].get('bait')
        ax_arcs.scatter(mb(position), 0, marker='s' if bait else 'o', s=16 if bait else 7, zorder=3,
                        color=INK if bait else 'white', edgecolor=INK, lw=0.5)
    labels = {}
    for node, (_, position) in positions.items():
        genes = [g for g in str(graph.nodes[node].get('genes', '')).split(',') if g and not g.startswith(('Gm', 'LOC'))]
        if graph.nodes[node].get('bait') and genes:
            labels[node] = (position, genes[0])
    last = -np.inf
    for node, (position, name) in sorted(labels.items(), key=lambda x: x[1][0]):
        # one label per cluster of neighbouring promoters
        if position - last < span * 0.025:
            continue
        last = position
        ax_arcs.text(mb(position), 0.03, name, rotation=90, fontsize=5.5, ha='center', va='bottom', color=INK,
                     style='italic')
    ax_arcs.axhline(0, color='#d9d8d2', lw=0.6, zorder=0)
    ax_arcs.set_yticks([])
    for side in ('left', 'top', 'right', 'bottom'):
        ax_arcs.spines[side].set_visible(False)
    ax_arcs.tick_params(axis='x', labelbottom=False, bottom=False)
    ax_arcs.set_ylim(-1.05, 1.25)
    ax_arcs.set_title(title or '{}:{:,}-{:,}'.format(chrom, start, end), loc='left', fontsize=8)

    # triangle map: cell (i, j) drawn at x = midpoint of the two bins, y = their distance
    edges = np.arange(n_bins + 1)
    X, Y = np.meshgrid(edges, edges)
    x = mb(origin + (X + Y) / 2 * bin_size)
    y = mb((Y - X) / 2 * bin_size)
    masked = np.ma.masked_where(np.tril(np.ones_like(matrix, dtype=bool)) | np.isnan(matrix), matrix)
    if metric == 'npmi':
        vmin, vmax, cmap, label = 0, np.nanpercentile(matrix, 98), 'Reds', 'NPMI'
    else:
        limit = np.nanpercentile(np.abs(matrix), 98)
        vmin, vmax, cmap, label = -limit, limit, 'RdBu_r', 'log$_2$ observed/expected'
    mesh = ax_map.pcolormesh(x, y, masked.T, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    fig.colorbar(mesh, cax=ax_bar, label=label)

    def mark(nodes_a, nodes_b, **kwargs):
        xs, ys = [], []
        for u, v in zip(nodes_a, nodes_b):
            p, q = sorted((positions[u][1], positions[v][1]))
            i, j = (p - origin) // bin_size, (q - origin) // bin_size
            if i == j:
                continue
            xs.append(mb(origin + ((i + j) / 2 + 0.5) * bin_size))
            ys.append(mb((j - i) / 2 * bin_size))
        ax_map.scatter(xs, ys, **kwargs)

    edges_list = list(sub.edges)
    if edges_list:
        mark([e[0] for e in edges_list], [e[1] for e in edges_list], s=26, facecolor='none', edgecolor=EDGE_COLOR,
             lw=1.2, label='significant interaction')
    if pairs is not None and len(pairs):
        mark(pairs['node_a'], pairs['node_b'], s=24, marker='x', color=INDIRECT_COLOR, lw=1.2,
             label='indirect contact')
    ax_map.set_ylim(0, mb(span / 2))
    ax_map.invert_yaxis()
    ax_map.set_xlim(mb(start), mb(end))
    ax_map.set_xlabel('{} (Mb)'.format(chrom))
    ax_map.set_ylabel('distance (Mb)')
    ax_map.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda value, _: '{:g}'.format(2 * value)))
    ax_map.legend(loc='lower left', frameon=False, fontsize=6)
    for side in ('top', 'right'):
        ax_map.spines[side].set_visible(False)
    fig.savefig(ensure_parent(output), dpi=300, bbox_inches='tight')
    plt.close(fig)


def community_bins(graph, data, community, liftover):
    members = [n for n, c in graph.nodes(data='community') if c is not None and int(c) == community]
    positions = node_positions(graph, members, liftover)
    chroms = {p[0] for p in positions.values()}
    if len(chroms) != 1:
        raise ValueError('community {} spans {} chromosomes after liftover'.format(community, len(chroms)))
    bins = {n: data.bins.bin_of(c, p) for n, (c, p) in positions.items()}
    return members, positions, {n: b for n, b in bins.items() if b >= 0 and data.coverage_raw[b] > 0}


def groups_with_loci(data, bins, min_loci):
    """Groups (clusters / tubes) containing at least `min_loci` of `bins`: dict group -> sorted bins."""
    bins = np.unique(bins)
    sub = data.incidence[bins].tocsc()
    hits = np.diff(sub.indptr)
    result = {}
    for group in np.nonzero(hits >= min_loci)[0]:
        rows = sub.indices[sub.indptr[group]:sub.indptr[group + 1]]
        result[int(group)] = bins[np.sort(rows)]
    return result


def plot_events(graph, data, community, output, liftover=None, min_loci=3, max_rows=60, shifts=50, seed=0,
                title=None):
    """Community graph above the SPRITE clusters / GAM tubes that hold several of its loci at the same time."""
    rng = np.random.default_rng(seed)
    members, positions, node_bins = community_bins(graph, data, community, liftover)
    bins = np.unique(list(node_bins.values()))
    chrom = next(iter(positions.values()))[0]
    first = data.bins.offsets[chrom]
    n_chrom = data.bins.n_chrom_bins(chrom)

    observed = groups_with_loci(data, bins, 2)
    covered = data.coverage_raw > 0
    all_bins = np.array(sorted({b for b in _graph_bins(graph, data, liftover, chrom)}))
    locus = np.zeros(len(covered), dtype=bool)
    locus[all_bins] = True
    null_counts = []
    while len(null_counts) < shifts:
        shifted = (bins - first + int(rng.integers(1, n_chrom))) % n_chrom + first
        if covered[shifted].mean() < 0.8 or locus[shifted].mean() < 0.5:
            continue
        found = groups_with_loci(data, shifted, 2)
        null_counts.append(np.bincount([len(v) for v in found.values()], minlength=len(bins) + 1))
    null_counts = np.array([np.pad(c, (0, max(0, len(bins) + 1 - len(c))))[:len(bins) + 1] for c in null_counts])
    obs_counts = np.bincount([len(v) for v in observed.values()], minlength=len(bins) + 1)[:len(bins) + 1]

    multi = {g: v for g, v in observed.items() if len(v) >= min_loci}
    order = sorted(multi, key=lambda g: (-len(multi[g]), multi[g][0]))
    if len(order) > max_rows:
        order = order[:max_rows]
    order = sorted(order, key=lambda g: (multi[g][0], -len(multi[g])))

    def mb(b):
        b = np.atleast_1d(b)
        return (data.bins.bin_start(b) + data.bins.bin_end(b)) / 2e6

    fig = plt.figure(figsize=(7.1, 5.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.7], width_ratios=[1.6, 1], hspace=0.12, wspace=0.3)
    ax_arcs = fig.add_subplot(grid[0, 0])
    ax_rows = fig.add_subplot(grid[1, 0], sharex=ax_arcs)
    ax_bar = fig.add_subplot(grid[1, 1])

    sub = graph.subgraph(members)
    x_of = {n: mb(node_bins[n])[0] for n in node_bins}
    span = max(x_of.values()) - min(x_of.values()) or 1.0
    for u, v in sub.edges:
        if u not in x_of or v not in x_of:
            continue
        x1, x2 = sorted((x_of[u], x_of[v]))
        xs = np.linspace(x1, x2, 60)
        c, h = (x1 + x2) / 2, max((x2 - x1) / 2, 1e-9)
        ax_arcs.plot(xs, (x2 - x1) / span * np.sqrt(np.clip(1 - ((xs - c) / h) ** 2, 0, 1)), color=EDGE_COLOR, lw=1.0)
    names = []
    for node, x in x_of.items():
        bait = graph.nodes[node].get('bait')
        ax_arcs.scatter(x, 0, marker='s' if bait else 'o', s=22 if bait else 9, zorder=3, color=INK if bait else 'white',
                        edgecolor=INK, lw=0.6)
        genes = [g for g in str(graph.nodes[node].get('genes', '')).split(',') if g and not g.startswith(('Gm', 'LOC',
                                                                                                         'RP', 'AC0'))]
        if bait and genes:
            names.append((x, genes[0]))
    last = -np.inf
    for x, name in sorted(names):
        if x - last < span * 0.04:
            continue
        last = x
        ax_arcs.text(x, 0.05, name, rotation=90, fontsize=6, ha='center', va='bottom', style='italic', color=INK)
    ax_arcs.axhline(0, color='#d9d8d2', lw=0.6, zorder=0)
    ax_arcs.set_ylim(-0.08, 1.35)
    ax_arcs.set_yticks([])
    for side in ('left', 'top', 'right', 'bottom'):
        ax_arcs.spines[side].set_visible(False)
    ax_arcs.tick_params(axis='x', labelbottom=False, bottom=False)
    group_name = 'tube' if data.kind == 'gam' else 'cluster'
    ax_arcs.set_title((title or 'Community {}'.format(community)) + ': capture Hi-C graph', loc='left', fontsize=8)

    for bin_value in bins:
        ax_rows.axvline(mb(bin_value)[0], color='#ecebe6', lw=0.6, zorder=0)
    for row, group in enumerate(order):
        xs = mb(multi[group])
        ax_rows.plot([xs.min(), xs.max()], [row, row], color='#9a9993', lw=0.6, zorder=1)
        ax_rows.scatter(xs, [row] * len(xs), s=10, color=INDIRECT_COLOR, zorder=2, edgecolor='none')
    ax_rows.set_ylim(len(order) - 0.5 + 0.5, -1)
    ax_rows.set_yticks([])
    assay = 'GAM' if data.kind == 'gam' else 'SPRITE'
    ax_rows.set_ylabel('{} {}s'.format(assay, group_name))
    ax_rows.set_title('{}: {} {}s containing \u2265 {} community loci (one row each)'.format(
        assay, len(multi), group_name, min_loci), loc='left', fontsize=7.5, color=INK)
    ax_rows.set_xlabel('{} (Mb)'.format(chrom))
    pad = span * 0.03
    ax_rows.set_xlim(min(x_of.values()) - pad, max(x_of.values()) + pad)
    for side in ('top', 'right', 'left'):
        ax_rows.spines[side].set_visible(False)

    k = np.arange(2, min(len(bins), 10) + 1)
    ax_bar.bar(k - 0.2, obs_counts[k], width=0.4, color=EDGE_COLOR, label='community loci')
    ax_bar.bar(k + 0.2, null_counts[:, k].mean(axis=0), width=0.4, facecolor='white', edgecolor='#6b6a64',
               hatch='////', lw=0.6, label='shifted captured loci')
    ax_bar.set_yscale('log')
    ax_bar.set_ylim(bottom=0.8)
    ax_bar.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax_bar.set_xlabel('community loci in one {}'.format(group_name))
    ax_bar.set_ylabel('{} {}s'.format('GAM' if data.kind == 'gam' else 'SPRITE', group_name))
    ax_bar.set_title('{} {}s per number of loci'.format('GAM' if data.kind == 'gam' else 'SPRITE', group_name),
                     loc='left', fontsize=7.5, color=INK)
    ax_bar.legend(frameon=False, fontsize=6)
    for side in ('top', 'right'):
        ax_bar.spines[side].set_visible(False)
    fig.savefig(ensure_parent(output), dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'groups_ge_min': len(multi), 'null_ge_min_mean': float(null_counts[:, min_loci:].sum(axis=1).mean()),
            'n_bins': int(len(bins)), 'chrom': chrom}


def event_purity(events, label_of):
    """Mean fraction of an event's labelled loci that belong to its most frequent community."""
    values = []
    for loci in events:
        labels = [label_of[b] for b in loci if b in label_of]
        if len(labels) < 2:
            continue
        values.append(np.bincount(labels).max() / len(labels))
    return float(np.mean(values)) if values else np.nan


def plot_communities(graph, data, region, output, liftover=None, min_loci=3, max_rows=60, max_communities=8,
                     title=None, seed=0, panels='all'):
    """The communities of a region as arcs, as a network and in the multi-way events that contain their loci."""
    import networkx as nx
    chrom, start, end = region
    nodes = [n for n, d in graph.nodes(data=True) if d.get('chrom') == chrom and
             start <= (int(d['start']) + int(d['end'])) // 2 < end]
    positions = node_positions(graph, nodes, liftover)
    data_chrom = next(iter(positions.values()))[0]
    node_bins = {n: data.bins.bin_of(c, p) for n, (c, p) in positions.items() if c == data_chrom}
    node_bins = {n: b for n, b in node_bins.items() if b >= 0 and data.coverage_raw[b] > 0}
    sub = graph.subgraph(node_bins)

    # same colour assignment as egaPlotRegion: communities ranked by nodes inside the region, in genomic order
    from .egaPlotRegion import community_colors, region_subgraph
    color_of = community_colors(region_subgraph(graph, region), max_communities)
    ranked = list(color_of)
    index_of = {c: k for k, c in enumerate(ranked)}
    palette = [color_of[c] for c in ranked]

    # bin labels: the colored community with most nodes in the bin
    votes = {}
    for n, b in node_bins.items():
        c = graph.nodes[n].get('community')
        if c in index_of:
            votes.setdefault(b, []).append(index_of[c])
    label_of = {b: int(np.bincount(v).argmax()) for b, v in votes.items()}
    labelled_bins = np.array(sorted(label_of))

    events = list(groups_with_loci(data, labelled_bins, min_loci).values())
    observed = event_purity(events, label_of)
    order_labels = [label_of[b] for b in labelled_bins]
    null = []
    for shift in range(1, len(labelled_bins)):
        rotated = {b: order_labels[(k + shift) % len(labelled_bins)] for k, b in enumerate(labelled_bins)}
        null.append(event_purity(events, rotated))
    null = np.array([v for v in null if not np.isnan(v)])

    def mb(b):
        b = np.atleast_1d(b)
        return ((data.bins.bin_start(b) + data.bins.bin_end(b)) / 2)[0] / 1e6

    x_of = {n: mb(b) for n, b in node_bins.items()}
    lo, hi = min(x_of.values()), max(x_of.values())
    span = hi - lo or 1.0

    if panels == 'evidence':
        fig = plt.figure(figsize=(7.1, 3.4))
        grid = fig.add_gridspec(1, 2, width_ratios=[3.2, 1], wspace=0.25)
        ax_arcs = ax_net = None
        ax_rows = fig.add_subplot(grid[0, 0])
        ax_bar = fig.add_subplot(grid[0, 1])
    else:
        fig = plt.figure(figsize=(7.1, 8.2))
        grid = fig.add_gridspec(3, 2, height_ratios=[0.75, 1.35, 1.25], width_ratios=[1.6, 1], hspace=0.32,
                                wspace=0.3)
        ax_arcs = fig.add_subplot(grid[0, :])
        ax_rows = fig.add_subplot(grid[1, :], sharex=ax_arcs)
        ax_net = fig.add_subplot(grid[2, 0])
        ax_bar = fig.add_subplot(grid[2, 1])

    def community_color(u, v=None):
        cu = graph.nodes[u].get('community')
        if v is None:
            return color_of.get(cu, OTHER_COLOR)
        cv = graph.nodes[v].get('community')
        return color_of.get(cu, OTHER_COLOR) if cu == cv else OTHER_COLOR

    if ax_arcs is not None:
        for u, v in sub.edges:
            x1, x2 = sorted((x_of[u], x_of[v]))
            xs = np.linspace(x1, x2, 60)
            c, h = (x1 + x2) / 2, max((x2 - x1) / 2, 1e-9)
            ax_arcs.plot(xs, (x2 - x1) / span * np.sqrt(np.clip(1 - ((xs - c) / h) ** 2, 0, 1)),
                         color=community_color(u, v), lw=1.0)
        for n in sub.nodes:
            bait = graph.nodes[n].get('bait')
            ax_arcs.scatter(x_of[n], 0, marker='s' if bait else 'o', s=20 if bait else 9, color=community_color(n),
                            edgecolor=INK, lw=0.4, zorder=3)
        last = -np.inf
        for n in sorted(sub.nodes, key=lambda n: x_of[n]):
            genes = [g for g in str(graph.nodes[n].get('genes', '')).split(',')
                     if g and not g.startswith(('Gm', 'LOC', 'RP', 'AC0', 'CTD', 'CTC'))]
            if graph.nodes[n].get('bait') and genes and x_of[n] - last > span * 0.035:
                ax_arcs.text(x_of[n], 0.05, genes[0], rotation=90, fontsize=5.5, style='italic', ha='center',
                             va='bottom', color=INK)
                last = x_of[n]
        ax_arcs.axhline(0, color='#d9d8d2', lw=0.6, zorder=0)
        ax_arcs.set_ylim(-0.08, 1.4)
        ax_arcs.set_yticks([])
        for side in ('left', 'top', 'right', 'bottom'):
            ax_arcs.spines[side].set_visible(False)
        ax_arcs.tick_params(axis='x', labelbottom=False, bottom=False)
    assay = 'GAM' if data.kind == 'gam' else 'SPRITE'
    group_name = 'tube' if data.kind == 'gam' else 'cluster'
    if ax_arcs is not None:
        ax_arcs.set_title((title or '{}:{:,}-{:,}'.format(chrom, start, end)) + ': communities as arcs', loc='left',
                      fontsize=8)

    rng = np.random.default_rng(seed)
    shown = sorted(events, key=lambda e: -len(e))[:max_rows * 3]
    if len(shown) > max_rows:
        shown = [shown[i] for i in sorted(rng.choice(len(shown), max_rows, replace=False))]

    def dominant(loci):
        labels = [label_of[b] for b in loci if b in label_of]
        return int(np.bincount(labels).argmax())

    shown = sorted(shown, key=lambda e: (dominant(e), e[0]))
    for row, loci in enumerate(shown):
        xs = [mb(b) for b in loci]
        ax_rows.plot([min(xs), max(xs)], [row, row], color='#c9c8c2', lw=0.5, zorder=1)
        ax_rows.scatter(xs, [row] * len(xs), s=9, zorder=2, edgecolor='none',
                        color=[palette[label_of[b]] for b in loci])
    ax_rows.set_ylim(len(shown), -1)
    ax_rows.set_yticks([])
    ax_rows.set_ylabel('{} {}s'.format(assay, group_name))
    ax_rows.set_xlabel('{} (Mb)'.format(data_chrom))
    ax_rows.set_xlim(lo - span * 0.02, hi + span * 0.02)
    if panels == 'evidence':
        ax_rows.set_title('{}, {}: {} {}s with \u2265 {} community loci ({} shown)'.format(
            title or chrom, assay, len(events), group_name, min_loci, len(shown)), loc='left', fontsize=7.5)
    else:
        ax_rows.set_title('{}: {} {}s containing \u2265 {} loci of these communities ({} shown, coloured by '
                          'community)'.format(assay, len(events), group_name, min_loci, len(shown)), loc='left',
                          fontsize=7.5)
    for side in ('top', 'right', 'left'):
        ax_rows.spines[side].set_visible(False)

    if ax_net is not None:
        layout = nx.spring_layout(sub, seed=seed, weight=None, k=1.2 / np.sqrt(max(1, sub.number_of_nodes())))
        for u, v in sub.edges:
            ax_net.plot([layout[u][0], layout[v][0]], [layout[u][1], layout[v][1]], color=community_color(u, v),
                        lw=0.6, zorder=1)
        for n in sub.nodes:
            bait = graph.nodes[n].get('bait')
            ax_net.scatter(layout[n][0], layout[n][1], marker='s' if bait else 'o', s=22 if bait else 10,
                           color=community_color(n), edgecolor=INK, lw=0.4, zorder=2)
        ax_net.set_axis_off()
        ax_net.set_title('The same communities as a graph', loc='left', fontsize=8)
    ax_bar.bar([0], [observed], color=EDGE_COLOR, width=0.6)
    if len(null):
        ax_bar.bar([1], [null.mean()], color='white', edgecolor='#6b6a64', hatch='////', width=0.6, lw=0.6)
        ax_bar.errorbar([1], [null.mean()], yerr=[null.std()], color=INK, lw=0.8, capsize=0)
    ax_bar.set_xticks([0, 1])
    ax_bar.set_xticklabels(['communities', 'rotated\nlabels'], fontsize=6.5)
    ax_bar.set_ylabel('purity')
    floor = min(observed, null.mean() if len(null) else observed)
    ax_bar.set_ylim(max(0, floor - 0.15), 1.0)
    ax_bar.set_title('Purity' if panels == 'evidence' else 'Community purity of {} {}s'.format(assay, group_name),
                     loc='left', fontsize=8)
    for side in ('top', 'right'):
        ax_bar.spines[side].set_visible(False)
    fig.savefig(ensure_parent(output), dpi=300, bbox_inches='tight')
    plt.close(fig)
    p = float((np.sum(null >= observed) + 1) / (len(null) + 1)) if len(null) else np.nan
    return {'communities': len(ranked), 'nodes': sub.number_of_nodes(), 'events': len(events),
            'purity': observed, 'rotated_mean': float(null.mean()) if len(null) else np.nan,
            'rotated_sd': float(null.std()) if len(null) else np.nan, 'empirical_p': p}


def _graph_bins(graph, data, liftover, chrom):
    nodes = [n for n in graph.nodes]
    for node, (c, p) in node_positions(graph, nodes, liftover).items():
        if c == chrom:
            b = data.bins.bin_of(c, p)
            if b >= 0:
                yield b


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    data = MultiwayData.load(args.multiway)
    liftover = LiftOver(args.chain) if args.chain else None
    if args.style == 'communities':
        if not args.region:
            raise ValueError('--style communities needs --region')
        result = plot_communities(graph, data, parse_region(args.region), args.output, liftover, args.minLoci,
                                  args.maxRows, title=args.title, panels=args.panels)
        log.info('%s', result)
    elif args.style == 'events':
        if args.highlightCommunity is None:
            raise ValueError('--style events needs --highlightCommunity')
        result = plot_events(graph, data, args.highlightCommunity, args.output, liftover, args.minLoci, args.maxRows,
                             args.shifts, title=args.title)
        log.info('%s', result)
    else:
        plot_locus(graph, data, parse_region(args.region), args.output, liftover, args.metric,
                   args.highlightCommunity, args.title)


if __name__ == '__main__':
    main()
