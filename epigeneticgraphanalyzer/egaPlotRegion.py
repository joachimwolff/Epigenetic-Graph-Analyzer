"""egaPlotRegion: a region of the graph as pyGenomeTracks arcs and as a promoter-centred graph."""

import logging
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ._cli import base_parser, ensure_parent, setup_logging  # noqa: E402
from .graph import read_graph  # noqa: E402
from .indirect import indirect_contacts  # noqa: E402

log = logging.getLogger(__name__)

COMMUNITY_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
OTHER_COLOR = '#b9b8b2'
INK = '#1f1f1e'


def parse_arguments(args=None):
    parser = base_parser('Write pyGenomeTracks input (interaction arcs coloured by community, indirect contacts) '
                         'for a region and draw the same region as a graph with the promoters in the centre.')
    parser.add_argument('--graph', '-g', required=True)
    parser.add_argument('--region', required=True, help='chrom:start-end')
    parser.add_argument('--genes', help='Gene annotation for the gene track (BED12 or GTF, uncompressed)')
    parser.add_argument('--maxCommunities', type=int, default=8)
    parser.add_argument('--noIndirect', action='store_true', help='Do not draw indirect contacts')
    parser.add_argument('--title', default=None)
    parser.add_argument('--geneFontSize', type=float, default=9, help='Font size of gene names in the gene track')
    parser.add_argument('--arcHeight', type=float, default=7, help='Height of the interaction track (cm)')
    parser.add_argument('--multiway', help='SPRITE / GAM .npz on the assembly of the graph; adds a track with the '
                                         'clusters or tubes that contain several loci of the region')
    parser.add_argument('--chain', help='Chain file from the graph assembly to the assembly of --multiway; the '
                                        'multi-way track of the pyGenomeTracks output is omitted in this case')
    parser.add_argument('--minLoci', type=int, default=3, help='Minimum region loci per cluster / tube')
    parser.add_argument('--maxEvents', type=int, default=60, help='Clusters / tubes drawn')
    parser.add_argument('--eventsHeight', type=float, default=6, help='Height of the multi-way track (cm)')
    parser.add_argument('--graphSize', nargs=2, type=float, default=[7.0, 5.0], metavar=('WIDTH', 'HEIGHT'),
                        help='Size of the graph figure in inches')
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>_tracks.ini, <prefix>_links_*.tsv and <prefix>_graph.pdf')
    return parser.parse_args(args)


def parse_region(region):
    chrom, span = region.split(':')
    start, end = span.replace(',', '').split('-')
    return chrom, int(start), int(end)


def region_subgraph(graph, region):
    chrom, start, end = region
    nodes = [n for n, d in graph.nodes(data=True)
             if d.get('chrom') == chrom and start <= int(d['start']) and int(d['end']) <= end]
    return graph.subgraph(nodes)


def community_colors(sub, max_communities):
    sizes = {}
    for n, c in sub.nodes(data='community'):
        sizes.setdefault(c, []).append(n)
    ranked = sorted((c for c in sizes if len(sizes[c]) >= 3), key=lambda c: -len(sizes[c]))[:max_communities]
    ranked = sorted(ranked, key=lambda c: min(int(sub.nodes[n]['start']) for n in sizes[c]))
    return {c: COMMUNITY_COLORS[k] for k, c in enumerate(ranked)}


def node_color(sub, colors, node):
    return colors.get(sub.nodes[node].get('community'), OTHER_COLOR)


def write_tracks(sub, colors, prefix, genes, title, indirect=True, gene_font_size=9, arc_height=7, events=None,
                 events_height=6, min_loci=3, max_events=60):
    """pyGenomeTracks configuration with one links file per community, the rest grey, and indirect contacts."""
    directory = os.path.dirname(os.path.abspath(prefix))
    os.makedirs(directory, exist_ok=True)
    groups = {}
    for u, v, data in sub.edges(data=True):
        cu, cv = sub.nodes[u].get('community'), sub.nodes[v].get('community')
        key = cu if (cu == cv and cu in colors) else 'other'
        groups.setdefault(key, []).append((u, v, data.get('weight', 1.0)))

    def line(u, v, score):
        a, b = sub.nodes[u], sub.nodes[v]
        if int(a['start']) > int(b['start']):
            a, b = b, a
        return '{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n'.format(a['chrom'], a['start'], a['end'], b['chrom'], b['start'],
                                                          b['end'], score)

    def distance(u, v):
        a, b = sub.nodes[u], sub.nodes[v]
        return abs(max(int(a['end']), int(b['end'])) - min(int(a['start']), int(b['start'])))

    # one shared y scale for all overlaid links tracks: the longest arc of the region
    ylim = int(max((distance(u, v) for u, v in sub.edges), default=1) * 1.1) + 1
    sections = ['[x-axis]\nwhere = top\n', '[spacer]\nheight = 0.3\n']
    if genes:
        sections.append('[genes]\nfile = {}\nfile_type = gtf\ntitle = genes\nheight = 4\nmerge_transcripts = true\n'
                        'prefered_name = gene_name\nlabels = true\nall_labels_inside = true\nmax_labels = 400\nfontsize = {}\n'.format(os.path.abspath(genes), gene_font_size))
        sections.append('[spacer]\nheight = 0.3\n')
    first = True
    for key in ['other'] + [c for c in colors if c in groups]:
        if key not in groups:
            continue
        path = '{}_links_{}.tsv'.format(prefix, key)
        with open(path, 'w') as handle:
            for u, v, score in groups[key]:
                handle.write(line(u, v, score))
        color = OTHER_COLOR if key == 'other' else colors[key]
        sections.append('[links {}]\nfile = {}\nfile_type = links\nlinks_type = arcs\ncolor = {}\nline_width = 0.8\nalpha = 0.7\n'
                        'height = {}\ncompact_arcs_level = 1\nylim = {}\n{}{}'.format(key, os.path.abspath(path), color, arc_height, ylim,
                                                  'title = {}\n'.format(title or 'capture Hi-C interactions')
                                                  if first else '', '' if first else 'overlay_previous = share-y\n'))
        first = False
    if indirect:
        pairs = indirect_contacts(sub)
        if len(pairs):
            path = prefix + '_links_indirect.tsv'
            with open(path, 'w') as handle:
                for row in pairs.itertuples():
                    handle.write(line(row.node_a, row.node_b, row.path_score))
            indirect_ylim = int(pairs['distance'].max() * 1.2) + 1
            sections.append('[spacer]\nheight = 0.2\n')
            sections.append('[links indirect]\nfile = {}\nfile_type = links\nlinks_type = arcs\ncolor = {}\nline_width = 0.6\n'
                            'line_style = dashed\nalpha = 0.6\norientation = inverted\ncompact_arcs_level = 1\nheight = 3\n'
                            'ylim = {}\ntitle = indirect contacts\n'.format(os.path.abspath(path), '#eb6834', indirect_ylim))
    if events is not None:
        path, count, assay = events
        sections.append('[spacer]\nheight = 0.3\n')
        group = 'tube' if assay == 'GAM' else 'cluster'
        sections.append('[events]\nfile = {}\nfile_type = bed\ndisplay = stacked\ncolor = bed_rgb\nborder_color = none\n'
                        'color_backbone = #cfcfcf\nlabels = false\nheight = {}\nline_width = 0.3\n'
                        'title = {} {}s, one per row; blocks: bins found together in the {} ({} of {})\n'
                        .format(os.path.abspath(path), events_height, assay, group, group, min(count, max_events),
                                count))
    ini = prefix + '_tracks.ini'
    with open(ini, 'w') as handle:
        handle.write('\n'.join(sections))
    return ini


def _hex_to_rgb(color):
    color = color.lstrip('#')
    return ','.join(str(int(color[i:i + 2], 16)) for i in (0, 2, 4))


def write_events(sub, colors, data, prefix, min_loci=3, max_events=60):
    """BED12 of the multi-way groups that contain several loci of the region, one feature per group.

    Blocks are the bins of the region's loci in the group; the colour is the community holding most of them.
    Returns (bed path, number of groups with at least `min_loci` loci).
    """
    from .egaPlotLocus import groups_with_loci
    node_bins = {}
    for n, d in sub.nodes(data=True):
        if sub.degree[n] == 0:
            continue
        b = data.bins.bin_of(d['chrom'], (int(d['start']) + int(d['end'])) // 2)
        if b >= 0 and data.coverage_raw[b] > 0:
            node_bins[n] = b
    community_of_bin = {}
    for n, b in node_bins.items():
        community_of_bin.setdefault(b, []).append(sub.nodes[n].get('community'))
    bins = np.array(sorted(community_of_bin))
    groups = groups_with_loci(data, bins, min_loci)

    def dominant(loci):
        votes = [c for b in loci for c in community_of_bin[b]]
        return max(set(votes), key=votes.count)

    ranked = sorted(groups.values(), key=lambda loci: (-len(loci), int(loci[0])))[:max_events]
    ranked = sorted(ranked, key=lambda loci: (str(dominant(loci)), int(loci[0])))
    chrom = sub.nodes[next(iter(node_bins))]['chrom']
    path = prefix + '_events.bed'
    with open(path, 'w') as handle:
        for k, loci in enumerate(ranked):
            starts = data.bins.bin_start(loci)
            ends = data.bins.bin_end(loci)
            start, end = int(starts.min()), int(ends.max())
            color = _hex_to_rgb(colors.get(dominant(loci), OTHER_COLOR))
            handle.write('\t'.join(map(str, [chrom, start, end, 'event{}'.format(k), 0, '.', start, end, color,
                                              len(loci), ','.join(str(int(e - b)) for b, e in zip(starts, ends)),
                                              ','.join(str(int(x) - start) for x in starts)])) + '\n')
    return path, len(groups)


def promoter_centred_layout(sub, seed=0, aspect=1.4):
    """Hub layout: every connected part of the region is laid out on its own, which places each promoter in the
    centre of its partner fragments and pulls promoters with shared partners together; the parts are arranged in
    genomic order on a grid."""
    import networkx as nx
    components = sorted(nx.connected_components(sub), key=lambda c: min(int(sub.nodes[n]['start']) for n in c))
    components = [c for c in components if len(c) >= 2]
    columns = int(np.ceil(np.sqrt(len(components) * aspect))) or 1
    positions = {}
    for k, nodes in enumerate(components):
        part = sub.subgraph(nodes)
        if len(nodes) == 2:
            layout = {n: np.array([i - 0.5, 0.0]) * 0.6 for i, n in enumerate(nodes)}
        else:
            layout = nx.kamada_kawai_layout(part)
        coordinates = np.array(list(layout.values()))
        coordinates -= coordinates.mean(axis=0)
        scale = np.abs(coordinates).max() or 1.0
        size = 0.42 * np.sqrt(len(nodes) / max(len(c) for c in components)) + 0.18
        row, column = divmod(k, columns)
        for n, xy in zip(layout, coordinates):
            positions[n] = xy / scale * size + np.array([column, -row])
    return positions


def plot_graph(sub, colors, output, title=None, size=(7.0, 5.0)):
    positions = promoter_centred_layout(sub, aspect=size[0] / size[1])
    fig, ax = plt.subplots(figsize=tuple(size))
    for u, v in sub.edges:
        if u not in positions or v not in positions:
            continue
        cu, cv = sub.nodes[u].get('community'), sub.nodes[v].get('community')
        color = colors.get(cu, OTHER_COLOR) if cu == cv else OTHER_COLOR
        ax.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]], color=color, lw=0.6,
                alpha=0.8, zorder=1)
    for n in sub.nodes:
        bait = sub.nodes[n].get('bait')
        if n not in positions:
            continue
        ax.scatter(*positions[n], marker='s' if bait else 'o', s=60 if bait else 10, color=node_color(sub, colors, n),
                   edgecolor=INK, lw=0.5, zorder=2 if not bait else 3)
        if bait:
            names = [g for g in str(sub.nodes[n].get('genes', '')).split(',')
                     if g and not g.startswith(('Gm', 'LOC', 'RP', 'AC0', 'CTD', 'CTC'))]
            if names:
                ax.text(positions[n][0], positions[n][1] + 0.05, names[0], fontsize=6, style='italic',
                        ha='center', va='bottom', color=INK, zorder=4,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.2))
    ax.set_axis_off()
    ax.set_title(title or 'promoter-centred graph', loc='left', fontsize=8)
    fig.savefig(ensure_parent(output), bbox_inches='tight')
    plt.close(fig)


def _locus_of(data, attrs, liftover=None):
    chrom, position = attrs['chrom'], (int(attrs['start']) + int(attrs['end'])) // 2
    if liftover is not None:
        lifted = liftover.position(chrom, position)
        if lifted is None:
            return -1
        chrom, position = lifted
    return data.bins.bin_of(chrom, position)


def plot_graph_multiway(sub, colors, data, output, min_loci=3, n_examples=4, max_overlap=0.5, title=None,
                        example_bins=(4, 8), liftover=None):
    """Promoter-centred graph with multi-way evidence.

    Panel a: interactions coloured by the number of groups (SPRITE clusters / GAM tubes) holding both fragments,
    nodes by the number of groups with at least `min_loci` region loci that hold them. Further panels: the same
    graph with the fragments of single groups highlighted; examples are the largest groups whose bin sets overlap
    each other by at most `max_overlap` (Jaccard).
    """
    from .egaPlotLocus import groups_with_loci
    positions = promoter_centred_layout(sub, aspect=1.3)
    node_bin = {}
    for n, d in sub.nodes(data=True):
        if n in positions:
            b = _locus_of(data, d, liftover)
            if b >= 0 and data.coverage_raw[b] > 0:
                node_bin[n] = b
    bins = np.array(sorted(set(node_bin.values())))
    groups = groups_with_loci(data, bins, min_loci)
    per_bin = {}
    for loci in groups.values():
        for b in loci:
            per_bin[b] = per_bin.get(b, 0) + 1
    edges = [(u, v) for u, v in sub.edges if u in node_bin and v in node_bin]
    # NPMI of the two ends relative to all detected locus pairs at the same offset (the metric of the validation)
    max_offset = max([abs(node_bin[u] - node_bin[v]) for u, v in edges] + [1])
    expected = data.expected_by_offset(max_offset, 'npmi')
    a = np.array([node_bin[u] for u, _ in edges], dtype=np.int64)
    b = np.array([node_bin[v] for _, v in edges], dtype=np.int64)
    if len(edges):
        pair = data.npmi(a, b) - expected[np.clip(np.abs(b - a), 1, max_offset)]
        pair[a == b] = np.nan
    else:
        pair = np.zeros(0)
    examples, chosen = [], []
    candidates = [g for g in groups if example_bins[0] <= len(groups[g]) <= example_bins[1]]
    for g in sorted(candidates, key=lambda g: (-len(groups[g]), g)):
        current = set(groups[g].tolist())
        if all(len(current & c) / len(current | c) <= max_overlap for c in chosen):
            chosen.append(current)
            examples.append(current)
        if len(examples) == n_examples:
            break
    group = 'tube' if data.kind == 'gam' else 'cluster'
    assay = 'GAM' if data.kind == 'gam' else 'SPRITE'

    fig = plt.figure(figsize=(7.1, 4.6))
    grid = fig.add_gridspec(2, 4, wspace=0.08, hspace=0.25)
    ax_main = fig.add_subplot(grid[:, :2])
    cmap = plt.get_cmap('RdBu_r')
    finite = pair[np.isfinite(pair)] if len(pair) else np.zeros(1)
    limit = max(np.nanpercentile(np.abs(finite), 95), 0.1)
    norm = plt.Normalize(-limit, limit)
    order = np.argsort(np.nan_to_num(pair, nan=-np.inf))
    for k in order:
        u, v = edges[k]
        if not np.isfinite(pair[k]):
            continue
        ax_main.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]], color=cmap(norm(pair[k])),
                     lw=0.5 + 1.2 * min(1.0, abs(pair[k]) / limit), zorder=1)
    nodes = [n for n in sub.nodes if n in positions]
    for n in nodes:
        bait = sub.nodes[n].get('bait')
        ax_main.scatter(*positions[n], marker='s' if bait else 'o', s=34 if bait else 8, color=node_color(sub, colors, n),
                        edgecolor=INK, lw=0.4, zorder=2)
    for n in nodes:
        names = [g for g in str(sub.nodes[n].get('genes', '')).split(',')
                 if g and not g.startswith(('Gm', 'LOC', 'RP', 'AC0', 'CTD', 'CTC'))]
        if sub.nodes[n].get('bait') and names:
            ax_main.text(positions[n][0], positions[n][1] + 0.04, names[0], fontsize=5.5, style='italic',
                         ha='center', va='bottom', color=INK, zorder=4,
                         bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.1))
    bar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax_main, orientation='horizontal',
                       fraction=0.05, pad=0.02, shrink=0.7)
    bar.set_label('{} NPMI of the two ends minus NPMI at the same offset'.format(assay), fontsize=6.5)
    bar.ax.tick_params(labelsize=6)
    ax_main.set_axis_off()
    ax_main.set_title('a  interactions by {} co-occurrence'.format(assay), loc='left', fontsize=7.5)
    for k, members in enumerate(examples):
        ax = fig.add_subplot(grid[k // 2, 2 + k % 2])
        inside = {n for n, b in node_bin.items() if b in members}
        for u, v in sub.edges:
            if u in positions and v in positions:
                both = u in inside and v in inside
                ax.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]],
                        color='#2a78d6' if both else '#e2e1dc', lw=0.7 if both else 0.3, zorder=2 if both else 1)
        for n in nodes:
            hit = n in inside
            ax.scatter(*positions[n], marker='s' if sub.nodes[n].get('bait') else 'o',
                       s=(26 if sub.nodes[n].get('bait') else 7), color='#eb6834' if hit else '#e2e1dc',
                       edgecolor=INK if hit else 'none', lw=0.3, zorder=3 if hit else 1)
        ax.set_axis_off()
        unit = 'fragments' if not data.bins.fixed else '{} kb windows'.format(data.bins.bin_size // 1000)
        ax.set_title('{}  one {}: {} {}'.format('bcde'[k], group, len(members), unit), loc='left', fontsize=6.5)
    fig.savefig(ensure_parent(output), bbox_inches='tight')
    plt.close(fig)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.graph)
    sub = region_subgraph(graph, parse_region(args.region))
    colors = community_colors(sub, args.maxCommunities)
    events = None
    if args.multiway:
        from .multiway import MultiwayData
        data = MultiwayData.load(args.multiway)
        liftover = None
        if args.chain:
            from .genome import LiftOver
            liftover = LiftOver(args.chain)
        else:
            path, count = write_events(sub, colors, data, args.outputPrefix, args.minLoci, args.maxEvents)
            events = (path, count, 'GAM' if data.kind == 'gam' else 'SPRITE')
            log.info('%d multi-way groups with >= %d loci of the region', count, args.minLoci)
        plot_graph_multiway(sub, colors, data, args.outputPrefix + '_graph_multiway.pdf', args.minLoci,
                            title=args.title, liftover=liftover)
    ini = write_tracks(sub, colors, args.outputPrefix, args.genes, args.title, indirect=not args.noIndirect,
                       gene_font_size=args.geneFontSize, arc_height=args.arcHeight, events=events,
                       events_height=args.eventsHeight, min_loci=args.minLoci, max_events=args.maxEvents)
    plot_graph(sub, colors, args.outputPrefix + '_graph.pdf', args.title, args.graphSize)
    # region spanned by the fragments with contacts, for plotting without empty flanks
    connected = [n for n in sub.nodes if sub.degree[n] > 0]
    if connected:
        first = min(int(sub.nodes[n]['start']) for n in connected)
        last = max(int(sub.nodes[n]['end']) for n in connected)
        margin = max(1000, int((last - first) * 0.01))
        with open(args.outputPrefix + '_extent.txt', 'w') as handle:
            handle.write('{}:{}-{}\n'.format(sub.nodes[connected[0]]['chrom'], max(0, first - margin), last + margin))
    log.info('%d nodes, %d edges, %d coloured communities; tracks: %s', sub.number_of_nodes(),
             sub.number_of_edges(), len(colors), ini)


if __name__ == '__main__':
    main()
