"""Differential graph and multi-way co-occurrence of the explaining paths of one gene, one PDF per panel.

<out>/<name>_step1_graph.pdf   promoter and its direct contacts in the union graph, nodes coloured by the log2 change
                               of one mark (B/A, shared scale -2..2), thick outline where |change| >= 1; edges by
                               status (shared grey, gained red dashed, lost blue dashed); labels give the position (Mb)
<out>/<name>_paths_graph.pdf   all shortest paths from the promoter through step 1 and step 2 to the nodes at step 3
                               whose mark changes at least two-fold in the direction of the expression change, one
                               column per step, every node drawn individually
<out>/<name>_multiway.pdf      number of multi-way groups (GAM tubes, SPRITE clusters, Pore-C concatemers) holding the
                               promoter together with the path constellations below, observed against the same
                               constellation shifted along the chromosome (controls)
<out>/<name>_gam_tubes.pdf     GAM tubes holding the promoter window (rows) by the windows of the path (columns)
<out>/<name>_multiway_counts.tsv   observed count, control mean, sd, quantiles and empirical P per data set and
                                   constellation; <name>_multiway_units.tsv the units (windows or fragments) per node

Constellations (a group "holds" a unit if it contains at least one bin of it):
  promoter_enhancer            promoter and at least one step-3 unit
  promoter_mediator_enhancer   promoter, the step-1 or the step-2 unit of the path, and at least one step-3 unit
  promoter_step1_enhancer / promoter_step2_enhancer   the same with one mediator fixed (TSV only)
  all_units                    promoter, both mediators and every step-3 unit

Units: in binned data (GAM, --bin data) every node is the bin of its midpoint (lifted with --chain), and nodes that
fall into one bin form one unit; mediator and step-3 units equal to the promoter unit (or a step-3 unit equal to a
mediator unit) are removed from their role, so that no constellation is satisfied by the promoter bin alone. In
fragment data (SPRITE, Pore-C, graph assembly) a node is the set of fragments it overlaps. Controls shift all units
by one random offset along the same chromosome (bins: an integer number of bins; fragments: base pairs, the shifted
intervals mapped to fragments again), keep shifts at least the span of the constellation away from it and with
every shifted unit covered by at least one group; --shifts controls per data set, fixed --seed. Empirical
P = (1 + controls >= observed) / (1 + controls).

Usage: python figure_gene_multiway.py <cmp prefix> <promoter node> <mark> <output dir>
           [--gam GAM.npz] [--chain hg19ToHg38.over.chain.gz] [--fragments LABEL=FILE.npz ...]
           [--conditions hESC hNEC] [--step 3] [--name PRRX2] [--shifts 1000] [--seed 1]
<cmp prefix> is the output prefix of egaDifferential (<prefix>_union.graphml, <prefix>_explanations_all_levels.tsv).
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_gene_report import log2_change, middle  # noqa: E402

INK, MUTED = '#222222', '#888888'
EDGE_COLOUR = {'gained': '#b2182b', 'lost': '#2166ac', 'shared': '#9a9a96'}
OBSERVED = '#b2182b'
CONSTELLATIONS = ['promoter_enhancer', 'promoter_mediator_enhancer', 'all_units']
CONSTELLATION_LABELS = {'promoter_enhancer': 'P+E', 'promoter_mediator_enhancer': 'P+M+E',
                        'all_units': 'all'}
COVERAGE_TOLERANCE = 0.2
NORM = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
CMAP = plt.get_cmap('RdBu_r')
plt.rcParams.update({'font.size': 7, 'pdf.fonttype': 42})


# ----------------------------------------------------------------------------------------- constellation counting
def union_groups(group_lists):
    """Sorted unique group ids of a list of group id arrays (a group holding any of the bins)."""
    if not group_lists:
        return np.zeros(0, dtype=np.int64)
    return np.unique(np.concatenate([np.asarray(g, dtype=np.int64) for g in group_lists]))


def intersect(*arrays):
    result = arrays[0]
    for array in arrays[1:]:
        result = np.intersect1d(result, array, assume_unique=True)
    return result


def count_constellations(promoter, mediators, enhancers):
    """Number of groups holding each constellation.

    promoter: sorted unique group ids holding the promoter unit; mediators: list of such arrays (step 1, step 2);
    enhancers: list of such arrays (step-3 units). Returns a dict of counts (see the module doc)."""
    any_enhancer = union_groups(enhancers)
    any_mediator = union_groups(mediators)
    with_enhancer = intersect(promoter, any_enhancer)
    counts = dict(promoter=len(promoter), promoter_enhancer=len(with_enhancer),
                  promoter_mediator_enhancer=len(intersect(with_enhancer, any_mediator)))
    for k, mediator in enumerate(mediators):
        counts['promoter_step{}_enhancer'.format(k + 1)] = len(intersect(with_enhancer, mediator))
    counts['all_units'] = len(intersect(promoter, *mediators, *enhancers)) if enhancers else 0
    return counts


class ChromGroups:
    """Groups holding each bin of one chromosome (rows of the incidence matrix as sorted group id arrays)."""

    def __init__(self, data, chrom):
        first = data.bins.offsets[chrom]
        n = data.bins.n_chrom_bins(chrom)
        block = data.incidence[first:first + n].tocsr()
        block.sort_indices()
        self.indptr, self.indices = block.indptr.astype(np.int64), block.indices.astype(np.int64)
        self.first, self.n = first, n
        self.fixed = data.bins.fixed
        self.n_groups = data.n_groups
        if self.fixed:
            self.bin_size = data.bins.bin_size
            self.size = data.bins.chrom_sizes[chrom]
        else:
            self.starts = data.bins.starts[first:first + n]
            self.ends = data.bins.ends[first:first + n]
            self.size = int(self.ends[-1])

    def groups(self, local_bins):
        return union_groups([self.indices[self.indptr[b]:self.indptr[b + 1]] for b in local_bins])

    def overlapping(self, start, end):
        """Local fragment indices overlapping [start, end) (fragment data)."""
        lo = np.searchsorted(self.ends, start, side='right')
        hi = np.searchsorted(self.starts, end, side='left')
        return list(range(lo, hi))


def unit_groups(chrom_groups, units):
    return [chrom_groups.groups(u) for u in units]


def shifted_controls(chrom_groups, units, roles, n_shifts, seed, span, max_tries=200):
    """Counts of the constellation shifted by random offsets.

    units: for bins, a list of lists of local bin indices; for fragments, a list of (start, end) intervals per node of
    the unit. roles: 'promoter', 'mediator' or 'enhancer' per unit. span: extent of the constellation (bins or bp)."""
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    if chrom_groups.fixed:
        low = -min(min(u) for u in units)
        high = chrom_groups.n - 1 - max(max(u) for u in units)
        candidates = np.array([k for k in range(low, high + 1) if abs(k) > span])
        rng.shuffle(candidates)
        offsets = iter(candidates)
    while len(rows) < n_shifts:
        tried += 1
        if chrom_groups.fixed:
            offset = next(offsets, None)
            if offset is None:
                break
            shifted = [[b + offset for b in u] for u in units]
            sets = unit_groups(chrom_groups, shifted)
        else:
            if tried > n_shifts * max_tries:
                break
            low = -min(s for u in units for s, _ in u)
            high = chrom_groups.size - max(e for u in units for _, e in u)
            offset = int(rng.integers(low, high))
            if abs(offset) <= span:
                continue
            sets = [chrom_groups.groups(sorted({f for s, e in u for f in chrom_groups.overlapping(s + offset,
                                                                                                    e + offset)}))
                    for u in units]
        if any(len(s) == 0 for s in sets):
            continue
        counts = split_count(sets, roles)
        counts['offset'] = int(offset)
        counts['coverage_sum'] = int(sum(len(s) for s in sets))
        rows.append(counts)
    return pd.DataFrame(rows)


def split_count(sets, roles):
    promoter = [s for s, r in zip(sets, roles) if r == 'promoter'][0]
    mediators = [s for s, r in zip(sets, roles) if r == 'mediator']
    enhancers = [s for s, r in zip(sets, roles) if r == 'enhancer']
    return count_constellations(promoter, mediators, enhancers)


def empirical_p(observed, controls):
    controls = np.asarray(controls)
    return (1 + int(np.sum(controls >= observed))) / (1 + len(controls))


# ---------------------------------------------------------------------------------------------- units per data set
def bin_units(graph, path_nodes, chrom_groups, data, liftover):
    """Units of binned data: midpoint bin per node, nodes in one bin merged, roles de-duplicated."""
    rows, unit_of_bin = [], {}
    for entry, (role, nodes) in enumerate(path_nodes):
        for node in nodes:
            attrs = graph.nodes[node]
            chrom, position = attrs['chrom'], middle(attrs)
            if liftover is not None:
                lifted = liftover.position(chrom, position)
                chrom, position = lifted if lifted is not None else (None, None)
            b = data.bins.bin_of(chrom, position) if chrom is not None else -1
            rows.append(dict(node=node, role=role, entry=entry, chrom=chrom, position=position, bin=int(b),
                             local_bin=int(b - chrom_groups.first) if b >= 0 else -1))
    table = pd.DataFrame(rows)
    units, roles, labels = [], [], []
    for entry, (role, _) in enumerate(path_nodes):
        for b, sub in table[table.entry == entry].groupby('local_bin', sort=True):
            if b < 0 or b in unit_of_bin:
                table.loc[sub.index, 'unit'] = unit_of_bin.get(b, -1)
                table.loc[sub.index, 'note'] = 'unmapped' if b < 0 else 'same bin as a {} unit'.format(
                    roles[unit_of_bin[b]])
                continue
            unit_of_bin[b] = len(units)
            table.loc[sub.index, 'unit'] = len(units)
            table.loc[sub.index, 'note'] = '' if len(sub) == 1 else '{} nodes in one bin'.format(len(sub))
            units.append([int(b)])
            roles.append(role)
            labels.append(sub.node.tolist())
    span = max(u[0] for u in units) - min(u[0] for u in units)
    return table, units, roles, labels, span


def fragment_units(graph, path_nodes, chrom_groups):
    """Units of fragment data: the fragments overlapping each node; nodes with the same fragments merged."""
    rows, seen = [], {}
    units, intervals, roles, labels = [], [], [], []
    for role, nodes in path_nodes:
        for node in nodes:
            attrs = graph.nodes[node]
            start, end = int(attrs['start']), int(attrs['end'])
            fragments = tuple(chrom_groups.overlapping(start + 1, end - 1))
            note = ''
            if not fragments:
                note = 'no fragment'
            elif fragments in seen:
                note = 'same fragments as a {} unit'.format(roles[seen[fragments]])
            else:
                seen[fragments] = len(units)
                units.append(list(fragments))
                intervals.append([(start + 1, end - 1)])
                roles.append(role)
                labels.append([node])
            rows.append(dict(node=node, role=role, chrom=attrs['chrom'], position=middle(attrs),
                             fragments=len(fragments), unit=seen.get(fragments, -1), note=note))
    starts = [s for u in intervals for s, _ in u]
    ends = [e for u in intervals for _, e in u]
    return pd.DataFrame(rows), units, intervals, roles, labels, max(ends) - min(starts)


def analyse(label, kind, data, graph, path_nodes, chrom, n_shifts, seed, liftover=None):
    chrom_groups = ChromGroups(data, chrom)
    if chrom_groups.fixed:
        table, units, roles, labels, span = bin_units(graph, path_nodes, chrom_groups, data, liftover)
        control_units = units
    else:
        table, units, control_units, roles, labels, span = fragment_units(graph, path_nodes, chrom_groups)
    sets = unit_groups(chrom_groups, units)
    for i, s in enumerate(sets):
        table.loc[table.unit == i, 'groups_holding_unit'] = len(s)
    observed = split_count(sets, roles)
    coverage = sum(len(s) for s in sets)
    controls = shifted_controls(chrom_groups, control_units, roles, n_shifts, seed, span)
    # coverage-matched subset: shifts whose units are held by about as many groups in total (+-20 %)
    matched = controls[np.abs(controls['coverage_sum'] - coverage) <= COVERAGE_TOLERANCE * coverage]
    rows = []
    for name in ['promoter_enhancer', 'promoter_mediator_enhancer', 'promoter_step1_enhancer',
                 'promoter_step2_enhancer', 'all_units']:
        if name not in observed:
            continue
        values = controls[name].to_numpy()
        # share of the groups holding the promoter, against the same share of the shifted promoter
        share = observed[name] / max(observed['promoter'], 1)
        control_share = (controls[name] / controls['promoter'].clip(lower=1)).to_numpy()
        rows.append(dict(data_set=label, kind=kind, constellation=name, observed=observed[name],
                         control_mean=float(values.mean()), control_sd=float(values.std(ddof=1)) if len(values) > 1 else np.nan,
                         control_median=float(np.median(values)), control_q95=float(np.quantile(values, 0.95)),
                         p_empirical=empirical_p(observed[name], values), shifts=len(values),
                         share_of_promoter_groups=share, control_share_mean=float(control_share.mean()),
                         p_share=empirical_p(share - 1e-12, control_share),
                         control_mean_coverage_matched=float(matched[name].mean()) if len(matched) else np.nan,
                         p_coverage_matched=empirical_p(observed[name], matched[name].to_numpy())
                         if len(matched) else np.nan, shifts_coverage_matched=len(matched),
                         unit_coverage_sum=coverage, control_unit_coverage_sum_mean=float(
                             controls['coverage_sum'].mean()),
                         groups_total=chrom_groups.n_groups, groups_with_promoter=observed['promoter'],
                         control_groups_with_promoter_mean=float(controls['promoter'].mean()),
                         units_promoter=roles.count('promoter'), units_mediator=roles.count('mediator'),
                         units_enhancer=roles.count('enhancer'), span=span))
    table.insert(0, 'data_set', label)
    print('{}: units {} (mediators {}, enhancers {}); {}'.format(
        label, len(units), roles.count('mediator'), roles.count('enhancer'),
        {r['constellation']: (r['observed'], round(r['control_mean'], 2), round(r['p_empirical'], 4))
         for r in rows}), flush=True)
    detail = None
    if chrom_groups.fixed:
        detail = dict(tubes=sets[roles.index('promoter')], sets=sets, roles=roles, labels=labels, units=units)
    return pd.DataFrame(rows), table, controls.assign(data_set=label), detail


# --------------------------------------------------------------------------------------------------- graph panels
def explaining_paths(graph, source, targets):
    paths = []
    for t in targets:
        paths.extend(nx.all_shortest_paths(graph, source, t))
    return paths


def node_style(graph, node, mark, a, b):
    change = log2_change(graph.nodes[node], mark, a, b)
    return change, CMAP(NORM(change)), 1.6 if abs(change) >= 1 else 0.4


def draw_edge(ax, p, q, status, lw=1.0):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=EDGE_COLOUR.get(status, EDGE_COLOUR['shared']),
            ls='--' if status in ('gained', 'lost') else '-', lw=lw, zorder=1)


def gene_label(graph, node, max_genes=2):
    genes = [g for g in str(graph.nodes[node].get('genes') or '').split(',') if g and not g.startswith('LOC')]
    return ','.join(genes[:max_genes])


def legends(fig, ax, mark, a, b, anchor_y=-0.02):
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), ax=ax, orientation='horizontal', fraction=0.05,
                       pad=0.02, shrink=0.55, aspect=22)
    bar.set_label('{} change {}/{}, log$_2$'.format(mark, b, a), fontsize=7)
    bar.set_ticks([-2, -1, 0, 1, 2])
    bar.ax.tick_params(labelsize=7)
    handles = [Line2D([0], [0], color=EDGE_COLOUR['shared'], label='shared'),
               Line2D([0], [0], color=EDGE_COLOUR['gained'], ls='--', label='gained in ' + b),
               Line2D([0], [0], color=EDGE_COLOUR['lost'], ls='--', label='lost in ' + b),
               Line2D([0], [0], marker='o', ls='', mfc='white', mec=INK, mew=0.4, ms=6, label='|change| < 1'),
               Line2D([0], [0], marker='o', ls='', mfc='white', mec=INK, mew=1.6, ms=6, label='|change| ≥ 1')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=7, frameon=False,
              handlelength=1.8, columnspacing=1.0)


def panel_step1(graph, source, mark, a, b, name, path):
    """Promoter (top) and its direct contacts on an arc below, ordered by position."""
    contacts = sorted(graph.neighbors(source), key=lambda n: middle(graph.nodes[n]))
    angles = np.linspace(np.pi * 1.12, np.pi * 1.88, len(contacts)) if len(contacts) > 1 else [1.5 * np.pi]
    pos = {source: (0.0, 0.55)}
    for n, t in zip(contacts, angles):
        pos[n] = (1.15 * np.cos(t), 0.55 + 1.15 * np.sin(t))
    fig, ax = plt.subplots(figsize=(3.3, 2.7))
    for n in contacts:
        draw_edge(ax, pos[source], pos[n], graph.edges[source, n].get('status', 'shared'), lw=1.2)
    for n in [source] + contacts:
        change, colour, lw = node_style(graph, n, mark, a, b)
        bait = int(graph.nodes[n].get('bait', 0) or 0) == 1
        ax.scatter(*pos[n], s=230 if n == source else 150, marker='s' if bait else 'o', c=[colour],
                   edgecolors=INK, linewidths=lw, zorder=3)
        mb = '{:.2f} Mb'.format(middle(graph.nodes[n]) / 1e6)
        if n == source:
            ax.text(pos[n][0], pos[n][1] + 0.17, '{}  {}  ({:+.2f})'.format(name, mb, change), ha='center',
                    va='bottom', fontsize=7.5, fontweight='bold')
        else:
            genes = gene_label(graph, n)
            text = '{}\n{:+.2f}'.format(mb, change) + ('\n' + genes if genes else '')
            ax.text(pos[n][0], pos[n][1] - 0.15, text, ha='center', va='top', fontsize=7, color=INK,
                    linespacing=1.1)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.15, 0.95)
    ax.set_aspect('equal')
    ax.axis('off')
    legends(fig, ax, mark, a, b)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return contacts


def panel_paths(graph, source, targets, step, mark, a, b, name, path):
    """Columns promoter, step 1 .. step: all nodes on the shortest paths to the targets, individually."""
    paths = explaining_paths(graph, source, targets)
    columns = {0: [source]}
    for k in range(1, step + 1):
        columns[k] = sorted({p[k] for p in paths}, key=lambda n: middle(graph.nodes[n]))
    height = max(2.6, 0.36 * max(len(v) for v in columns.values()) + 0.6)
    pos = {}
    for k, nodes in columns.items():
        ys = np.linspace(1, 0, len(nodes)) if len(nodes) > 1 else [0.5]
        for n, y in zip(nodes, ys):
            pos[n] = (float(k), float(y))
    fig, ax = plt.subplots(figsize=(3.4, height))
    drawn = set()
    for p in paths:
        for u, v in zip(p[:-1], p[1:]):
            if (u, v) not in drawn:
                drawn.add((u, v))
                draw_edge(ax, pos[u], pos[v], graph.edges[u, v].get('status', 'shared'))
    for n, (x, y) in pos.items():
        change, colour, lw = node_style(graph, n, mark, a, b)
        bait = int(graph.nodes[n].get('bait', 0) or 0) == 1
        ax.scatter([x], [y], s=150 if n == source else 80, marker='s' if bait else 'o', c=[colour],
                   edgecolors=INK, linewidths=lw, zorder=3)
        mb = '{:.2f} Mb'.format(middle(graph.nodes[n]) / 1e6)
        genes = gene_label(graph, n)
        if x == step:
            text = '{}  {:+.2f}'.format(mb, change) + ('  ' + genes if genes else '')
            ax.text(x + 0.16, y, text, ha='left', va='center', fontsize=7, color=INK)
        else:
            text = '{}\n{}\n{:+.2f}'.format(name if n == source else genes.replace(',', '\n'), mb,
                                            change).lstrip('\n')
            above = int(x) % 2 == 1 and len(columns[int(x)]) == 1
            ax.text(x, y + 0.06 if above else y - 0.06, text, ha='center', va='bottom' if above else 'top',
                    fontsize=7, color=INK, linespacing=1.1, fontweight='bold' if n == source else 'normal')
    for k in range(step + 1):
        ax.text(k, 1.1, 'promoter' if k == 0 else 'step {}'.format(k), ha='center', va='bottom', fontsize=7.5)
    ax.set_xlim(-0.45, step + 1.25)
    ax.set_ylim(-0.12, 1.18)
    ax.axis('off')
    legends(fig, ax, mark, a, b)
    ax.get_legend().set_bbox_to_anchor((0.5, 1.06))
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return columns, paths


# ------------------------------------------------------------------------------------------------ multi-way panels
def panel_controls(counts, controls, path):
    sets = list(dict.fromkeys(counts.data_set))
    ncols = 2 if len(sets) > 1 else 1
    nrows = int(np.ceil(len(sets) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4, 1.55 * nrows + 0.2), squeeze=False)
    for ax, label in zip(axes.ravel(), sets):
        sub = counts[counts.data_set == label].set_index('constellation')
        ctrl = controls[controls.data_set == label]
        # counts of fragment data have long control tails: densities and points on log10(1 + count)
        top = 1.0
        for x, name in enumerate(CONSTELLATIONS):
            values = np.log10(1 + ctrl[name].to_numpy().astype(float))
            observed = np.log10(1 + sub.loc[name, 'observed'])
            top = max(top, values.max(), observed)
            if values.max() > values.min():
                parts = ax.violinplot([values], positions=[x], widths=0.75, showextrema=False)
                for body in parts['bodies']:
                    body.set_facecolor('#c8c8c4')
                    body.set_edgecolor('none')
                    body.set_alpha(1)
            mean = np.log10(1 + sub.loc[name, 'control_mean'])
            ax.plot([x - 0.25, x + 0.25], [mean] * 2, color=MUTED, lw=1.0)
            ax.scatter([x], [observed], marker='D', s=18, color=OBSERVED, zorder=3)
        ticks = [t for t in (0, 1, 10, 100, 1000, 10000) if np.log10(1 + t) <= top * 1.1]
        ax.set_yticks(np.log10(1 + np.array(ticks)))
        ax.set_yticklabels([str(t) for t in ticks])
        ax.set_ylim(-0.05, top * 1.1)
        ax.set_xticks(range(len(CONSTELLATIONS)))
        ax.set_xticklabels(['{}\nP={:.2g}'.format(CONSTELLATION_LABELS[c], sub.loc[c, 'p_empirical'])
                            for c in CONSTELLATIONS], fontsize=7)
        ax.set_title(label, fontsize=7.5, pad=2)
        ax.tick_params(axis='y', labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    for ax in axes.ravel()[len(sets):]:
        ax.axis('off')
    for ax in axes[:, 0]:
        ax.set_ylabel('groups', fontsize=7)
    handles = [Line2D([0], [0], marker='D', ls='', color=OBSERVED, ms=4, label='observed'),
               plt.Rectangle((0, 0), 1, 1, color='#c8c8c4', label='shifted controls'),
               Line2D([0], [0], color=MUTED, label='control mean')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def panel_tubes(detail, graph, path):
    """Tubes holding the promoter window (rows) by the path windows (columns), ordered by position."""
    order = sorted(range(len(detail['units'])), key=lambda i: detail['units'][i][0])
    tubes = detail['tubes']
    matrix = np.array([[np.isin(t, detail['sets'][i]) for i in order] for t in tubes], dtype=float)
    weight = matrix.sum(axis=1)
    matrix = matrix[sorted(range(len(tubes)), key=lambda r: (-weight[r], tuple(-matrix[r])))]
    role_tag = {'promoter': 'P', 'mediator': 'M', 'enhancer': 'E'}
    labels = []
    for i in order:
        mbs = '/'.join('{:.2f}'.format(middle(graph.nodes[n]) / 1e6) for n in detail['labels'][i])
        labels.append('{} {}'.format(role_tag[detail['roles'][i]], mbs))
    fig, ax = plt.subplots(figsize=(3.3, max(2.0, 0.035 * len(tubes) + 1.3)))
    ax.imshow(matrix, aspect='auto', cmap=plt.matplotlib.colors.ListedColormap(['#f2f2ef', INK]),
              interpolation='nearest', vmin=0, vmax=1)
    for x in range(len(order)):
        ax.text(x, -0.8, '{:d}'.format(int(matrix[:, x].sum())), ha='center', va='bottom', fontsize=7)
    ax.text(len(order) - 0.4, -0.8, 'tubes', ha='left', va='bottom', fontsize=7, color=MUTED)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('H1 GAM tubes with the promoter window (n = {})'.format(len(tubes)), fontsize=7)
    ax.set_yticks([])
    ax.tick_params(axis='x', length=0)
    for x in range(len(order) + 1):
        ax.axvline(x - 0.5, color='white', lw=0.6)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('prefix')
    parser.add_argument('promoter')
    parser.add_argument('mark')
    parser.add_argument('output')
    parser.add_argument('--gam', help='GAM incidence (MultiwayData npz, binned, assembly of --chain target)')
    parser.add_argument('--gamLabel', default='H1-hESC GAM')
    parser.add_argument('--chain', help='Chain from the graph assembly to the GAM assembly')
    parser.add_argument('--fragments', nargs='*', default=[],
                        help='LABEL=FILE.npz multi-way data on fragments of the graph assembly')
    parser.add_argument('--conditions', nargs=2, default=['hESC', 'hNEC'])
    parser.add_argument('--step', type=int, default=3)
    parser.add_argument('--name', default=None)
    parser.add_argument('--shifts', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()
    a, b = args.conditions
    os.makedirs(args.output, exist_ok=True)
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    rows = levels[levels.node == args.promoter]
    graph = nx.read_graphml(args.prefix + '_union.graphml')
    name = args.name or (graph.nodes[args.promoter].get('genes') or args.promoter).split(',')[0]
    out = os.path.join(args.output, name)
    sign = np.sign(float(rows['rna_log2FC'].iloc[0]))
    distance = nx.single_source_shortest_path_length(graph, args.promoter, cutoff=args.step)
    targets = [n for n, k in distance.items() if k == args.step
               and sign * log2_change(graph.nodes[n], args.mark, a, b) >= 1]
    panel_step1(graph, args.promoter, args.mark, a, b, name, out + '_step1_graph.pdf')
    columns, paths = panel_paths(graph, args.promoter, targets, args.step, args.mark, a, b, name,
                                 out + '_paths_graph.pdf')
    print('{}: {} concordant nodes at step {}; path nodes per step {}'.format(
        name, len(targets), args.step, {k: len(v) for k, v in columns.items()}), flush=True)
    path_nodes = [('promoter', [args.promoter])] + [('mediator', [n]) for k in range(1, args.step)
                                                     for n in columns[k]] + [('enhancer', columns[args.step])]
    chrom = graph.nodes[args.promoter]['chrom']
    from epigeneticgraphanalyzer.genome import LiftOver
    from epigeneticgraphanalyzer.multiway import MultiwayData
    results, units, controls, gam_detail = [], [], [], None
    data_sets = ([(args.gamLabel, 'GAM', args.gam)] if args.gam else []) + \
        [tuple([item.split('=', 1)[0], 'fragments', item.split('=', 1)[1]]) for item in args.fragments]
    for label, kind, path in data_sets:
        data = MultiwayData.load(path)
        liftover = LiftOver(args.chain) if (kind == 'GAM' and args.chain) else None
        target_chrom = chrom
        if liftover is not None:
            target_chrom = liftover.position(chrom, middle(graph.nodes[args.promoter]))[0]
        result, table, ctrl, detail = analyse(label, kind, data, graph, path_nodes, target_chrom, args.shifts,
                                              args.seed, liftover)
        del data
        results.append(result)
        units.append(table)
        controls.append(ctrl)
        if detail is not None and gam_detail is None:
            gam_detail = detail
    if not results:
        return
    counts = pd.concat(results, ignore_index=True)
    counts.to_csv(out + '_multiway_counts.tsv', sep='\t', index=False)
    pd.concat(units, ignore_index=True).to_csv(out + '_multiway_units.tsv', sep='\t', index=False)
    controls = pd.concat(controls, ignore_index=True)
    controls.to_csv(out + '_multiway_controls.tsv.gz', sep='\t', index=False)
    panel_controls(counts, controls, out + '_multiway.pdf')
    if gam_detail is not None:
        panel_tubes(gam_detail, graph, out + '_gam_tubes.pdf')
    print(counts[['data_set', 'constellation', 'observed', 'control_mean', 'p_empirical']].to_string(index=False))


if __name__ == '__main__':
    main()
