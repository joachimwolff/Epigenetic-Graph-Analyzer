"""Null models and pair sets of the graph for the multi-way step curves.

* graph_distance_pairs: node pairs of one chromosome at graph distance 1..max_steps within a genomic distance limit
* degree_bin_classes: per multi-way bin, the quantile class of the graph degree of the node(s) in the bin, computed
  separately for bait and other-end loci, so that controls can match the degree of both ends of a tested pair
* rewire_by_distance: degree-preserving double-edge swaps that stay on one chromosome, keep the bait / other-end
  structure of every edge and keep each edge length within its logarithmic length bin
"""

import bisect
import logging
import math

import networkx as nx
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def node_midpoints(graph):
    return {n: (int(d['start']) + int(d['end'])) // 2 for n, d in graph.nodes(data=True)}


def graph_distance_pairs(graph, max_distance=2000000, max_steps=5):
    """Pairs (node_a < node_b) on one chromosome at shortest-path distance 1..max_steps with genomic distance (of the
    node midpoints) <= max_distance; max_distance None keeps all pairs of one chromosome.
    Returns a DataFrame node_a, node_b, chrom, distance, steps."""
    middle = node_midpoints(graph)
    rows = []
    for source in graph.nodes:
        lengths = nx.single_source_shortest_path_length(graph, source, cutoff=max_steps)
        chrom = graph.nodes[source]['chrom']
        for target, steps in lengths.items():
            if target <= source or graph.nodes[target]['chrom'] != chrom:
                continue
            distance = abs(middle[target] - middle[source])
            if max_distance is not None and distance > max_distance:
                continue
            rows.append((source, target, chrom, distance, steps))
    return pd.DataFrame(rows, columns=['node_a', 'node_b', 'chrom', 'distance', 'steps'])


def degree_bin_classes(graph, node_bins, n_bins, n_classes=5):
    """Degree class of every multi-way bin that hosts a graph node, -1 for the other bins.

    The degree of a bin is the largest degree of its nodes; a bin hosting a bait is a bait locus. Classes are the
    quantile classes of the degree among bait loci (n_classes .. 2 n_classes - 1) and, separately, among other-end
    loci (0 .. n_classes - 1), so that the class also encodes the kind of the locus. Tied degrees (most other ends
    have degree 1) merge quantile classes, so fewer than n_classes classes per kind can occur."""
    degree = np.zeros(n_bins, dtype=np.int64)
    is_bait = np.zeros(n_bins, dtype=bool)
    hosted = np.zeros(n_bins, dtype=bool)
    for node, bin_id in node_bins.items():
        degree[bin_id] = max(degree[bin_id], graph.degree(node))
        is_bait[bin_id] |= bool(graph.nodes[node].get('bait'))
        hosted[bin_id] = True
    classes = np.full(n_bins, -1, dtype=np.int64)
    for kind, mask in ((0, hosted & ~is_bait), (1, hosted & is_bait)):
        if not mask.any():
            continue
        values = degree[mask]
        edges = np.unique(np.quantile(values, np.linspace(0, 1, n_classes + 1)[1:-1]))
        classes[mask] = kind * n_classes + np.searchsorted(edges, values, side='left')
    return classes


def _length_bin(distance, width):
    return np.floor(np.log10(np.maximum(np.asarray(distance, dtype=float), 1.0)) / width).astype(np.int64)


def rewire_by_distance(graph, swaps_per_edge=10, log_bin_width=0.1, seed=0, max_attempts_factor=20, max_length=None):
    """Degree-preserving rewiring by double-edge swaps restricted to one chromosome and similar edge length.

    Every edge is oriented bait -> other end (bait-bait edges from the upstream node); two edges (a1, b1), (a2, b2)
    of the same chromosome, the same kind (bait-other or bait-bait) and the same length bin (log10 of the midpoint
    distance in steps of log_bin_width) become (a1, b2), (a2, b1) if both new edges fall in the same length bin, are no
    self loops and do not exist yet. Node degrees, the bait / other-end kind of every edge end, the number of edges per
    chromosome and kind, and the number of edges per length bin are kept. For an edge (a1, b1), the partner edge is
    drawn among the edges of its group whose other end b2 lies at a distance from a1 within the length bin, which
    keeps the acceptance rate high. With max_length, new edges longer than max_length bp are rejected as well (for
    wide length bins). Trans edges are left unchanged. Returns (rewired graph, statistics dict)."""
    rng = np.random.default_rng(seed)
    middle = node_midpoints(graph)
    nodes = list(graph.nodes)
    index = {n: i for i, n in enumerate(nodes)}
    position = np.array([middle[n] for n in nodes], dtype=np.int64)
    bait = np.array([bool(graph.nodes[n].get('bait')) for n in nodes])
    chrom = [graph.nodes[n]['chrom'] for n in nodes]
    anchor, partner, groups = [], [], {}
    fixed = []
    for u, v in graph.edges():
        i, j = index[u], index[v]
        if chrom[i] != chrom[j]:
            fixed.append((i, j))
            continue
        # orientation: bait end first; between two baits (or two other ends) the upstream node first
        if bait[j] and not bait[i] or (bait[i] == bait[j] and position[j] < position[i]):
            i, j = j, i
        kind = int(bait[i]) + int(bait[j])
        key = (chrom[i], kind, int(_length_bin(abs(position[j] - position[i]), log_bin_width)))
        groups.setdefault(key, []).append(len(anchor))
        anchor.append(i)
        partner.append(j)
    anchor = np.array(anchor, dtype=np.int64)
    partner = np.array(partner, dtype=np.int64)
    original = set(zip(np.minimum(anchor, partner).tolist(), np.maximum(anchor, partner).tolist()))
    present = set(original)
    # per group, the partner positions in sorted slots; swaps exchange the owners of two slots, the positions stay.
    # The loop uses plain Python containers and pre-drawn random numbers (millions of attempts).
    n_edges = len(anchor)
    pos = position.tolist()
    anchor_list, partner_list = anchor.tolist(), partner.tolist()
    group_of = [0] * n_edges
    edge_slot = [0] * n_edges
    slot_position, slot_edge, group_bin = [], [], []
    for g, (key, edge_ids) in enumerate(groups.items()):
        edge_ids = sorted(edge_ids, key=lambda e: pos[partner_list[e]])
        slot_position.append([pos[partner_list[e]] for e in edge_ids])
        slot_edge.append(list(edge_ids))
        for k, e in enumerate(edge_ids):
            edge_slot[e] = k
            group_of[e] = g
        group_bin.append(key[2])
    bounds = [(10 ** (b * log_bin_width), 10 ** ((b + 1) * log_bin_width)) for b in group_bin]

    def length_bin(d):
        return math.floor(math.log10(max(d, 1)) / log_bin_width)

    target = int(swaps_per_edge * n_edges)
    attempts = accepted = 0
    max_attempts = max_attempts_factor * target
    batch = 1 << 20
    draws_e, draws_u, k_draw = [], [], batch
    while accepted < target and attempts < max_attempts and n_edges > 1:
        if k_draw >= len(draws_e):
            draws_e = rng.integers(0, n_edges, batch).tolist()
            draws_u = rng.random(batch).tolist()
            k_draw = 0
        e1, u = draws_e[k_draw], draws_u[k_draw]
        k_draw += 1
        attempts += 1
        g = group_of[e1]
        slots = slot_position[g]
        if len(slots) < 2:
            continue
        a1, b1 = anchor_list[e1], partner_list[e1]
        low, high = bounds[g]
        p1 = pos[a1]
        # slots whose position p satisfies low <= |p - a1| < high (downstream and upstream ranges)
        d0, d1 = bisect.bisect_left(slots, p1 + low), bisect.bisect_left(slots, p1 + high)
        u0, u1 = bisect.bisect_right(slots, p1 - high), bisect.bisect_right(slots, p1 - low)
        n_down = d1 - d0
        n_candidates = n_down + (u1 - u0)
        if n_candidates <= 0:
            continue
        pick = min(int(u * n_candidates), n_candidates - 1)
        slot = d0 + pick if pick < n_down else u0 + pick - n_down
        e2 = slot_edge[g][slot]
        if e2 == e1:
            continue
        a2, b2 = anchor_list[e2], partner_list[e2]
        if a1 == a2 or b1 == b2 or a1 == b2 or a2 == b1:
            continue
        g_bin = group_bin[g]
        l1, l2 = abs(pos[b2] - p1), abs(pos[b1] - pos[a2])
        if length_bin(l1) != g_bin or length_bin(l2) != g_bin:
            continue
        if max_length is not None and (l1 > max_length or l2 > max_length):
            continue
        new1 = (a1, b2) if a1 < b2 else (b2, a1)
        new2 = (a2, b1) if a2 < b1 else (b1, a2)
        if new1 in present or new2 in present:
            continue
        present.discard((a1, b1) if a1 < b1 else (b1, a1))
        present.discard((a2, b2) if a2 < b2 else (b2, a2))
        present.add(new1)
        present.add(new2)
        partner_list[e1], partner_list[e2] = b2, b1
        s1, s2 = edge_slot[e1], edge_slot[e2]
        edge_slot[e1], edge_slot[e2] = s2, s1
        slot_edge[g][s1], slot_edge[g][s2] = e2, e1
        accepted += 1
    partner = np.array(partner_list, dtype=np.int64)
    rewired = nx.Graph()
    rewired.add_nodes_from(graph.nodes(data=True))
    rewired.add_edges_from((nodes[i], nodes[j]) for i, j in zip(anchor, partner))
    rewired.add_edges_from((nodes[i], nodes[j]) for i, j in fixed)
    final = set(zip(np.minimum(anchor, partner).tolist(), np.maximum(anchor, partner).tolist()))
    statistics = dict(edges=n_edges + len(fixed), swappable_edges=n_edges, attempts=attempts, accepted=accepted,
                      fraction_edges_changed=float(len(final - original) / max(1, n_edges)))
    log.info('rewiring: %s', statistics)
    return rewired, statistics


def graph_statistics(graph, max_distance=2000000):
    """Size and shape of an interaction graph: nodes, baits, edges, cis edges (and those within max_distance),
    cis edge length (node midpoints: median, quartiles), mean degree of baits and other ends, largest component."""
    middle = node_midpoints(graph)
    lengths = np.array([abs(middle[u] - middle[v]) for u, v in graph.edges()
                        if graph.nodes[u]['chrom'] == graph.nodes[v]['chrom']], dtype=float)
    bait = {n: bool(d.get('bait')) for n, d in graph.nodes(data=True)}
    bait_degree = [graph.degree(n) for n in graph.nodes if bait[n]]
    other_degree = [graph.degree(n) for n in graph.nodes if not bait[n]]
    largest = max((len(c) for c in nx.connected_components(graph)), default=0)
    quartiles = np.percentile(lengths, [25, 50, 75]) if len(lengths) else [np.nan] * 3
    return dict(nodes=graph.number_of_nodes(), baits=len(bait_degree), other_ends=len(other_degree),
                edges=graph.number_of_edges(), cis_edges=len(lengths),
                cis_edges_within_limit=int((lengths <= max_distance).sum()),
                bait_bait_edges=sum(1 for u, v in graph.edges() if bait[u] and bait[v]),
                median_edge_length=float(quartiles[1]), edge_length_q1=float(quartiles[0]),
                edge_length_q3=float(quartiles[2]),
                mean_bait_degree=float(np.mean(bait_degree)) if bait_degree else np.nan,
                median_bait_degree=float(np.median(bait_degree)) if bait_degree else np.nan,
                mean_other_end_degree=float(np.mean(other_degree)) if other_degree else np.nan,
                largest_component_fraction=largest / max(1, graph.number_of_nodes()))
