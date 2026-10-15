"""Inference of indirect contacts from the capture Hi-C graph.

An indirect contact is a pair of nodes (u, w) on the same chromosome that are
not connected by a significant interaction but share at least one neighbour v
(the mediator). Pairs are classified by the capture design:

* ``oe-oe``: neither end is baited. Promoter capture Hi-C cannot measure this
  contact at all, so the graph is the only source of evidence for it.
* ``bait-oe``: the pair was measurable but not called significant.
* ``bait-bait``: both ends baited, measurable but not significant.

Each pair is scored by the number of mediators, the Adamic-Adar index and a
weighted path score (sum over mediators of the weaker of the two edges).
"""

import logging
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def pair_class(is_bait_u, is_bait_w):
    if is_bait_u and is_bait_w:
        return 'bait-bait'
    if is_bait_u or is_bait_w:
        return 'bait-oe'
    return 'oe-oe'


def _node_midpoint(graph, node):
    data = graph.nodes[node]
    return (int(data['start']) + int(data['end'])) // 2


def indirect_contacts(graph, max_distance=2000000, min_distance=0, max_mediator_degree=None,
                      weight='weight'):
    """Enumerate all indirect pairs (two-hop paths without a direct edge).

    Returns a DataFrame with one row per pair and the column 'mediators'
    holding the comma separated mediator nodes.
    """
    pairs = defaultdict(lambda: [0, 0.0, 0.0, []])
    for mediator in graph.nodes:
        degree = graph.degree[mediator]
        if degree < 2:
            continue
        if max_mediator_degree is not None and degree > max_mediator_degree:
            continue
        adamic_adar = 1.0 / np.log(degree)
        # genomic order, so that node_a is always the upstream end of a pair
        neighbours = sorted(graph.neighbors(mediator), key=lambda n: (graph.nodes[n]['chrom'], _node_midpoint(graph, n)))
        chrom = graph.nodes[mediator]['chrom']
        positions = {n: _node_midpoint(graph, n) for n in neighbours}
        for u, w in combinations(neighbours, 2):
            if graph.nodes[u]['chrom'] != chrom or graph.nodes[w]['chrom'] != chrom:
                continue
            distance = abs(positions[u] - positions[w])
            if distance > max_distance or distance < min_distance:
                continue
            if graph.has_edge(u, w):
                continue
            entry = pairs[(u, w)]
            entry[0] += 1
            entry[1] += adamic_adar
            entry[2] += min(graph[u][mediator].get(weight, 1.0), graph[mediator][w].get(weight, 1.0))
            entry[3].append(mediator)
    rows = []
    for (u, w), (n_mediators, adamic_adar, path_score, mediators) in pairs.items():
        bait_u = bool(graph.nodes[u].get('bait'))
        bait_w = bool(graph.nodes[w].get('bait'))
        rows.append({
            'node_a': u, 'node_b': w, 'chrom': graph.nodes[u]['chrom'],
            'distance': abs(_node_midpoint(graph, u) - _node_midpoint(graph, w)),
            'pair_class': pair_class(bait_u, bait_w),
            'n_mediators': n_mediators, 'adamic_adar': adamic_adar, 'path_score': path_score,
            'mediator_is_bait': ','.join('1' if graph.nodes[m].get('bait') else '0' for m in mediators),
            'mediators': ','.join(mediators),
            'same_community': graph.nodes[u].get('community', -1) == graph.nodes[w].get('community', -2),
        })
    columns = ['node_a', 'node_b', 'chrom', 'distance', 'pair_class', 'n_mediators', 'adamic_adar',
               'path_score', 'mediator_is_bait', 'mediators', 'same_community']
    result = pd.DataFrame(rows, columns=columns)
    log.info('%d indirect pairs (%s)', len(result),
             ', '.join('{}={}'.format(k, v) for k, v in result['pair_class'].value_counts().items()))
    return result


def direct_contacts(graph, max_distance=2000000):
    """All edges of the graph as a pair table (positive reference set)."""
    rows = []
    for u, w, data in graph.edges(data=True):
        if graph.nodes[u]['chrom'] != graph.nodes[w]['chrom']:
            continue
        if _node_midpoint(graph, u) > _node_midpoint(graph, w):
            u, w = w, u
        distance = abs(_node_midpoint(graph, u) - _node_midpoint(graph, w))
        if distance > max_distance:
            continue
        rows.append({'node_a': u, 'node_b': w, 'chrom': graph.nodes[u]['chrom'], 'distance': distance,
                     'pair_class': pair_class(bool(graph.nodes[u].get('bait')), bool(graph.nodes[w].get('bait'))),
                     'weight': data.get('weight', 1.0)})
    return pd.DataFrame(rows, columns=['node_a', 'node_b', 'chrom', 'distance', 'pair_class', 'weight'])


def two_hop_triplets(graph, max_distance=2000000, require_open=True, max_mediator_degree=None):
    """Paths u - v - w (v the mediator). With `require_open` only paths whose ends are not adjacent."""
    rows = []
    for mediator in graph.nodes:
        degree = graph.degree[mediator]
        if degree < 2 or (max_mediator_degree is not None and degree > max_mediator_degree):
            continue
        neighbours = sorted(graph.neighbors(mediator), key=lambda n: (graph.nodes[n]['chrom'], _node_midpoint(graph, n)))
        chrom = graph.nodes[mediator]['chrom']
        for u, w in combinations(neighbours, 2):
            if graph.nodes[u]['chrom'] != chrom or graph.nodes[w]['chrom'] != chrom:
                continue
            span = max(_node_midpoint(graph, n) for n in (u, mediator, w)) - \
                min(_node_midpoint(graph, n) for n in (u, mediator, w))
            if span > max_distance:
                continue
            closed = graph.has_edge(u, w)
            if require_open and closed:
                continue
            rows.append({'node_a': u, 'mediator': mediator, 'node_b': w, 'chrom': chrom, 'span': span,
                         'closed': closed, 'mediator_is_bait': bool(graph.nodes[mediator].get('bait')),
                         'pair_class': pair_class(bool(graph.nodes[u].get('bait')),
                                                  bool(graph.nodes[w].get('bait')))})
    return pd.DataFrame(rows, columns=['node_a', 'mediator', 'node_b', 'chrom', 'span', 'closed',
                                       'mediator_is_bait', 'pair_class'])
