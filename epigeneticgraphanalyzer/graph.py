"""Interaction graph of promoter capture Hi-C data.

Nodes are genomic fragments identified by 'chrom:start-end'. A fragment that
is baited in the capture design is a bait node, every other fragment an
other-end (OE) node. Because nodes are keyed by coordinates, a baited fragment
that appears as the other end of a promoter-promoter interaction is the same
node as the bait, so paths through promoters are preserved and graphs of
different conditions share their node space.
"""

import logging
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd

from .genome import parse_node_id

log = logging.getLogger(__name__)

WEIGHT_CHOICES = ('score', 'count', 'unit')


def edge_weights(df, weight):
    """Positive edge weights from an interaction table."""
    if weight == 'unit':
        return np.ones(len(df))
    if weight == 'count':
        values = df['count'].to_numpy(dtype=float)
    elif weight == 'score':
        values = df['score'].to_numpy(dtype=float)
    else:
        raise ValueError('Unknown weight {}'.format(weight))
    values = np.nan_to_num(values, nan=0.0)
    return np.clip(values, 1e-6, None)


def annotate_genes(fragments, tss):
    """Map transcription start sites onto fragments.

    `fragments` is a DataFrame with chrom, start, end, node; returns a dict
    node -> (sorted gene names, sorted gene ids).
    """
    result = {}
    if tss is None or len(tss) == 0:
        return result
    for chrom, frag in fragments.groupby('chrom'):
        sites = tss[tss['chrom'] == chrom].sort_values('tss')
        if len(sites) == 0:
            continue
        positions = sites['tss'].to_numpy()
        names = sites['gene_name'].to_numpy()
        ids = sites['gene_id'].to_numpy()
        lo = np.searchsorted(positions, frag['start'].to_numpy(), side='left')
        hi = np.searchsorted(positions, frag['end'].to_numpy(), side='left')
        for node, a, b in zip(frag['node'].to_numpy(), lo, hi):
            if b > a:
                result[node] = (sorted({str(x) for x in names[a:b] if str(x)}),
                                sorted({str(x) for x in ids[a:b] if str(x)}))
    return result


def build_graph(interactions, tss=None, weight='score', significant_only=True, bait_fragments=None):
    """Create an undirected weighted graph from a (significant) interaction table.

    Edge attributes: weight, score, count, expected, distance.
    Node attributes: chrom, start, end, bait (bool), genes, gene_ids, position.
    """
    df = interactions
    if significant_only and 'significant' in df.columns:
        df = df[df['significant'].astype(bool)]
    df = df.reset_index(drop=True)
    baits = df['chrom'] + ':' + df['bait_start'].astype(str) + '-' + df['bait_end'].astype(str)
    others = df['chrom'] + ':' + df['other_start'].astype(str) + '-' + df['other_end'].astype(str)
    weights = edge_weights(df, weight)

    graph = nx.Graph()
    bait_nodes = set(baits)
    bait_nodes.update(others[df['other_is_bait'].to_numpy()])
    if bait_fragments is not None:
        bait_nodes.update(bait_fragments)

    count = df['count'].to_numpy(dtype=float)
    expected = df['expected'].to_numpy(dtype=float)
    score = df['score'].to_numpy(dtype=float)
    distance = df['distance'].to_numpy()
    for i, (u, v) in enumerate(zip(baits.to_numpy(), others.to_numpy())):
        if u == v:
            continue
        if graph.has_edge(u, v):
            # the same pair reported twice (e.g. promoter-promoter both ways): keep the stronger
            if graph[u][v]['weight'] >= weights[i]:
                continue
        graph.add_edge(u, v, weight=float(weights[i]), score=float(score[i]), count=float(count[i]),
                       expected=float(expected[i]), distance=int(distance[i]))

    fragments = []
    for node in graph.nodes:
        chrom, start, end = parse_node_id(node)
        fragments.append((node, chrom, start, end))
    fragments = pd.DataFrame(fragments, columns=['node', 'chrom', 'start', 'end'])
    genes = annotate_genes(fragments, tss)
    for node, chrom, start, end in fragments.itertuples(index=False):
        attributes = graph.nodes[node]
        attributes['chrom'] = chrom
        attributes['start'] = int(start)
        attributes['end'] = int(end)
        attributes['position'] = node
        attributes['bait'] = node in bait_nodes
        names, ids = genes.get(node, ([], []))
        attributes['genes'] = ','.join(names)
        attributes['gene_ids'] = ','.join(ids)
    log.info('graph: %d nodes (%d baits), %d edges', graph.number_of_nodes(),
             sum(1 for _, b in graph.nodes(data='bait') if b), graph.number_of_edges())
    return graph


def detect_communities(graph, method='leiden', resolution=1.0, seed=42, weight='weight'):
    """Partition the graph into communities. Returns a list of node sets, largest first.

    'leiden' uses igraph/leidenalg (modularity with resolution parameter),
    'louvain' the networkx implementation.
    """
    if graph.number_of_edges() == 0:
        return [{node} for node in graph.nodes]
    if method == 'louvain':
        communities = nx.community.louvain_communities(graph, weight=weight, resolution=resolution, seed=seed)
    elif method == 'leiden':
        import igraph as ig
        import leidenalg
        nodes = list(graph.nodes)
        index = {node: i for i, node in enumerate(nodes)}
        edges = [(index[u], index[v]) for u, v in graph.edges]
        weights = [graph[u][v].get(weight, 1.0) for u, v in graph.edges]
        ig_graph = ig.Graph(n=len(nodes), edges=edges)
        partition = leidenalg.find_partition(
            ig_graph, leidenalg.RBConfigurationVertexPartition, weights=weights,
            resolution_parameter=resolution, seed=seed)
        communities = [{nodes[i] for i in members} for members in partition]
    else:
        raise ValueError('Unknown community method {}'.format(method))
    communities = sorted(communities, key=len, reverse=True)
    for community_id, members in enumerate(communities):
        for node in members:
            graph.nodes[node]['community'] = community_id
    return communities


def community_table(graph, communities):
    """One row per community: size, number of baits, genes, genomic span."""
    rows = []
    for community_id, members in enumerate(communities):
        starts = [graph.nodes[n]['start'] for n in members]
        ends = [graph.nodes[n]['end'] for n in members]
        genes = sorted({g for n in members for g in graph.nodes[n].get('genes', '').split(',') if g})
        rows.append({
            'community': community_id,
            'chrom': graph.nodes[next(iter(members))]['chrom'],
            'start': min(starts), 'end': max(ends),
            'n_nodes': len(members),
            'n_baits': sum(1 for n in members if graph.nodes[n].get('bait')),
            'n_edges': graph.subgraph(members).number_of_edges(),
            'genes': ','.join(genes),
            'nodes': ','.join(sorted(members)),
        })
    return pd.DataFrame(rows)


def node_table(graph):
    rows = []
    for node, data in graph.nodes(data=True):
        rows.append({'node': node, 'chrom': data.get('chrom'), 'start': data.get('start'),
                     'end': data.get('end'), 'bait': bool(data.get('bait')), 'degree': graph.degree[node],
                     'community': data.get('community', -1), 'genes': data.get('genes', '')})
    return pd.DataFrame(rows)


def edge_table(graph):
    rows = []
    for u, v, data in graph.edges(data=True):
        row = {'node_a': u, 'node_b': v}
        row.update(data)
        rows.append(row)
    return pd.DataFrame(rows)


def write_graph(graph, path):
    """GraphML export; booleans are written as integers for tool compatibility."""
    out = graph.copy()
    for _, data in out.nodes(data=True):
        for key, value in list(data.items()):
            if isinstance(value, (bool, np.bool_)):
                data[key] = int(value)
            elif isinstance(value, np.generic):
                data[key] = value.item()
    for _, _, data in out.edges(data=True):
        for key, value in list(data.items()):
            if isinstance(value, np.generic):
                data[key] = value.item()
    nx.write_graphml(out, path)


def read_graph(path):
    graph = nx.read_graphml(path)
    for _, data in graph.nodes(data=True):
        if 'bait' in data:
            data['bait'] = bool(int(data['bait']))
    return graph


def neighbourhood_layers(graph, source, max_steps):
    """Nodes at exactly 1..max_steps hops from `source` (list of sets)."""
    lengths = nx.single_source_shortest_path_length(graph, source, cutoff=max_steps)
    layers = defaultdict(set)
    for node, length in lengths.items():
        if length > 0:
            layers[length].add(node)
    return [layers[k] for k in range(1, max_steps + 1)]
