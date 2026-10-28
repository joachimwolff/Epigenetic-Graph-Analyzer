"""Robustness of the communities and of the community-level results to the choices of community detection.

The reference partitions are those of the manuscript: Leiden (RB configuration modularity), resolution 1, edge weight
= CHiCAGO score; seed 42 for the GM12878 graph (02_build_graphs.sh, score >= 5, <= 2 Mb) and seed 0 for the hESC and
hNEC graphs of the differential analysis (04_differential.sh, score >= 12). Each variant changes one choice:
resolution 0.5 or 2, ten other seeds, Louvain instead of Leiden, edge weight log(score) or unweighted, and the score
threshold (7 for GM12878, 15 for hESC/hNEC). Every variant is applied to both conditions of the comparison.

GM12878: similarity of each partition to the reference (adjusted Rand index, normalised mutual information over the
nodes present in both, fraction of reference communities with a best match of node Jaccard >= 0.5) and the number and
size of communities.

hESC vs hNEC: the same similarities for both condition graphs; the transitions (community_transitions), the stability
of the transition label of each reference community (label of its best matching variant community); the community
expression test (community_expression_changes, same permutations and seed as egaDifferential) with the number of
coherent communities (BH q < 0.1) and the Spearman correlation of the rewired edge fraction with |mean RNA log2FC|;
and whether the coherent communities of the reference are recovered: a reference community counts as recovered when
the tested variant community sharing most genes with it (gene Jaccard) is coherent (q < 0.1), and as contained when
a coherent variant community holds at least half of its genes.

Ties in community_transitions (main B community) and community_expression_changes are broken by the iteration order
of node sets, which follows Python's string hashing; the script pins PYTHONHASHSEED=2 (reproduces the stored
egaDifferential results) and reruns the hESC/hNEC reference under PYTHONHASHSEED 0..9 to measure this effect.

Usage: python community_robustness.py <data_dir> <results_dir> <output_dir> [--permutations 10000] [--hashSeeds 10]
"""

import argparse
import json
import logging
import os
import resource
import subprocess
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from epigeneticgraphanalyzer.differential import (community_expression_changes, community_transitions,
                                                  consensus_scores, differential_edges_from_scores)
from epigeneticgraphanalyzer.graph import build_graph, compare_partitions, detect_communities, read_graph
from epigeneticgraphanalyzer.io import read_baits, read_interactions

log = logging.getLogger('community_robustness')
TRANSITIONS = ['stable', 'split', 'merged', 'remodelled', 'dissolved', 'formed']
HASH_COLUMNS = ['communities_tested', 'coherent_q_lt_0_1', 'coherent_as_stored', 'rewiring_abs_expression_rho',
                'rewiring_abs_expression_p'] + ['n_' + t for t in TRANSITIONS]


def variants(reference_seed, threshold, alternative_threshold, seeds):
    base = {'method': 'leiden', 'resolution': 1.0, 'seed': reference_seed, 'weight': 'score', 'threshold': threshold}
    out = [dict(base, name='reference', group='reference')]
    for resolution in (0.5, 2.0):
        out.append(dict(base, name='resolution_{:g}'.format(resolution), group='resolution', resolution=resolution))
    for seed in seeds:
        out.append(dict(base, name='seed_{}'.format(seed), group='seed', seed=seed))
    out.append(dict(base, name='louvain', group='method', method='louvain'))
    out.append(dict(base, name='weight_log_score', group='weight', weight='log_score'))
    out.append(dict(base, name='weight_unweighted', group='weight', weight='unit'))
    out.append(dict(base, name='threshold_{:g}'.format(alternative_threshold), group='threshold',
                    threshold=alternative_threshold))
    return out


def add_weights(graph):
    for _, _, data in graph.edges(data=True):
        data['log_score'] = float(np.log(data['score']))
        data['unit'] = 1.0
    return graph


def thresholded(graph, threshold):
    return graph.edge_subgraph([(u, v) for u, v, s in graph.edges(data='score') if s >= threshold]).copy()


def size_stats(communities, graph):
    """Number and size of communities, and how many of them are whole connected components of the graph."""
    sizes = np.array([len(c) for c in communities])
    components = {frozenset(c) for c in nx.connected_components(graph)}
    whole = [c for c in communities if frozenset(c) in components]
    return {'n_communities': int(len(sizes)), 'n_communities_ge5': int((sizes >= 5).sum()),
            'median_size': float(np.median(sizes)), 'mean_size': float(sizes.mean()), 'max_size': int(sizes.max()),
            'n_components': len(components), 'n_communities_whole_component': len(whole),
            'fraction_nodes_in_whole_component_communities': sum(len(c) for c in whole) / graph.number_of_nodes()}


def communities_of(graph, variant):
    return detect_communities(graph, method=variant['method'], resolution=variant['resolution'],
                              seed=variant['seed'], weight=variant['weight'])


def gm12878(results, output):
    graph = add_weights(read_graph(results + '/graphs/GM_chicago.graphml'))
    stored = pd.Series(dict(graph.nodes(data='community')))
    reference = [set(group.index) for _, group in stored.groupby(stored)]
    rows = []
    for variant in variants(42, 5.0, 7.0, range(10)):
        start = time.time()
        g = graph if variant['threshold'] == 5.0 else thresholded(graph, variant['threshold'])
        communities = communities_of(g, variant)
        similarity, _ = compare_partitions(reference, communities)
        row = dict(variant, n_nodes=g.number_of_nodes(), n_edges=g.number_of_edges(), **size_stats(communities, g),
                   **similarity, seconds=round(time.time() - start, 1))
        rows.append(row)
        log.info('GM %s: %s', variant['name'], {k: row[k] for k in ('n_communities', 'ari', 'nmi', 'fraction_matched')})
    table = pd.DataFrame(rows)
    if table.loc[0, 'ari'] < 1.0:
        log.warning('GM12878: recomputed reference differs from the stored communities (ARI %.4f)', table.loc[0, 'ari'])
    table.to_csv(output + '/gm12878_variants.tsv', sep='\t', index=False)
    return table


def unit_genes(row, communities_a, communities_b, expression, node_genes):
    nodes = set(communities_a[row.community_a])
    if row.community_b >= 0:
        nodes |= communities_b[row.community_b]
    return {g for n in nodes if n in expression for g in node_genes.get(n, '').split(',') if g}


def hesc_hnec(data, results, output, permutations, reference_only=False):
    """All variants of the hESC vs hNEC comparison; with `reference_only` the reference row only, nothing written."""
    study = data + '/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC'
    baits = read_baits(data + '/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap')
    raw_a = read_interactions([study + '/GSE86821_hESC12_merged_washU.txt.gz'], 'washu', baits=baits)
    raw_b = read_interactions([study + '/GSE86821_hNPC12_merged_washU.txt.gz'], 'washu', baits=baits)
    # bait-level RNA log2FC and gene names of the union graph written by egaDifferential
    union = read_graph(results + '/differential/hESC_vs_hNEC_union.graphml')
    expression_all = {n: float(v) for n, v in union.nodes(data='rna_log2FC') if v is not None}
    node_genes = {n: g or '' for n, g in union.nodes(data='genes')}
    del union
    log.info('%d bait nodes with expression change', len(expression_all))

    inputs = {}
    for threshold in (12.0, 15.0):
        table_a = consensus_scores([raw_a[raw_a['score'] >= threshold]], 'union')
        table_b = consensus_scores([raw_b[raw_b['score'] >= threshold]], 'union')
        edges = differential_edges_from_scores(table_a, table_b, 'hESC', 'hNEC', 12.0)
        graph_a = add_weights(build_graph(table_a, tss=None, weight='score'))
        graph_b = add_weights(build_graph(table_b, tss=None, weight='score'))
        nodes = set(graph_a.nodes) | set(graph_b.nodes)
        expression = {n: v for n, v in expression_all.items() if n in nodes}
        inputs[threshold] = (edges, graph_a, graph_b, expression)
        log.info('threshold %g: %d union edges, hESC %d nodes, hNEC %d nodes, %d genes with expression', threshold,
                 len(edges), graph_a.number_of_nodes(), graph_b.number_of_nodes(), len(expression))

    stored = pd.read_csv(results + '/differential/hESC_vs_hNEC_community_expression.tsv', sep='\t')
    rows, label_rows, recovery_rows = [], [], []
    reference = None
    for variant in variants(0, 12.0, 15.0, range(1, 11))[:1 if reference_only else None]:
        start = time.time()
        edges, graph_a, graph_b, expression = inputs[variant['threshold']]
        communities_a = communities_of(graph_a, variant)
        communities_b = communities_of(graph_b, variant)
        transitions, main = community_transitions(communities_a, communities_b, set(graph_a.nodes),
                                                  set(graph_b.nodes))
        changes = community_expression_changes(communities_a, communities_b, main, transitions, edges, expression,
                                               n_permutations=permutations, seed=0)
        if not reference_only:
            changes.to_csv('{}/hesc_hnec_community_expression_{}.tsv'.format(output, variant['name']), sep='\t',
                           index=False)
        rho, p = stats.spearmanr(changes['rewired_fraction'], changes['mean_rna_log2FC'].abs(), nan_policy='omit')
        coherent = changes[changes['coherence_q'] < 0.1]
        labels_a = transitions[transitions['side'] == 'A'].set_index('community_a')['transition']
        formed_b = set(transitions.loc[transitions['side'] == 'B', 'community_b'])
        genes = [unit_genes(r, communities_a, communities_b, expression, node_genes) for r in changes.itertuples()]
        if reference is None:
            reference = {'a': communities_a, 'b': communities_b, 'labels_a': labels_a, 'formed_b': formed_b,
                         'coherent': coherent.assign(genes=[genes[i] for i in coherent.index])}
            check = stored[stored['coherence_q'] < 0.1]
            same_coherent = sorted(check['communities_a']) == sorted(coherent['communities_a'])
            log.info('reference reproduces the stored coherent communities: %s', same_coherent)
        row = dict(variant)
        if variant['name'] == 'reference':
            row['coherent_as_stored'] = same_coherent
        for side, communities, ref, g in (('hESC', communities_a, reference['a'], graph_a),
                                          ('hNEC', communities_b, reference['b'], graph_b)):
            similarity, matching = compare_partitions(ref, communities)
            row.update({side + '_' + k: v for k, v in size_stats(communities, g).items()})
            row.update({side + '_' + k: v for k, v in similarity.items()})
            if side == 'hESC':
                matching_a = matching
            else:
                matching_b = matching
        # transition label of every reference community against that of its best matching variant community
        kept = {}
        for label in TRANSITIONS:
            if label == 'formed':
                ref_ids = sorted(reference['formed_b'])
                match = matching_b.set_index('community').loc[ref_ids]
                same = match['best_match'].isin(formed_b)
            else:
                ref_ids = list(reference['labels_a'][reference['labels_a'] == label].index)
                match = matching_a.set_index('community').loc[ref_ids]
                same = match['best_match'].map(labels_a).eq(label)
            matched = match['jaccard'] >= 0.5
            kept[label] = (int(same.sum()), len(ref_ids))
            label_rows.append({'variant': variant['name'], 'group': variant['group'], 'transition': label,
                               'n_reference': len(ref_ids), 'n_same_label': int(same.sum()),
                               'fraction_same_label': float(same.mean()) if len(ref_ids) else np.nan,
                               'n_matched_jaccard_ge_0.5': int(matched.sum()),
                               'fraction_same_label_among_matched': float(same[matched].mean())
                               if matched.any() else np.nan})
        row['fraction_same_transition'] = sum(k for k, _ in kept.values()) / sum(n for _, n in kept.values())
        row.update({'n_' + label: int((transitions['transition'] == label).sum()) for label in TRANSITIONS})
        row.update(communities_tested=int(len(changes)), coherent_q_lt_0_1=int(len(coherent)),
                   rewiring_abs_expression_rho=float(rho), rewiring_abs_expression_p=float(p))
        # recovery of the coherent reference communities, matched by their genes
        coherent_genes = [genes[i] for i in coherent.index]
        recovered = contained = 0
        for ref_row in reference['coherent'].itertuples():
            ref_genes = ref_row.genes
            jaccards = [len(ref_genes & g) / len(ref_genes | g) if g else 0.0 for g in genes]
            best = int(np.argmax(jaccards))
            best_row = changes.iloc[best]
            held = max([len(ref_genes & g) / len(ref_genes) for g in coherent_genes], default=0.0)
            is_recovered = bool(best_row['coherence_q'] < 0.1)
            recovered += is_recovered
            contained += held >= 0.5
            recovery_rows.append({
                'variant': variant['name'], 'group': variant['group'],
                'reference_communities_a': ref_row.communities_a, 'reference_n_genes': ref_row.n_genes,
                'reference_mean_rna_log2FC': ref_row.mean_rna_log2FC, 'reference_q': ref_row.coherence_q,
                'reference_genes': ','.join(sorted(ref_genes)),
                'best_match_gene_jaccard': jaccards[best], 'best_match_n_genes': int(best_row['n_genes']),
                'best_match_mean_rna_log2FC': float(best_row['mean_rna_log2FC']),
                'best_match_q': float(best_row['coherence_q']), 'best_match_transition': best_row['transition'],
                'recovered': is_recovered,
                'same_direction': bool(np.sign(best_row['mean_rna_log2FC']) == np.sign(ref_row.mean_rna_log2FC)),
                'max_fraction_in_coherent': held, 'contained': bool(held >= 0.5)})
        row.update(reference_coherent_recovered=recovered, reference_coherent_contained=contained,
                   n_reference_coherent=int(len(reference['coherent'])), seconds=round(time.time() - start, 1))
        rows.append(row)
        log.info('hESC/hNEC %s: %s', variant['name'],
                 {k: row[k] for k in ('hESC_ari', 'hNEC_ari', 'fraction_same_transition', 'coherent_q_lt_0_1',
                                      'reference_coherent_recovered', 'rewiring_abs_expression_rho',
                                      'rewiring_abs_expression_p', 'seconds')})
    table = pd.DataFrame(rows)
    if reference_only:
        return table, None
    table.to_csv(output + '/hesc_hnec_variants.tsv', sep='\t', index=False)
    pd.DataFrame(label_rows).to_csv(output + '/hesc_hnec_transition_stability.tsv', sep='\t', index=False)
    pd.DataFrame(recovery_rows).to_csv(output + '/hesc_hnec_coherent_recovery.tsv', sep='\t', index=False)
    ref = reference['coherent'].assign(genes=lambda d: d['genes'].map(lambda g: ','.join(sorted(g))))
    ref.to_csv(output + '/hesc_hnec_reference_coherent.tsv', sep='\t', index=False)
    return table, pd.DataFrame(recovery_rows)


def ranges(table, columns):
    """min / median / max of each column for the seed variants and the value for every other variant."""
    out = {}
    seeds = table[table['group'] == 'seed']
    for column in columns:
        out[column] = {'seeds_min': float(seeds[column].min()), 'seeds_median': float(seeds[column].median()),
                       'seeds_max': float(seeds[column].max())}
        out[column].update({name: float(v) for name, v in zip(table['name'], table[column])
                            if name.split('_')[0] != 'seed'})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('data')
    parser.add_argument('results')
    parser.add_argument('output')
    parser.add_argument('--permutations', type=int, default=10000,
                        help='Random gene sets per community size (egaDifferential --communityPermutations)')
    parser.add_argument('--hashSeeds', type=int, default=10,
                        help='Rerun the hESC/hNEC reference under PYTHONHASHSEED 0..N-1 (tie-breaking check)')
    parser.add_argument('--referenceOnly', action='store_true',
                        help='Print the hESC/hNEC reference statistics as JSON and exit (used by the check)')
    args = parser.parse_args()
    if 'PYTHONHASHSEED' not in os.environ:
        # string hashing orders the node sets, which breaks ties in community_transitions and
        # community_expression_changes; pin it so that the run is reproducible (2 reproduces the stored results)
        os.execve(sys.executable, [sys.executable] + sys.argv, dict(os.environ, PYTHONHASHSEED='2'))
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    if args.referenceOnly:
        table, _ = hesc_hnec(args.data, args.results, args.output, args.permutations, reference_only=True)
        print('REFERENCE ' + table.iloc[0][HASH_COLUMNS].to_json())
        return
    start = time.time()
    log.info('PYTHONHASHSEED=%s', os.environ['PYTHONHASHSEED'])
    gm = gm12878(args.results, args.output)
    hesc, recovery = hesc_hnec(args.data, args.results, args.output, args.permutations)
    hash_rows = []
    for hash_seed in range(args.hashSeeds):
        command = [sys.executable, __file__, args.data, args.results, args.output, '--referenceOnly',
                   '--permutations', str(args.permutations)]
        done = subprocess.run(command, env=dict(os.environ, PYTHONHASHSEED=str(hash_seed)), check=True,
                              capture_output=True, text=True)
        line = [x for x in done.stdout.splitlines() if x.startswith('REFERENCE ')][0]
        hash_rows.append(dict(python_hash_seed=hash_seed, **json.loads(line[len('REFERENCE '):])))
        log.info('reference with PYTHONHASHSEED=%d: %s', hash_seed, hash_rows[-1])
    hashes = pd.DataFrame(hash_rows)
    hashes.to_csv(args.output + '/hesc_hnec_reference_hash_seeds.tsv', sep='\t', index=False)
    summary = {
        'python_hash_seed': os.environ['PYTHONHASHSEED'],
        'reference_under_hash_seeds': {c: [float(hashes[c].min()), float(hashes[c].max())]
                                       for c in HASH_COLUMNS if c != 'coherent_as_stored'}
        if len(hashes) else {},
        'reference_hash_seeds_reproducing_stored_coherent': int(hashes['coherent_as_stored'].sum())
        if len(hashes) else None,
        'gm12878': ranges(gm, ['ari', 'nmi', 'fraction_matched', 'n_communities', 'median_size',
                               'n_communities_whole_component', 'fraction_nodes_in_whole_component_communities']),
        'hesc_hnec': ranges(hesc, ['hESC_ari', 'hESC_nmi', 'hESC_fraction_matched', 'hNEC_ari', 'hNEC_nmi',
                                   'hNEC_fraction_matched', 'hESC_n_communities', 'fraction_same_transition',
                                   'coherent_q_lt_0_1', 'reference_coherent_recovered',
                                   'reference_coherent_contained', 'rewiring_abs_expression_rho',
                                   'rewiring_abs_expression_p']),
        'coherent_recovery_per_reference_community': recovery[recovery['variant'] != 'reference'].groupby(
            'reference_communities_a')[['recovered', 'contained']].mean().to_dict(orient='index'),
        'runtime_seconds': round(time.time() - start, 1),
        'peak_memory_gb': round(max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                    resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024 ** 2, 2)}
    with open(args.output + '/summary.json', 'w') as handle:
        json.dump(summary, handle, indent=2)
    log.info('done in %.0f s, peak memory %.2f GB', summary['runtime_seconds'], summary['peak_memory_gb'])


if __name__ == '__main__':
    main()
