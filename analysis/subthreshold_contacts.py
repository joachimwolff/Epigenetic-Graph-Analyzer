"""Sub-threshold direct contacts between candidate promoters and their explaining two-hop nodes.

Candidates are the differential genes explained at indirect neighbours only, beyond the mirrored loci, with a per-gene
permutation P < --maxP (the _indirect_only_genes.tsv of analysis/indirect_gene_test.py). Their explaining nodes are
the two-hop nodes with a positive concordant contribution to an explaining mark (as in
differential.indirect_only_genes; the three largest per mark are the listed drivers). For every explaining node the
largest score of a direct pair with the candidate promoter is looked up in each replicate of the interaction files
given here, which may hold pairs below the calling threshold (e.g. CHiCAGO exports at score >= 1).

Reclassification: the pairs of the candidate promoters that pass the replicate rule of the run at --scoreThreshold
(the sub-threshold level) and are not in the union graph are added as direct edges with strength score /
--referenceScore per condition; new fragments get their mark signal from the bigWig files. Promoter, direct and
indirect features and the mirrored control of the candidates are recomputed on this graph and the candidates are
classified again (indirect only, beyond the mirror); the per-gene test is rerun on the same graph with
--permutations. A candidate remains if it is still indirect only, beyond the mirror and has P < --maxP.

Writes <out>_explaining_nodes.tsv (one row per candidate and explaining node), <out>_candidates.tsv (one row per
candidate) and <out>_summary.json.
"""

import argparse
import json
import logging
import os
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.annotation import bigwig_node_signal, parse_track_arguments
from epigeneticgraphanalyzer.differential import (ActivityOperators, consensus_scores, differential_edges_from_scores,
                                                   gene_features, gene_level_indirect_test, indirect_only_genes,
                                                   mirrored_features)
from epigeneticgraphanalyzer.io import read_baits, read_interactions

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indirect_gene_test import load  # noqa: E402

log = logging.getLogger('subthreshold_contacts')


def node_ids(table):
    a = table['chrom'] + ':' + table['bait_start'].astype(np.int64).astype(str) + '-' + \
        table['bait_end'].astype(np.int64).astype(str)
    b = table['chrom'] + ':' + table['other_start'].astype(np.int64).astype(str) + '-' + \
        table['other_end'].astype(np.int64).astype(str)
    return a, b


def explaining_nodes(graph, operators, vectors, candidate, marks, signs, top=3):
    """Two-hop nodes with positive concordant contribution to an explaining mark of a candidate."""
    i = operators.sources.index(candidate['node'])
    direction = np.sign(candidate['rna_log2FC'])
    rows = {}
    for mark in str(candidate['explaining_marks']).split(','):
        vector_a, vector_b = vectors[mark]
        row_a = operators.indirect_a.getrow(i)
        row_b = operators.indirect_b.getrow(i)
        columns = np.union1d(row_a.indices, row_b.indices)
        weight_a = np.asarray(row_a[:, columns].todense()).ravel()
        weight_b = np.asarray(row_b[:, columns].todense()).ravel()
        contribution = (weight_b * vector_b[columns] - weight_a * vector_a[columns]) * signs.get(mark, 1) * direction
        positive = contribution > 0
        total = contribution[positive].sum()
        order = np.argsort(contribution)[::-1]
        rank = np.empty(len(order), dtype=int)
        rank[order] = np.arange(1, len(order) + 1)
        for k in np.nonzero(positive)[0]:
            node = operators.nodes[columns[k]]
            entry = rows.setdefault(node, {'node': node, 'marks': [], 'share': 0.0, 'driver': False, 'best_rank': 10 ** 9})
            entry['marks'].append(mark)
            entry['share'] = max(entry['share'], contribution[k] / total)
            entry['best_rank'] = min(entry['best_rank'], int(rank[k]))
            entry['driver'] |= rank[k] <= top
    for entry in rows.values():
        mediators = set(graph.neighbors(candidate['node'])) & set(graph.neighbors(entry['node']))
        entry['mediators'] = len(mediators)
        entry['genes'] = graph.nodes[entry['node']].get('genes', '')
        entry['bait'] = bool(graph.nodes[entry['node']].get('bait'))
        entry['marks'] = ','.join(entry['marks'])
    return list(rows.values())


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prefix', required=True, help='Output prefix of the egaDifferential run')
    parser.add_argument('--candidates', required=True, help='_indirect_only_genes.tsv of indirect_gene_test.py')
    parser.add_argument('--names', nargs=2, required=True)
    parser.add_argument('--markSign', nargs='*', default=['H3K27me3=-1', 'H3K9me3=-1'])
    parser.add_argument('--rnaThreshold', type=float, default=1.0)
    parser.add_argument('--chromatinThreshold', type=float, default=1.0)
    parser.add_argument('--maxP', type=float, default=0.05)
    parser.add_argument('--interactionsA', nargs='+', required=True)
    parser.add_argument('--interactionsB', nargs='+', required=True)
    parser.add_argument('--format', default='ibed')
    parser.add_argument('--baits')
    parser.add_argument('--replicates', default='all', help='Replicate rule of the sub-threshold edges')
    parser.add_argument('--scoreThreshold', type=float, default=None, help='Sub-threshold level (e.g. 2)')
    parser.add_argument('--referenceScore', type=float, default=5.0)
    parser.add_argument('--bins', nargs='+', type=float, default=[1, 2, 5],
                        help='Score bin edges of the report (pairs >= the last edge counted together)')
    parser.add_argument('--tracksA', nargs='+', required=True)
    parser.add_argument('--tracksB', nargs='+', required=True)
    parser.add_argument('--permutations', type=int, default=20000)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--output', '-o', required=True)
    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
    start = time.time()

    graph, marks, signals, change = load(args.prefix, args.names)
    signs = {m: 1 for m in marks}
    for entry in args.markSign:
        mark, sign = entry.split('=')
        if mark in signs:
            signs[mark] = int(sign)
    table = pd.read_csv(args.candidates, sep='\t')
    candidates = table[table['beyond_mirror'] & (table['indirect_p'] < args.maxP)].reset_index(drop=True)
    candidate_nodes = list(candidates['node'])
    log.info('%d candidates', len(candidates))
    operators = ActivityOperators(graph, list(change))
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in signals.items()}

    # replicate tables, restricted to pairs of a candidate promoter
    baits = read_baits(args.baits) if args.baits else None
    wanted = set(candidate_nodes)
    replicates, score_range = {}, {}
    for condition, paths in (('A', args.interactionsA), ('B', args.interactionsB)):
        replicates[condition] = []
        for path in paths:
            t = read_interactions([path], args.format, baits=baits)
            score_range[path] = [float(t['score'].min()), float(t['score'].max()), int(len(t))]
            a, b = node_ids(t)
            hit = (a.isin(wanted) | b.isin(wanted)).to_numpy()
            t = t[hit].assign(node_a=a[hit].to_numpy(), node_b=b[hit].to_numpy())
            replicates[condition].append(t)
            log.info('%s: score %.2f to %.2f, %d pairs, %d at candidate promoters', path, *score_range[path], len(t))

    def pair_scores(t):
        hit = t[t['node_a'].isin(wanted) | t['node_b'].isin(wanted)]
        keys = [tuple(sorted(p)) for p in zip(hit['node_a'], hit['node_b'])]
        return pd.Series(hit['score'].to_numpy(), index=pd.MultiIndex.from_tuples(keys)).groupby(level=[0, 1]).max() \
            if len(keys) else pd.Series(dtype=float)

    scores = {c: [pair_scores(t) for t in tables] for c, tables in replicates.items()}

    # explaining nodes and their direct scores with the candidate promoter
    rows = []
    for _, candidate in candidates.iterrows():
        for entry in explaining_nodes(graph, operators, vectors, candidate, marks, signs):
            key = tuple(sorted((candidate['node'], entry['node'])))
            row = {'candidate': candidate['node'], 'candidate_genes': candidate['genes'], **entry}
            best = 0.0
            for condition, name in (('A', args.names[0]), ('B', args.names[1])):
                values = [float(s.get(key, 0.0)) for s in scores[condition]]
                row['score_' + name] = ','.join('{:.2f}'.format(v) for v in values)
                row['max_score_' + name] = max(values)
                best = max(best, max(values))
            row['max_score'] = best
            rows.append(row)
    nodes_table = pd.DataFrame(rows)
    edges = list(args.bins) + [np.inf]

    def bin_of(score):
        if score <= 0:
            return 'none'
        if score < edges[0]:
            return '<{:g}'.format(edges[0])
        for lo, hi in zip(edges[:-1], edges[1:]):
            if lo <= score < hi:
                return '>={:g}'.format(lo) if np.isinf(hi) else '{:g}-{:g}'.format(lo, hi)
        return 'none'

    if len(nodes_table):
        nodes_table['score_bin'] = nodes_table['max_score'].map(bin_of)
    nodes_table.to_csv(args.output + '_explaining_nodes.tsv', sep='\t', index=False)

    # sub-threshold edges of the candidate promoters added to the union graph
    def condition(tables):
        if len(tables) == 1 and args.scoreThreshold is not None:
            return consensus_scores([tables[0][tables[0]['score'] >= args.scoreThreshold]], 'union')
        return consensus_scores(tables, args.replicates, args.scoreThreshold)

    keep = ['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end', 'other_is_bait', 'distance', 'score']
    table_a = condition([t[keep] for t in replicates['A']])
    table_b = condition([t[keep] for t in replicates['B']])
    sub = differential_edges_from_scores(table_a, table_b, args.names[0], args.names[1], args.referenceScore)
    sub = sub[sub['node_a'].isin(wanted) | sub['node_b'].isin(wanted)]
    new = sub[[not graph.has_edge(u, v) for u, v in zip(sub['node_a'], sub['node_b'])]]
    augmented = graph.copy()
    added_nodes = []
    for row in new.itertuples(index=False):
        row = row._asdict()
        for node in (row['node_a'], row['node_b']):
            if node not in augmented:
                chrom, rest = node.split(':')
                s, e = rest.split('-')
                augmented.add_node(node, chrom=chrom, start=int(s), end=int(e), bait=False, genes='')
                added_nodes.append(node)
        augmented.add_edge(row['node_a'], row['node_b'], r_a=float(row['enrichment_' + args.names[0]]),
                           r_b=float(row['enrichment_' + args.names[1]]), status=row['status'])
    tracks_a = parse_track_arguments(args.tracksA)
    tracks_b = parse_track_arguments(args.tracksB)
    if added_nodes:
        sub_graph = augmented.subgraph(added_nodes)
        for mark in marks:
            for tracks, name in ((tracks_a, args.names[0]), (tracks_b, args.names[1])):
                values = [bigwig_node_signal(sub_graph, p) for p in tracks[mark]]
                for node in added_nodes:
                    signals[mark][0 if name == args.names[0] else 1][node] = float(np.mean([v[node] for v in values]))
    log.info('%d sub-threshold edges of candidate promoters added (%d new fragments)', len(new), len(added_nodes))

    candidate_change = {n: change[n] for n in candidate_nodes}

    def classify(g):
        ops = ActivityOperators(g, candidate_nodes)
        vec = {m: (ops.vector(a), ops.vector(b)) for m, (a, b) in signals.items()}
        feats = gene_features(g, signals, candidate_change, operators=ops, mark_vectors=vec)
        feats = feats.join(mirrored_features(g, ops, tracks_a, tracks_b), on='node')
        only = indirect_only_genes(g, ops, vec, feats, marks, signs, args.rnaThreshold, args.chromatinThreshold)
        return feats, only

    # the original graph reproduces the features of the run
    original_features, _ = classify(graph)
    stored = pd.read_csv(args.prefix + '_gene_features.tsv', sep='\t').set_index('node').loc[candidate_nodes]
    columns = [c for c in original_features.columns if c in stored.columns and c not in ('node', 'genes')]
    deviation = float(np.nanmax(np.abs(original_features.set_index('node')[columns].to_numpy(dtype=float) -
                                       stored[columns].to_numpy(dtype=float))))
    log.info('largest deviation of recomputed candidate features from the run: %g', deviation)

    features, only = classify(augmented)
    operators_all = ActivityOperators(augmented, list(change))
    vectors_all = {m: (operators_all.vector(a), operators_all.vector(b)) for m, (a, b) in signals.items()}
    tested = pd.read_csv(args.prefix + '_gene_classes.tsv', sep='\t', usecols=['node', 'rna_log2FC'])
    test = gene_level_indirect_test(operators_all, vectors_all, tested, marks, signs, args.rnaThreshold,
                                    args.permutations, args.seed).set_index('node')

    out = candidates[['node', 'genes', 'rna_log2FC', 'explaining_marks', 'indirect_p', 'indirect_q']].copy()
    by_candidate = nodes_table.groupby('candidate') if len(nodes_table) else None
    lo_sub = args.scoreThreshold if args.scoreThreshold is not None else edges[0]
    records = []
    for _, c in out.iterrows():
        n = nodes_table[nodes_table['candidate'] == c['node']] if by_candidate is not None else nodes_table
        f = features.set_index('node').loc[c['node']]
        o = only[only['node'] == c['node']] if len(only) else only
        record = {
            'explaining_nodes': int(len(n)),
            'explaining_with_any_pair': int((n['max_score'] > 0).sum()),
            'explaining_with_pair_ge_sub': int((n['max_score'] >= lo_sub).sum()),
            'drivers': int(n['driver'].sum()),
            'drivers_with_pair_ge_sub': int((n['driver'] & (n['max_score'] >= lo_sub)).sum()),
            'max_pair_score': float(n['max_score'].max()) if len(n) else 0.0,
            'contribution_share_ge_sub': float(n.loc[n['max_score'] >= lo_sub, 'share'].sum()),
            'added_direct_edges': int(((new['node_a'] == c['node']) | (new['node_b'] == c['node'])).sum()),
            'n_direct_after': int(f['n_direct']), 'n_indirect_after': int(f['n_indirect']),
            'indirect_only_after': bool(len(o)),
            'beyond_mirror_after': bool(len(o) and o['beyond_mirror'].iloc[0]),
            'explaining_marks_after': o['explaining_marks'].iloc[0] if len(o) else '',
            'indirect_p_after': float(test.loc[c['node'], 'indirect_p']),
            'indirect_q_after': float(test.loc[c['node'], 'indirect_q']),
        }
        for mark in marks:
            for level in ('promoter', 'direct', 'indirect', 'indirect_mirror'):
                record['{}_{}_after'.format(mark, level)] = float(f['{}_{}'.format(mark, level)])
        record['remains'] = record['indirect_only_after'] and record['beyond_mirror_after'] and \
            record['indirect_p_after'] < args.maxP
        records.append(record)
    out = pd.concat([out.reset_index(drop=True), pd.DataFrame(records)], axis=1)
    out.to_csv(args.output + '_candidates.tsv', sep='\t', index=False)
    summary = {
        'prefix': args.prefix, 'interactionsA': args.interactionsA, 'interactionsB': args.interactionsB,
        'sub_threshold': args.scoreThreshold, 'replicates': args.replicates, 'reference_score': args.referenceScore,
        'score_range_min_max_pairs': score_range,
        'candidates': int(len(out)),
        'candidates_with_explaining_pair_ge_sub': int((out['explaining_with_pair_ge_sub'] > 0).sum()),
        'candidates_with_driver_pair_ge_sub': int((out['drivers_with_pair_ge_sub'] > 0).sum()),
        'candidates_with_any_explaining_pair': int((out['explaining_with_any_pair'] > 0).sum()),
        'explaining_nodes_by_score_bin': nodes_table['score_bin'].value_counts().to_dict() if len(nodes_table) else {},
        'added_edges': int(len(new)), 'added_fragments': int(len(added_nodes)),
        'remain_indirect_only': int(out['indirect_only_after'].sum()),
        'remain_beyond_mirror': int(out['beyond_mirror_after'].sum()),
        'remain_all': int(out['remains'].sum()),
        'feature_reproduction_max_deviation': deviation,
        'permutations': args.permutations, 'seconds': round(time.time() - start, 1),
    }
    with open(args.output + '_summary.json', 'w') as handle:
        json.dump(summary, handle, indent=2)
    log.info('%s', summary)


if __name__ == '__main__':
    main()
