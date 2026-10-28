"""Recompute the per-gene permutation test of the indirect level from the outputs of an egaDifferential run.

The union graph (<prefix>_union.graphml) holds the edge strengths of both conditions and the mark signal of every
node, so the contact operators and mark vectors of egaDifferential are rebuilt without reading interactions or
bigWig files again. The tested genes are those of <prefix>_gene_classes.tsv (the differential genes with
|log2FC| >= --rnaThreshold), the test is differential.gene_level_indirect_test with the given number of permutations
and seed (egaDifferential uses --genePermutations 1000, --seed 0). Writes

    <out>_indirect_gene_test.tsv     node, genes, statistic, P, BH q (across the tested genes)
    <out>_indirect_only_genes.tsv    <prefix>_indirect_only_genes.tsv with the new P and q
    <out>_indirect_gene_test.json    counts at P < 0.05, q < 0.1, q < 0.05 and of the candidates
"""

import argparse
import json
import logging
import time

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.differential import ActivityOperators, default_threads, gene_level_indirect_test

log = logging.getLogger('indirect_gene_test')


def load(prefix, names):
    """Union graph, marks (name -> (signal A, signal B)) and expression change of an egaDifferential run."""
    graph = nx.read_graphml(prefix + '_union.graphml')
    name_a, name_b = names
    attributes = set()
    for _, data in graph.nodes(data=True):
        attributes.update(data)
    marks = [a[:-len(name_a) - 1] for a in sorted(attributes)
             if a.endswith('_' + name_a) and a[:-len(name_a) - 1] + '_' + name_b in attributes]
    features = pd.read_csv(prefix + '_gene_features.tsv', sep='\t', usecols=['node', 'rna_log2FC'])
    change = dict(zip(features['node'], features['rna_log2FC']))
    signals = {}
    for mark in marks:
        signals[mark] = ({n: float(d.get(mark + '_' + name_a, 0.0)) for n, d in graph.nodes(data=True)},
                         {n: float(d.get(mark + '_' + name_b, 0.0)) for n, d in graph.nodes(data=True)})
    return graph, marks, signals, change


def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prefix', required=True, help='Output prefix of the egaDifferential run')
    parser.add_argument('--names', nargs=2, required=True, help='Condition names of the run')
    parser.add_argument('--markSign', nargs='*', default=['H3K27me3=-1', 'H3K9me3=-1'])
    parser.add_argument('--rnaThreshold', type=float, default=1.0)
    parser.add_argument('--permutations', type=int, default=20000)
    parser.add_argument('--threads', type=int, default=default_threads(),
                        help='Worker processes (default: the available cores, at most 32)')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--output', '-o', required=True, help='Output prefix')
    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')

    start = time.time()
    graph, marks, signals, change = load(args.prefix, args.names)
    signs = {m: 1 for m in marks}
    for entry in args.markSign:
        mark, sign = entry.split('=')
        if mark in signs:
            signs[mark] = int(sign)
    log.info('%d nodes, %d edges, marks %s, signs %s', graph.number_of_nodes(), graph.number_of_edges(), marks, signs)
    operators = ActivityOperators(graph, list(change))
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in signals.items()}
    tested = pd.read_csv(args.prefix + '_gene_classes.tsv', sep='\t', usecols=['node', 'genes', 'rna_log2FC'])
    loaded = time.time()
    log.info('%d differential genes, %d permutations', len(tested), args.permutations)
    test = gene_level_indirect_test(operators, vectors, tested, marks, signs, args.rnaThreshold,
                                    args.permutations, args.seed, threads=args.threads)
    finished = time.time()
    test = test.merge(tested[['node', 'genes']], on='node', how='left')
    test = test[['node', 'genes', 'indirect_statistic', 'indirect_p', 'indirect_q']]
    test.to_csv(args.output + '_indirect_gene_test.tsv', sep='\t', index=False)

    indirect_only = pd.read_csv(args.prefix + '_indirect_only_genes.tsv', sep='\t')
    previous = indirect_only[['node', 'indirect_p', 'indirect_q']].rename(
        columns={'indirect_p': 'indirect_p_previous', 'indirect_q': 'indirect_q_previous'})
    indirect_only = indirect_only.drop(columns=['indirect_statistic', 'indirect_p', 'indirect_q'])
    indirect_only = indirect_only.merge(test.drop(columns=['genes']), on='node', how='left')
    indirect_only = indirect_only.merge(previous, on='node', how='left')
    indirect_only.to_csv(args.output + '_indirect_only_genes.tsv', sep='\t', index=False)

    candidates = indirect_only[indirect_only['beyond_mirror'] & (indirect_only['indirect_p'] < 0.05)]
    summary = {
        'prefix': args.prefix, 'marks': marks, 'mark_signs': signs, 'permutations': args.permutations,
        'seed': args.seed, 'p_floor': 1 / (args.permutations + 1),
        'tested_genes': int(len(test)),
        'genes_p_lt_0.05': int((test['indirect_p'] < 0.05).sum()),
        'genes_q_lt_0.1': int((test['indirect_q'] < 0.1).sum()),
        'genes_q_lt_0.05': int((test['indirect_q'] < 0.05).sum()),
        'genes_at_p_floor': int((test['indirect_p'] <= 1 / (args.permutations + 1) + 1e-15).sum()),
        'indirect_only': int(len(indirect_only)),
        'indirect_only_beyond_mirror': int(indirect_only['beyond_mirror'].sum()),
        'candidates_beyond_mirror_p_lt_0.05': int(len(candidates)),
        'candidates_q_lt_0.1': int((candidates['indirect_q'] < 0.1).sum()),
        'candidates_q_lt_0.05': int((candidates['indirect_q'] < 0.05).sum()),
        'candidates': candidates[['node', 'genes', 'rna_log2FC', 'explaining_marks', 'indirect_statistic',
                                  'indirect_p', 'indirect_q', 'indirect_p_previous',
                                  'indirect_q_previous']].to_dict(orient='records'),
        'seconds_loading': round(loaded - start, 1), 'seconds_permutations': round(finished - loaded, 1),
    }
    with open(args.output + '_indirect_gene_test.json', 'w') as handle:
        json.dump(summary, handle, indent=2)
    log.info('%s', {k: v for k, v in summary.items() if k != 'candidates'})


if __name__ == '__main__':
    main()
