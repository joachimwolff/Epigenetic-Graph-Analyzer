"""egaDifferential: differential interaction graph with chromatin and expression changes."""

import json
import logging

import numpy as np
import pandas as pd
from scipy import stats

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .annotation import bigwig_node_signal, gene_values_to_nodes, parse_track_arguments
from .differential import (REPLICATE_MODES, ActivityOperators, annotate_chicdiff, classify_genes, consensus_scores,
                           differential_edges_from_scores, feature_correlations, gene_features,
                           community_expression_changes, community_transitions, gene_community_changes,
                           gene_level_indirect_test,
                           indirect_only_genes, match_communities, mirror_comparison, mirrored_features,
                           nested_models, permutation_test, union_graph)
from .egaCreateGraph import load_tss
from .graph import annotate_genes, build_graph, detect_communities, write_graph
from .io import INTERACTION_FORMATS, read_baits, read_chicdiff, read_expression, read_interactions

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Compare the interaction graphs of two conditions: gained, lost and shared interactions '
                         '(optionally tested by Chicdiff), contact-weighted chromatin changes at promoter, direct '
                         'and indirect (two-hop) level, and their relation to expression changes.')
    parser.add_argument('--interactionsA', nargs='+', required=True,
                        help='Called interactions of condition A, one file per replicate (e.g. CHiCAGO ibed)')
    parser.add_argument('--interactionsB', nargs='+', required=True,
                        help='Called interactions of condition B, one file per replicate')
    parser.add_argument('--format', '-f', choices=INTERACTION_FORMATS, default='ibed')
    parser.add_argument('--baits',
                        help="Baited fragments of the capture design (BED or CHiCAGO .baitmap); required for "
                             "--format cool, for other formats it sets the bait status and orientation of each pair")
    parser.add_argument('--coolField', default='count', help="Pixel column read from cool matrices")
    parser.add_argument('--maxFdr', type=float, default=None,
                        help='Keep pairs with an FDR at most this value (inputs with an FDR column, e.g. HOMER)')
    parser.add_argument('--scoreThreshold', type=float, default=None,
                        help='Keep pairs with score >= threshold in both conditions (e.g. CHiCAGO 5)')
    parser.add_argument('--replicates', choices=REPLICATE_MODES, default='all',
                        help="'all': pairs present in every replicate and reaching --scoreThreshold in at least one; "
                             "'intersect': reaching it in every replicate; 'union': in any replicate")
    parser.add_argument('--referenceScore', type=float, default=5.0,
                        help='Score that maps to edge strength 1 (the calling threshold, CHiCAGO 5)')
    parser.add_argument('--chicdiff', help='Chicdiff results table (fwrite of <prefix>_results.Rds) testing the '
                                           'change of each interaction')
    parser.add_argument('--chicdiffPadj', type=float, default=0.05)
    parser.add_argument('--names', nargs=2, default=['A', 'B'], help='Condition names')
    parser.add_argument('--tss', required=True, help='GTF or TSS BED for bait gene annotation')
    parser.add_argument('--geneTypes', nargs='*', default=['protein_coding', 'lincRNA', 'antisense'])
    parser.add_argument('--tracksA', nargs='+', required=True, help='mark=file.bw[,rep.bw] for condition A')
    parser.add_argument('--tracksB', nargs='+', required=True, help='mark=file.bw[,rep.bw] for condition B')
    parser.add_argument('--markSign', nargs='*', default=['H3K27me3=-1', 'H3K9me3=-1'],
                        help='mark=+1|-1, the expected direction of a mark on expression (default +1)')
    parser.add_argument('--expression', required=True,
                        help='Gene table: either a log2 fold change column (--log2FCColumn) or two expression '
                             'columns (--expressionColumns)')
    parser.add_argument('--log2FCColumn', help='Column holding log2(B/A), e.g. DESeq2 log2FoldChange')
    parser.add_argument('--padjColumn', help='Adjusted p-value column; restricts differential genes to padj <= --alpha')
    parser.add_argument('--alpha', type=float, default=0.05)
    parser.add_argument('--expressionColumns', nargs=2, help='Expression columns of A and B (e.g. RPKM)')
    parser.add_argument('--minExpression', type=float, default=1.0,
                        help='With --expressionColumns: keep genes expressed at least this much in A or B')
    parser.add_argument('--expressionPseudoCount', type=float, default=1.0)
    parser.add_argument('--geneKey', choices=('gene_ids', 'genes'), default='gene_ids')
    parser.add_argument('--rnaThreshold', type=float, default=1.0, help='|log2FC| of a differential gene')
    parser.add_argument('--chromatinThreshold', type=float, default=1.0,
                        help='|log2FC| of a chromatin feature counted as explaining a gene')
    parser.add_argument('--permutations', type=int, default=100)
    parser.add_argument('--genePermutations', type=int, default=1000,
                        help='Permutations of the per-gene test of the indirect level (0 skips it)')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Bootstrap resamples of the mirror comparison')
    parser.add_argument('--resolution', type=float, default=1.0, help='Community resolution parameter')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--outputPrefix', '-o', required=True)
    return parser.parse_args(args)


def expression_change(table, args):
    if args.log2FCColumn:
        return table[args.log2FCColumn].astype(float)
    if not args.expressionColumns:
        raise ValueError('Give --log2FCColumn or --expressionColumns')
    a, b = args.expressionColumns
    values = table[[a, b]].astype(float)
    values = values[(values[a] >= args.minExpression) | (values[b] >= args.minExpression)]
    return np.log2((values[b] + args.expressionPseudoCount) / (values[a] + args.expressionPseudoCount))


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    name_a, name_b = args.names
    prefix = args.outputPrefix
    baits = read_baits(args.baits) if args.baits else None

    def condition(paths):
        tables = []
        for path in paths:
            table = read_interactions([path], args.format, baits=baits, cool_field=args.coolField)
            if args.maxFdr is not None and 'fdr' in table:
                table = table[table['fdr'] <= args.maxFdr]
            tables.append(table)
        if len(tables) == 1 and args.scoreThreshold is not None:
            return consensus_scores([tables[0][tables[0]['score'] >= args.scoreThreshold]], 'union')
        return consensus_scores(tables, args.replicates, args.scoreThreshold)

    table_a = condition(args.interactionsA)
    table_b = condition(args.interactionsB)
    edges = differential_edges_from_scores(table_a, table_b, name_a, name_b, args.referenceScore)
    if args.chicdiff:
        edges = annotate_chicdiff(edges, read_chicdiff(args.chicdiff), args.chicdiffPadj)
    community_weight = 'score'
    log.info('%d union edges: %s', len(edges), edges['status'].value_counts().to_dict())
    write_table(edges, prefix + '_edges.tsv.gz')

    tss = load_tss(args.tss, args.geneTypes)
    graph = union_graph(edges, name_a, name_b)
    nodes = pd.DataFrame([(n, d['chrom'], d['start'], d['end']) for n, d in graph.nodes(data=True)],
                         columns=['node', 'chrom', 'start', 'end'])
    genes = annotate_genes(nodes, tss)
    for node in graph.nodes:
        names, ids = genes.get(node, ([], []))
        graph.nodes[node]['genes'] = ','.join(names)
        graph.nodes[node]['gene_ids'] = ','.join(ids)

    tracks_a = parse_track_arguments(args.tracksA)
    tracks_b = parse_track_arguments(args.tracksB)
    if set(tracks_a) != set(tracks_b):
        raise ValueError('Conditions need the same marks: {} vs {}'.format(sorted(tracks_a), sorted(tracks_b)))
    marks = {}
    for mark in tracks_a:
        signal_a = [bigwig_node_signal(graph, path) for path in tracks_a[mark]]
        signal_b = [bigwig_node_signal(graph, path) for path in tracks_b[mark]]
        mean_a = {n: float(np.mean([s[n] for s in signal_a])) for n in graph.nodes}
        mean_b = {n: float(np.mean([s[n] for s in signal_b])) for n in graph.nodes}
        marks[mark] = (mean_a, mean_b)
        for n in graph.nodes:
            graph.nodes[n][mark + '_' + name_a] = mean_a[n]
            graph.nodes[n][mark + '_' + name_b] = mean_b[n]
    mark_signs = {m: 1 for m in marks}
    for entry in args.markSign:
        mark, sign = entry.split('=')
        if mark in mark_signs:
            mark_signs[mark] = int(sign)

    table = read_expression(args.expression)
    differential_ids = None
    if args.padjColumn:
        differential_ids = set(table.index[table[args.padjColumn].astype(float) <= args.alpha])
    change = gene_values_to_nodes(graph, expression_change(table, args), key=args.geneKey)
    change = {n: v for n, v in change.items() if graph.nodes[n].get('bait')}
    for node, value in change.items():
        graph.nodes[node]['rna_log2FC'] = value
    log.info('%d bait nodes with expression change', len(change))

    operators = ActivityOperators(graph, list(change))
    mark_vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=mark_vectors)
    mirror = mirrored_features(graph, operators, tracks_a, tracks_b)
    features = features.join(mirror, on='node')
    write_table(features, prefix + '_gene_features.tsv')
    # genes with direct and indirect neighbours carry information on every level
    informative = features[(features['n_direct'] > 0) & (features['n_indirect'] > 0)].reset_index(drop=True)
    write_table(feature_correlations(informative, list(marks)), prefix + '_correlations.tsv')
    write_table(nested_models(informative, list(marks)), prefix + '_nested_models.tsv')
    comparison = mirror_comparison(informative, list(marks), n_bootstrap=args.bootstrap, seed=args.seed)
    differential_nodes = None
    if differential_ids is not None:
        differential_nodes = {n for n in features['node']
                              if set(str(graph.nodes[n].get(args.geneKey, '')).split(',')) & differential_ids}
    classified = features if differential_nodes is None else features[features['node'].isin(differential_nodes)]
    classes = classify_genes(classified, list(marks), mark_signs, args.rnaThreshold, args.chromatinThreshold)
    write_table(classes, prefix + '_gene_classes.tsv')
    indirect_only = indirect_only_genes(graph, operators, mark_vectors, features, list(marks), mark_signs,
                                        args.rnaThreshold, args.chromatinThreshold, genes=differential_nodes)
    if args.genePermutations:
        tested = classified if differential_nodes is not None else features
        gene_test = gene_level_indirect_test(operators, mark_vectors, tested, list(marks), mark_signs,
                                             args.rnaThreshold, args.genePermutations, args.seed)
        write_table(gene_test, prefix + '_indirect_gene_test.tsv')
        if len(indirect_only):
            indirect_only = indirect_only.merge(gene_test, on='node', how='left')
    write_table(indirect_only, prefix + '_indirect_only_genes.tsv')
    informative_change = {n: change[n] for n in informative['node']}
    result = permutation_test(graph, marks, informative_change, mark_signs, n_permutations=args.permutations,
                              seed=args.seed,
                              rna_threshold=args.rnaThreshold, chromatin_threshold=args.chromatinThreshold)
    result['class_counts'] = classes['explanation'].value_counts().to_dict()
    result['informative_genes'] = int(len(informative))
    result['mirror_comparison'] = comparison
    with open(ensure_parent(prefix + '_permutation.json'), 'w') as handle:
        json.dump(result, handle, indent=2)

    # communities per condition and their correspondence
    graph_a = build_graph(table_a, tss=None, weight=community_weight)
    graph_b = build_graph(table_b, tss=None, weight=community_weight)
    communities_a = detect_communities(graph_a, resolution=args.resolution, seed=args.seed)
    communities_b = detect_communities(graph_b, resolution=args.resolution, seed=args.seed)
    baits = {n for n, b in graph.nodes(data='bait') if b}
    write_table(match_communities(communities_a, communities_b, baits), prefix + '_community_matching.tsv')
    transitions, main_b = community_transitions(communities_a, communities_b, set(graph_a.nodes), set(graph_b.nodes))
    write_table(transitions, prefix + '_community_transitions.tsv')
    node_genes = {n: graph.nodes[n].get('genes', '') for n in graph.nodes}
    de_nodes = differential_nodes if differential_nodes is not None else \
        {n for n, v in change.items() if abs(v) >= args.rnaThreshold}
    gene_changes, switch_test = gene_community_changes(communities_a, communities_b, transitions, node_genes, baits,
                                                       expression=change, differential=de_nodes)
    write_table(gene_changes, prefix + '_gene_community_changes.tsv')
    community_changes = community_expression_changes(communities_a, communities_b, main_b, transitions, edges, change,
                                                     n_permutations=args.genePermutations or 1000, seed=args.seed)
    write_table(community_changes, prefix + '_community_expression.tsv')
    summary = {'transitions': transitions['transition'].value_counts().to_dict(),
               'gene_status': gene_changes['status'].value_counts().to_dict() if len(gene_changes) else {},
               'switching_vs_differential': switch_test,
               'communities_tested': int(len(community_changes)),
               'coherent_communities_q_lt_0.1': int((community_changes['coherence_q'] < 0.1).sum())
               if len(community_changes) else 0}
    if len(community_changes) > 3:
        rho, p = stats.spearmanr(community_changes['rewired_fraction'], community_changes['mean_rna_log2FC'].abs(),
                                 nan_policy='omit')
        summary['rewiring_vs_abs_expression_spearman'] = {'rho': float(rho), 'p': float(p)}
    with open(ensure_parent(prefix + '_community_summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=2)
    write_graph(graph, ensure_parent(prefix + '_union.graphml'))


if __name__ == '__main__':
    main()
