"""egaDifferential: differential interaction graph with chromatin and expression changes."""

import json
import logging
import time

import numpy as np
import pandas as pd
from scipy import stats

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .annotation import gene_values_to_nodes, parse_track_arguments
from .differential import (REPLICATE_MODES, ActivityOperators, _benjamini_hochberg, annotate_chicdiff,
                           classify_genes, consensus_scores, differential_edges_from_scores, feature_correlations,
                           gene_features, community_expression_changes, community_transitions,
                           gene_community_changes, gene_level_hop_test, indirect_only_genes, level_name,
                           match_communities, mirror_comparison, mirrored_features, nested_models, permutation_test,
                           _interval_signals, cascade_table, default_threads, explanation_cascade, read_proximity_support,
                           union_graph, explanation_details, explanation_types)
from .states import (MARK_ROLES, bigwig_genome_mean, gene_states, learn_mark_roles, learned_roles,
                     learned_state_transitions, parse_mark_roles, resolve_mark_roles, role_signs,
                     state_transition_table, states_callable)
from .egaCreateGraph import load_tss
from .graph import annotate_genes, build_graph, detect_communities, write_graph
from .io import INTERACTION_FORMATS, read_baits, read_chicdiff, read_expression, read_interactions

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Compare the interaction graphs of two conditions: gained, lost and shared interactions '
                         '(optionally tested by Chicdiff), contact-weighted chromatin changes at promoter, direct '
                         'and indirect (two-hop) level and up to --maxHops graph steps, and their relation to '
                         'expression changes.')
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
                        help='mark=+1|-1|0, the expected direction of a mark on expression (default +1); an alias of '
                             '--markRoles with +1 = active, -1 = repressive, 0 = neutral')
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
    parser.add_argument('--communityPermutations', type=int, default=10000,
                        help='Random gene sets per community size for the coherence test of community expression')
    parser.add_argument('--genePermutations', type=int, default=1000,
                        help='Permutations of the per-gene test of the indirect level (0 skips it)')
    parser.add_argument('--maxHops', type=int, default=3,
                        help='Graph levels explaining a gene: promoter, 1 hop (direct), 2 hops (indirect), ... up to '
                             'this many hops (at least 2; level k >= 3 writes columns <mark>_hop<k>, n_hop<k>, ...)')
    parser.add_argument('--proximitySupport',
                        help='Step curve of analysis/hop_distance.py (TSV with steps, concordance, ci_low, ci_high): '
                             'adds to every explained gene the probability that graph neighbours at its explaining '
                             'level are closer in the nucleus than distance-matched pairs (proximity_support, '
                             'proximity_support_ci; 1.0 / NA at the promoter)')
    parser.add_argument('--stateThreshold', type=float, default=2.0,
                        help='Chromatin states: a mark is present on a fragment if its mean signal over the fragment '
                             'is at least this value (fold enrichment; see --stateScale). Bivalent promoter (bait): an '
                             'active-role and a repressive-role mark present; poised enhancer (other fragment): an '
                             'active- or enhancer-role and a repressive-role mark present')
    parser.add_argument('--stateThresholds', type=float, nargs='*', default=[1.5, 2.0, 3.0],
                        help='Presence thresholds of the state transition report (<prefix>_bivalent_transitions.tsv, '
                             '<prefix>_state_thresholds.tsv); the cascade uses --stateThreshold')
    parser.add_argument('--stateScale', choices=('fold', 'genomeMean'), default='fold',
                        help="'fold': the bigWig signal is a fold enrichment over input (e.g. Roadmap .fc.signal); "
                             "'genomeMean': coverage tracks without input, divided by their genome-wide mean "
                             "(sumData / nBasesCovered of the bigWig header) before the threshold")
    parser.add_argument('--stateMinFragments', type=int, default=50,
                        help='Learned states: mark combinations seen on at least this many fragments (both '
                             'conditions) are kept, the rest pooled')
    parser.add_argument('--markRoles', nargs='*', default=[],
                        help='Declared mark roles, mark=role entries and / or TSV files with the columns mark and role; '
                             'roles: ' + ', '.join(MARK_ROLES) + ' (enhancer: active-type at fragments that are not '
                             'baited promoters; architectural, neutral, unknown: no expected direction). They take '
                             'precedence over --markSign; marks declared nowhere are active')
    parser.add_argument('--roles', choices=('declared', 'learned', 'both'), default='declared',
                        help='Roles used by the explanation cascade and the states: declared (--markRoles / '
                             '--markSign), learned from the data (<prefix>_mark_roles.tsv), or both (declared for the '
                             'main outputs, learned in addition as <prefix>_explanation_cascade_learned.tsv and '
                             'columns *_learned)')
    parser.add_argument('--roleBootstrap', type=int, default=1000,
                        help='Chromosome block-bootstrap resamples of the learned roles and learned states')
    parser.add_argument('--stateTracksA', nargs='+',
                        help='mark=file.bw for the chromatin states of condition A (e.g. fold enrichment over input), '
                             'when the --tracksA signal is not on a fold-enrichment scale; marks not listed use --tracksA')
    parser.add_argument('--stateTracksB', nargs='+', help='as --stateTracksA for condition B')
    parser.add_argument('--stateMinChange', type=float, default=0.5,
                        help='A state transition needs every mark crossing the presence threshold to change by at '
                             'least this log2 fold (pseudo count 1) in the direction of the crossing')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Bootstrap resamples of the mirror comparison')
    parser.add_argument('--resolution', type=float, default=1.0, help='Community resolution parameter')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threads', type=int, default=default_threads(),
                        help='Worker processes for the bigWig reads of the mirrored loci, the per-gene permutation '
                             'test and the community coherence test (default: the available cores, at most 32); '
                             'results do not depend on it')
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
    # mean signal over every node fragment (as annotation.bigwig_node_signal), read in parallel
    node_list = list(graph.nodes)
    state_tracks_a = dict(tracks_a, **parse_track_arguments(args.stateTracksA)) if args.stateTracksA else tracks_a
    state_tracks_b = dict(tracks_b, **parse_track_arguments(args.stateTracksB)) if args.stateTracksB else tracks_b
    paths = list(dict.fromkeys(p for tracks in (tracks_a, tracks_b, state_tracks_a, state_tracks_b)
                               for mark in tracks for p in tracks[mark]))
    started = time.time()
    read = _interval_signals(paths, [graph.nodes[n]['chrom'] for n in node_list],
                             np.array([graph.nodes[n]['start'] for n in node_list], dtype=np.int64),
                             np.array([graph.nodes[n]['end'] for n in node_list], dtype=np.int64), args.threads)
    read = {p: dict(zip(node_list, values.tolist())) for p, values in read.items()}
    log.info('node signals of %d bigWig files: %.0f s', len(paths), time.time() - started)
    for mark in tracks_a:
        signal_a = [read[path] for path in tracks_a[mark]]
        signal_b = [read[path] for path in tracks_b[mark]]
        mean_a = {n: float(np.mean([s[n] for s in signal_a])) for n in graph.nodes}
        mean_b = {n: float(np.mean([s[n] for s in signal_b])) for n in graph.nodes}
        marks[mark] = (mean_a, mean_b)
        for n in graph.nodes:
            graph.nodes[n][mark + '_' + name_a] = mean_a[n]
            graph.nodes[n][mark + '_' + name_b] = mean_b[n]
    # mark roles (declared) and the mark signal on the scale of the presence threshold of the states
    declared = resolve_mark_roles(list(marks), args.markSign, parse_mark_roles(args.markRoles))
    scale = {p: 1.0 for p in paths}
    if args.stateScale == 'genomeMean':
        scale = {p: 1.0 / bigwig_genome_mean(p) for p in paths}
    thresholds = sorted(set(args.stateThresholds) | {args.stateThreshold})
    state_signals = {}
    state_rows = []
    for mark in marks:
        values = []
        for condition, tracks in ((name_a, state_tracks_a), (name_b, state_tracks_b)):
            signal = np.mean([np.array([read[p][n] for n in node_list]) * scale[p] for p in tracks[mark]], axis=0)
            values.append(dict(zip(node_list, signal.tolist())))
            row = {'mark': mark, 'declared_role': declared[mark], 'condition': condition,
                   'scale': ';'.join('{:.4g}'.format(scale[p]) for p in tracks[mark]), 'fragments': len(signal),
                   'median': float(np.median(signal))}
            for threshold in thresholds:
                row['fraction_ge_{:g}'.format(threshold)] = float(np.mean(signal >= threshold))
            state_rows.append(row)
        state_signals[mark] = tuple(values)
    write_table(pd.DataFrame(state_rows), prefix + '_state_marks.tsv')

    table = read_expression(args.expression)
    differential_ids = None
    if args.padjColumn:
        differential_ids = set(table.index[table[args.padjColumn].astype(float) <= args.alpha])
    change = gene_values_to_nodes(graph, expression_change(table, args), key=args.geneKey)
    change = {n: v for n, v in change.items() if graph.nodes[n].get('bait')}
    for node, value in change.items():
        graph.nodes[node]['rna_log2FC'] = value
    log.info('%d bait nodes with expression change', len(change))

    if args.maxHops < 2:
        raise ValueError('--maxHops must be at least 2')
    proximity = read_proximity_support(args.proximitySupport) if args.proximitySupport else None
    operators = ActivityOperators(graph, list(change), max_hops=args.maxHops)
    log.info('nodes per level (median over genes): %s',
             {level: float(np.median(operators.level(level)[2])) for level in operators.levels})
    mark_vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, change, operators=operators, mark_vectors=mark_vectors)
    started = time.time()
    mirror = mirrored_features(graph, operators, tracks_a, tracks_b, threads=args.threads)
    log.info('mirrored loci: %.0f s with %d threads', time.time() - started, args.threads)
    features = features.join(mirror, on='node')
    write_table(features, prefix + '_gene_features.tsv')
    # genes with direct and indirect neighbours carry information on the promoter, direct and indirect levels (the
    # same genes for every --maxHops; genes without nodes at k >= 3 hops have a k-hop change of 0)
    informative = features[(features['n_direct'] > 0) & (features['n_indirect'] > 0)].reset_index(drop=True)
    write_table(feature_correlations(informative, list(marks)), prefix + '_correlations.tsv')
    write_table(nested_models(informative, list(marks)), prefix + '_nested_models.tsv')
    comparison = mirror_comparison(informative, list(marks), n_bootstrap=args.bootstrap, seed=args.seed)
    differential_nodes = None
    if differential_ids is not None:
        differential_nodes = {n for n in features['node']
                              if set(str(graph.nodes[n].get(args.geneKey, '')).split(',')) & differential_ids}
    classified = features if differential_nodes is None else features[features['node'].isin(differential_nodes)]
    # learned roles: association of every mark change with the expression change (step 0, step 1) and with the
    # change of the promoter's own interactions, chromosome block-bootstrap intervals
    learned_table = learn_mark_roles(features, list(marks), declared, args.chromatinThreshold, args.rnaThreshold,
                                     n_boot=args.roleBootstrap, seed=args.seed)
    write_table(learned_table, prefix + '_mark_roles.tsv')
    learned = learned_roles(learned_table, declared)
    log.info('mark roles declared %s, learned %s', declared, learned)
    roles = learned if args.roles == 'learned' else declared
    mark_signs = role_signs(roles)
    classes = classify_genes(classified, list(marks), mark_signs, args.rnaThreshold, args.chromatinThreshold,
                             proximity=proximity)
    # learned states: recurring mark combinations and their transitions, rated by the expression changes
    write_table(learned_state_transitions(features, graph, state_signals,
                                          dict(zip(learned_table['mark'], learned_table['learned_sign'])),
                                          args.stateThreshold, args.rnaThreshold, args.stateMinFragments,
                                          n_boot=args.roleBootstrap, seed=args.seed),
                prefix + '_state_transitions.tsv')

    def cascade_with(roles, suffix):
        """Explanation cascade with the states of `roles` (gene: marks / interactions, then bivalent promoter states;
        step 1: marks, then poised states of contacting fragments; then the first step with a concordant mark
        change); transitions and cascade counts at every presence threshold."""
        signs = role_signs(roles)
        callable_ = states_callable(roles)
        log.info('roles %s: bivalent promoters %s, poised contacting fragments %s', roles,
                 'callable' if callable_['promoter'] else 'not callable',
                 'callable' if callable_['step1'] else 'not callable')
        tables = {t: gene_states(classes, graph, state_signals, roles, t, args.stateMinChange, names=(name_a, name_b))
                  for t in thresholds}
        transitions_report, thresholds_report = [], []
        for t in thresholds:
            for level in ('promoter', 'step1'):
                transitions_report.append(state_transition_table(classes, tables[t], level).assign(threshold=t))
            variant = explanation_cascade(classes, list(marks), signs, args.chromatinThreshold, states=tables[t])
            thresholds_report.append(cascade_table(variant, max_hops=args.maxHops, states_callable=callable_)
                                     .assign(threshold=t))
        for name, report in (('_bivalent_transitions', transitions_report), ('_state_thresholds', thresholds_report)):
            report = pd.concat(report, ignore_index=True)
            write_table(report[['threshold'] + [c for c in report if c != 'threshold']], prefix + name + suffix + '.tsv')
        states = tables[args.stateThreshold]
        cascade = explanation_cascade(classes, list(marks), signs, args.chromatinThreshold, states=states)
        write_table(explanation_types(classes, cascade, list(marks), signs, args.chromatinThreshold, callable_),
                    prefix + '_explanation_types' + suffix + '.tsv')
        return cascade, states, callable_

    cascade, states, callable_ = cascade_with(roles, '')
    details = explanation_details(classes, cascade, graph, marks, mark_signs, args.chromatinThreshold, states=states)
    classes = classes.assign(explained_at=cascade['explained_at'].to_numpy(),
                             explained_type=cascade['explained_type'].to_numpy(),
                             explained_by=cascade['explained_by'].to_numpy())
    classes = classes.join(details, on='node').join(states, on='node')
    if args.roles == 'both':
        other, _, other_callable = cascade_with(learned, '_learned')
        write_table(cascade_table(other, max_hops=args.maxHops, states_callable=other_callable),
                    prefix + '_explanation_cascade_learned.tsv')
        classes = classes.assign(explained_at_learned=other['explained_at'].to_numpy(),
                                 explained_type_learned=other['explained_type'].to_numpy(),
                                 explained_by_learned=other['explained_by'].to_numpy())
    write_table(classes, prefix + '_gene_classes.tsv')
    indirect_only = indirect_only_genes(graph, operators, mark_vectors, features, list(marks), mark_signs,
                                        args.rnaThreshold, args.chromatinThreshold, genes=differential_nodes,
                                        proximity=proximity)
    new_columns = [c for c in indirect_only.columns
                   if c in ('explaining_hops', 'n_hop_nodes_tested', 'proximity_support', 'proximity_support_ci')]
    if args.genePermutations:
        tested = classified if differential_nodes is not None else features
        # the indirect level for every differential gene (as before), every level k >= 3 for the genes first
        # explained there; one set of permutations for all levels
        # and every cascade step for the genes first explained there (p-values do not depend on the tested set)
        by_marks = cascade['explained_type'] == 'marks'
        levels = {'direct': cascade.index[(cascade['explained_at'] == 'step1') & by_marks].tolist(), 'indirect': None}
        for k in range(3, args.maxHops + 1):
            levels[level_name(k)] = sorted(set(classes.loc[classes['explanation'] == level_name(k), 'node']) |
                                           set(cascade.index[cascade['explained_at'] == 'step{}'.format(k)]))
        started = time.time()
        tests = gene_level_hop_test(operators, mark_vectors, tested, list(marks), mark_signs, args.rnaThreshold,
                                    args.genePermutations, args.seed, genes=levels, threads=args.threads)
        log.info('per-gene test (%s): %.0f s with %d threads', ', '.join(levels), time.time() - started, args.threads)
        gene_test = tests['indirect']
        write_table(gene_test, prefix + '_indirect_gene_test.tsv')
        for level, test in tests.items():
            if level not in ('direct', 'indirect'):
                # genes first explained at this level by classify_genes, BH across them
                test = test[test['node'].isin(set(classes.loc[classes['explanation'] == level, 'node']))].copy()
                test[level + '_q'] = _benjamini_hochberg(test[level + '_p'].to_numpy())
                write_table(test, prefix + '_{}_gene_test.tsv'.format(level))
        cascade_p = {}
        for k in range(1, args.maxHops + 1):
            step = set(cascade.index[(cascade['explained_at'] == 'step{}'.format(k)) & by_marks])
            test = tests[level_name(k)]
            cascade_p.update({n: p for n, p in zip(test['node'], test[level_name(k) + '_p']) if n in step})
        if len(indirect_only):
            indirect_only = indirect_only.merge(gene_test, on='node', how='left')
            # per explaining level: statistic, P and BH q across the genes first explained at that level
            parts = []
            for hops, part in indirect_only.groupby('explaining_hops', sort=False):
                level = level_name(int(hops))
                values = tests[level].set_index('node').loc[part['node']]
                part = part.assign(hops_statistic=values[level + '_statistic'].to_numpy(),
                                   hops_p=values[level + '_p'].to_numpy())
                parts.append(part.assign(hops_q=_benjamini_hochberg(part['hops_p'].to_numpy())))
            indirect_only = pd.concat(parts).loc[indirect_only.index]
            new_columns += ['hops_statistic', 'hops_p', 'hops_q']
    if len(indirect_only):
        indirect_only = indirect_only[[c for c in indirect_only.columns if c not in new_columns] + new_columns]
    write_table(indirect_only, prefix + '_indirect_only_genes.tsv')
    names = dict(zip(classes['node'], classes['genes']))
    write_table(cascade_table(cascade, names, cascade_p if args.genePermutations else None, proximity, args.maxHops,
                              callable_), prefix + '_explanation_cascade.tsv')
    informative_change = {n: change[n] for n in informative['node']}
    result = permutation_test(graph, marks, informative_change, mark_signs, n_permutations=args.permutations,
                              seed=args.seed,
                              rna_threshold=args.rnaThreshold, chromatin_threshold=args.chromatinThreshold,
                              max_hops=args.maxHops)
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
                                                     n_permutations=args.communityPermutations, seed=args.seed,
                                                     threads=args.threads)
    write_table(community_changes, prefix + '_community_expression.tsv')
    summary = {'transitions': transitions['transition'].value_counts().to_dict(),
               'gene_status': gene_changes['status'].value_counts().to_dict() if len(gene_changes) else {},
               'switching_vs_differential': switch_test,
               'communities_tested': int(len(community_changes)),
               'coherent_community_pairs_q_lt_0.1': int((community_changes['coherence_q'] < 0.1).sum())
               if len(community_changes) else 0,
               'coherent_communities_b_q_lt_0.1': int(community_changes.loc[community_changes['coherence_q'] < 0.1,
                                                                            'community_b'].nunique())
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
