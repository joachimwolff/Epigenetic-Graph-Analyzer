"""egaDifferential: differential interaction graph with chromatin and expression changes."""

import json
import logging
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .annotation import gene_values_to_nodes, parse_track_arguments
from .differential import (REPLICATE_MODES, STEP_FEATURES, ActivityOperators, _benjamini_hochberg, annotate_chicdiff,
                           classify_genes, consensus_scores, differential_edges_from_scores, feature_correlations,
                           gene_features, community_expression_changes, community_transitions,
                           gene_community_changes, gene_level_hop_test, indirect_only_genes, level_name,
                           match_communities, mirror_comparison, mirrored_features, nested_models, permutation_test,
                           _interval_signals, cascade_table, default_threads, explanation_cascade, read_proximity_support,
                           union_graph, explanation_details, explanation_types, explanation_levels,
                           explanation_level_summary, explanation_level_counts, explanation_combinations,
                           PRESENCE_STATISTICS)
from .states import (MARK_ROLES, PeakPresence, bigwig_genome_mean, gene_states, learn_mark_roles,
                     learned_level_signs, learned_roles, learned_state_transitions, level_states_callable,
                     parse_mark_roles, peak_nodes, read_peaks, resolve_mark_roles, role_signs,
                     state_transition_table)
from .elements import (assign_elements, element_agreement, element_class_counts, element_transitions,
                       node_elements, node_types, read_element_bed, roles_per_element_type, stratum_masks, tss_status)
from .egaCreateGraph import load_tss
from .graph import annotate_genes, build_graph, detect_communities, write_graph
from .structural import (add_sv_edges, breakpoint_calls, flag_sv_edges, read_copy_number, read_structural_variants,
                         read_sv_contacts, read_trans_calls, structural_genes, sv_edge_table, sv_status)
from .egaTadCrossing import add_tad_arguments, tad_kwargs
from .support import all_levels_provenance, default_floor, edge_support, input_detection, support_summary
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
    parser.add_argument('--supportThreshold', type=float, default=None,
                        help='Calling score of the support classes of gained and lost edges (default: '
                             '--scoreThreshold, else --referenceScore)')
    parser.add_argument('--supportFloor', type=float, default=None,
                        help='Gained / lost edges whose pair is listed in the other condition at a score of at least '
                             'this value (below the calling score) are threshold-sensitive, below it stable; absence '
                             'from a file whose lowest score is above the floor is not testable (default: 0.6 x the '
                             'calling score, 3 at CHiCAGO 5). Columns support, support_reason, other_max_score, '
                             'replicates_floor_<name> of <prefix>_edges.tsv.gz, edge attribute support, '
                             '<prefix>_support_inputs.tsv (lowest score per input file), <prefix>_support.tsv '
                             '(edges per status and class) and the support columns of the all-level explanations')
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
    parser.add_argument('--stepFeature', choices=STEP_FEATURES, default='combined',
                        help="What the step levels (1 .. --maxHops) measure. 'combined': log2 (S_B + 1) / (S_A + 1) of "
                             "the contact-weighted mark signal S with the contact strengths r and marks x of each "
                             "condition (changes of r and of x are confounded); 'marks': the mark term, S with the "
                             "marks of B versus A under the fixed strengths r_bar = (r_A + r_B) / 2 on the union "
                             "graph (an edge called in one condition enters with half its strength); 'contacts': the "
                             "contact term, S with the strengths of B versus A under the fixed marks x_bar = "
                             "(x_A + x_B) / 2; 'decomposed': the mark term for the cascade, the learned roles, the "
                             "mirror control and the per-gene permutation test, and in addition the columns "
                             "<mark>_<level>_combined / _marks / _contacts (and _combined_mirror, _contacts_mirror) "
                             "and in the all-level table combined_changes, mark_term_changes, contact_term_changes, "
                             "concordant_combined_R, concordant_mark_term_R, concordant_contact_term_R and driver_R "
                             "(marks, contacts, both, neither: which term changes concordantly at the level). The "
                             "promoter levels and contact_<level> are the same under every choice")
    parser.add_argument('--proximitySupport',
                        help='Step curve of analysis/hop_distance.py (TSV with steps, concordance, ci_low, ci_high): '
                             'adds to every explained gene the probability that graph neighbours at its explaining '
                             'level are closer in the nucleus than distance-matched pairs (proximity_support, '
                             'proximity_support_ci; 1.0 / NA at the promoter)')
    parser.add_argument('--stateThreshold', type=float, default=2.0,
                        help='Chromatin states: a mark is present on a fragment if its mean signal over the fragment '
                             'is at least this value (fold enrichment; see --stateScale). Bivalent promoter (bait): an '
                             'active-role and a repressive-role mark present; poised enhancer (other fragment): an '
                             'active- or enhancer-role and a repressive-role mark present. Marks with peak calls '
                             '(--statePeaksA / --statePeaksB) are present where a peak overlaps the fragment instead')
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
                             '--markSign), learned from the data separately at every level (<prefix>_mark_roles.tsv: '
                             'promoter, step 1 .. --maxHops and the own interactions; at each level the cascade uses '
                             'the role learned there, a mark without a consistent effect at a level explains no gene '
                             'there), or both (declared for the main outputs, learned in addition as '
                             '<prefix>_explanation_cascade_learned.tsv and columns *_learned)')
    parser.add_argument('--roleBootstrap', type=int, default=1000,
                        help='Chromosome block-bootstrap resamples of the learned roles and learned states')
    parser.add_argument('--stateTracksA', nargs='+',
                        help='mark=file.bw for the chromatin states of condition A (e.g. fold enrichment over input), '
                             'when the --tracksA signal is not on a fold-enrichment scale; marks not listed use --tracksA')
    parser.add_argument('--stateTracksB', nargs='+', help='as --stateTracksA for condition B')
    parser.add_argument('--statePeaksA', nargs='+',
                        help='mark=peaks.bed[,rep.bed] for the chromatin states of condition A: peak calls of any '
                             'caller (MACS2 / MACS3 narrowPeak or broadPeak, SEACR, generic BED; columns chrom, start, '
                             'end are read, several files of a mark are pooled). A mark with peaks is present on a '
                             'fragment when at least one peak overlaps it by at least --statePeakMinOverlap bp, under '
                             'both presence statistics and every --stateThresholds value (bivalent promoters, poised '
                             'enhancers, element classes, learned states); marks without peaks keep the threshold '
                             'rule. The signal of the mark (--stateTracksA, else --tracksA) still gives the change '
                             'required by --stateMinChange and the values in the state summaries; the mark changes '
                             'of the explanation levels always come from the signal tracks. A mark given here needs '
                             'peaks in --statePeaksB too')
    parser.add_argument('--statePeaksB', nargs='+', help='as --statePeaksA for condition B')
    parser.add_argument('--statePeakMinOverlap', type=int, default=1,
                        help='Overlap in bp of a peak with a fragment that makes the mark present (--statePeaksA / '
                             '--statePeaksB)')
    parser.add_argument('--stateWindow', type=int, default=1000,
                        help='Chromatin states are called under two presence statistics of a mark on a fragment, both '
                             'reported: the mean signal over the fragment (outputs without suffix, columns *_mean) '
                             'and the maximum of the mean signal over windows of at most this many bp within the '
                             'fragment (outputs *_max1kb, columns *_max1kb), which keeps narrow peaks on long '
                             'fragments')
    parser.add_argument('--stateMinChange', type=float, default=0.5,
                        help='A state transition needs every mark crossing the presence threshold to change by at '
                             'least this log2 fold (pseudo count 1) in the direction of the crossing')
    parser.add_argument('--annotation',
                        help='Element classes: GTF or TSS BED whose transcription start sites make a fragment a '
                             'promoter, baited or not (gene types --geneTypes; default: --tss)')
    parser.add_argument('--elements', nargs='+',
                        help='External element annotation, BED with a name column (e.g. ENCODE cCREs, ChromHMM '
                             'states): one file for both conditions or one per condition (A B); the agreement of the '
                             'element classes with it is written to <prefix>_element_agreement.tsv')
    parser.add_argument('--primaryRepressive', nargs='*', default=None,
                        help='Element classes: primary repressive marks (e.g. the Polycomb-type mark); a mark that is '
                             'repressive without being primary makes a heterochromatin-like element (default: the '
                             'marks declared repressive)')
    parser.add_argument('--elementMinGenes', type=int, default=50,
                        help='Roles learned per chromatin class of the nodes: classes with at least this many genes '
                             'at a level')
    parser.add_argument('--structuralVariants', nargs='+',
                        help='Breakpoints of structural variants: BEDPE with the conditions in which every join is '
                             'present (column conditions, present_<condition>, pairs_<condition> / '
                             'obs_exp_<condition>; analysis/capture_sv_calls.py writes <prefix>_sv.bedpe), '
                             'condition=file.bedpe (every join of the file present in that condition) or '
                             'condition=<prefix>_trans_hotspots.tsv of analysis/capture_sv_cnv.py. Adds the level '
                             "'structural variant' after the promoter levels, flags the graph edges that join the two "
                             'sides of a breakpoint (sv=True) and writes <prefix>_sv_genes.tsv, _sv_breakpoints.tsv')
    parser.add_argument('--svContacts',
                        help='Contacts of baited fragments across the junctions (<prefix>_sv_contacts.tsv of '
                             'analysis/capture_sv_calls.py): the promoter\'s own contacts with the partner side, '
                             'written as edges sv=True to <prefix>_sv_edges.tsv and <prefix>_union_sv.graphml. With '
                             '--format ibed the called inter-chromosomal and long-range pairs across a breakpoint '
                             '(dropped from the cis graph) are written to <prefix>_sv_called_edges.tsv and added '
                             'to <prefix>_union_sv.graphml as edges sv=True')
    parser.add_argument('--copyNumber', nargs='+',
                        help="condition=file.bedGraph of the relative copy number per condition: adds the level "
                             "'copy number' after the structural variant level")
    parser.add_argument('--svWindow', type=float, default=1.0,
                        help='Mb: promoters within this distance of a breakend, partner nodes within it of the '
                             'partner breakend')
    parser.add_argument('--svMinPairs', type=int, default=50,
                        help='capture_sv_cnv.py hotspot tables: read pairs of a breakpoint bin pair')
    parser.add_argument('--svMinObsExp', type=float, default=20.0,
                        help='capture_sv_cnv.py hotspot tables: observed / expected of a breakpoint bin pair')
    parser.add_argument('--svThreshold', type=float, default=1.0,
                        help='log2 change of the support (obs/exp) of a join present in both conditions counted as '
                             'gained or lost')
    parser.add_argument('--copyNumberThreshold', type=float, default=0.5,
                        help='log2 B/A of the relative copy number at the promoter counted as a change')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Bootstrap resamples of the mirror comparison')
    parser.add_argument('--resolution', type=float, default=1.0, help='Community resolution parameter')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threads', type=int, default=default_threads(),
                        help='Worker processes for the bigWig reads of the mirrored loci, the per-gene permutation '
                             'test and the community coherence test (default: the available cores, at most 32); '
                             'results do not depend on it')
    add_tad_arguments(parser)
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Output prefix; with --tads also <prefix>_tad_crossing_* (egaTadCrossing)')
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


def state_peak_presence(args, graph, node_list, marks, names):
    """mark -> (set of nodes with a peak in A, set in B) of the marks given in --statePeaksA / --statePeaksB."""
    peaks_a = parse_track_arguments(args.statePeaksA)
    peaks_b = parse_track_arguments(args.statePeaksB)
    if set(peaks_a) != set(peaks_b):
        raise ValueError('--statePeaksA and --statePeaksB need the same marks: {} vs {}'.format(sorted(peaks_a),
                                                                                               sorted(peaks_b)))
    unknown = sorted(set(peaks_a) - set(marks))
    if unknown:
        raise ValueError('Peaks of marks without tracks (--tracksA / --tracksB): {}'.format(', '.join(unknown)))
    chroms = [graph.nodes[n]['chrom'] for n in node_list]
    starts = [int(graph.nodes[n]['start']) for n in node_list]
    ends = [int(graph.nodes[n]['end']) for n in node_list]
    out = {}
    for mark in peaks_a:
        out[mark] = tuple(peak_nodes(read_peaks(paths[mark]), node_list, chroms, starts, ends,
                                     args.statePeakMinOverlap) for paths in (peaks_a, peaks_b))
        log.info('%s: presence from peaks (overlap >= %d bp): %d fragments in %s, %d in %s', mark,
                 args.statePeakMinOverlap, len(out[mark][0]), names[0], len(out[mark][1]), names[1])
    return out


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    name_a, name_b = args.names
    prefix = args.outputPrefix
    baits = read_baits(args.baits) if args.baits else None

    replicate_tables = {}

    def condition(paths, name):
        tables = []
        for path in paths:
            table = read_interactions([path], args.format, baits=baits, cool_field=args.coolField)
            if args.maxFdr is not None and 'fdr' in table:
                table = table[table['fdr'] <= args.maxFdr]
            tables.append(table)
        replicate_tables[name] = tables
        if len(tables) == 1 and args.scoreThreshold is not None:
            return consensus_scores([tables[0][tables[0]['score'] >= args.scoreThreshold]], 'union')
        return consensus_scores(tables, args.replicates, args.scoreThreshold)

    table_a = condition(args.interactionsA, name_a)
    table_b = condition(args.interactionsB, name_b)
    edges = differential_edges_from_scores(table_a, table_b, name_a, name_b, args.referenceScore)
    if args.chicdiff:
        edges = annotate_chicdiff(edges, read_chicdiff(args.chicdiff), args.chicdiffPadj)
    community_weight = 'score'
    log.info('%d union edges: %s', len(edges), edges['status'].value_counts().to_dict())
    # support of the gained and lost edges: the scores of the other condition below the calling threshold, and
    # whether its input files hold sub-threshold scores at all (lowest score per file)
    support_threshold = args.supportThreshold if args.supportThreshold is not None else \
        args.scoreThreshold if args.scoreThreshold is not None else args.referenceScore
    support_floor = args.supportFloor if args.supportFloor is not None else default_floor(support_threshold)
    detection = input_detection(replicate_tables, {name_a: args.interactionsA, name_b: args.interactionsB},
                                (name_a, name_b), support_threshold, support_floor)
    write_table(detection, prefix + '_support_inputs.tsv')
    for row in detection.itertuples():
        log.info('input %s replicate %d: %d pairs, lowest score %.3g: absence %s', row.condition, row.replicate,
                 row.pairs, row.min_score, row.absence)
    edges = pd.concat([edges, edge_support(edges, replicate_tables, detection, (name_a, name_b), support_threshold,
                                           support_floor)], axis=1)
    replicate_tables.clear()
    support_table = support_summary(edges)
    write_table(support_table.assign(threshold=support_threshold, floor=support_floor), prefix + '_support.tsv')
    # inputs of the run for the displays (egaReport --overlay draws the interactions and mark tracks of a region)
    with open(ensure_parent(prefix + '_support.json'), 'w') as handle:
        json.dump({'names': [name_a, name_b], 'threshold': support_threshold, 'floor': support_floor,
                   'format': args.format, 'replicates': args.replicates, 'scoreThreshold': args.scoreThreshold,
                   'interactions': {name_a: [os.path.abspath(p) for p in args.interactionsA],
                                    name_b: [os.path.abspath(p) for p in args.interactionsB]},
                   'tracks': {name_a: {m: [os.path.abspath(p) for p in v]
                                       for m, v in parse_track_arguments(args.tracksA).items()},
                              name_b: {m: [os.path.abspath(p) for p in v]
                                       for m, v in parse_track_arguments(args.tracksB).items()}},
                   'inputs': json.loads(detection.to_json(orient='records'))}, handle, indent=2)
    log.info('edge support (threshold %g, floor %g): %s', support_threshold, support_floor,
             {s: dict(zip(p['support'], p['edges'])) for s, p in support_table.groupby('status', sort=False)})
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
    # presence statistic max1kb: largest window mean within every fragment, for the tracks of the states
    state_paths = list(dict.fromkeys(p for tracks in (state_tracks_a, state_tracks_b) for mark in tracks
                                     for p in tracks[mark]))
    started = time.time()
    read_max = _interval_signals(state_paths, [graph.nodes[n]['chrom'] for n in node_list],
                                 np.array([graph.nodes[n]['start'] for n in node_list], dtype=np.int64),
                                 np.array([graph.nodes[n]['end'] for n in node_list], dtype=np.int64), args.threads,
                                 window=args.stateWindow)
    read_max = {p: dict(zip(node_list, values.tolist())) for p, values in read_max.items()}
    log.info('maximum %d bp window signals of %d bigWig files: %.0f s', args.stateWindow, len(state_paths),
             time.time() - started)
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
    # peak calls defining the presence of a mark (instead of the threshold rule)
    peak_present = state_peak_presence(args, graph, node_list, list(marks), (name_a, name_b))
    state_signals_by = {}
    state_rows = []
    for presence, source in (('mean', read), ('max1kb', read_max)):
        state_signals_by[presence] = {}
        for mark in marks:
            values = []
            for c, (condition, tracks) in enumerate(((name_a, state_tracks_a), (name_b, state_tracks_b))):
                signal = np.mean([np.array([source[p][n] for n in node_list]) * scale[p] for p in tracks[mark]],
                                 axis=0)
                values.append(dict(zip(node_list, signal.tolist())))
                row = {'mark': mark, 'declared_role': declared[mark], 'condition': condition, 'presence': presence,
                       'scale': ';'.join('{:.4g}'.format(scale[p]) for p in tracks[mark]), 'fragments': len(signal),
                       'median': float(np.median(signal))}
                if mark in peak_present:
                    # presence from the peaks: the same fraction at every threshold
                    row['scale'] = 'peaks'
                    fraction = float(np.mean([n in peak_present[mark][c] for n in node_list]))
                    for threshold in thresholds:
                        row['fraction_ge_{:g}'.format(threshold)] = fraction
                else:
                    for threshold in thresholds:
                        row['fraction_ge_{:g}'.format(threshold)] = float(np.mean(signal >= threshold))
                state_rows.append(row)
            state_signals_by[presence][mark] = PeakPresence(values, peak_present[mark]) if mark in peak_present \
                else tuple(values)
    state_signals = state_signals_by['mean']
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
    operators = ActivityOperators(graph, list(change), max_hops=args.maxHops, step_feature=args.stepFeature)
    log.info('step feature: %s', args.stepFeature)
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
                                     n_boot=args.roleBootstrap, seed=args.seed, levels=operators.levels)
    # learned roles per level: no level decides for another
    learned = learned_roles(learned_table, declared)
    learned_signs = learned_level_signs(learned_table)
    log.info('mark roles declared %s, learned per level %s', declared, learned)
    # element classes of every fragment per condition (declared and learned roles, both presence statistics)
    started = time.time()
    if args.annotation:
        status = tss_status(graph, annotate_genes(nodes, load_tss(args.annotation, args.geneTypes)))
    else:
        status = tss_status(graph)
    primary = args.primaryRepressive if args.primaryRepressive is not None else \
        [m for m, r in declared.items() if r == 'repressive']
    learned_other = {m: learned.get('direct', {}).get(m, 'neutral') for m in marks}
    learned_bait = {m: learned.get('promoter', {}).get(m, 'neutral') for m in marks}
    elements = node_elements(graph, state_signals_by, args.stateThreshold,
                             {'declared': declared, 'learned': learned_other}, status, names=(name_a, name_b),
                             primary=primary, bait_roles_by={'declared': declared, 'learned': learned_bait})
    write_table(elements, prefix + '_elements.tsv')
    other_ends = set(edges['node_b'])
    write_table(element_class_counts(elements, other_ends), prefix + '_element_classes.tsv')
    if args.elements:
        external = [assign_elements(graph, read_element_bed(path)) for path in args.elements]
        external = {name_a: external[0], name_b: external[-1]}
        write_table(element_agreement(elements, external), prefix + '_element_agreement.tsv')
    log.info('element classes of %d fragments: %.0f s', graph.number_of_nodes(), time.time() - started)
    # roles learned per element type of the nodes (promoter / non-promoter element, chromatin class)
    started = time.time()
    masks = stratum_masks(elements, operators)
    per_type = roles_per_element_type(features, operators, mark_vectors, list(marks), declared, masks,
                                      args.chromatinThreshold, args.rnaThreshold, n_boot=args.roleBootstrap,
                                      seed=args.seed, min_genes=args.elementMinGenes, threads=args.threads)
    write_table(pd.concat([learned_table.assign(element_type='all'), per_type], ignore_index=True),
                prefix + '_mark_roles.tsv')
    log.info('roles per element type (%d strata): %.0f s', len(masks), time.time() - started)
    declared_setting = (declared, declared, role_signs(declared))
    learned_setting = (learned.get('promoter', {}), learned.get('direct', {}), learned_signs)
    promoter_roles, step1_roles, mark_signs = learned_setting if args.roles == 'learned' else declared_setting
    classes = classify_genes(classified, list(marks), mark_signs, args.rnaThreshold, args.chromatinThreshold,
                             proximity=proximity)
    # learned states: recurring mark combinations and their transitions, rated by the expression changes
    write_table(learned_state_transitions(features, graph, state_signals, learned_signs,
                                          args.stateThreshold, args.rnaThreshold, args.stateMinFragments,
                                          n_boot=args.roleBootstrap, seed=args.seed),
                prefix + '_state_transitions.tsv')

    # structural variants and copy number near the promoters (levels after the promoter levels)
    structural, structural_levels, sv_edges = None, [], None
    if args.structuralVariants or args.copyNumber:
        window = int(args.svWindow * 1e6)
        svs = read_structural_variants(args.structuralVariants or [], (name_a, name_b), args.svMinPairs,
                                       args.svMinObsExp)
        contacts = read_sv_contacts(args.svContacts, (name_a, name_b)) if args.svContacts else None
        copy_number = read_copy_number(args.copyNumber, (name_a, name_b)) if args.copyNumber else None
        if args.structuralVariants:
            structural_levels.append('structural_variant')
            flagged = flag_sv_edges(graph, svs, window)
            breakpoints = pd.concat([svs, sv_status(svs, (name_a, name_b), args.svThreshold)], axis=1) \
                if len(svs) else svs
            write_table(breakpoints, prefix + '_sv_breakpoints.tsv')
            log.info('%d structural variant joins, %d graph edges join the two sides of a breakpoint', len(svs),
                     flagged)
        if copy_number:
            structural_levels.append('copy_number')
        # called inter-chromosomal and long-range pairs across the junctions (CHiCAGO ibed lists them; the cis
        # graph leaves them out): written as edges sv=True, never mixed into the cis levels
        called = None
        if args.structuralVariants and args.format == 'ibed':
            calls = {name_a: read_trans_calls(args.interactionsA, args.scoreThreshold),
                     name_b: read_trans_calls(args.interactionsB, args.scoreThreshold)}
            if sum(len(t) for t in calls.values()):
                called = breakpoint_calls(calls, svs, window, (name_a, name_b))
                write_table(called, prefix + '_sv_called_edges.tsv')
                log.info('%d called pairs join the two sides of a breakpoint (%s)', len(called),
                         called['status'].value_counts().to_dict() if len(called) else {})
            else:
                log.info('no inter-chromosomal or long-range pairs among the called interactions')
        sv_genes, structural = structural_genes(
            graph, change, svs, (name_a, name_b), window, elements=elements, marks=marks, contacts=contacts,
            copy_number=copy_number, differential=set(classes['node']), sv_threshold=args.svThreshold,
            cn_threshold=args.copyNumberThreshold, called=called, tss=tss)
        if not args.structuralVariants:
            structural = structural.drop(columns=[c for c in structural if c.startswith('sv_')])
        if args.structuralVariants:
            write_table(sv_genes, prefix + '_sv_genes.tsv')
            log.info('%d genes with the promoter within %.1f Mb of a breakend (%d differential, %d enhancer-hijacking '
                     'candidates)', sv_genes['node'].nunique() if len(sv_genes) else 0, args.svWindow,
                     sv_genes.loc[sv_genes['differential'].astype(bool), 'node'].nunique() if len(sv_genes) else 0,
                     sv_genes.loc[sv_genes['hijacking_candidate'].astype(bool), 'node'].nunique()
                     if len(sv_genes) else 0)
        if contacts is not None:
            sv_edges = sv_edge_table(graph, contacts, svs, (name_a, name_b))
            write_table(sv_edges, prefix + '_sv_edges.tsv')
        if called is not None and len(called):
            sv_edges = called if sv_edges is None or not len(sv_edges) else \
                pd.concat([sv_edges.assign(source='read pairs'), called], ignore_index=True)

    def cascade_with(setting, suffix, presence='mean', write=True):
        """Explanation cascade with the roles of `setting` (promoter roles, step-1 roles, mark signs per level):
        gene: marks / interactions, then bivalent promoter states; step 1: marks, then poised states of contacting
        fragments; then the first step with a concordant mark change. Transitions and cascade counts at every
        presence threshold, with the presence statistic `presence` (mean, max1kb) of the marks on a fragment."""
        roles, roles_step1, signs = setting
        state_signals = state_signals_by[presence]
        callable_ = level_states_callable(roles, roles_step1)
        log.info('roles %s (step 1: %s): bivalent promoters %s, poised contacting fragments %s', roles, roles_step1,
                 'callable' if callable_['promoter'] else 'not callable',
                 'callable' if callable_['step1'] else 'not callable')
        tables = {t: gene_states(classes, graph, state_signals, roles, t, args.stateMinChange, names=(name_a, name_b),
                                 step1_roles=roles_step1)
                  for t in thresholds}
        transitions_report, thresholds_report = [], []
        for t in thresholds:
            for level in ('promoter', 'step1'):
                transitions_report.append(state_transition_table(classes, tables[t], level).assign(threshold=t))
            variant = explanation_cascade(classes, list(marks), signs, args.chromatinThreshold, states=tables[t],
                                          structural=structural)
            thresholds_report.append(cascade_table(variant, max_hops=args.maxHops, states_callable=callable_,
                                                   structural=structural is not None)
                                     .assign(threshold=t))
        states = tables[args.stateThreshold]
        cascade = explanation_cascade(classes, list(marks), signs, args.chromatinThreshold, states=states,
                                      structural=structural)
        if write:
            for name, report in (('_bivalent_transitions', transitions_report),
                                 ('_state_thresholds', thresholds_report)):
                report = pd.concat(report, ignore_index=True)
                write_table(report[['threshold'] + [c for c in report if c != 'threshold']],
                            prefix + name + suffix + '.tsv')
            write_table(explanation_types(classes, cascade, list(marks), signs, args.chromatinThreshold, callable_),
                        prefix + '_explanation_types' + suffix + '.tsv')
        return cascade, states, callable_

    # every role set under every presence statistic; the main role set (--roles) writes the outputs without suffix
    # (mean) and *_max1kb, the other role set *_learned and *_learned_max1kb with --roles both
    main_roles = 'learned' if args.roles == 'learned' else 'declared'
    settings = {'declared': declared_setting, 'learned': learned_setting}
    variants = {}
    for roles in ('declared', 'learned'):
        for presence in PRESENCE_STATISTICS:
            suffix = ('' if roles == main_roles else '_learned') + ('' if presence == 'mean' else '_' + presence)
            write = roles == main_roles or args.roles == 'both'
            variants[(roles, presence)] = cascade_with(settings[roles], suffix, presence, write) + (suffix, write)
    cascade, states, callable_ = variants[(main_roles, 'mean')][:3]
    # class transitions of the fragments and the expression change of the genes in contact (step 0 .. maxHops)
    write_table(element_transitions(elements, operators, features, args.rnaThreshold, n_boot=args.roleBootstrap,
                                    seed=args.seed, differential=set(classes['node'])),
                prefix + '_element_transitions.tsv')
    details = explanation_details(classes, cascade, graph, marks, mark_signs, args.chromatinThreshold, states=states,
                                  proximity=proximity, structural=structural)
    classes = classes.assign(explained_at=cascade['explained_at'].to_numpy(),
                             explained_type=cascade['explained_type'].to_numpy(),
                             explained_by=cascade['explained_by'].to_numpy())
    classes = classes.join(details, on='node').join(states, on='node')
    if structural is not None:
        classes = classes.join(structural, on='node')
    if args.roles == 'both':
        other, other_states, other_callable = variants[('learned', 'mean')][:3]
        other_details = explanation_details(classes, other, graph, marks, learned_signs, args.chromatinThreshold,
                                            states=other_states, proximity=proximity, structural=structural)
        classes = classes.assign(explained_at_learned=other['explained_at'].to_numpy(),
                                 explained_type_learned=other['explained_type'].to_numpy(),
                                 explained_by_learned=other['explained_by'].to_numpy(),
                                 explanation_detail_learned=other_details['explanation_detail']
                                 .reindex(classes['node']).to_numpy())
    indirect_only = indirect_only_genes(graph, operators, mark_vectors, features, list(marks), mark_signs,
                                        args.rnaThreshold, args.chromatinThreshold, genes=differential_nodes,
                                        proximity=proximity)
    new_columns = [c for c in indirect_only.columns
                   if c in ('explaining_hops', 'n_hop_nodes_tested', 'proximity_support', 'proximity_support_ci')]

    def first_step_p(which, roles):
        """node -> per-gene P of the step at which the cascade `which` first explains the gene by marks"""
        out = {}
        marks_type = which['explained_type'] == 'marks'
        for k in range(1, args.maxHops + 1):
            step = set(which.index[(which['explained_at'] == 'step{}'.format(k)) & marks_type])
            test = tests_by[roles][level_name(k)]
            out.update({n: p for n, p in zip(test['node'], test[level_name(k) + '_p']) if n in step})
        return out

    tests_by = {}
    if args.genePermutations:
        tested = classified if differential_nodes is not None else features
        # every step for every differential gene, under the declared and under the learned signs (the P of a gene
        # does not depend on the other genes tested)
        levels = {level: None for level in operators.levels}
        tests_by = {}
        for roles, roles_signs in (('declared', declared_setting[2]), ('learned', learned_signs)):
            started = time.time()
            tests_by[roles] = gene_level_hop_test(operators, mark_vectors, tested, list(marks), roles_signs,
                                                  args.rnaThreshold, args.genePermutations, args.seed, genes=levels,
                                                  threads=args.threads)
            log.info('per-gene test, %s signs (%s): %.0f s with %d threads', roles, ', '.join(levels),
                     time.time() - started, args.threads)
        tests = tests_by[main_roles]
        gene_test = tests['indirect']
        write_table(gene_test, prefix + '_indirect_gene_test.tsv')
        for level, test in tests.items():
            if level not in ('direct', 'indirect'):
                # genes first explained at this level by classify_genes, BH across them
                test = test[test['node'].isin(set(classes.loc[classes['explanation'] == level, 'node']))].copy()
                test[level + '_q'] = _benjamini_hochberg(test[level + '_p'].to_numpy())
                write_table(test, prefix + '_{}_gene_test.tsv'.format(level))
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
    # first-level cascade of every role set and presence statistic (one summary each)
    for (roles, presence), (which, _, which_callable, suffix, write) in variants.items():
        if write:
            write_table(cascade_table(which, names, first_step_p(which, roles) if args.genePermutations else None,
                                      proximity, args.maxHops, which_callable, structural=structural is not None),
                        prefix + '_explanation_cascade' + suffix + '.tsv')
    # every level for every differential gene, under both role sets and both presence statistics
    started = time.time()
    all_levels = explanation_levels(classes, graph, operators, mark_vectors, list(marks),
                                    {'declared': declared_setting[2], 'learned': learned_signs},
                                    states={key: value[1] for key, value in variants.items()},
                                    tests=tests_by if args.genePermutations else None, proximity=proximity,
                                    chromatin_threshold=args.chromatinThreshold, names=(name_a, name_b),
                                    node_types=node_types(elements), structural=structural,
                                    structural_levels=structural_levels)
    write_table(all_levels, prefix + '_explanations_all_levels.tsv')
    write_table(all_levels_provenance(all_levels.columns), prefix + '_explanations_provenance.tsv')
    write_table(explanation_level_counts(all_levels), prefix + '_explanation_levels.tsv')
    write_table(explanation_combinations(all_levels, top=20), prefix + '_explanation_combinations.tsv')
    classes = classes.join(explanation_level_summary(all_levels), on='node')
    write_table(classes, prefix + '_gene_classes.tsv')
    log.info('all-level explanations: %.0f s', time.time() - started)
    informative_change = {n: change[n] for n in informative['node']}
    result = permutation_test(graph, marks, informative_change, mark_signs, n_permutations=args.permutations,
                              seed=args.seed,
                              rna_threshold=args.rnaThreshold, chromatin_threshold=args.chromatinThreshold,
                              max_hops=args.maxHops, step_feature=args.stepFeature)
    result['step_feature'] = args.stepFeature
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
    if sv_edges is not None and len(sv_edges):
        # the junction contacts as edges sv=True, kept out of the cis graph used by the other levels
        write_graph(add_sv_edges(graph, sv_edges), ensure_parent(prefix + '_union_sv.graphml'))
    if args.tads:
        # boundaries crossed by the union edges: written next to the run, the other outputs are not changed
        from .tads import run_for_prefix
        run_for_prefix(prefix, args.tads, prefix + '_tad_crossing', graph=graph, names=[name_a, name_b],
                       **tad_kwargs(args))


if __name__ == '__main__':
    main()
