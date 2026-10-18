"""egaDesign: perturbation experiments (CRISPRi / CRISPRa / deletion) that test the explanations of egaDifferential.

Differential mode (-p <egaDifferential prefix>): for every explained gene (or the genes of --genes) every node of
every explaining level (promoter, own contacts, step-1 state, step 1 .. step k) is perturbed in silico and the level
statistic recomputed with the egaDifferential formulas (design module): flatten (the node's mark change set to 0),
silence (its marks set to 0 in both conditions) and, for the strongest candidates, remove (node and edges deleted,
distances recomputed). The per-gene permutation P of a step level is recomputed before and after with a reduced number
of permutations (--permutations). Candidates are ranked per gene by the predicted effect (--effect), the node's own
evidence (own mark change, level P, mirror control, 3D support of the step, contact in both conditions) and practical
criteria (non-promoter element, not pleiotropic); the gene's own promoter is listed last. Outputs (prefix -o):

  <out>_design.tsv            per gene the --top best candidates: rank, node, element class in plain words, nearest
                              gene, coordinates, predicted effect (statistic before / after, fraction removed, level
                              lost, P before / after), evidence, pleiotropy (other genes for which the node removes
                              >= half of a level), expected direction of the gene's change and one sentence.
  <out>_design.txt            compact text per gene (the three best candidates).
  <out>_design_top.bed        the --bedTop best non-promoter-of-the-gene nodes per gene as the 1 kb window with the
                              highest signal of the explaining mark in the condition where it is higher (--tracksA,
                              --tracksB; without tracks the central 1 kb), for primer design.
  <out>_design_candidates.tsv.gz  every candidate of every gene and level before ranking.
  <out>_design_levels.tsv     recomputed statistic of every level against the egaDifferential call (agreement).
  <out>_design_summary.json   counts.

Single-condition mode (--graph and --track): no differential run; the nodes within --maxHops steps of every bait
(or of the --genes) are ranked by their contribution to the bait's contact-weighted mark activity in that condition
(design.single_condition_contributions): <out>_single_condition.tsv.gz.

Candidates are reported, not decided.
"""

import json
import logging
import os
import re
import time

import numpy as np
import pandas as pd

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .design import (EVIDENCE, GraphArrays, TssIndex, _gene_task, _map, candidate_sentence, expected_direction,
                     level_hops_of, level_nulls, level_signs, mark_level_of, null_p, plain_class, pleiotropy,
                     rank_candidates, read_gene_strands, relative_position, single_condition_contributions)
from .graph import read_graph

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser(__doc__.split('\n\n')[0], experimental=True)
    parser.add_argument('--prefix', '-p', help='Output prefix of egaDifferential (union graph, '
                        'explanations_all_levels, gene_features, elements)')
    parser.add_argument('--names', nargs=2, default=None,
                        help='Condition names A B (default: from the prefix <A>_vs_<B>)')
    parser.add_argument('--genes', nargs='*', default=None,
                        help='Gene names to report (default: every explained gene); a file with one name per line '
                             'is read when the single argument is a path')
    parser.add_argument('--roles', default='declared', choices=['declared', 'learned'])
    parser.add_argument('--presence', default='mean', choices=['mean', 'max1kb'])
    parser.add_argument('--chromatinThreshold', type=float, default=1.0)
    parser.add_argument('--effect', default='flatten', choices=['flatten', 'silence', 'remove'],
                        help='Perturbation that ranks the candidates (remove is computed for the --nRemove '
                             'strongest candidates of a gene only)')
    parser.add_argument('--nRemove', type=int, default=10,
                        help='Candidates per gene with the remove perturbation (node and edges deleted)')
    parser.add_argument('--permutations', type=int, default=1000,
                        help='Permutations of the recomputed per-gene P (0: none)')
    parser.add_argument('--strong', type=float, default=0.5,
                        help='Pleiotropy: fraction of a level statistic a node must remove to count for a gene')
    parser.add_argument('--top', type=int, default=10, help='Candidates per gene in <out>_design.tsv')
    parser.add_argument('--bedTop', type=int, default=3, help='Nodes per gene in the BED file')
    parser.add_argument('--window', type=int, default=1000, help='Window of the BED file (bp)')
    parser.add_argument('--tracksA', nargs='*', default=[],
                        help='mark=bigWig[,bigWig] of condition A (BED windows; single-condition mode: --track)')
    parser.add_argument('--tracksB', nargs='*', default=[], help='mark=bigWig[,bigWig] of condition B')
    parser.add_argument('--tss', default=None, help='GTF for the gene strands (upstream / downstream in the text)')
    parser.add_argument('--graph', default=None, help='Single-condition mode: graph (graphml)')
    parser.add_argument('--track', default=None, help='Single-condition mode: mark=bigWig')
    parser.add_argument('--referenceScore', type=float, default=5.0,
                        help='Single-condition mode: edge strength = score / reference score')
    parser.add_argument('--maxHops', type=int, default=3, help='Single-condition mode: steps')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--outPrefix', '-o', required=True)
    return parser.parse_args(args)


def _tracks(entries):
    out = {}
    for entry in entries or []:
        mark, paths = entry.split('=', 1)
        out[mark] = [p for p in paths.split(',') if p]
    return out


def _gene_list(genes):
    if genes and len(genes) == 1 and os.path.exists(genes[0]):
        with open(genes[0]) as handle:
            return [line.strip() for line in handle if line.strip()]
    return genes


def _matches(label, wanted):
    names = {g.strip().upper() for g in re.split(r'[,;]', str(label)) if g.strip()}
    return bool(names & wanted)


def element_table(path, roles, presence, names):
    """node -> element_type, tss_genes, class in A and in B (elements.tsv of egaDifferential)."""
    if not os.path.exists(path):
        return pd.DataFrame(columns=['element_type', 'tss_genes', 'class_A', 'class_B']).rename_axis('node')
    table = pd.read_csv(path, sep='\t', dtype={'tss_genes': str})
    column = 'class_{}_{}'.format(roles, presence)
    if column not in table:
        column = [c for c in table if c.startswith('class_')][0]
    a = table[table['condition'] == names[0]].set_index('node')
    b = table[table['condition'] == names[1]].set_index('node')
    out = pd.DataFrame(index=a.index.union(b.index))
    base = pd.concat([a, b[~b.index.isin(a.index)]])
    out['element_type'] = base['element_type'].reindex(out.index)
    out['tss_genes'] = base['tss_genes'].reindex(out.index).fillna('')
    out['class_A'] = a[column].reindex(out.index)
    out['class_B'] = b[column].reindex(out.index)
    return out


def bed_windows(rows, arrays, tracks_a, tracks_b, window, threads):
    """1 kb window of every BED row: the highest signal of its window mark in its window condition."""
    from .triggers import node_windows
    rows = rows.copy()
    rows['win_start'] = ((rows['start'] + rows['end']) // 2 - window // 2).clip(lower=0)
    rows['win_end'] = rows['win_start'] + window
    rows['win_source'] = 'central window (no track)'
    for (mark, condition), part in rows.groupby(['window_mark', 'window_condition']):
        tracks = (tracks_a if condition == 'A' else tracks_b).get(mark)
        if not tracks:
            continue
        nodes = part[['node', 'chrom', 'start', 'end']].drop_duplicates('node')
        win = node_windows(nodes, tracks, window=window, threads=threads).set_index('node')
        rows.loc[part.index, 'win_start'] = win.loc[part['node'], 'win_start'].to_numpy()
        rows.loc[part.index, 'win_end'] = win.loc[part['node'], 'win_end'].to_numpy()
        rows.loc[part.index, 'win_source'] = '{} {}'.format(mark, condition)
    return rows


def single_condition(args):
    graph = read_graph(args.graph)
    mark, path = args.track.split('=', 1)
    from .annotation import bigwig_interval_signal
    nodes = list(graph.nodes)
    chroms = [graph.nodes[n]['chrom'] for n in nodes]
    starts = [int(graph.nodes[n]['start']) for n in nodes]
    ends = [int(graph.nodes[n]['end']) for n in nodes]
    signal = np.mean([bigwig_interval_signal(p, chroms, starts, ends) for p in path.split(',')], axis=0)
    for n, v in zip(nodes, signal):
        graph.nodes[n]['_signal'] = float(v)
    weight = 'score' if any('score' in d for _, _, d in graph.edges(data=True)) else 'weight'
    arrays = GraphArrays(graph, {mark: ('_signal', '_signal')}, weights=(weight, weight),
                         reference_score=args.referenceScore)
    sources = [n for n, d in graph.nodes(data=True) if d.get('bait')]
    genes = _gene_list(args.genes)
    if genes:
        wanted = {g.upper() for g in genes}
        sources = [n for n in sources if _matches(graph.nodes[n].get('genes', ''), wanted)]
    table = single_condition_contributions(arrays, sources, args.maxHops, mark)
    table.insert(1, 'genes', [graph.nodes[n].get('genes', '') for n in table['gene_node']])
    write_table(table, ensure_parent(args.outPrefix + '_single_condition.tsv.gz'))
    log.info('single-condition contributions of %d genes, %d rows', len(set(table['gene_node'])), len(table))


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    if args.graph:
        if not args.track:
            raise SystemExit('--graph needs --track')
        return single_condition(args)
    if not args.prefix:
        raise SystemExit('--prefix (differential mode) or --graph and --track (single-condition mode) required')
    start = time.time()
    names = args.names or os.path.basename(args.prefix).split('_vs_')
    roles, presence, threshold = args.roles, args.presence, args.chromatinThreshold
    explains_col = 'explains_{}_{}'.format(roles, presence)
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    levels[explains_col] = levels[explains_col].astype(str).str.lower() == 'true'
    features = pd.read_csv(args.prefix + '_gene_features.tsv', sep='\t', usecols=['node'])
    graph = read_graph(args.prefix + '_union.graphml')
    from .triggers import graph_marks
    marks = graph_marks(graph, names)
    arrays = GraphArrays(graph, {m: ('{}_{}'.format(m, names[0]), '{}_{}'.format(m, names[1])) for m in marks},
                         sources=features['node'])
    del graph
    max_hops = max(level_hops_of(lv) for lv in levels['level'].unique())
    signs = level_signs(levels, roles)
    log.info('%d marks (%s), %d steps, signs %s; graph read in %.0f s', len(marks), ','.join(marks), max_hops,
             dict(signs), time.time() - start)

    explained = levels.groupby('node')[explains_col].any()
    genes_wanted = _gene_list(args.genes)
    nodes = list(explained[explained].index)
    if genes_wanted:
        wanted = {g.upper() for g in genes_wanted}
        requested = [n for n in levels['node'].unique()
                     if _matches(levels.loc[levels['node'] == n, 'genes'].iloc[0], wanted)]
        nodes = list(dict.fromkeys(nodes + requested))
    else:
        requested = None
    by_gene = {n: part.set_index('level') for n, part in levels[levels['node'].isin(set(nodes))].groupby('node')}
    shared = dict(arrays=arrays, by_gene=by_gene, max_hops=max_hops, signs=signs, threshold=threshold,
                  n_remove=args.nRemove, explains_col=explains_col, concordant_col='concordant_marks_' + roles,
                  detail_col='detail_{}_{}'.format(roles, presence),
                  state_col='state_fragment_{}_{}'.format(roles, presence))
    t = time.time()
    results = _map(_gene_task, nodes, args.threads, shared)
    log.info('perturbation of %d genes: %.0f s', len(nodes), time.time() - t)
    candidates = pd.DataFrame([c for r in results for c in r[0]])
    level_check = pd.DataFrame([row for r in results for row in r[1]])
    hops_by_level = {}
    for node, r in zip(nodes, results):
        for k, h in r[2].items():
            hops_by_level.setdefault(k, {})[node] = h
    write_table(level_check, ensure_parent(args.outPrefix + '_design_levels.tsv'))
    mark_levels = level_check[level_check['level'] != 'promoter_contacts']
    agreement = float((mark_levels['explains_table'] == mark_levels['explains_recomputed']).mean()) \
        if len(mark_levels) else np.nan
    log.info('recomputed level calls agree with egaDifferential for %.4f of gene levels', agreement)
    if candidates.empty:
        log.warning('no candidates')
        return

    # permutation P before and after (step levels)
    for column in ('p_before_reduced', 'p_after_reduced', 'p_remove_reduced'):
        candidates[column] = np.nan
    if args.permutations > 0:
        directions = {n: float(np.sign(by_gene[n].iloc[0]['rna_log2FC'])) for n in nodes}
        for k, by_node in hops_by_level.items():
            t = time.time()
            nulls = level_nulls(arrays, by_node, directions, _signs_at(signs, k), args.permutations, args.seed + k,
                                args.threads)
            rows = candidates.index[(candidates['hops'] == k) & (candidates['level'] == mark_level_of(k))]
            for i in rows:
                null = nulls.get(candidates.at[i, 'gene_node'])
                candidates.at[i, 'p_before_reduced'] = null_p(null, candidates.at[i, 'statistic_before'])
                candidates.at[i, 'p_after_reduced'] = null_p(null, candidates.at[i, 'statistic_' + (
                    args.effect if args.effect != 'remove' else 'flatten')])
                candidates.at[i, 'p_remove_reduced'] = null_p(null, candidates.at[i, 'statistic_remove'])
            log.info('null of step %d (%d genes, %d permutations): %.0f s', k, len(by_node), args.permutations,
                     time.time() - t)

    # annotation
    elements = element_table(args.prefix + '_elements.tsv', roles, presence, names)
    idx = candidates['node_index'].to_numpy()
    candidates['chrom'] = arrays.chrom[idx]
    candidates['start'] = arrays.start[idx]
    candidates['end'] = arrays.end[idx]
    gidx = np.array([arrays.index[n] for n in candidates['gene_node']])
    mid = (candidates['start'] + candidates['end']) // 2
    gmid = (arrays.start[gidx] + arrays.end[gidx]) // 2
    same = candidates['chrom'].to_numpy() == arrays.chrom[gidx]
    candidates['distance_to_gene'] = np.where(same, np.abs(mid - gmid), np.nan)
    candidates['element_type'] = elements['element_type'].reindex(candidates['node']).fillna('').to_numpy()
    tss_genes = elements['tss_genes'].reindex(candidates['node']).fillna('').to_numpy()
    class_a = elements['class_A'].reindex(candidates['node']).to_numpy()
    class_b = elements['class_B'].reindex(candidates['node']).to_numpy()
    candidates['element_class'] = [plain_class(a, b, names) for a, b in zip(class_a, class_b)]
    tss_nodes = elements.index[elements['tss_genes'].astype(str).str.len() > 0]
    tss_idx = np.array([arrays.index[n] for n in tss_nodes if n in arrays.index], dtype=np.int64)
    tss_index = TssIndex(np.array(arrays.nodes, dtype=object)[tss_idx], arrays.chrom[tss_idx], arrays.start[tss_idx],
                         arrays.end[tss_idx], elements['tss_genes'].reindex(
                             np.array(arrays.nodes, dtype=object)[tss_idx]).to_numpy())
    nearest, nearest_distance = [], []
    for node, own, chrom, m in zip(candidates['node'], tss_genes, candidates['chrom'], mid):
        if own:
            nearest.append(own)
            nearest_distance.append(0)
        else:
            g, dist = tss_index.nearest(chrom, m, exclude={node})
            nearest.append(g)
            nearest_distance.append(dist)
    candidates['nearest_gene'] = nearest
    candidates['nearest_gene_distance'] = nearest_distance
    candidates['gene_label'] = candidates['genes'].astype(str)
    # own change of the explaining mark and its role sign
    own_change, mark_sign, window_marks = [], [], []
    for r in candidates.itertuples():
        mark = r.explaining_mark
        if mark == 'state' and not r.own_promoter:
            # state fragment: the signed mark whose own change at the fragment is most concordant with the gene
            fl = _feature(max(1, r.hops))
            d_gene = np.sign(r.rna_log2FC)
            options = [(signs.get(fl, {}).get(m, 0) * d_gene *
                        float(np.log2((arrays.vectors[m][1][r.node_index] + 1) /
                                      (arrays.vectors[m][0][r.node_index] + 1))), m)
                       for m in arrays.marks if signs.get(fl, {}).get(m, 0)]
            # an active mark that changes concordantly is preferred (silencing the element is the experiment);
            # otherwise the most concordant mark (a repressive mark: deletion)
            active = [o for o in options if signs.get(fl, {}).get(o[1], 0) > 0 and o[0] > 0]
            if active:
                mark = max(active)[1]
            elif options:
                mark = max(options)[1]
        window_marks.append(mark)
        if mark in arrays.vectors:
            va, vb = arrays.vectors[mark]
            own_change.append(float(np.log2((vb[r.node_index] + 1) / (va[r.node_index] + 1))))
            fl = 'promoter' if r.hops == 0 else _feature(r.hops)
            mark_sign.append(signs.get(fl, {}).get(mark, 0))
        elif mark == 'contact':
            own_change.append(float(np.log2((r.path_weight_B + 1) / (r.path_weight_A + 1))))
            mark_sign.append(signs.get('promoter', {}).get('contact', 1))
        else:
            own_change.append(np.nan)
            mark_sign.append(1)
    candidates['own_change'] = own_change
    candidates['mark_sign'] = mark_sign
    candidates['window_mark'] = window_marks
    # level evidence from the egaDifferential call
    level_keys = list(zip(levels['node'], levels['level']))
    keys = list(zip(candidates['gene_node'], candidates['level'].str.split('+').str[0]))

    def level_value(column):
        if column not in levels:
            return np.full(len(keys), np.nan, dtype=object)
        values = dict(zip(level_keys, levels[column]))
        return np.array([values.get(k, np.nan) for k in keys], dtype=object)

    candidates['p_level'] = pd.to_numeric(pd.Series(level_value('p_' + roles)), errors='coerce').to_numpy()
    candidates['q_level'] = pd.to_numeric(pd.Series(level_value('q_' + roles)), errors='coerce').to_numpy()
    candidates['beyond_mirror'] = level_value('beyond_mirror_' + roles)
    candidates['proximity_support'] = level_value('proximity_support')
    candidates['proximity_support_ci'] = level_value('proximity_support_ci')
    d = np.sign(candidates['rna_log2FC'].to_numpy())
    candidates['ev_own_change'] = (np.sign(candidates['own_change']) * candidates['mark_sign'] == d) & \
        (candidates['own_change'].abs() >= threshold)
    candidates['ev_level_p'] = candidates['p_level'] < 0.05
    candidates['ev_beyond_mirror'] = candidates['beyond_mirror'].astype(str).str.lower() == 'true'
    ci_low = pd.to_numeric(candidates['proximity_support_ci'].astype(str).str.split('-').str[0], errors='coerce')
    candidates['ev_support_3d'] = (ci_low > 0.5) & (candidates['hops'] > 0)
    candidates['ev_contact_both'] = (candidates['path_weight_A'] > 0) & (candidates['path_weight_B'] > 0)
    candidates['ev_non_promoter'] = ~candidates['own_promoter'].astype(bool) & \
        (candidates['element_type'] != 'promoter')
    strong = pleiotropy(candidates, args.effect if args.effect != 'remove' else 'flatten', args.strong)
    candidates['n_genes_strong'] = [len(strong.get(n, [])) for n in candidates['node']]
    candidates['also_explains'] = [', '.join(g for g in strong.get(n, []) if g != label)
                                   for n, label in zip(candidates['node'], candidates['gene_label'])]
    candidates['pleiotropic'] = candidates['also_explains'].str.len() > 0
    candidates['ev_not_pleiotropic'] = ~candidates['pleiotropic']
    write_table(candidates.drop(columns=['node_index']), ensure_parent(args.outPrefix + '_design_candidates.tsv.gz'))

    # rank, position, expected direction, sentence
    ranked = rank_candidates(candidates, args.effect)
    strands = read_gene_strands(args.tss) if args.tss else {}
    positions, directions, shorts = [], [], []
    for r in ranked.to_dict('records'):
        gene = str(r['gene_label']).split(',')[0]
        gnode = arrays.index[r['gene_node']]
        chrom, tss, strand = strands.get(gene, (arrays.chrom[gnode], (arrays.start[gnode] + arrays.end[gnode]) // 2,
                                                None))
        positions.append(relative_position(chrom, tss, strand, r['chrom'], (r['start'] + r['end']) // 2) + ' '
                         + gene)
        label, _ = expected_direction(r['mark_sign'], r['own_change'] if not pd.isna(r['own_change']) else 0.0,
                                      names, gene)
        directions.append(label if not r['own_promoter'] else 'down (gene silenced directly)')
        shorts.append(_short_class(r, names))
    ranked['position'] = positions
    ranked['expected_direction'] = directions
    ranked['element_short'] = shorts
    ranked['sentence'] = [candidate_sentence(r, args.effect) for r in ranked.to_dict('records')]
    if requested is not None:
        ranked = ranked[ranked['gene_node'].isin(set(requested))]
    columns = ['genes', 'gene_node', 'rna_log2FC', 'rank', 'node', 'chrom', 'start', 'end', 'element_class',
               'element_type', 'nearest_gene', 'nearest_gene_distance', 'distance_to_gene', 'position', 'level',
               'other_levels', 'hops', 'explaining_mark', 'own_change', 'path_weight_A', 'path_weight_B',
               'contribution', 'statistic_before', 'statistic_flatten', 'statistic_silence', 'statistic_remove',
               'effect_fraction', 'level_lost_flatten', 'level_lost_silence', 'level_lost_remove',
               'still_explained_flatten', 'still_explained_remove', 'p_level', 'q_level', 'p_before_reduced',
               'p_after_reduced', 'p_remove_reduced', 'beyond_mirror', 'proximity_support', 'proximity_support_ci',
               'state_fragments'] + ['ev_' + e for e in EVIDENCE] + \
              ['evidence_points', 'n_genes_strong', 'also_explains', 'pleiotropic', 'expected_direction',
               'design_score', 'sentence']
    columns = [c for c in columns if c in ranked]
    top = ranked[ranked['rank'] <= args.top][columns].rename(columns={'genes': 'gene'})
    write_table(top, ensure_parent(args.outPrefix + '_design.tsv'))

    # text
    lines = []
    for gene_node, part in top.groupby('gene_node', sort=False):
        first = part.iloc[0]
        explaining = by_gene[gene_node]
        levels_text = ', '.join(lv for lv in explaining.index if bool(explaining.at[lv, explains_col]))
        lines.append('{} ({}, log2FC {:+.1f}; explained at {})'.format(first['gene'], gene_node, first['rna_log2FC'],
                                                                     levels_text))
        for r in part.head(3).itertuples():
            lines.append('  {}. {}; expected: {}'.format(r.rank, r.sentence, r.expected_direction))
        lines.append('')
    with open(ensure_parent(args.outPrefix + '_design.txt'), 'w') as handle:
        handle.write('\n'.join(lines))

    # BED of the best nodes per gene (not the gene's own promoter)
    bed = ranked[~ranked['own_promoter'].astype(bool)].copy()
    bed['bed_rank'] = bed.groupby('gene_node').cumcount() + 1
    bed = bed[bed['bed_rank'] <= args.bedTop].copy()
    tracks_a, tracks_b = _tracks(args.tracksA), _tracks(args.tracksB)
    active = [m for m in marks if signs.get('promoter', {}).get(m, 0) > 0] or list(marks)
    bed['window_mark'] = [m if m in arrays.vectors else active[0] for m in bed['window_mark']]
    change = []
    for r in bed.itertuples():
        va, vb = arrays.vectors[r.window_mark]
        change.append(vb[r.node_index] - va[r.node_index])
    bed['window_condition'] = np.where(np.array(change) > 0, 'B', 'A')
    if len(bed):
        bed = bed_windows(bed, arrays, tracks_a, tracks_b, args.window, args.threads)
        out = pd.DataFrame({'chrom': bed['chrom'], 'start': bed['win_start'].astype(int),
                            'end': bed['win_end'].astype(int),
                            'name': [('{}|rank{}|{}|{}'.format(str(g).replace(',', '+'), r, n, w)).replace(' ', '_')
                                     for g, r, n, w in zip(bed['genes'], bed['rank'], bed['node'],
                                                           bed['win_source'])],
                            'score': np.clip((1000 * bed['effect_fraction']).round(), 0, 1000).astype(int),
                            'strand': '.'})
        out.to_csv(ensure_parent(args.outPrefix + '_design_top.bed'), sep='\t', header=False, index=False)

    # summary
    first_rank = ranked[(ranked['rank'] == 1) & ~ranked['own_promoter'].astype(bool)]
    beyond = candidates[~candidates['own_promoter'].astype(bool)]
    by_level = {}
    for level, part in beyond.groupby('level'):
        lost = part.groupby('gene_node')['level_lost_' + ('flatten' if args.effect == 'remove' else args.effect)] \
            .apply(lambda s: bool(s.astype('boolean').fillna(False).any()))
        by_level[level] = dict(genes=int(len(lost)), genes_single_node_removes_level=int(lost.sum()))
    summary = dict(comparison=os.path.basename(args.prefix), names=names, roles=roles, presence=presence,
                   effect=args.effect, permutations=args.permutations, genes=len(nodes),
                   genes_reported=int(ranked['gene_node'].nunique()), candidates=int(len(candidates)),
                   level_agreement=agreement, by_level=by_level,
                   pleiotropic_nodes=int(sum(len(v) > 1 for v in strong.values())),
                   genes_with_non_promoter_top=int(len(first_rank)), seconds=round(time.time() - start, 1))
    with open(ensure_parent(args.outPrefix + '_design_summary.json'), 'w') as handle:
        json.dump(summary, handle, indent=2)
    log.info('done: %s', json.dumps(summary))


def _feature(hops):
    from .differential import level_name
    return level_name(int(hops))


def _signs_at(signs, k):
    from .differential import level_name
    return signs.get(level_name(k), {})


def _short_class(r, names):
    """Element class of the condition where the explaining mark is higher, in short words for the sentence."""
    text = str(r['element_class'])
    if ' in both' in text:
        return text.replace(' in both', '')
    parts = re.match(r'(.+) in (\S+), (.+) in (\S+)$', text)
    if not parts:
        return text
    change = r['own_change'] if not pd.isna(r['own_change']) else 0.0
    return parts.group(3) if change > 0 else parts.group(1)


if __name__ == '__main__':
    main()
