"""egaTimecourse: trajectories and ordering of contact, mark and expression changes over ordered conditions
(e.g. days of differentiation) from pairwise egaDifferential runs, and knockdown / UMI-4C checks of candidate
trigger factors (egaTriggers).

Inputs: the egaDifferential prefixes <differential>/<a>_vs_<b> of every consecutive pair of --days and of the first
and last day; DESeq2-like tables <expressionDir>/deseq2_<b>_vs_<a>.tsv (gene_name, log2FoldChange = b relative to
a, padj) and <expressionDir>/normalised_counts.tsv (gene_name and one column per sample whose name contains the day,
e.g. B1_D0); the called interactions of every day (--interactions day=rep1,rep2, consensus as in egaDifferential);
the mark bigWigs per day (--markTracks [mark:]day=file[,file]).

Outputs (<out> = --outputPrefix):
  <out>_explaining_nodes.tsv.gz   explaining nodes of every pairwise run (triggers.explaining_nodes, declared roles)
  <out>_trajectories.tsv          per differential gene node and day: expression, promoter marks, promoter contacts,
                                  contact-weighted marks at step 1..maxHops (log2(value + 1)); interval changes,
                                  expression pattern and the levels explaining the gene in every pairwise run
  <out>_units.tsv                 per gene and explaining node (step 0 = promoter): values per day, interval changes
                                  of contact, mark and expression and their event intervals
  <out>_ordering.tsv              order categories against the within-gene permutation null, per unit subset, layer
                                  set and chromatin threshold
  <out>_ordering_patterns.tsv     joint event-interval patterns
  <out>_knockdown_nodes.tsv       knockdown versus control change at every explaining node of the --triggers run
  <out>_knockdown_test.tsv        candidate-trigger genes versus other explained genes per knocked-down factor
  <out>_umi4c.tsv                 UMI-4C viewpoints overlapping explained genes: contact frequency with the
                                  explaining nodes per condition
"""

import logging
import os
import time

import numpy as np
import pandas as pd

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .timecourse import (LAYER_SETS, bigwig_scale, bin_totals, day_operators, day_tables, event_patterns,
                         expression_pattern, find_viewpoints, first_interval, intervals_of, knockdown_test,
                         log2_change, node_gene_choice, node_signals, order_test, profile_over, read_umi4c,
                         timecourse_graph, viewpoint_profile, UMI4C_NAME)

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser(__doc__.split('\n\n')[0])
    parser.add_argument('--days', nargs='+', required=True, help='Ordered condition names (at least three)')
    parser.add_argument('--differential', required=True, help='Directory of the egaDifferential runs <a>_vs_<b>')
    parser.add_argument('--expressionDir', required=True)
    parser.add_argument('--interactions', nargs='+', required=True, help='day=rep1[,rep2...] called interactions')
    parser.add_argument('--format', default='ibed')
    parser.add_argument('--baits', help='Baited fragments (as egaDifferential --baits)')
    parser.add_argument('--replicates', default='all', choices=('all', 'intersect', 'union'))
    parser.add_argument('--scoreThreshold', type=float, default=5.0)
    parser.add_argument('--referenceScore', type=float, default=5.0)
    parser.add_argument('--markTracks', nargs='+', required=True,
                        help='[mark:]day=file[,file] (mark default H3K27ac)')
    parser.add_argument('--orderMark', default='H3K27ac', help='Mark of the ordering test')
    parser.add_argument('--maxHops', type=int, default=3)
    parser.add_argument('--rnaThreshold', type=float, default=1.0)
    parser.add_argument('--alpha', type=float, default=0.05)
    parser.add_argument('--chromatinThreshold', type=float, default=1.0)
    parser.add_argument('--contactThreshold', type=float, default=1.0)
    parser.add_argument('--sensitivityThresholds', type=float, nargs='*', default=[0.5, 1.5],
                        help='Further chromatin / contact thresholds of the ordering test (all units)')
    parser.add_argument('--triggers', help='egaTriggers output prefix (gene_triggers, explaining_nodes) of the '
                                           'day pair of the knockdowns')
    parser.add_argument('--knockdownControl', nargs='*', default=[], help='Control bigWigs (replicates)')
    parser.add_argument('--knockdown', nargs='*', default=[], help='FACTOR=rep1.bw[,rep2.bw] knockdown bigWigs')
    parser.add_argument('--motifs', help='JASPAR file: reports whether a knocked-down factor has a motif')
    parser.add_argument('--umi4c', nargs='*', default=[], help='UMI-4C pair count files (names with '
                                                               'UMI4C_KD_<cond>i_Day<d>_Rep<r>)')
    parser.add_argument('--umi4cWindow', type=int, default=1000000)
    parser.add_argument('--permutations', type=int, default=10000)
    parser.add_argument('--bootstrap', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--outputPrefix', '-o', required=True)
    parsed = parser.parse_args(args)
    if len(parsed.days) < 3:
        parser.error('--days needs at least three conditions')
    return parsed


def parse_day_files(entries, default_mark=None):
    result = {}
    for entry in entries:
        key, files = entry.split('=', 1)
        if default_mark is not None:
            mark, day = key.split(':', 1) if ':' in key else (default_mark, key)
            result.setdefault(mark, {})[day] = [f for f in files.split(',') if f]
        else:
            result[key] = [f for f in files.split(',') if f]
    return result


def read_deseq(path):
    table = pd.read_csv(path, sep='\t')
    table = table.dropna(subset=['gene_name']).drop_duplicates('gene_name').set_index('gene_name')
    return table[['log2FoldChange', 'padj', 'baseMean']]


def pair_name(a, b):
    return '{}_vs_{}'.format(a, b)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    from .graph import read_graph
    from .io import read_baits
    from .triggers import explaining_nodes
    out = args.outputPrefix
    ensure_parent(out)
    days = args.days
    intervals = intervals_of(days)
    pairs = intervals + [(days[0], days[-1])]
    marks_by_day = parse_day_files(args.markTracks, default_mark='H3K27ac')
    if args.orderMark not in marks_by_day:
        raise ValueError('--orderMark {} has no --markTracks'.format(args.orderMark))

    # expression
    expression = {p: read_deseq(os.path.join(args.expressionDir, 'deseq2_{}_vs_{}.tsv'.format(p[1], p[0])))
                  for p in pairs}
    counts = pd.read_csv(os.path.join(args.expressionDir, 'normalised_counts.tsv'), sep='\t')
    counts = counts.dropna(subset=['gene_name']).drop_duplicates('gene_name').set_index('gene_name')
    day_columns = {d: [c for c in counts.columns if c != 'gene_id' and d in c.split('_')] for d in days}

    # pairwise runs: differential gene nodes, explanation flags, explaining nodes
    levels_by_pair, explaining = {}, []
    node_genes = {}
    for a, b in pairs:
        prefix = os.path.join(args.differential, pair_name(a, b))
        levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t')
        levels_by_pair[(a, b)] = levels
        node_genes.update(dict(zip(levels['node'], levels['genes'])))
        t = time.time()
        graph = read_graph(prefix + '_union.graphml')
        nodes = explaining_nodes(graph, levels, [a, b], max_hops=args.maxHops,
                                 chromatin_threshold=args.chromatinThreshold, roles='declared', presence='mean')
        nodes.insert(0, 'pair', pair_name(a, b))
        explaining.append(nodes[nodes['explaining']])
        log.info('%s: %d differential gene nodes, %d explaining nodes (%.0f s)', pair_name(a, b),
                 levels['node'].nunique(), int(nodes['explaining'].sum()), time.time() - t)
        del graph
    explaining = pd.concat(explaining, ignore_index=True)
    explaining.to_csv(out + '_explaining_nodes.tsv.gz', sep='\t', index=False)
    gene_nodes = list(node_genes)

    choice = node_gene_choice(node_genes, expression, args.rnaThreshold, args.alpha)
    gene_nodes = [n for n in gene_nodes if n in choice]

    # graph of all days, operators per day, node signals
    t = time.time()
    baits = read_baits(args.baits) if args.baits else None
    tables = day_tables(parse_day_files(args.interactions), args.format, baits, args.replicates, args.scoreThreshold)
    tables = {d: tables[d] for d in days}
    graph = timecourse_graph(tables, args.referenceScore)
    sources = [n for n in gene_nodes if n in graph]
    ops, hops = day_operators(graph, sources, days, max_hops=args.maxHops)
    node_table = pd.DataFrame([(n, d['chrom'], d['start'], d['end']) for n, d in graph.nodes(data=True)],
                              columns=['node', 'chrom', 'start', 'end']).set_index('node')
    signals = node_signals(ops.nodes, node_table, marks_by_day, threads=args.threads)
    log.info('graph of %d days: %d nodes, %d edges, %d gene nodes as sources (%.0f s)', len(days),
             graph.number_of_nodes(), graph.number_of_edges(), len(sources), time.time() - t)
    for d in days:
        log.info('day %s: %d called pairs', d, len(tables[d]))
    row_of = {s: i for i, s in enumerate(ops.sources)}
    col_of = ops.index

    # trajectories
    interval_lfc = np.array([[expression[p]['log2FoldChange'].get(choice[s], np.nan) for p in intervals]
                             for s in ops.sources])
    interval_padj = np.array([[expression[p]['padj'].get(choice[s], np.nan) for p in intervals]
                              for s in ops.sources])
    span_lfc = np.array([expression[pairs[-1]]['log2FoldChange'].get(choice[s], np.nan) for s in ops.sources])
    span_padj = np.array([expression[pairs[-1]]['padj'].get(choice[s], np.nan) for s in ops.sources])
    patterns = expression_pattern(np.nan_to_num(interval_lfc), np.nan_to_num(interval_padj, nan=1.0),
                                  args.rnaThreshold, args.alpha)
    flags = {}
    for (a, b), levels in levels_by_pair.items():
        sub = levels.drop_duplicates(['node', 'level']).set_index(['node', 'level'])['explains_declared_mean']
        flags[pair_name(a, b)] = sub
    level_names = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state'] + \
        ['step{}'.format(k) for k in range(2, args.maxHops + 1)]
    rows = []
    for s in ops.sources:
        i = row_of[s]
        gene = choice[s]
        rec = {'node': s, 'genes': node_genes[s], 'gene': gene, 'pattern': patterns[i]}
        for j, (a, b) in enumerate(intervals):
            rec['lfc_{}'.format(pair_name(a, b))] = interval_lfc[i, j]
            rec['padj_{}'.format(pair_name(a, b))] = interval_padj[i, j]
        rec['lfc_' + pair_name(*pairs[-1])] = span_lfc[i]
        rec['padj_' + pair_name(*pairs[-1])] = span_padj[i]
        for d in days:
            values = counts.loc[gene, day_columns[d]] if gene in counts.index else pd.Series(dtype=float)
            rec['expression_' + d] = float(np.log2(values.mean() + 1)) if len(values) else np.nan
            rec['contacts_' + d] = float(np.log2(hops[d][1].getrow(i).sum() + 1))
            for mark in marks_by_day:
                vector = signals[mark][d]
                rec['{}_promoter_{}'.format(mark, d)] = float(np.log2(vector[col_of[s]] + 1))
                for k in range(1, args.maxHops + 1):
                    rec['{}_step{}_{}'.format(mark, k, d)] = float(np.log2(hops[d][k].getrow(i).dot(vector)[0] + 1))
            rec['n_step_nodes_' + d] = int(sum(hops[d][k].getrow(i).nnz for k in range(1, args.maxHops + 1)))
        for pair, sub in flags.items():
            explained = [lv for lv in level_names if bool(sub.get((s, lv), False))]
            rec['explained_' + pair] = ','.join(explained)
            rec['differential_' + pair] = (s, 'promoter_marks') in sub.index
        rows.append(rec)
    trajectories = pd.DataFrame(rows)
    write_table(trajectories, out + '_trajectories.tsv')
    log.info('trajectories of %d gene nodes; patterns %s', len(trajectories),
             trajectories['pattern'].value_counts().to_dict())

    # units: gene node x explaining node (step 0 = promoter)
    mark = args.orderMark
    units = explaining[explaining['gene_node'].isin(row_of)].copy()
    units['pair_step'] = units['level'].map(lambda lv: 0 if lv == 'promoter' else int(lv.replace('step', '')))
    units = units.groupby(['gene_node', 'node'], as_index=False).agg(
        pair_step=('pair_step', 'min'), pairs=('pair', lambda s: ','.join(sorted(set(s)))))
    records = []
    for u in units.itertuples(index=False):
        i = row_of[u.gene_node]
        rec = {'gene_node': u.gene_node, 'gene': choice[u.gene_node], 'node': u.node, 'pairs': u.pairs,
               'pair_step': u.pair_step}
        if u.node == u.gene_node:
            step = 0
            contact = [hops[d][1].getrow(i).sum() for d in days]
        else:
            if u.node not in col_of:
                continue
            c = col_of[u.node]
            step = next((k for k in range(1, args.maxHops + 1) if any(hops[d][k][i, c] > 0 for d in days)), -1)
            if step < 0:
                continue
            contact = [hops[d][step][i, c] for d in days]
        rec['step'] = step
        for d, v in zip(days, contact):
            rec['contact_' + d] = float(v)
            rec['mark_' + d] = float(signals[mark][d][col_of[u.node]])
        records.append(rec)
    units = pd.DataFrame(records)
    for j, (a, b) in enumerate(intervals):
        units['d_contact_' + pair_name(a, b)] = log2_change(units['contact_' + a], units['contact_' + b])
        units['d_mark_' + pair_name(a, b)] = log2_change(units['mark_' + a], units['mark_' + b])
        units['d_expression_' + pair_name(a, b)] = units['gene_node'].map(
            lambda n, j=j: interval_lfc[row_of[n], j])
        units['padj_expression_' + pair_name(a, b)] = units['gene_node'].map(
            lambda n, j=j: interval_padj[row_of[n], j])
    names = [pair_name(a, b) for a, b in intervals]

    def event_times(table, chromatin, contact_threshold):
        result = table.copy()
        result['t_contact'] = first_interval(table[['d_contact_' + n for n in names]].to_numpy(), contact_threshold)
        result['t_mark'] = first_interval(table[['d_mark_' + n for n in names]].to_numpy(), chromatin)
        result['t_expression'] = first_interval(
            table[['d_expression_' + n for n in names]].to_numpy(), args.rnaThreshold,
            significant=(table[['padj_expression_' + n for n in names]].fillna(1.0).to_numpy() <= args.alpha))
        return result

    timed = event_times(units, args.chromatinThreshold, args.contactThreshold)
    write_table(timed, out + '_units.tsv')
    log.info('%d units (%s); events contact %d, mark %d, expression %d', len(timed),
             timed['step'].value_counts().sort_index().to_dict(), int((timed['t_contact'] >= 0).sum()),
             int((timed['t_mark'] >= 0).sum()), int((timed['t_expression'] >= 0).sum()))

    subsets = {'all': timed, 'promoter': timed[timed['step'] == 0], 'elements': timed[timed['step'] >= 1]}
    for k in range(1, args.maxHops + 1):
        subsets['step{}'.format(k)] = timed[timed['step'] == k]
    results, patterns_out = [], []
    seed = args.seed
    for subset, table in subsets.items():
        for set_name, layers in LAYER_SETS.items():
            res = order_test(table, layers, args.permutations, args.bootstrap, seed, gene_column='gene_node')
            results.append(res.assign(subset=subset, layers=set_name, chromatin_threshold=args.chromatinThreshold,
                                      contact_threshold=args.contactThreshold))
            patterns_out.append(event_patterns(table, layers).assign(subset=subset, layers=set_name))
    for threshold in args.sensitivityThresholds:
        table = event_times(units, threshold, threshold)
        for set_name, layers in LAYER_SETS.items():
            res = order_test(table, layers, args.permutations, args.bootstrap, seed, gene_column='gene_node')
            results.append(res.assign(subset='all', layers=set_name, chromatin_threshold=threshold,
                                      contact_threshold=threshold))
    ordering = pd.concat(results, ignore_index=True)
    lead = ['subset', 'layers', 'chromatin_threshold', 'contact_threshold']
    ordering = ordering[lead + [c for c in ordering.columns if c not in lead]]
    write_table(ordering, out + '_ordering.tsv')
    write_table(pd.concat(patterns_out, ignore_index=True), out + '_ordering_patterns.tsv')
    main_rows = ordering[(ordering['subset'] == 'all') & (ordering['layers'] == 'contact_mark_expression') &
                         (ordering['chromatin_threshold'] == args.chromatinThreshold)]
    log.info('ordering (all units, contact / mark / expression):\n%s', main_rows.to_string(index=False))

    # knockdown at the explaining nodes of the trigger pair
    if args.triggers and args.knockdown:
        knockdown(args, out)
    if args.umi4c and args.triggers:
        umi4c(args, out, node_table)


def knockdown(args, out):
    from .differential import _interval_signals
    from .triggers import motif_genes, read_jaspar
    triggers = pd.read_csv(args.triggers + '_gene_triggers.tsv', sep='\t')
    nodes = pd.read_csv(args.triggers + '_explaining_nodes.tsv.gz', sep='\t')
    nodes = nodes[nodes['explaining']].drop_duplicates(['gene_node', 'node'])
    factors = parse_day_files(args.knockdown)
    motif_factors = set()
    motif_nodes, scanned = None, set()
    if args.motifs:
        motifs = read_jaspar(args.motifs)
        motif_factors = {g for m in motifs for g in motif_genes(m['name'])}
        scan_path, window_path = args.triggers + '_scan.npz', args.triggers + '_windows.tsv.gz'
        if os.path.exists(scan_path) and os.path.exists(window_path):
            scan = np.load(scan_path, allow_pickle=False)
            windows = pd.read_csv(window_path, sep='\t', usecols=['node'])
            scanned = set(windows['node'])
            ids = list(scan['motif_ids'])
            present = scan['counts'] > 0
            motif_nodes = {}
            for factor in parse_day_files(args.knockdown):
                columns = [ids.index(m['id']) for m in motifs if factor.upper() in motif_genes(m['name'])
                           and m['id'] in ids]
                if columns:
                    motif_nodes[factor.upper()] = set(windows['node'][present[:, columns].any(axis=1)])
    unique = nodes.drop_duplicates('node')
    coords = unique['node'].str.extract(r'^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$')
    paths = list(args.knockdownControl) + [p for files in factors.values() for p in files]
    read = _interval_signals(paths, coords['chrom'].tolist(), coords['start'].astype(np.int64).to_numpy(),
                             coords['end'].astype(np.int64).to_numpy(), args.threads)
    scaled = {p: read[p] * bigwig_scale(p) for p in paths}
    table = pd.DataFrame({'node': unique['node'].to_numpy()})
    control = np.mean([scaled[p] for p in args.knockdownControl], axis=0)
    table['control'] = control
    for r, p in enumerate(args.knockdownControl):
        table['control_rep{}'.format(r + 1)] = scaled[p]
    gene_lfc = nodes.drop_duplicates('gene_node').set_index('gene_node')['rna_log2FC']
    rows = []
    for factor, files in factors.items():
        kd = np.mean([scaled[p] for p in files], axis=0)
        table[factor] = kd
        table['log2_' + factor] = log2_change(control, kd)
        for r, p in enumerate(files):
            if r < len(args.knockdownControl):
                table['log2_{}_rep{}'.format(factor, r + 1)] = log2_change(scaled[args.knockdownControl[r]], scaled[p])
        hits = triggers[triggers['tf'].astype(str).str.upper() == factor.upper()]
        groupings = [('any', set(hits['gene_node']), None),
                     ('named factor', set(hits.loc[hits['evidence'] == 'named factor', 'gene_node']), None)]
        table['candidate_' + factor] = table['node'].isin(nodes.loc[nodes['gene_node'].isin(groupings[0][1]),
                                                                        'node'])
        if motif_nodes is not None and factor.upper() in motif_nodes:
            groupings.append(('motif at node', None, motif_nodes[factor.upper()]))
            table['motif_' + factor] = table['node'].isin(motif_nodes[factor.upper()])
            table['scanned'] = table['node'].isin(scanned)
        for evidence, cand_genes, with_motif in groupings:
            for direction in ('all', 'up', 'down'):
                if direction == 'all':
                    keep = gene_lfc.index
                else:
                    keep = gene_lfc.index[(gene_lfc > 0) if direction == 'up' else (gene_lfc < 0)]
                pool = nodes[nodes['gene_node'].isin(keep)]
                if with_motif is None:
                    cand = pool[pool['gene_node'].isin(cand_genes)]
                    ctrl = pool[~pool['gene_node'].isin(cand_genes)]
                else:
                    # node level: explaining nodes whose scanned window carries a motif of the factor
                    cand = pool[pool['node'].isin(with_motif)]
                    ctrl = pool[pool['node'].isin(scanned) & ~pool['node'].isin(with_motif)]
                    cand_genes = set(cand['gene_node'])
                for column in ['log2_' + factor] + ['log2_{}_rep{}'.format(factor, r + 1) for r in range(len(files))
                                                     if r < len(args.knockdownControl)]:
                    values = dict(zip(table['node'], table[column]))
                    res = knockdown_test(values, cand, ctrl, args.bootstrap, args.seed)
                    rows.append(dict(factor=factor, motif_in_jaspar=(factor.upper() in motif_factors)
                                     if args.motifs else np.nan, evidence=evidence, genes=direction,
                                     measure='mean' if column == 'log2_' + factor
                                     else column.replace('log2_{}_'.format(factor), ''),
                                     candidate_genes_total=len(cand_genes), **res))
    write_table(table, out + '_knockdown_nodes.tsv')
    result = pd.DataFrame(rows)
    write_table(result, out + '_knockdown_test.tsv')
    log.info('knockdown test:\n%s', result[(result['measure'] == 'mean') & (result['evidence'] == 'any')][
        ['factor', 'genes', 'genes_candidate', 'nodes_candidate', 'nodes_control', 'median_difference', 'md_low',
         'md_high', 'cliffs_delta', 'p_two_sided']].to_string(index=False))


def umi4c(args, out, node_table):
    nodes = pd.read_csv(args.triggers + '_explaining_nodes.tsv.gz', sep='\t')
    samples = {}
    for path in args.umi4c:
        match = UMI4C_NAME.search(os.path.basename(path))
        if not match:
            log.warning('UMI-4C file name not recognised: %s', path)
            continue
        key = '{}_D{}'.format(match.group('condition'), match.group('day'))
        samples.setdefault(key, []).append(read_umi4c(path))
    pooled = {k: pd.concat(v, ignore_index=True) for k, v in samples.items()}
    reference = next((k for k in pooled if k.startswith('CTR') and not k.endswith('_D0')), next(iter(pooled)))
    viewpoints = find_viewpoints(bin_totals(pooled[reference]), min_fraction=0.001)
    log.info('UMI-4C viewpoints (%s): %s', reference, viewpoints)
    gene_nodes = nodes.drop_duplicates('gene_node')
    coords = gene_nodes['gene_node'].str.extract(r'^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$')
    coords[['start', 'end']] = coords[['start', 'end']].astype(np.int64)
    rows = []
    for chrom, start, end in viewpoints:
        hit = gene_nodes[(coords['chrom'] == chrom).to_numpy() & (coords['start'] < end + 5000).to_numpy() &
                         (coords['end'] > start - 5000).to_numpy()]
        profiles = {k: viewpoint_profile(t, (chrom, start, end), args.umi4cWindow) for k, t in pooled.items()}
        totals = {k: float(p.sum()) for k, p in profiles.items()}
        if not len(hit):
            rows.append({'viewpoint': '{}:{}-{}'.format(chrom, start, end), 'gene_node': '', 'genes': '',
                         'explained': False, **{'total_' + k: v for k, v in totals.items()}})
            continue
        for g in hit.itertuples(index=False):
            sub = nodes[(nodes['gene_node'] == g.gene_node) & (nodes['level'] != 'promoter')]
            for explaining_flag, part in sub.groupby('explaining'):
                rec = {'viewpoint': '{}:{}-{}'.format(chrom, start, end), 'gene_node': g.gene_node,
                       'genes': g.genes, 'rna_log2FC': g.rna_log2FC, 'explained': True,
                       'node_set': 'explaining' if explaining_flag else 'other step nodes', 'nodes': len(part)}
                intervals = []
                for n in part['node'].unique():
                    c, s_e = n.split(':')
                    s, e = (int(x) for x in s_e.split('-'))
                    if c == chrom and abs(s - start) <= args.umi4cWindow:
                        intervals.append((s, e))
                rec['nodes_in_window'] = len(intervals)
                for k, p in profiles.items():
                    reads = profile_over(p, intervals)
                    rec['reads_' + k] = reads
                    rec['total_' + k] = totals[k]
                    rec['freq_' + k] = reads / totals[k] if totals[k] else np.nan
                rows.append(rec)
    table = pd.DataFrame(rows)
    control = reference
    for k in pooled:
        if k == control or 'freq_' + k not in table:
            continue
        # knockdowns against the day-3 control; the day-0 control as day 3 against day 0
        num, den = (control, k) if k.endswith('_D0') else (k, control)
        table['log2_{}_vs_{}'.format(num, den)] = np.log2((table['freq_' + num] + 1e-5) /
                                                          (table['freq_' + den] + 1e-5))
    write_table(table, out + '_umi4c.tsv')
    log.info('UMI-4C: %d viewpoints, %d rows', len(viewpoints), len(table))


if __name__ == '__main__':
    main()
