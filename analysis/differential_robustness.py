"""Robustness of the gained, lost and shared edges of the differential graph and of its communities.

Subcommands

statuses     Recompute the union edges of egaDifferential (consensus_scores + differential_edges_from_scores, the same
             reading of one file per replicate) for a reference setting and for variants (another score threshold,
             replicate rule or input files, e.g. recalled from downsampled counts). For every edge status of the
             reference the fraction of edges that keep it, that become shared / gained / lost, or that are absent from
             the variant union graph is reported, and the reverse (the reference status of the variant edges).
             For the gained and lost edges of the reference the largest score of the pair in the replicate files of
             the condition where it is absent (files given by --rawA/--rawB, e.g. CHiCAGO exports at score >= 1) is
             binned, and the edges whose bait fragment is missing from replicates of that condition (CHiCAGO
             removes baits with fewer than minNPerBait reads) are counted; the status transfer is repeated for the
             edges whose bait is listed in every replicate file. The community part repeats egaDifferential's community analysis for each setting: Leiden
             communities of both condition graphs (resolution 1, seed 0, weight = score), community_transitions and
             community_expression_changes with the bait-level RNA log2FC mapped as in egaDifferential. Communities of
             a variant are compared with the recomputed reference (ARI, transition label of the best matching
             community) and the coherent units (q < 0.1) of the reference and of the stored egaDifferential output
             are looked up: a unit counts as recovered when the tested variant unit sharing most expressed bait
             fragments with it (Jaccard) is coherent, with the direction of its mean RNA log2FC.

downsample   Binomial thinning of a CHiCAGO input (chinput) file: every read count N of a bait / other-end pair is
             replaced by Binomial(N, p); pairs with 0 reads are dropped.

chicdiff     Overlap of the gained / lost edges of an egaDifferential edge table with Chicdiff calls (padj < --padj)
             by direction (log2FC of B over A), restricted to the interactions Chicdiff tested.

Usage: python differential_robustness.py <subcommand> --help
"""

import argparse
import json
import logging
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

from epigeneticgraphanalyzer.annotation import gene_values_to_nodes
from epigeneticgraphanalyzer.differential import (_key, community_expression_changes, community_transitions,
                                                  consensus_scores, differential_edges_from_scores, union_graph)
from epigeneticgraphanalyzer.egaCreateGraph import load_tss
from epigeneticgraphanalyzer.graph import annotate_genes, build_graph, compare_partitions, detect_communities
from epigeneticgraphanalyzer.io import read_expression, read_interactions

log = logging.getLogger('differential_robustness')
STATUSES = ['gained', 'lost', 'shared']
TRANSITIONS = ['stable', 'split', 'merged', 'remodelled', 'dissolved', 'formed']


# ----------------------------------------------------------------------------------------------------------- statuses

def parse_variant(text):
    """name:mode:threshold:A1,A2:B1,B2 (threshold 'none' for no threshold)."""
    name, mode, threshold, files_a, files_b = text.split(':')
    return {'name': name, 'mode': mode, 'threshold': None if threshold == 'none' else float(threshold),
            'files_a': files_a.split(','), 'files_b': files_b.split(',')}


_cache = {}


def read_one(path, file_format):
    if path not in _cache:
        _cache[path] = read_interactions([path], file_format)
    return _cache[path]


def condition_table(paths, file_format, mode, threshold):
    """As egaDifferential.main.condition (no --maxFdr)."""
    tables = [read_one(p, file_format) for p in paths]
    if len(tables) == 1 and threshold is not None:
        return consensus_scores([tables[0][tables[0]['score'] >= threshold]], 'union')
    return consensus_scores(tables, mode, threshold)


def edge_keys(edges):
    return edges['node_a'] + '|' + edges['node_b']


def max_raw_scores(paths, file_format):
    """Largest score per pair over the replicate files and the number of replicates listing the pair."""
    frames = []
    for path in paths:
        t = read_one(path, file_format)
        frames.append(pd.DataFrame({'key': _key(t).to_numpy(), 'score': t['score'].to_numpy()})
                      .sort_values('score', ascending=False).drop_duplicates('key'))
    joined = pd.concat(frames)
    return joined.groupby('key')['score'].agg(['max', 'count'])


def bait_sets(paths, file_format):
    """Bait fragments listed in each replicate file (CHiCAGO drops baits with too few reads, minNPerBait)."""
    out = []
    for path in paths:
        t = read_one(path, file_format)
        out.append(set(t['chrom'] + ':' + t['bait_start'].astype(str) + '-' + t['bait_end'].astype(str)))
    return out


def status_transfer(reference, variant, name, baits=None):
    """Rows per reference status: where its edges go in the variant (and the reverse).

    With `baits` (one set per replicate file of both conditions), the transfer is repeated for the edges whose bait
    is listed in every replicate (subset 'bait_in_all_replicates'), i.e. not removed by CHiCAGO's bait filter.
    """
    ref = pd.Series(reference['status'].to_numpy(), index=edge_keys(reference))
    var = pd.Series(variant['status'].to_numpy(), index=edge_keys(variant))
    rows = []
    subsets = [('all', None)]
    if baits is not None:
        subsets.append(('bait_in_all_replicates', baits))
    for subset, sets in subsets:
        for direction, source, target in (('reference_to_variant', ref, var), ('variant_to_reference', var, ref)):
            if sets is not None:
                bait = source.index.str.split('|').str[0]
                keep = np.ones(len(source), dtype=bool)
                for s in sets:
                    keep &= bait.isin(s)
                src = source[keep]
            else:
                src = source
            rows += _transfer_rows(src, target, name, direction, subset)
    return rows


def _transfer_rows(source, target, name, direction, subset):
    rows = []
    mapped = target.reindex(source.index).fillna('absent')
    for status in STATUSES:
        selected = mapped[source == status]
        n = len(selected)
        row = {'variant': name, 'subset': subset, 'direction': direction, 'status': status, 'n': n}
        for target_status in STATUSES + ['absent']:
            count = int((selected == target_status).sum())
            row['n_' + target_status] = count
            row['fraction_' + target_status] = count / n if n else np.nan
        rows.append(row)
    return rows


def retention_by_score(reference, variant, name, names, bins):
    """Status of the reference gained / lost edges in the variant by their score in the condition where present."""
    var = pd.Series(variant['status'].to_numpy(), index=edge_keys(variant))
    mapped = var.reindex(edge_keys(reference)).fillna('absent').to_numpy()
    edges = [-np.inf] + list(bins) + [np.inf]
    rows = []
    for status, present in (('gained', names[1]), ('lost', names[0])):
        selected = (reference['status'] == status).to_numpy()
        score = reference['score_' + present].to_numpy()[selected]
        target = mapped[selected]
        for lo, hi in zip(edges[:-1], edges[1:]):
            mask = (score >= lo) & (score < hi)
            n = int(mask.sum())
            if not n:
                continue
            row = {'variant': name, 'status': status, 'score_present_condition': '[{:g},{:g})'.format(lo, hi), 'n': n}
            for t in STATUSES + ['absent']:
                row['fraction_' + t] = float((target[mask] == t).mean())
            rows.append(row)
    return rows


def other_condition_scores(reference, raw_a, raw_b, bins):
    """For gained (absent in A) and lost (absent in B) reference edges: largest raw score in the absent condition."""
    keys = edge_keys(reference)
    rows = []
    labels = ['not listed'] + ['[{:g},{:g})'.format(lo, hi) for lo, hi in zip(bins[:-1], bins[1:])] + \
             ['>={:g}'.format(bins[-1])]
    for status, raw, absent in (('gained', raw_a, 'A'), ('lost', raw_b, 'B')):
        selected = keys[(reference['status'] == status).to_numpy()]
        best = raw['max'].reindex(selected)
        listed_in = raw['count'].reindex(selected).fillna(0).astype(int)
        values = best.to_numpy()
        index = np.where(np.isnan(values), 0, np.searchsorted(bins, np.nan_to_num(values, nan=-1), side='right'))
        for i, label in enumerate(labels):
            mask = index == i
            rows.append({'status': status, 'absent_in': absent, 'max_score_other_condition': label,
                         'n': int(mask.sum()), 'fraction': float(mask.mean()) if len(mask) else np.nan,
                         'n_listed_in_all_replicates': int(((listed_in == raw['count'].max()) & mask).sum())})
    return rows


def expression_per_node(edges, names, tss, expression, args):
    """Bait-level RNA log2FC of the union graph, as in egaDifferential."""
    graph = union_graph(edges, *names)
    nodes = pd.DataFrame([(n, d['chrom'], d['start'], d['end']) for n, d in graph.nodes(data=True)],
                         columns=['node', 'chrom', 'start', 'end'])
    genes = annotate_genes(nodes, tss)
    for node in graph.nodes:
        gene_names, ids = genes.get(node, ([], []))
        graph.nodes[node]['genes'] = ','.join(gene_names)
        graph.nodes[node]['gene_ids'] = ','.join(ids)
    change = gene_values_to_nodes(graph, expression, key=args.geneKey)
    return {n: v for n, v in change.items() if graph.nodes[n].get('bait')}


def unit_nodes(changes, communities_a, communities_b, expression):
    out = []
    for row in changes.itertuples():
        nodes = set(communities_a[row.community_a])
        if row.community_b >= 0:
            nodes |= communities_b[row.community_b]
        out.append({n for n in nodes if n in expression})
    return out


def stored_units(prefix):
    """Coherent units (q < 0.1) of a stored egaDifferential run as sets of expressed bait fragments."""
    changes = pd.read_csv(prefix + '_community_expression.tsv', sep='\t')
    genes = pd.read_csv(prefix + '_gene_community_changes.tsv', sep='\t')
    genes = genes[genes['rna_log2FC'].notna()]
    units = []
    for row in changes[changes['coherence_q'] < 0.1].itertuples():
        nodes = set(genes.loc[(genes['community_a'] == row.community_a) |
                              ((genes['community_b'] == row.community_b) & (row.community_b >= 0)), 'node'])
        units.append({'communities_a': str(row.communities_a), 'n_genes': int(row.n_genes),
                      'n_nodes_found': len(nodes), 'mean_rna_log2FC': float(row.mean_rna_log2FC),
                      'q': float(row.coherence_q), 'transition': row.transition, 'nodes': nodes})
    return units


def recover(units, changes, variant_units, label, variant_name):
    rows = []
    coherent = [i for i in range(len(changes)) if changes['coherence_q'].iloc[i] < 0.1]
    for unit in units:
        ref_nodes = unit['nodes']
        jaccards = [len(ref_nodes & g) / len(ref_nodes | g) if (ref_nodes | g) else 0.0 for g in variant_units]
        best = int(np.argmax(jaccards)) if jaccards else -1
        held = max([len(ref_nodes & variant_units[i]) / len(ref_nodes) for i in coherent], default=0.0) \
            if ref_nodes else 0.0
        best_row = changes.iloc[best] if best >= 0 else None
        rows.append({'variant': variant_name, 'reference': label, 'reference_communities_a': unit['communities_a'],
                     'reference_n_genes': unit['n_genes'], 'reference_mean_rna_log2FC': unit['mean_rna_log2FC'],
                     'reference_q': unit['q'], 'reference_transition': unit['transition'],
                     'best_match_jaccard': jaccards[best] if best >= 0 else 0.0,
                     'best_match_n_genes': int(best_row['n_genes']) if best_row is not None else 0,
                     'best_match_mean_rna_log2FC': float(best_row['mean_rna_log2FC']) if best_row is not None else np.nan,
                     'best_match_q': float(best_row['coherence_q']) if best_row is not None else np.nan,
                     'best_match_transition': best_row['transition'] if best_row is not None else '',
                     'recovered': bool(best_row is not None and best_row['coherence_q'] < 0.1),
                     'same_direction': bool(best_row is not None and np.sign(best_row['mean_rna_log2FC']) ==
                                            np.sign(unit['mean_rna_log2FC'])),
                     'max_fraction_in_coherent': held, 'contained': bool(held >= 0.5)})
    return rows


def statuses(args):
    names = args.names
    variants = [parse_variant(v) for v in args.variant]
    tss = load_tss(args.tss, args.geneTypes)
    table = read_expression(args.expression)
    expression = table[args.log2FCColumn].astype(float)
    stored = stored_units(args.storedPrefix) if args.storedPrefix else []
    for unit in stored:
        if unit['n_nodes_found'] != unit['n_genes']:
            log.warning('stored unit %s: %d nodes found for %d genes', unit['communities_a'], unit['n_nodes_found'],
                        unit['n_genes'])
    reference = None
    transfer_rows, community_rows, recovery_rows, label_rows, retention_rows = [], [], [], [], []
    for variant in variants:
        start = time.time()
        table_a = condition_table(variant['files_a'], args.format, variant['mode'], variant['threshold'])
        table_b = condition_table(variant['files_b'], args.format, variant['mode'], variant['threshold'])
        edges = differential_edges_from_scores(table_a, table_b, names[0], names[1], args.referenceScore)
        counts = edges['status'].value_counts().to_dict()
        log.info('%s: %d union edges %s', variant['name'], len(edges), counts)
        if reference is None:
            reference = {'edges': edges}
            if args.storedPrefix:
                stored_edges = pd.read_csv(args.storedPrefix + '_edges.tsv.gz', sep='\t',
                                           usecols=['node_a', 'node_b', 'status'])
                same = pd.Series(stored_edges['status'].to_numpy(), index=edge_keys(stored_edges)).reindex(
                    edge_keys(edges)).eq(edges['status'].to_numpy())
                reference['reproduces_stored_edges'] = bool(len(stored_edges) == len(edges) and same.all())
                log.info('reference reproduces the stored edge table: %s', reference['reproduces_stored_edges'])
            raw_a = max_raw_scores(args.rawA, args.format)
            raw_b = max_raw_scores(args.rawB, args.format)
            baits_a, baits_b = bait_sets(args.rawA, args.format), bait_sets(args.rawB, args.format)
            reference['baits'] = baits_a + baits_b
            # gained / lost edges whose bait is missing from replicates of the condition where the edge is absent
            bait = edges['node_a']
            reference['bait_filter'] = {}
            for status, sets in (('gained', baits_a), ('lost', baits_b)):
                selected = bait[(edges['status'] == status).to_numpy()]
                missing = np.column_stack([~selected.isin(s).to_numpy() for s in sets])
                reference['bait_filter'][status] = {
                    'n': int(len(selected)), 'bait_missing_in_any_replicate_of_absent_condition': int(missing.any(1).sum()),
                    'bait_missing_in_all_replicates_of_absent_condition': int(missing.all(1).sum())}
            log.info('bait filter: %s', reference['bait_filter'])
            pd.DataFrame(other_condition_scores(edges, raw_a, raw_b, args.scoreBins)).to_csv(
                args.output + '_other_condition_scores.tsv', sep='\t', index=False)
        transfer_rows += status_transfer(reference['edges'], edges, variant['name'], reference['baits'])
        retention_rows += retention_by_score(reference['edges'], edges, variant['name'], names, args.retentionBins)
        row = {'variant': variant['name'], 'mode': variant['mode'], 'threshold': variant['threshold'],
               'files_a': ','.join(map(os.path.basename, variant['files_a'])),
               'files_b': ','.join(map(os.path.basename, variant['files_b'])),
               'n_edges': len(edges), **{'n_' + s: int(counts.get(s, 0)) for s in STATUSES}}
        if args.skipCommunities:
            community_rows.append(row)
            continue
        change = expression_per_node(edges, names, tss, expression, args)
        graph_a = build_graph(table_a, tss=None, weight='score')
        graph_b = build_graph(table_b, tss=None, weight='score')
        communities_a = detect_communities(graph_a, resolution=1.0, seed=args.seed)
        communities_b = detect_communities(graph_b, resolution=1.0, seed=args.seed)
        transitions, main_b = community_transitions(communities_a, communities_b, set(graph_a.nodes),
                                                    set(graph_b.nodes))
        changes = community_expression_changes(communities_a, communities_b, main_b, transitions, edges, change,
                                               n_permutations=args.communityPermutations, seed=args.seed)
        changes.to_csv('{}_community_expression_{}.tsv'.format(args.output, variant['name']), sep='\t', index=False)
        units = unit_nodes(changes, communities_a, communities_b, change)
        coherent = changes[changes['coherence_q'] < 0.1]
        rho, p = stats.spearmanr(changes['rewired_fraction'], changes['mean_rna_log2FC'].abs(), nan_policy='omit')
        labels_a = transitions[transitions['side'] == 'A'].set_index('community_a')['transition']
        formed_b = set(transitions.loc[transitions['side'] == 'B', 'community_b'])
        if 'communities' not in reference:
            reference.update(communities={'a': communities_a, 'b': communities_b}, labels_a=labels_a,
                             formed_b=formed_b,
                             units=[{'communities_a': str(r.communities_a), 'n_genes': int(r.n_genes),
                                     'mean_rna_log2FC': float(r.mean_rna_log2FC), 'q': float(r.coherence_q),
                                     'transition': r.transition, 'nodes': units[i]}
                                    for i, r in zip(coherent.index, coherent.itertuples())])
        row.update(n_bait_nodes_expression=len(change), n_nodes_a=graph_a.number_of_nodes(),
                   n_nodes_b=graph_b.number_of_nodes(), n_communities_a=len(communities_a),
                   n_communities_b=len(communities_b))
        for side, communities, ref in (('a', communities_a, reference['communities']['a']),
                                       ('b', communities_b, reference['communities']['b'])):
            similarity, matching = compare_partitions(ref, communities)
            row.update({'{}_{}'.format(k, side): v for k, v in similarity.items()})
            if side == 'a':
                matching_a = matching.set_index('community')
            else:
                matching_b = matching.set_index('community')
        same_total = n_total = 0
        for label in TRANSITIONS:
            if label == 'formed':
                ref_ids = sorted(reference['formed_b'])
                same = matching_b.loc[ref_ids, 'best_match'].isin(formed_b)
            else:
                ref_ids = list(reference['labels_a'][reference['labels_a'] == label].index)
                same = matching_a.loc[ref_ids, 'best_match'].map(labels_a).eq(label)
            same_total += int(same.sum())
            n_total += len(ref_ids)
            label_rows.append({'variant': variant['name'], 'transition': label, 'n_reference': len(ref_ids),
                               'n_variant': int((transitions['transition'] == label).sum()),
                               'n_same_label': int(same.sum()),
                               'fraction_same_label': float(same.mean()) if len(ref_ids) else np.nan})
            row['n_' + label] = int((transitions['transition'] == label).sum())
        row['fraction_same_transition_label'] = same_total / n_total if n_total else np.nan
        row.update(communities_tested=int(len(changes)), coherent_q_lt_0_1=int(len(coherent)),
                   rewiring_abs_expression_rho=float(rho), rewiring_abs_expression_p=float(p))
        for label, ref_units in (('recomputed', reference['units']), ('stored', stored)):
            if not ref_units:
                continue
            rec = recover(ref_units, changes, units, label, variant['name'])
            recovery_rows += rec
            row['{}_coherent'.format(label)] = len(ref_units)
            row['{}_coherent_recovered'.format(label)] = sum(r['recovered'] for r in rec)
            row['{}_coherent_recovered_same_direction'.format(label)] = sum(r['recovered'] and r['same_direction']
                                                                            for r in rec)
            row['{}_coherent_contained'.format(label)] = sum(r['contained'] for r in rec)
        row['seconds'] = round(time.time() - start, 1)
        community_rows.append(row)
        log.info('%s: %s', variant['name'], {k: row[k] for k in row if k.startswith(('coherent', 'recomputed',
                                                                                      'stored', 'fraction', 'ari'))})
    pd.DataFrame(transfer_rows).to_csv(args.output + '_status_transfer.tsv', sep='\t', index=False)
    pd.DataFrame(retention_rows).to_csv(args.output + '_retention_by_score.tsv', sep='\t', index=False)
    pd.DataFrame(community_rows).to_csv(args.output + '_variants.tsv', sep='\t', index=False)
    if label_rows:
        pd.DataFrame(label_rows).to_csv(args.output + '_transition_labels.tsv', sep='\t', index=False)
    if recovery_rows:
        pd.DataFrame(recovery_rows).to_csv(args.output + '_coherent_recovery.tsv', sep='\t', index=False)
    summary = {'reproduces_stored_edges': reference.get('reproduces_stored_edges'),
               'bait_filter': reference['bait_filter'],
               'stored_coherent_units': [{k: v for k, v in u.items() if k != 'nodes'} for u in stored],
               'python_hash_seed': os.environ.get('PYTHONHASHSEED')}
    with open(args.output + '_summary.json', 'w') as handle:
        json.dump(summary, handle, indent=2)


# --------------------------------------------------------------------------------------------------------- downsample

def downsample(args):
    rng = np.random.default_rng(args.seed)
    with open(args.chinput) as handle:
        header = [handle.readline(), handle.readline()]
    total_in = total_out = 0
    with open(args.output, 'w') as out:
        out.writelines(header)
        reader = pd.read_csv(args.chinput, sep='\t', skiprows=2, header=None, chunksize=5_000_000,
                             dtype={0: np.int64, 1: np.int64, 2: np.int64, 3: np.int64, 4: str}, keep_default_na=False)
        for chunk in reader:
            counts = chunk[2].to_numpy()
            thinned = rng.binomial(counts, args.fraction)
            total_in += int(counts.sum())
            total_out += int(thinned.sum())
            chunk[2] = thinned
            chunk[thinned > 0].to_csv(out, sep='\t', header=False, index=False)
    log.info('%s: %d -> %d reads (%.4f, target %.4f)', args.chinput, total_in, total_out, total_out / total_in,
             args.fraction)


# ----------------------------------------------------------------------------------------------------------- chicdiff

def chicdiff(args):
    from epigeneticgraphanalyzer.differential import annotate_chicdiff
    from epigeneticgraphanalyzer.io import read_chicdiff
    edges = pd.read_csv(args.edges, sep='\t', usecols=['node_a', 'node_b', 'status', 'other_is_bait', 'distance',
                                                       'log2FC'])
    edges = edges.rename(columns={'log2FC': 'score_log2FC'})
    edges['log2FC'] = np.nan
    edges['pvalue'] = np.nan
    edges['qvalue'] = np.nan
    calls = read_chicdiff(args.chicdiff)
    if args.chromosomes:
        edges = edges[edges['node_a'].str.split(':').str[0].isin(args.chromosomes)].reset_index(drop=True)
    annotated = annotate_chicdiff(edges, calls, args.padj)
    tested = annotated[annotated['chicdiff_tested']].copy()
    # orientation of the Chicdiff fold change (B over A expected): checked against the score log2 ratio of EGA
    rho = float(stats.spearmanr(tested['log2FC'], tested['score_log2FC'], nan_policy='omit')[0])
    flipped = rho < 0
    if flipped:
        tested['log2FC'] = -tested['log2FC']
    log.info('Spearman rho of Chicdiff log2FC and EGA score log2 ratio: %.3f (sign flipped: %s)', rho, flipped)
    # direction of the Chicdiff fold change by EGA status, and nominal (weighted) p < 0.05 by direction
    by_status = {}
    for status in STATUSES:
        sel = tested[tested['status'] == status]
        by_status[status] = {
            'n_tested': int(len(sel)), 'median_log2FC': float(sel['log2FC'].median()) if len(sel) else None,
            'fraction_log2FC_positive': float((sel['log2FC'] > 0).mean()) if len(sel) else None,
            'n_p_lt_0.05_up': int(((sel['pvalue'] < 0.05) & (sel['log2FC'] > 0)).sum()),
            'n_p_lt_0.05_down': int(((sel['pvalue'] < 0.05) & (sel['log2FC'] < 0)).sum()),
            'n_padj_lt_{:g}'.format(args.padj): int(sel['chicdiff_significant'].sum())}
    with open(args.output + '_chicdiff_summary.json', 'w') as handle:
        json.dump({'union_edges': int(len(annotated)), 'tested': int(len(tested)),
                   'significant': int(tested['chicdiff_significant'].sum()), 'padj': args.padj,
                   'chicdiff_regions_padj_lt': int((calls['padj'] < args.padj).sum()),
                   'chicdiff_regions_min_padj': float(calls['padj'].min()),
                   'by_status': by_status,
                   'chicdiff_regions': int(len(calls)), 'spearman_log2FC_vs_score_ratio': rho,
                   'sign_flipped': bool(flipped),
                   'tested_by_status': tested['status'].value_counts().to_dict(),
                   'union_by_status': annotated['status'].value_counts().to_dict()}, handle, indent=2)
    up = tested['chicdiff_significant'] & (tested['log2FC'] > 0)
    down = tested['chicdiff_significant'] & (tested['log2FC'] < 0)
    rows = []
    for status, call, label in (('gained', up, 'up in B'), ('lost', down, 'down in B')):
        ega = tested['status'] == status
        tp = int((ega & call).sum())
        rows.append({'ega_status': status, 'chicdiff_call': label, 'n_union_edges': int((annotated['status'] ==
                                                                                        status).sum()),
                     'n_tested': int(ega.sum()), 'n_chicdiff': int(call.sum()), 'n_both': tp,
                     'precision_ega_vs_chicdiff': tp / ega.sum() if ega.sum() else np.nan,
                     'recall_ega_vs_chicdiff': tp / call.sum() if call.sum() else np.nan,
                     'n_opposite_direction': int((ega & (down if status == 'gained' else up)).sum()),
                     'n_shared_with_call': int(((tested['status'] == 'shared') & call).sum())})
    table = pd.DataFrame(rows)
    table.to_csv(args.output + '_chicdiff_overlap.tsv', sep='\t', index=False)
    crosstab = pd.crosstab(tested['status'], np.where(up, 'up', np.where(down, 'down', 'ns')))
    crosstab.to_csv(args.output + '_chicdiff_crosstab.tsv', sep='\t')
    log.info('Chicdiff tested %d of %d union edges\n%s\n%s', len(tested), len(annotated), table.to_string(),
             crosstab.to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    s = sub.add_parser('statuses')
    s.add_argument('--variant', action='append', required=True,
                   help='name:mode:threshold:A1,A2:B1,B2; the first is the reference')
    s.add_argument('--format', default='ibed')
    s.add_argument('--names', nargs=2, required=True)
    s.add_argument('--rawA', nargs='+', required=True, help='Replicate files of A with the lowest exported score')
    s.add_argument('--rawB', nargs='+', required=True)
    s.add_argument('--scoreBins', nargs='+', type=float, default=[1, 3, 5])
    s.add_argument('--retentionBins', nargs='+', type=float, default=[8, 12, 20],
                   help='Score bins (condition where the edge is present) of the retention table')
    s.add_argument('--referenceScore', type=float, default=5.0)
    s.add_argument('--tss', required=True)
    s.add_argument('--geneTypes', nargs='*', default=['protein_coding', 'lincRNA', 'antisense'])
    s.add_argument('--expression', required=True)
    s.add_argument('--log2FCColumn', default='log2FoldChange')
    s.add_argument('--geneKey', default='genes')
    s.add_argument('--communityPermutations', type=int, default=10000)
    s.add_argument('--seed', type=int, default=0)
    s.add_argument('--storedPrefix', help='Output prefix of the stored egaDifferential run')
    s.add_argument('--skipCommunities', action='store_true')
    s.add_argument('-o', '--output', required=True)
    d = sub.add_parser('downsample')
    d.add_argument('chinput')
    d.add_argument('fraction', type=float)
    d.add_argument('output')
    d.add_argument('--seed', type=int, default=0)
    c = sub.add_parser('chicdiff')
    c.add_argument('--edges', required=True)
    c.add_argument('--chicdiff', required=True)
    c.add_argument('--padj', type=float, default=0.05)
    c.add_argument('--chromosomes', nargs='*', help='Restrict the union edges to these bait chromosomes')
    c.add_argument('-o', '--output', required=True)
    args = parser.parse_args()
    if 'PYTHONHASHSEED' not in os.environ:
        os.execve(sys.executable, [sys.executable] + sys.argv, dict(os.environ, PYTHONHASHSEED='0'))
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    {'statuses': statuses, 'downsample': downsample, 'chicdiff': chicdiff}[args.command](args)


if __name__ == '__main__':
    main()
