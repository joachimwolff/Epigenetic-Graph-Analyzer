"""Baselines of the expression models: do graph neighbours explain expression changes better than features that use
the linear genome only?

For every comparison of an egaDifferential run (<run>/<cmp>_gene_features.tsv, _nested_models.tsv, _union.graphml,
_support.json, _gene_classes.tsv, _mark_roles.tsv) the nested models of differential.nested_models (OLS of
rna_log2FC on the promoter features, then + direct, + indirect, + hop3, on the 'informative' genes with direct and
indirect neighbours) are reproduced from the feature table and from the union graph. On the same genes and with the
same promoter model as base, the following feature blocks are added (one feature per mark, log2 B/A with
pseudo-count 1, as the existing features):

window   distance-weighted mark signal of the linear neighbourhood, without graph edges: the signal x_c(m) of every
         element m on the chromosome of the promoter fragment (the nodes of the union graph, or bins of the bigWig
         tracks) whose centre lies at a distance lo < d <= hi from the centre of the promoter fragment (elements
         overlapping the promoter fragment excluded), weighted by d^-1 (power-law decay with exponent 1):
             S_c = sum_m x_c(m) / d(m) / sum_m 1 / d(m),   feature = log2 (S_B + 1) / (S_A + 1).
         The weights are the same in both conditions (as the mark term of the step features). Dividing by the sum
         of the weights keeps S_c on the scale of the signal of one fragment, so that the pseudo-count acts as on
         the promoter feature; the ratio S_B / S_A is that of the weighted sums. Ranges: (0, largest distance of the
         gene's direct partners], (0, largest distance of its nodes up to three steps], and fixed windows of
         radius 50 kb, 250 kb, 1 Mb and 2 Mb, and the ring 50 kb - 2 Mb;
degree   number of interactions of the promoter in condition A and in B (r > 0) and log2 (deg_B + 1) / (deg_A + 1);
contact  contact_direct of the feature table (log2 change of the summed strength of the direct interactions);
direct   the existing direct features (promoter + direct);
graph    the existing direct + indirect (graph2) and direct + indirect + hop3 (graph3) features.

Per block: R^2 of base and base + block, delta R^2, partial F test. For graph versus window: the graph features
added on top of base + window and the window features added on top of base + graph. Among the differentially
expressed genes of the explanation cascade (<cmp>_gene_classes.tsv): genes with a concordant change (|log2| >= 1 in
the direction of the declared role of the mark and of the expression change) at the promoter, at step 1, 2, 3 and
in every window.
"""

import argparse
import json
import logging
import multiprocessing
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

from epigeneticgraphanalyzer.differential import ActivityOperators, gene_features, nested_models, sign_of
from epigeneticgraphanalyzer.states import role_signs

log = logging.getLogger('expression_baselines')

FIXED_WINDOWS = (('w50k', 0, 50000), ('w250k', 0, 250000), ('w1M', 0, 1000000), ('w2M', 0, 2000000),
                 ('ring50k_2M', 50000, 2000000))


# ---------------------------------------------------------------------------------------------------------------
# window feature
# ---------------------------------------------------------------------------------------------------------------

def window_signals(sources, elements, values, lo, hi, exponent=1.0, min_distance=1.0):
    """Distance-weighted mean signal of the elements around every source, without graph edges.

    sources   DataFrame with chrom, start, end (the promoter fragments);
    elements  DataFrame with chrom, start, end (fragments or bins);
    values    array (n_elements, n_signals) of element signals;
    lo, hi    distance range per source (arrays or scalars, bp, centre to centre): lo < d <= hi;
    Elements overlapping the source interval are left out; the weight of an element is max(d, min_distance) ^
    -exponent. Returns (S, n): S (n_sources, n_signals) = sum w x / sum w (0 without elements in range) and the
    number of elements n per source.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    n_sources = len(sources)
    lo = np.broadcast_to(np.asarray(lo, dtype=float), (n_sources,))
    hi = np.broadcast_to(np.asarray(hi, dtype=float), (n_sources,))
    out = np.zeros((n_sources, values.shape[1]))
    count = np.zeros(n_sources, dtype=np.int64)
    e_chrom = elements['chrom'].to_numpy()
    e_start = elements['start'].to_numpy(dtype=np.int64)
    e_end = elements['end'].to_numpy(dtype=np.int64)
    e_centre = (e_start + e_end) / 2.0
    s_chrom = sources['chrom'].to_numpy()
    s_start = sources['start'].to_numpy(dtype=np.int64)
    s_end = sources['end'].to_numpy(dtype=np.int64)
    s_centre = (s_start + s_end) / 2.0
    for chrom in pd.unique(s_chrom):
        on = np.nonzero(e_chrom == chrom)[0]
        if not len(on):
            continue
        order = on[np.argsort(e_centre[on], kind='stable')]
        centre = e_centre[order]
        for i in np.nonzero(s_chrom == chrom)[0]:
            if hi[i] <= lo[i]:
                continue
            left = np.searchsorted(centre, s_centre[i] - hi[i], side='left')
            right = np.searchsorted(centre, s_centre[i] + hi[i], side='right')
            idx = order[left:right]
            d = np.abs(e_centre[idx] - s_centre[i])
            keep = (d > lo[i]) & (d <= hi[i]) & ((e_end[idx] <= s_start[i]) | (e_start[idx] >= s_end[i]))
            if not keep.any():
                continue
            idx, d = idx[keep], d[keep]
            w = np.maximum(d, min_distance) ** -exponent
            out[i] = w @ values[idx] / w.sum()
            count[i] = len(idx)
    return out, count


def log2_change(a, b, pseudo_count=1.0):
    return np.log2((np.asarray(b, dtype=float) + pseudo_count) / (np.asarray(a, dtype=float) + pseudo_count))


# ---------------------------------------------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------------------------------------------

def ols(y, X):
    """(R^2, RSS, number of parameters with intercept) of an OLS fit (as differential._ols_r2)."""
    X = np.column_stack([np.ones(len(y)), X]) if X.size else np.ones((len(y), 1))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    residual = y - X @ beta
    rss = float(residual @ residual)
    tss = float(((y - y.mean()) ** 2).sum())
    return (1 - rss / tss if tss > 0 else 0.0), rss, X.shape[1]


def design(data, columns):
    return np.nan_to_num(data[columns].to_numpy(dtype=float)) if columns else np.zeros((len(data), 0))


def nested_test(data, base, added):
    """R^2 of base and base + added columns, delta R^2 and the partial F test (as differential.nested_models)."""
    y = data['rna_log2FC'].to_numpy(dtype=float)
    r2_0, rss_0, k_0 = ols(y, design(data, base))
    r2_1, rss_1, k_1 = ols(y, design(data, base + added))
    df1, df2 = k_1 - k_0, len(y) - k_1
    f = ((rss_0 - rss_1) / df1) / (rss_1 / df2) if rss_1 > 0 and df2 > 0 and df1 > 0 else np.nan
    return {'n_genes': len(y), 'r2_base': r2_0, 'r2_full': r2_1, 'delta_r2': r2_1 - r2_0, 'df1': df1, 'df2': df2,
            'f_statistic': f, 'f_pvalue': float(stats.f.sf(f, df1, df2)) if not np.isnan(f) else np.nan}


# ---------------------------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------------------------

def marks_of(features):
    return [c[:-len('_promoter')] for c in features.columns if c.endswith('_promoter')]


def declared_signs(path, marks):
    table = pd.read_csv(path, sep='\t')
    table = table[table['element_type'] == 'all'] if 'element_type' in table else table
    roles = table.drop_duplicates('mark').set_index('mark')['declared_role'].to_dict()
    return role_signs({m: roles.get(m, 'active') for m in marks})


def read_union_graph(path):
    import networkx as nx
    graph = nx.read_graphml(path)
    for _, data in graph.nodes(data=True):
        data['start'] = int(data['start'])
        data['end'] = int(data['end'])
    return graph


def node_table(graph, marks, names):
    rows = {'node': list(graph.nodes)}
    for key in ('chrom', 'start', 'end'):
        rows[key] = [graph.nodes[n][key] for n in rows['node']]
    for mark in marks:
        for name in names:
            rows['{}_{}'.format(mark, name)] = [float(graph.nodes[n].get('{}_{}'.format(mark, name), 0.0))
                                                for n in rows['node']]
    return pd.DataFrame(rows)


def _bigwig_bins(task):
    import pyBigWig
    path, chrom, bin_size = task
    handle = pyBigWig.open(path)
    try:
        size = handle.chroms().get(chrom)
        if not size:
            return None
        n = size // bin_size
        if n == 0:
            return np.zeros(0)
        stat = handle.stats(chrom, 0, n * bin_size, type='mean', nBins=n, exact=True)
        return np.array([float(s) if s is not None else 0.0 for s in stat])
    finally:
        handle.close()


def bigwig_bins(tracks, chroms, bin_size, threads):
    """Mean signal of consecutive bins of `bin_size` bp per mark and condition (mean over the files of a mark; 0
    without data, as annotation.bigwig_interval_signal). tracks: {condition: {mark: [paths]}}. Returns the bin
    table (chrom, start, end) and {(mark, condition): values}."""
    paths = sorted({p for by_mark in tracks.values() for ps in by_mark.values() for p in ps})
    tasks = [(p, c, bin_size) for p in paths for c in chroms]
    context = multiprocessing.get_context('fork')
    with context.Pool(max(1, min(threads, len(tasks)))) as pool:
        parts = pool.map(_bigwig_bins, tasks, chunksize=1)
    read = dict(zip(tasks, parts))
    lengths = {}
    for (p, c, _), v in read.items():
        if v is not None:
            lengths[c] = max(lengths.get(c, 0), len(v))
    bins = pd.concat([pd.DataFrame({'chrom': c, 'start': np.arange(n) * bin_size,
                                    'end': np.arange(1, n + 1) * bin_size}) for c, n in lengths.items()],
                     ignore_index=True)
    def padded(path, chrom, n):
        v = read[(path, chrom, bin_size)]
        v = np.zeros(0) if v is None else v
        return np.pad(v, (0, n - len(v)))

    signals = {}
    for condition, by_mark in tracks.items():
        for mark, ps in by_mark.items():
            signals[(mark, condition)] = np.mean([np.concatenate([padded(p, c, n) for c, n in lengths.items()])
                                                  for p in ps], axis=0)
    return bins, signals


# ---------------------------------------------------------------------------------------------------------------
# one comparison
# ---------------------------------------------------------------------------------------------------------------

def neighbour_ranges(graph, operators):
    """Largest distance (centre to centre, bp) of the nodes at 1 step and at 1 - 3 steps of every source (paths of
    either condition)."""
    centre = np.array([(graph.nodes[n]['start'] + graph.nodes[n]['end']) / 2.0 for n in operators.nodes])
    src_centre = centre[operators.source_index]
    ranges = {}
    reach = np.zeros(len(operators.sources))
    for k in sorted(operators.condition_a):
        pattern = ((abs(operators.condition_a[k]) > 0).astype(np.int8) +
                   (abs(operators.condition_b[k]) > 0).astype(np.int8)).tocoo()
        d = np.abs(centre[pattern.col] - src_centre[pattern.row])
        far = np.zeros(len(operators.sources))
        np.maximum.at(far, pattern.row, d)
        reach = np.maximum(reach, far)
        if k == 1:
            ranges['direct'] = far.copy()
        ranges['upto{}'.format(k)] = reach.copy()
    return ranges


def concordant(values, sign, direction, threshold=1.0):
    values = np.nan_to_num(np.asarray(values, dtype=float))
    return (np.sign(values) * sign == direction) & (np.abs(values) >= threshold)


def analyse(run, comparison, out, bin_size=5000, threads=8, pseudo_count=1.0):
    started = time.time()
    prefix = os.path.join(run, comparison)
    features = pd.read_csv(prefix + '_gene_features.tsv', sep='\t')
    marks = marks_of(features)
    signs = declared_signs(prefix + '_mark_roles.tsv', marks)
    support = json.load(open(prefix + '_support.json'))
    tracks = support['tracks']
    names = list(tracks)
    log.info('%s: %d genes, marks %s, conditions %s, declared signs %s', comparison, len(features), marks, names,
             signs)
    informative_mask = (features['n_direct'] > 0) & (features['n_indirect'] > 0)
    informative = features[informative_mask].reset_index(drop=True)

    # 1. reproduction of the nested models from the feature table
    reproduced = nested_models(informative, marks)
    stored = pd.read_csv(prefix + '_nested_models.tsv', sep='\t')
    check = reproduced.merge(stored, on='model', suffixes=('', '_stored'))
    check['r2_abs_difference'] = (check['r2'] - check['r2_stored']).abs()
    check['n_genes_equal'] = check['n_genes'] == check['n_genes_stored']
    log.info('nested models reproduced, largest |R2 difference| %.3g', check['r2_abs_difference'].max())

    # 2. union graph: node signals, recomputed features (promoter, direct, indirect, hop3; mark term)
    graph = read_union_graph(prefix + '_union.graphml')
    change = dict(zip(features['node'], features['rna_log2FC']))
    max_hops = 3 if any(c.endswith('_hop3') for c in features.columns) else 2
    operators = ActivityOperators(graph, list(change), max_hops=max_hops, step_feature='decomposed')
    signal = {m: ({n: float(graph.nodes[n].get('{}_{}'.format(m, names[0]), 0.0)) for n in graph.nodes},
                  {n: float(graph.nodes[n].get('{}_{}'.format(m, names[1]), 0.0)) for n in graph.nodes})
              for m in marks}
    recomputed = gene_features(graph, signal, change, operators=operators).set_index('node')
    stored_features = features.set_index('node').loc[recomputed.index]
    feature_check = []
    for level in ['promoter'] + operators.levels:
        for mark in marks:
            column = '{}_{}'.format(mark, level)
            diff = np.nanmax(np.abs(recomputed[column].to_numpy() - stored_features[column].to_numpy()))
            feature_check.append({'column': column, 'max_abs_difference': float(diff)})
    feature_check = pd.DataFrame(feature_check)
    log.info('features recomputed from the union graph, largest |difference| %.3g',
             feature_check['max_abs_difference'].max())

    # 3. baseline features of every gene of the feature table
    sources = pd.DataFrame({'node': operators.sources,
                            'chrom': [graph.nodes[n]['chrom'] for n in operators.sources],
                            'start': [graph.nodes[n]['start'] for n in operators.sources],
                            'end': [graph.nodes[n]['end'] for n in operators.sources]})
    ranges = neighbour_ranges(graph, operators)
    windows = [('direct', 0, ranges['direct']), ('upto{}'.format(max_hops), 0, ranges['upto{}'.format(max_hops)])]
    windows += list(FIXED_WINDOWS)
    baseline = pd.DataFrame({'node': operators.sources, 'range_direct_bp': ranges['direct'],
                             'range_upto{}_bp'.format(max_hops): ranges['upto{}'.format(max_hops)]})
    nodes = node_table(graph, marks, names)
    node_values = nodes[['{}_{}'.format(m, c) for c in names for m in marks]].to_numpy()
    element_sets = [('nodes', nodes, node_values)]
    chroms = sorted(set(sources['chrom']))
    t = time.time()
    bins, bin_signals = bigwig_bins(tracks, chroms, bin_size, threads)
    bin_values = np.column_stack([bin_signals[(m, c)] for c in names for m in marks])
    element_sets.append(('bigwig', bins, bin_values))
    log.info('bigWig bins of %d bp: %d bins, %.0f s', bin_size, len(bins), time.time() - t)
    n_marks = len(marks)
    for source_name, elements, values in element_sets:
        for window, lo, hi in windows:
            S, count = window_signals(sources, elements, values, lo, hi)
            block = '{}_{}'.format(window, source_name)
            baseline['n_' + block] = count
            for j, mark in enumerate(marks):
                baseline['{}_window_{}'.format(mark, block)] = log2_change(S[:, j], S[:, n_marks + j], pseudo_count)
    degree_a = np.zeros(len(operators.sources))
    degree_b = np.zeros(len(operators.sources))
    index = {n: i for i, n in enumerate(operators.sources)}
    for u, v, data in graph.edges(data=True):
        for node in (u, v):
            if node in index:
                degree_a[index[node]] += float(data['r_a']) > 0
                degree_b[index[node]] += float(data['r_b']) > 0
    baseline['degree_a'] = degree_a
    baseline['degree_b'] = degree_b
    baseline['degree_log2_change'] = log2_change(degree_a, degree_b)
    table = features.merge(baseline, on='node', how='left')
    table.drop(columns=[c for c in table if c.endswith('_mirror')]).to_csv(
        os.path.join(out, comparison + '_baseline_features.tsv'), sep='\t', index=False)
    data = table[informative_mask.to_numpy()].reset_index(drop=True)

    # 4. models on the informative genes, base = promoter
    promoter = [m + '_promoter' for m in marks]
    levels = operators.levels
    blocks = {'degree': ['degree_a', 'degree_b', 'degree_log2_change'], 'contact': ['contact_direct'],
              'direct': [m + '_direct' for m in marks],
              'graph2': [m + '_' + level for level in levels[:2] for m in marks],
              'graph{}'.format(max_hops): [m + '_' + level for level in levels for m in marks]}
    window_blocks = {}
    for source_name, _, _ in element_sets:
        for window, _, _ in windows:
            block = '{}_{}'.format(window, source_name)
            window_blocks['window_' + block] = ['{}_window_{}'.format(m, block) for m in marks]
    rows = []
    for name, columns in list(blocks.items()) + list(window_blocks.items()):
        rows.append(dict(comparison=comparison, test='block_over_promoter', block=name, added=name,
                         **nested_test(data, promoter, columns)))
    graph_blocks = ['direct', 'graph2', 'graph{}'.format(max_hops)]
    for w_name, w_columns in window_blocks.items():
        for g_name in graph_blocks:
            g_columns = blocks[g_name]
            rows.append(dict(comparison=comparison, test='graph_over_window', block=w_name, added=g_name,
                             **nested_test(data, promoter + w_columns, g_columns)))
            rows.append(dict(comparison=comparison, test='window_over_graph', block=g_name, added=w_name,
                             **nested_test(data, promoter + g_columns, w_columns)))
    models = pd.DataFrame(rows)
    models.to_csv(os.path.join(out, comparison + '_baseline_models.tsv'), sep='\t', index=False)

    # 5. concordant genes among the differential genes of the cascade
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t', usecols=['node', 'rna_log2FC', 'explained_at'])
    de = classes.merge(table.drop(columns=['rna_log2FC']), on='node', how='left')
    direction = np.sign(de['rna_log2FC'].to_numpy(dtype=float))

    def any_concordant(level_of_mark):
        hit = np.zeros(len(de), dtype=bool)
        for mark in marks:
            column, level = level_of_mark(mark)
            hit |= concordant(de[column], sign_of(signs, mark, level), direction)
        return hit

    hits = {'promoter': any_concordant(lambda m: (m + '_promoter', 'promoter'))}
    for k, level in enumerate(levels, start=1):
        hits['step{}'.format(k)] = any_concordant(lambda m, level=level: (m + '_' + level, level))
    hits['step1_or_step2'] = hits['step1'] | hits['step2']
    for name in window_blocks:
        block = name[len('window_'):]
        hits[name] = any_concordant(lambda m, block=block: ('{}_window_{}'.format(m, block), 'direct'))
    reaching = (de['explained_at'] != 'gene').to_numpy()
    without_promoter = ~hits['promoter']
    rows = []
    for name, hit in hits.items():
        row = {'comparison': comparison, 'level': name, 'differential_genes': len(de), 'concordant': int(hit.sum()),
               'concordant_without_promoter_marks': int((hit & without_promoter).sum()),
               'genes_reaching_step1': int(reaching.sum()), 'concordant_reaching_step1': int((hit & reaching).sum())}
        if name.startswith('window_'):
            for ref in ('step1', 'step2', 'step1_or_step2'):
                row['both_with_' + ref] = int((hit & hits[ref]).sum())
                row[ref + '_only'] = int((hits[ref] & ~hit).sum())
                row['window_only_vs_' + ref] = int((hit & ~hits[ref]).sum())
        rows.append(row)
    counts = pd.DataFrame(rows)
    counts.to_csv(os.path.join(out, comparison + '_baseline_concordance.tsv'), sep='\t', index=False)
    check.assign(comparison=comparison).to_csv(os.path.join(out, comparison + '_reproduction.tsv'), sep='\t',
                                               index=False)
    feature_check.assign(comparison=comparison).to_csv(os.path.join(out, comparison + '_feature_check.tsv'),
                                                       sep='\t', index=False)
    log.info('%s done: %.0f s', comparison, time.time() - started)
    return check.assign(comparison=comparison), feature_check.assign(comparison=comparison), models, counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('run', help='directory of the egaDifferential outputs (<cmp>_*)')
    parser.add_argument('out', help='output directory')
    parser.add_argument('--comparisons', nargs='+', required=True)
    parser.add_argument('--binSize', type=int, default=5000)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')
    os.makedirs(args.out, exist_ok=True)
    parts = [analyse(args.run, c, args.out, bin_size=args.binSize, threads=min(args.threads, 8))
             for c in args.comparisons]
    for i, name in enumerate(('reproduction', 'feature_check', 'models', 'concordance')):
        pd.concat([p[i] for p in parts], ignore_index=True).to_csv(
            os.path.join(args.out, 'baselines_{}.tsv'.format(name)), sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
