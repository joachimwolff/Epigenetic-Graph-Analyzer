"""Validation of the egaDesign ranking against K562 CRISPRi element-gene pairs (single-condition mode).

No differential run exists for K562 alone, so the nodes of a gene's neighbourhood (up to three steps on the K562 DMSO
pcHi-C graph, edge strength = CHiCAGO score / 5) are ranked by their contribution to the gene's contact-weighted
DNase activity (design.single_condition_contributions: rank by the contribution, path weight x DNase, whose share of
the gene's summed activity is the design score; in addition silence_drop, the drop of the log2 activity of the node's
step when its DNase signal is set to 0). Every CRISPRi pair of the ENCODE-rE2G K562 training set (the pairs and
filters of crispr_benchmark.py: valid connection, measured, |distance| <= 2 Mb) is mapped onto the graph: the gene to
the bait fragment holding its TSS, the element to the fragments it overlaps (best rank among them); pairs whose element lies
on the TSS fragment or whose TSS fragment is not a bait are not scored. Outputs (prefix):

  _pairs.tsv.gz        every pair with its design rank, share, silence_drop, steps and the number of ranked nodes
  _bins.tsv            fraction of pairs with a significant CRISPRi effect per design rank bin (Wilson 95% CI),
                       overall and within distance strata
  _metrics.tsv         AUPRC / AUROC of the design score against distance, element DNase and a direct edge
  _added_value.tsv     likelihood-ratio tests of the design score added to log distance + element DNase (logistic)
  _within_gene.tsv     per gene with positive and negative reached elements: AUROC of the design score and of
                       distance within the gene
  _summary.json

Usage: python design_crispr_validation.py <graph.graphml> <benchmark.tsv.gz> <DNase.bigWig> <output prefix>
       [threads=8] [max hops=3] [reference score=5]
"""

import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

from epigeneticgraphanalyzer.design import GraphArrays, single_condition_contributions
from epigeneticgraphanalyzer.differential import _interval_signals
from epigeneticgraphanalyzer.graph import read_graph

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crispr_benchmark import MAX_DISTANCE, NodeIndex  # noqa: E402

RANK_BINS = [(1, 1, '1'), (2, 3, '2-3'), (4, 10, '4-10'), (11, 50, '11-50'), (51, 10 ** 9, '>50')]
DISTANCE_BINS = [(0, 20000, '<20 kb'), (20000, 100000, '20-100 kb'), (100000, MAX_DISTANCE + 1, '>=100 kb')]


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return p, max(0.0, centre - half), min(1.0, centre + half)


def rank_bin(rank):
    if pd.isna(rank):
        return 'not reached'
    for low, high, label in RANK_BINS:
        if low <= rank <= high:
            return label
    return 'not reached'


def bins_table(pairs, label):
    rows = []
    order = [b[2] for b in RANK_BINS] + ['not reached']
    for name in order:
        part = pairs[pairs['rank_bin'] == name]
        k, n = int(part['Regulated'].sum()), len(part)
        p, low, high = wilson(k, n)
        rows.append(dict(set=label, rank_bin=name, pairs=n, positives=k, fraction=p, ci_low=low, ci_high=high,
                         median_distance=float(part['abs_distance'].median()) if n else np.nan))
    return rows


def metrics(table, features, label):
    rows = []
    y = table['Regulated'].astype(int).to_numpy()
    if y.sum() == 0 or y.sum() == len(y):
        return rows
    for feature in features:
        x = table[feature].astype(float).to_numpy()
        rows.append(dict(set=label, feature=feature, n=len(y), positives=int(y.sum()),
                         auprc=average_precision_score(y, x), auroc=roc_auc_score(y, x),
                         baseline_precision=float(y.mean())))
    return rows


def added_value(table, base, feature, label):
    y = table['Regulated'].astype(int).to_numpy()
    try:
        m0 = sm.Logit(y, sm.add_constant(table[base].astype(float).to_numpy())).fit(disp=0)
        m1 = sm.Logit(y, sm.add_constant(table[base + [feature]].astype(float).to_numpy())).fit(disp=0)
    except Exception as error:  # separation, singular design
        return dict(set=label, feature=feature, base='+'.join(base), error=str(error))
    lr = 2 * (m1.llf - m0.llf)
    return dict(set=label, feature=feature, base='+'.join(base), n=len(y), positives=int(y.sum()),
                coefficient=float(m1.params[-1]), z=float(m1.tvalues[-1]), lr_statistic=float(lr),
                lr_p=float(stats.chi2.sf(lr, 1)))


def main(argv):
    graph_path, bench_path, dnase_path, prefix = argv[1:5]
    threads = int(argv[5]) if len(argv) > 5 else 8
    max_hops = int(argv[6]) if len(argv) > 6 else 3
    reference = float(argv[7]) if len(argv) > 7 else 5.0
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    graph = read_graph(graph_path)
    nodes = list(graph.nodes)
    signal = _interval_signals([dnase_path], [graph.nodes[n]['chrom'] for n in nodes],
                               [int(graph.nodes[n]['start']) for n in nodes],
                               [int(graph.nodes[n]['end']) for n in nodes], threads)[dnase_path]
    for n, v in zip(nodes, signal):
        graph.nodes[n]['DNase'] = float(v)
    arrays = GraphArrays(graph, {'DNase': ('DNase', 'DNase')}, weights=('score', 'score'), reference_score=reference)
    index = NodeIndex(graph)

    bench = pd.read_csv(bench_path, sep='\t')
    valid = bench['ValidConnection'].astype(str).str.upper() == 'TRUE'
    bench = bench[(bench['CellType'] == 'K562') & valid & bench['Regulated'].notna()
                  & (bench['distanceToTSS'].abs() <= MAX_DISTANCE)].copy()
    bench['Regulated'] = bench['Regulated'].astype(str).str.upper() == 'TRUE'
    bench = bench.reset_index(drop=True)
    tss_node, element_nodes = [], []
    for row in bench.itertuples():
        t = index.overlapping(row.chrTSS, int(row.startTSS), int(row.endTSS))
        tss_node.append(t[0] if t and graph.nodes[t[0]].get('bait') else None)
        element_nodes.append(index.overlapping(row.chrom, int(row.chromStart), int(row.chromEnd)))
    bench['tss_node'] = tss_node
    bench['element_nodes'] = [','.join(e) for e in element_nodes]
    bench['same_fragment'] = [t is not None and t in e for t, e in zip(tss_node, element_nodes)]
    sources = sorted({t for t in tss_node if t is not None})
    ranking = single_condition_contributions(arrays, sources, max_hops, 'DNase')
    ranking.to_csv(prefix + '_single_condition.tsv.gz', sep='\t', index=False)
    lookup = {(g, n): (r, lr, sh, d, h) for g, n, r, lr, sh, d, h in
              zip(ranking['gene_node'], ranking['node'], ranking['rank'], ranking['level_rank'], ranking['share'],
                  ranking['silence_drop'], ranking['hops'])}
    n_ranked = ranking.groupby('gene_node')['node'].size().to_dict()
    best = {key: [] for key in ('design_rank', 'level_rank', 'share', 'silence_drop', 'steps')}
    direct = []
    for t, elements in zip(tss_node, element_nodes):
        hits = [lookup[(t, e)] for e in elements if (t, e) in lookup] if t is not None else []
        if hits:
            best['design_rank'].append(min(h[0] for h in hits))
            best['level_rank'].append(min(h[1] for h in hits))
            best['share'].append(max(h[2] for h in hits))
            best['silence_drop'].append(max(h[3] for h in hits))
            best['steps'].append(min(h[4] for h in hits))
        else:
            best['design_rank'].append(np.nan)
            best['level_rank'].append(np.nan)
            best['share'].append(0.0)
            best['silence_drop'].append(0.0)
            best['steps'].append(np.nan)
        direct.append(int(t is not None and any(graph.has_edge(t, e) for e in elements)))
    for key, values in best.items():
        bench[key] = values
    bench['ranked_nodes'] = [n_ranked.get(t, 0) if t is not None else 0 for t in tss_node]
    bench['direct_edge'] = direct
    bench['reached'] = bench['design_rank'].notna()
    bench['rank_bin'] = [rank_bin(r) for r in bench['design_rank']]
    bench['abs_distance'] = bench['distanceToTSS'].abs()
    bench['neg_distance'] = -bench['abs_distance']
    bench['log_distance'] = np.log10(bench['abs_distance'] + 1000)
    bench['log_dnase'] = np.log10(bench['DHS.RPM'].astype(float).fillna(0) + 0.1)
    bench['neg_rank'] = -bench['design_rank'].fillna(10 ** 6)
    bench['neg_level_rank'] = -bench['level_rank'].fillna(10 ** 6)
    bench['top3'] = (bench['design_rank'] <= 3).astype(int)
    bench['top10'] = (bench['design_rank'] <= 10).astype(int)
    keep = ['chrom', 'chromStart', 'chromEnd', 'name', 'measuredGeneSymbol', 'EffectSize', 'Regulated', 'Significant',
            'distanceToTSS', 'DHS.RPM', 'tss_node', 'element_nodes', 'same_fragment', 'direct_edge', 'reached',
            'design_rank', 'rank_bin', 'share', 'level_rank', 'silence_drop', 'steps', 'ranked_nodes']
    bench[keep].to_csv(prefix + '_pairs.tsv.gz', sep='\t', index=False)

    scorable = bench[bench['tss_node'].notna() & ~bench['same_fragment']].copy()
    reached = scorable[scorable['reached']]
    rows = bins_table(scorable, 'scorable')
    for low, high, label in DISTANCE_BINS:
        part = scorable[(scorable['abs_distance'] >= low) & (scorable['abs_distance'] < high)]
        rows += bins_table(part, 'scorable, distance ' + label)
    pd.DataFrame(rows).to_csv(prefix + '_bins.tsv', sep='\t', index=False)

    features = ['share', 'neg_rank', 'silence_drop', 'neg_level_rank', 'neg_distance', 'log_dnase', 'direct_edge']
    metric_rows = metrics(scorable, features, 'scorable') + metrics(reached, features, 'reached')
    pd.DataFrame(metric_rows).to_csv(prefix + '_metrics.tsv', sep='\t', index=False)

    base = ['log_distance', 'log_dnase']
    tests = []
    for label, table in (('scorable', scorable), ('reached', reached)):
        for feature in ('share', 'silence_drop', 'top3', 'top10', 'direct_edge'):
            tests.append(added_value(table, base, feature, label))
    for label, table in (('scorable_given_direct', scorable), ('reached_given_direct', reached)):
        for feature in ('share', 'top3', 'top10'):
            tests.append(added_value(table, base + ['direct_edge'], feature, label))
    pd.DataFrame(tests).to_csv(prefix + '_added_value.tsv', sep='\t', index=False)

    within = []
    for gene, part in reached.groupby('measuredGeneSymbol'):
        y = part['Regulated'].astype(int)
        if y.sum() == 0 or y.sum() == len(y):
            continue
        within.append(dict(gene=gene, pairs=len(part), positives=int(y.sum()),
                           auroc_design=roc_auc_score(y, part['share']),
                           auroc_level_drop=roc_auc_score(y, part['silence_drop']),
                           auroc_distance=roc_auc_score(y, part['neg_distance']),
                           auroc_dnase=roc_auc_score(y, part['log_dnase'])))
    within = pd.DataFrame(within)
    within.to_csv(prefix + '_within_gene.tsv', sep='\t', index=False)

    summary = dict(pairs=len(bench), positives=int(bench['Regulated'].sum()), scorable=len(scorable),
                   scorable_positives=int(scorable['Regulated'].sum()), reached=len(reached),
                   reached_positives=int(reached['Regulated'].sum()), genes_ranked=len(sources), max_hops=max_hops,
                   reference_score=reference)
    if len(within):
        summary.update(within_genes=len(within),
                       within_mean_auroc_design=float(within['auroc_design'].mean()),
                       within_mean_auroc_distance=float(within['auroc_distance'].mean()),
                       within_mean_auroc_dnase=float(within['auroc_dnase'].mean()),
                       within_wilcoxon_design_vs_distance_p=float(stats.wilcoxon(
                           within['auroc_design'], within['auroc_distance']).pvalue)
                       if len(within) > 5 else np.nan,
                       within_wilcoxon_design_vs_half_p=float(stats.wilcoxon(within['auroc_design'] - 0.5).pvalue)
                       if len(within) > 5 else np.nan)
    for name in ('1', '2-3', '4-10'):
        part = scorable[scorable['rank_bin'] == name]
        summary['fraction_positive_rank_' + name] = float(part['Regulated'].mean()) if len(part) else np.nan
    for name in ('>50', 'not reached'):
        part = scorable[scorable['rank_bin'] == name]
        summary['fraction_positive_' + name.replace(' ', '_')] = float(part['Regulated'].mean()) if len(part) \
            else np.nan
    summary['fraction_positive_top3'] = float(scorable.loc[scorable['top3'] == 1, 'Regulated'].mean())
    summary['fraction_positive_rest'] = float(scorable.loc[scorable['top3'] == 0, 'Regulated'].mean())
    table = [[int(((scorable['top3'] == 1) & scorable['Regulated']).sum()),
              int(((scorable['top3'] == 1) & ~scorable['Regulated']).sum())],
             [int(((scorable['top3'] == 0) & scorable['Regulated']).sum()),
              int(((scorable['top3'] == 0) & ~scorable['Regulated']).sum())]]
    odds, p = stats.fisher_exact(table)
    summary.update(top3_odds_ratio=float(odds), top3_fisher_p=float(p))
    with open(prefix + '_summary.json', 'w') as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(rows).to_string())
    print(pd.DataFrame(metric_rows).to_string())
    print(pd.DataFrame(tests).to_string())


if __name__ == '__main__':
    main(sys.argv)
