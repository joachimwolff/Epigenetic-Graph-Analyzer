"""Benchmark graph features of a pcHi-C graph against CRISPRi-tested enhancer-gene pairs.

Every tested element-gene pair is mapped onto the graph: the element to the fragments it overlaps, the gene to the
fragment holding its TSS. Pairs whose TSS fragment is not a bait with interactions, and pairs whose element lies on
the TSS fragment, cannot be measured by promoter capture and are reported but not scored. Scores per pair:

  distance        negative TSS distance
  direct          largest CHiCAGO score of an edge element - TSS fragment (0 without edge)
  two_hop         largest min(score) over paths element - mediator - TSS fragment (0 without such path)
  hops            1 / length of the shortest path (0 if longer than six edges or disconnected)
  community       element and TSS fragment in one community
  ppr             personalised PageRank of the element fragment with restart at the TSS fragment
  abc, re2g       ABC score and ENCODE-rE2G score of the overlapping candidate element (ENCODE predictions)
  powerlaw        3D contact estimate of ENCODE-rE2G (power law of distance)

AUPRC and AUROC are computed for all scorable pairs and for those without a direct edge; logistic models test
whether graph features add to log distance and element activity (degree in the graph).

Usage: python crispr_benchmark.py <graph.graphml> <benchmark.tsv.gz> <predictions.bed.gz> <output prefix>
       [cell type=K562]
"""

import sys

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

from epigeneticgraphanalyzer.graph import read_graph

MAX_DISTANCE = 2000000
MAX_HOPS = 6


class NodeIndex:
    def __init__(self, graph):
        self.nodes = list(graph.nodes)
        self.position = {n: i for i, n in enumerate(self.nodes)}
        table = pd.DataFrame([(n, d['chrom'], int(d['start']), int(d['end'])) for n, d in graph.nodes(data=True)],
                             columns=['node', 'chrom', 'start', 'end']).sort_values(['chrom', 'start'])
        self.by_chrom = {c: g.reset_index(drop=True) for c, g in table.groupby('chrom')}

    def overlapping(self, chrom, start, end):
        table = self.by_chrom.get(chrom)
        if table is None:
            return []
        hi = np.searchsorted(table['start'].to_numpy(), end, 'left')
        lo = max(0, np.searchsorted(table['start'].to_numpy(), start, 'right') - 50)
        sub = table.iloc[lo:hi]
        return list(sub.loc[(sub['end'] > start) & (sub['start'] < end), 'node'])


def personalised_pagerank(adjacency, source, restart=0.15, iterations=50):
    vector = np.zeros(adjacency.shape[0])
    vector[source] = 1.0
    start = vector.copy()
    for _ in range(iterations):
        vector = (1 - restart) * (adjacency @ vector) + restart * start
    return vector


def load_predictions(path, genes):
    columns = ['#chr', 'start', 'end', 'TargetGene', 'ABC.Score.Feature', 'Score', '3DContact.Feature']
    parts = []
    for chunk in pd.read_csv(path, sep='\t', usecols=columns, chunksize=2000000):
        parts.append(chunk[chunk['TargetGene'].isin(genes)])
    predictions = pd.concat(parts, ignore_index=True)
    return predictions.rename(columns={'#chr': 'chrom', 'ABC.Score.Feature': 'abc', 'Score': 're2g',
                                       '3DContact.Feature': 'powerlaw'})


def prediction_scores(pairs, predictions):
    out = np.zeros((len(pairs), 3))
    grouped = {k: g for k, g in predictions.groupby('TargetGene')}
    for i, row in enumerate(pairs.itertuples()):
        g = grouped.get(row.measuredGeneSymbol)
        if g is None:
            continue
        hit = g[(g['chrom'] == row.chrom) & (g['end'] > row.chromStart) & (g['start'] < row.chromEnd)]
        if len(hit):
            out[i] = hit[['abc', 're2g', 'powerlaw']].max().to_numpy()
    return out


def evaluate(table, features, label):
    rows = []
    y = table['Regulated'].to_numpy().astype(int)
    for feature in features:
        x = table[feature].to_numpy().astype(float)
        if y.sum() == 0 or y.sum() == len(y):
            continue
        rows.append(dict(set=label, feature=feature, n=len(y), positives=int(y.sum()),
                         auprc=average_precision_score(y, x), auroc=roc_auc_score(y, x),
                         nonzero=float(np.mean(x != 0)),
                         precision_nonzero=float(y[x != 0].mean()) if np.any(x != 0) else np.nan,
                         recall_nonzero=float(y[x != 0].sum() / y.sum())))
    return rows


def added_value(table, base, feature, label):
    """Likelihood-ratio test of `feature` added to a logistic model with the `base` covariates."""
    y = table['Regulated'].astype(int).to_numpy()
    x0 = sm.add_constant(table[base].astype(float).to_numpy())
    x1 = sm.add_constant(table[base + [feature]].astype(float).to_numpy())
    try:
        m0 = sm.Logit(y, x0).fit(disp=0)
        m1 = sm.Logit(y, x1).fit(disp=0)
    except Exception as error:  # perfect separation or singular design
        return dict(set=label, feature=feature, error=str(error))
    lr = 2 * (m1.llf - m0.llf)
    return dict(set=label, feature=feature, base='+'.join(base), coefficient=m1.params[-1], z=m1.tvalues[-1],
                lr_statistic=lr, lr_p=stats.chi2.sf(lr, 1), n=len(y), positives=int(y.sum()))


def main():
    graph_path, benchmark_path, predictions_path, prefix = sys.argv[1:5]
    cell_type = sys.argv[5] if len(sys.argv) > 5 else 'K562'
    graph = read_graph(graph_path)
    index = NodeIndex(graph)
    bench = pd.read_csv(benchmark_path, sep='\t')
    valid = bench['ValidConnection'].astype(str).str.upper() == 'TRUE'
    bench = bench[(bench['CellType'] == cell_type) & valid & bench['Regulated'].notna()
                  & (bench['distanceToTSS'].abs() <= MAX_DISTANCE)].copy()
    bench['Regulated'] = bench['Regulated'].astype(str).str.upper() == 'TRUE'
    bench = bench.reset_index(drop=True)

    # weighted adjacency (column-normalised) for personalised PageRank
    n = len(index.nodes)
    rows, cols, vals = [], [], []
    for u, v, d in graph.edges(data=True):
        w = float(d.get('score', d.get('weight', 1.0)))
        i, j = index.position[u], index.position[v]
        rows += [i, j]
        cols += [j, i]
        vals += [w, w]
    adjacency = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    degree = np.asarray(adjacency.sum(axis=0)).ravel()
    adjacency = adjacency @ sp.diags(np.where(degree > 0, 1.0 / np.maximum(degree, 1e-12), 0.0))

    records = []
    ppr_cache = {}
    for row in bench.itertuples():
        tss_nodes = index.overlapping(row.chrTSS, int(row.startTSS), int(row.endTSS))
        element_nodes = index.overlapping(row.chrom, int(row.chromStart), int(row.chromEnd))
        tss = tss_nodes[0] if tss_nodes else None
        record = dict(tss_in_graph=tss is not None and bool(graph.nodes[tss].get('bait')),
                      element_in_graph=bool(element_nodes), same_fragment=tss is not None and tss in element_nodes,
                      direct=0.0, two_hop=0.0, n_mediators=0, hops=0.0, community=0, ppr=0.0, element_degree=0)
        if element_nodes:
            record['element_degree'] = max(graph.degree[e] for e in element_nodes)
        if record['tss_in_graph'] and element_nodes and not record['same_fragment']:
            if tss not in ppr_cache:
                ppr_cache[tss] = (personalised_pagerank(adjacency, index.position[tss]),
                                  nx.single_source_shortest_path_length(graph, tss, cutoff=MAX_HOPS))
            ppr, lengths = ppr_cache[tss]
            tss_neighbours = set(graph.neighbors(tss))
            community = graph.nodes[tss].get('community')
            for e in element_nodes:
                if graph.has_edge(e, tss):
                    record['direct'] = max(record['direct'], float(graph.edges[e, tss].get('score', 1.0)))
                mediators = tss_neighbours & set(graph.neighbors(e))
                record['n_mediators'] = max(record['n_mediators'], len(mediators))
                for m in mediators:
                    path = min(float(graph.edges[e, m].get('score', 1.0)), float(graph.edges[m, tss].get('score', 1.0)))
                    record['two_hop'] = max(record['two_hop'], path)
                if e in lengths:
                    record['hops'] = max(record['hops'], 1.0 / lengths[e])
                if community is not None and graph.nodes[e].get('community') == community and int(community) >= 0:
                    record['community'] = 1
                record['ppr'] = max(record['ppr'], float(ppr[index.position[e]]))
        records.append(record)
    table = pd.concat([bench, pd.DataFrame(records)], axis=1)
    table['distance'] = -table['distanceToTSS'].abs()
    table['log_distance'] = np.log10(table['distanceToTSS'].abs() + 1000)
    predictions = load_predictions(predictions_path, set(table['measuredGeneSymbol']))
    table[['abc', 're2g', 'powerlaw']] = prediction_scores(table, predictions)
    table.to_csv(prefix + '_pairs.tsv.gz', sep='\t', index=False)

    features = ['distance', 'powerlaw', 'direct', 'two_hop', 'hops', 'community', 'ppr', 'abc', 're2g']
    scorable = table[table['tss_in_graph'] & ~table['same_fragment']]
    no_direct = scorable[scorable['direct'] == 0]
    coverage = dict(pairs=len(table), positives=int(table['Regulated'].sum()),
                    tss_bait_in_graph=int(table['tss_in_graph'].sum()),
                    scorable=len(scorable), scorable_positives=int(scorable['Regulated'].sum()),
                    same_fragment=int(table['same_fragment'].sum()))
    print(coverage, flush=True)
    rows = (evaluate(table, ['distance', 'powerlaw', 'abc', 're2g'], 'all_pairs') +
            evaluate(scorable, features, 'scorable') + evaluate(no_direct, features, 'scorable_without_direct_edge'))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(prefix + '_metrics.tsv', sep='\t', index=False)
    print(metrics.to_string(), flush=True)

    tests = []
    base = ['log_distance', 'element_degree']
    for feature in ['direct', 'two_hop', 'hops', 'community', 'ppr']:
        tests.append(added_value(scorable, base, feature, 'scorable'))
    for feature in ['two_hop', 'hops', 'community', 'ppr']:
        tests.append(added_value(no_direct, base, feature, 'scorable_without_direct_edge'))
        tests.append(added_value(scorable, base + ['direct'], feature, 'scorable_given_direct'))
        tests.append(added_value(scorable[scorable['abc'] > 0], ['abc', 'log_distance'], feature, 'scorable_given_abc'))
    tests = pd.DataFrame(tests)
    tests.to_csv(prefix + '_added_value.tsv', sep='\t', index=False)
    print(tests.to_string(), flush=True)
    pd.Series(coverage).to_csv(prefix + '_coverage.tsv', sep='\t', header=False)


if __name__ == '__main__':
    main()
