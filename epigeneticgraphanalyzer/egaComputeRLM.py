import argparse
import os

import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig


def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute per-node and neighborhood RLM values")
    parser.add_argument('-g', '--graph', type=str, required=True, help='Path to the input graph file (GraphML)')
    parser.add_argument(
        '-n', '--conditions',
        type=str,
        nargs='+',
        required=True,
        help='Condition names. Must match --bigwigFiles length (one entry per condition).',
    )
    parser.add_argument(
        '-b', '--bigwigFiles',
        type=str,
        nargs='+',
        required=True,
        help='Per-condition comma-separated bigWig replicates, e.g. cond1_rep1.bw,cond1_rep2.bw',
    )
    parser.add_argument(
        '-d', '--deseqCounts',
        type=str,
        required=False,
        default=None,
        help='Optional DESeq2 table with gene symbol and log2 fold-change columns',
    )
    parser.add_argument('-s', '--steps', type=int, default=2, help='Neighborhood radius for RLM aggregation (default: 2)')
    parser.add_argument('-o', '--outputFolder', type=str, required=True, help='Path to the output folder')
    return parser.parse_args()


def _parse_condition_replicates(conditions, bigwig_files):
    if len(conditions) != len(bigwig_files):
        raise ValueError('The number of --conditions entries must equal the number of --bigwigFiles entries.')

    parsed = {}
    for condition, files in zip(conditions, bigwig_files):
        reps = [x.strip() for x in files.split(',') if x.strip()]
        if len(reps) == 0:
            raise ValueError(f'No bigWig files provided for condition {condition}.')
        parsed[condition] = reps
    return parsed


def _node_interval(node_attrs):
    position = node_attrs.get('position')
    if not position or ':' not in position or '-' not in position:
        return None
    chromosome, value = position.split(':', 1)
    start, end = value.split('-', 1)
    return chromosome, int(start), int(end)


def _bw_region_mean(bw_handle, chromosome, start, end):
    try:
        stat = bw_handle.stats(chromosome, start, end, type='mean')[0]
        if stat is None:
            return 0.0
        return float(stat)
    except RuntimeError:
        return 0.0


def annotate_graph_with_bigwig(graph, condition_to_replicates):
    for condition, replicates in condition_to_replicates.items():
        handles = [pyBigWig.open(rep) for rep in replicates]
        try:
            for node in graph.nodes:
                interval = _node_interval(graph.nodes[node])
                if interval is None:
                    graph.nodes[node][f'{condition}_mean'] = 0.0
                    continue
                chrom, start, end = interval
                values = [_bw_region_mean(handle, chrom, start, end) for handle in handles]
                graph.nodes[node][f'{condition}_mean'] = float(np.mean(values))
        finally:
            for handle in handles:
                handle.close()


def _infer_deseq_columns(df):
    candidates_gene = ['Gene', 'gene', 'symbol', 'GeneName']
    candidates_l2fc = ['log2FoldChange', 'log2FC', 'LFC']

    gene_col = next((c for c in candidates_gene if c in df.columns), None)
    l2fc_col = next((c for c in candidates_l2fc if c in df.columns), None)
    return gene_col, l2fc_col


def annotate_graph_with_rna(graph, deseq_path):
    if deseq_path is None:
        return

    df = pd.read_csv(deseq_path, sep=None, engine='python')
    gene_col, l2fc_col = _infer_deseq_columns(df)
    if gene_col is None or l2fc_col is None:
        return

    mapping = dict(zip(df[gene_col].astype(str), df[l2fc_col].astype(float)))
    for node in graph.nodes:
        name = graph.nodes[node].get('name')
        if not name:
            continue
        symbol = str(name).split('_')[0]
        if symbol in mapping:
            graph.nodes[node]['RNA_log2FC'] = mapping[symbol]


def create_rlm_table(graph, steps, condition_names):
    rows = []
    for node in graph.nodes:
        path_lengths = nx.single_source_shortest_path_length(graph, node, cutoff=steps)
        neighborhood = [n for n, d in path_lengths.items() if d > 0]

        row = {
            'node_id': str(node),
            'name': graph.nodes[node].get('name', ''),
            'position': graph.nodes[node].get('position', ''),
            'neighborhood_size': len(neighborhood),
        }

        for condition in condition_names:
            node_value = float(graph.nodes[node].get(f'{condition}_mean', 0.0))
            neighborhood_values = [float(graph.nodes[n].get(f'{condition}_mean', 0.0)) for n in neighborhood]
            neighborhood_mean = float(np.mean(neighborhood_values)) if neighborhood_values else 0.0

            row[f'{condition}_node_mean'] = node_value
            row[f'{condition}_neighborhood_mean'] = neighborhood_mean
            row[f'{condition}_rlm'] = node_value - neighborhood_mean

        if 'RNA_log2FC' in graph.nodes[node]:
            row['RNA_log2FC'] = float(graph.nodes[node]['RNA_log2FC'])

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    args = parse_arguments()
    os.makedirs(args.outputFolder, exist_ok=True)

    condition_to_replicates = _parse_condition_replicates(args.conditions, args.bigwigFiles)
    graph = nx.read_graphml(args.graph)

    annotate_graph_with_bigwig(graph, condition_to_replicates)
    annotate_graph_with_rna(graph, args.deseqCounts)

    rlm_df = create_rlm_table(graph, args.steps, args.conditions)
    rlm_df.to_csv(os.path.join(args.outputFolder, 'rlm_table.csv'), index=False)
    nx.write_graphml(graph, os.path.join(args.outputFolder, 'graph_annotated.graphml'))


if __name__ == '__main__':
    main()
