import argparse
import os

import networkx as nx
import numpy as np
import pandas as pd
import pyBigWig


def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge condition-wise RLM signal summaries")
    parser.add_argument('-g', '--graph', type=str, required=True, help='Path to the input graph file (GraphML)')
    parser.add_argument(
        '-n', '--conditions',
        type=str,
        nargs='+',
        required=True,
        help='Condition names. Must match --bigwigFiles length.',
    )
    parser.add_argument(
        '-b', '--bigwigFiles',
        type=str,
        nargs='+',
        required=True,
        help='Per-condition comma-separated bigWig replicates',
    )
    parser.add_argument('-o', '--outputFolder', type=str, required=True, help='Path to the output folder')
    return parser.parse_args()


def _parse_condition_replicates(conditions, bigwig_files):
    if len(conditions) != len(bigwig_files):
        raise ValueError('The number of --conditions entries must equal the number of --bigwigFiles entries.')
    mapping = {}
    for condition, files in zip(conditions, bigwig_files):
        reps = [x.strip() for x in files.split(',') if x.strip()]
        if len(reps) == 0:
            raise ValueError(f'No bigWig files provided for condition {condition}.')
        mapping[condition] = reps
    return mapping


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


def annotate_graph(graph, condition_to_replicates):
    for condition, replicates in condition_to_replicates.items():
        handles = [pyBigWig.open(rep) for rep in replicates]
        try:
            for node in graph.nodes:
                interval = _node_interval(graph.nodes[node])
                if interval is None:
                    graph.nodes[node][f'{condition}_mean'] = 0.0
                    continue
                chrom, start, end = interval
                values = [_bw_region_mean(h, chrom, start, end) for h in handles]
                graph.nodes[node][f'{condition}_mean'] = float(np.mean(values))
        finally:
            for h in handles:
                h.close()


def create_merged_table(graph, conditions):
    rows = []
    c1 = conditions[0] if len(conditions) >= 1 else None
    c2 = conditions[1] if len(conditions) >= 2 else None

    for node in graph.nodes:
        row = {
            'node_id': str(node),
            'name': graph.nodes[node].get('name', ''),
            'position': graph.nodes[node].get('position', ''),
        }
        for condition in conditions:
            row[f'{condition}_mean'] = float(graph.nodes[node].get(f'{condition}_mean', 0.0))

        if c1 and c2:
            numerator = row[f'{c2}_mean'] + 1e-9
            denominator = row[f'{c1}_mean'] + 1e-9
            row[f'{c2}_vs_{c1}_log2FC'] = float(np.log2(numerator / denominator))

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    args = parse_arguments()
    os.makedirs(args.outputFolder, exist_ok=True)

    condition_to_replicates = _parse_condition_replicates(args.conditions, args.bigwigFiles)
    graph = nx.read_graphml(args.graph)
    annotate_graph(graph, condition_to_replicates)

    merged_df = create_merged_table(graph, args.conditions)
    merged_df.to_csv(os.path.join(args.outputFolder, 'merged_rlm_table.csv'), index=False)
    nx.write_graphml(graph, os.path.join(args.outputFolder, 'graph_merged_rlm.graphml'))


if __name__ == '__main__':
    main()
