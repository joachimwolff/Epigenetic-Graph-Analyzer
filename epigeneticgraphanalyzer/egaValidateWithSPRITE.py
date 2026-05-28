import argparse
from collections import Counter
from itertools import combinations
import gzip

import networkx as nx
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Validate indirect graph mappings with SPRITE triplets')
    parser.add_argument('--sprite', required=True, type=str, help='Path to the SPRITE cluster file (.gz)')
    parser.add_argument('--graph', required=True, type=str, help='Path to the graph input file (GraphML)')
    parser.add_argument('--genomeSize', required=True, type=str, help='Path to the genome size file (tab-separated)')
    parser.add_argument('--outputPrefix', required=True, type=str, help='Output path prefix for results')
    parser.add_argument('--binSize', required=False, type=int, default=10000, help='Bin size for genomic intervals (default: 10000)')
    parser.add_argument('--minSpriteCount', required=False, type=int, default=1, help='Minimum SPRITE triplet count to consider a graph triplet supported')
    parser.add_argument(
        '--memoryMode',
        required=False,
        choices=['targeted', 'full'],
        default='targeted',
        help='targeted: count only graph triplets (low memory), full: count all SPRITE triplets (high memory)',
    )
    parser.add_argument(
        '--maxClusterBins',
        required=False,
        type=int,
        default=0,
        help='Skip clusters with more than this many unique bins (0 disables this filter)',
    )
    return parser.parse_args()


def read_chrom_sizes(p_genome_size_file):
    chrom_sizes = {}
    genome_size_df = pd.read_csv(
        p_genome_size_file,
        sep='\t',
        header=None,
        names=['chromosome', 'size']
    )
    for _, row in genome_size_df.iterrows():
        chrom_sizes[str(row['chromosome'])] = int(row['size'])
    return chrom_sizes


def build_bin_offsets(p_chrom_sizes, p_bin_size):
    offsets = {}
    current = 0
    for chrom, size in p_chrom_sizes.items():
        offsets[chrom] = current
        current += (size + p_bin_size - 1) // p_bin_size
    return offsets


def locus_to_bin_id(p_chrom, p_position, p_chrom_sizes, p_bin_offsets, p_bin_size):
    if p_chrom not in p_chrom_sizes:
        return None
    if p_position < 0 or p_position >= p_chrom_sizes[p_chrom]:
        return None
    return p_bin_offsets[p_chrom] + (p_position // p_bin_size)


def parse_position_attribute(p_position):
    chrom, value = p_position.split(':', 1)
    start, end = value.split('-', 1)
    start_i = int(start)
    end_i = int(end)
    if end_i < start_i:
        start_i, end_i = end_i, start_i
    midpoint = start_i + ((end_i - start_i) // 2)
    return chrom, midpoint


def build_node_bin_mapping(p_graph, p_chrom_sizes, p_bin_offsets, p_bin_size):
    mapping = {}
    for node, attrs in p_graph.nodes(data=True):
        position = attrs.get('position')
        if position is None:
            continue
        try:
            chrom, midpoint = parse_position_attribute(position)
        except Exception:
            continue
        bin_id = locus_to_bin_id(chrom, midpoint, p_chrom_sizes, p_bin_offsets, p_bin_size)
        if bin_id is not None:
            mapping[node] = bin_id
    return mapping


def count_sprite_triplets(p_sprite_file, p_chrom_sizes, p_bin_offsets, p_bin_size):
    triplet_counts = Counter()
    pair_counts = Counter()
    with gzip.open(p_sprite_file, 'rt') as handle:
        for line in handle:
            tokens = line.strip().split('\t')
            if len(tokens) < 3:
                continue

            loci_tokens = tokens[1:]
            cluster_bins = set()
            for locus in loci_tokens:
                if ':' not in locus:
                    continue
                chrom, start = locus.split(':', 1)
                try:
                    start_i = int(start)
                except ValueError:
                    continue
                bin_id = locus_to_bin_id(chrom, start_i, p_chrom_sizes, p_bin_offsets, p_bin_size)
                if bin_id is not None:
                    cluster_bins.add(bin_id)

            if len(cluster_bins) < 2:
                continue

            sorted_bins = sorted(cluster_bins)
            for pair in combinations(sorted_bins, 2):
                pair_counts[pair] += 1

            if len(sorted_bins) < 3:
                continue

            for triplet in combinations(sorted_bins, 3):
                triplet_counts[triplet] += 1
    return triplet_counts, pair_counts


def _build_triplet_index(target_triplets):
    index = {}
    for triplet in target_triplets:
        for bin_id in triplet:
            if bin_id not in index:
                index[bin_id] = set()
            index[bin_id].add(triplet)
    return index


def count_target_triplets(
    p_sprite_file,
    p_chrom_sizes,
    p_bin_offsets,
    p_bin_size,
    p_target_triplets,
    p_max_cluster_bins=0,
):
    counts = Counter()
    if len(p_target_triplets) == 0:
        return counts

    triplet_index = _build_triplet_index(p_target_triplets)

    with gzip.open(p_sprite_file, 'rt') as handle:
        for line in handle:
            tokens = line.strip().split('\t')
            if len(tokens) < 4:
                continue

            cluster_bins = set()
            for locus in tokens[1:]:
                if ':' not in locus:
                    continue
                chrom, start = locus.split(':', 1)
                try:
                    start_i = int(start)
                except ValueError:
                    continue
                bin_id = locus_to_bin_id(chrom, start_i, p_chrom_sizes, p_bin_offsets, p_bin_size)
                if bin_id is not None:
                    cluster_bins.add(bin_id)

            if len(cluster_bins) < 3:
                continue
            if p_max_cluster_bins > 0 and len(cluster_bins) > p_max_cluster_bins:
                continue

            candidate_triplets = set()
            for bin_id in cluster_bins:
                candidate_triplets.update(triplet_index.get(bin_id, set()))

            for triplet in candidate_triplets:
                if triplet[0] in cluster_bins and triplet[1] in cluster_bins and triplet[2] in cluster_bins:
                    counts[triplet] += 1

    return counts


def build_graph_triplets(p_graph, p_node_bin_mapping):
    graph_triplets = []
    seen_node_triplets = set()
    for middle in p_graph.nodes:
        neighbors = list(p_graph.neighbors(middle))
        if len(neighbors) < 2:
            continue
        for left, right in combinations(neighbors, 2):
            node_triplet = tuple(sorted((left, middle, right), key=str))
            if node_triplet in seen_node_triplets:
                continue
            seen_node_triplets.add(node_triplet)

            if left not in p_node_bin_mapping or middle not in p_node_bin_mapping or right not in p_node_bin_mapping:
                continue

            bin_triplet = tuple(sorted((
                p_node_bin_mapping[left],
                p_node_bin_mapping[middle],
                p_node_bin_mapping[right],
            )))
            if len(set(bin_triplet)) < 3:
                continue

            graph_triplets.append((node_triplet, bin_triplet, middle, left, right))
    return graph_triplets


def validate_graph_triplets(p_graph_triplets, p_sprite_triplet_counts, p_min_sprite_count):
    rows = []
    supported = 0
    for node_triplet, bin_triplet, middle, left, right in p_graph_triplets:
        sprite_count = p_sprite_triplet_counts.get(bin_triplet, 0)
        is_supported = int(sprite_count >= p_min_sprite_count)
        supported += is_supported
        rows.append({
            'node_a': str(node_triplet[0]),
            'node_b': str(node_triplet[1]),
            'node_c': str(node_triplet[2]),
            'path_middle_node': str(middle),
            'path_left_node': str(left),
            'path_right_node': str(right),
            'bin_a': bin_triplet[0],
            'bin_b': bin_triplet[1],
            'bin_c': bin_triplet[2],
            'sprite_triplet_count': sprite_count,
            'supported': is_supported,
        })
    return rows, supported


def write_results(p_rows, p_total_triplets, p_supported_triplets, p_output_prefix, p_min_sprite_count):
    details_df = pd.DataFrame(p_rows)
    details_path = p_output_prefix + '_indirect_triplets.tsv'
    details_df.to_csv(details_path, sep='\t', index=False)

    support_rate = 0.0
    if p_total_triplets > 0:
        support_rate = p_supported_triplets / p_total_triplets

    summary_df = pd.DataFrame([{
        'total_indirect_triplets': p_total_triplets,
        'supported_indirect_triplets': p_supported_triplets,
        'unsupported_indirect_triplets': p_total_triplets - p_supported_triplets,
        'support_rate': support_rate,
        'min_sprite_count': p_min_sprite_count,
    }])
    summary_path = p_output_prefix + '_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)


def main():
    args = parse_arguments()
    graph = nx.read_graphml(args.graph)

    chrom_sizes = read_chrom_sizes(args.genomeSize)
    bin_offsets = build_bin_offsets(chrom_sizes, args.binSize)
    node_bin_mapping = build_node_bin_mapping(graph, chrom_sizes, bin_offsets, args.binSize)

    graph_triplets = build_graph_triplets(graph, node_bin_mapping)

    if args.memoryMode == 'full':
        sprite_triplet_counts, _ = count_sprite_triplets(
            args.sprite,
            chrom_sizes,
            bin_offsets,
            args.binSize,
        )
    else:
        target_triplets = {x[1] for x in graph_triplets}
        sprite_triplet_counts = count_target_triplets(
            args.sprite,
            chrom_sizes,
            bin_offsets,
            args.binSize,
            target_triplets,
            p_max_cluster_bins=args.maxClusterBins,
        )

    rows, supported = validate_graph_triplets(graph_triplets, sprite_triplet_counts, args.minSpriteCount)
    write_results(rows, len(graph_triplets), supported, args.outputPrefix, args.minSpriteCount)


if __name__ == '__main__':
    main()
