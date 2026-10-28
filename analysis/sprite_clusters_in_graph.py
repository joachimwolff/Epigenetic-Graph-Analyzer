"""Large SPRITE clusters seen from the graph: do the graph nodes of one cluster fall into one community?

For SPRITE clusters with at least 6 fragments (size classes 6-10, 11-100, 101-1000 fragments), the fragments on the
chromosome holding most of the cluster are mapped to the nodes of the interaction graph. For clusters with at least
three graph nodes, the script reports the co-membership (fraction of node pairs in the same community), the same
fraction for node pairs within 2 Mb (the range of the graph), the share of the largest community and the number of
communities, and the same measures for the cluster shifted circularly along
its chromosome (ten shifts, spacing of the fragments kept, shifted sets with at least three graph nodes).

Usage: python sprite_clusters_in_graph.py <sprite fragments .npz> <graph>_nodes.tsv <output prefix>
           [clusters per class, default 20000]
"""

import sys

import numpy as np
import pandas as pd
from scipy import stats

CLASSES = ((6, 10), (11, 100), (101, 1000))
N_SHIFTS = 10
MIN_NODES = 3
MAX_DISTANCE = 2000000


def measures(communities, positions):
    """Co-membership over all node pairs, largest community share, number of communities, number of nodes, and
    co-membership over node pairs within MAX_DISTANCE (the range of the graph) with the number of such pairs."""
    hit = communities >= 0
    members = communities[hit]
    if len(members) < MIN_NODES:
        return None
    counts = np.bincount(members)
    counts = counts[counts > 0]
    n = len(members)
    where = positions[hit]
    i, j = np.triu_indices(n, 1)
    close = np.abs(where[i] - where[j]) <= MAX_DISTANCE
    local = float((members[i][close] == members[j][close]).mean()) if close.any() else np.nan
    return (float((counts * (counts - 1)).sum() / (n * (n - 1))), float(counts.max() / n), len(counts), n, local,
            int(close.sum()))


def main():
    npz, nodes_path, prefix = sys.argv[1:4]
    per_class = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    data = np.load(npz)
    indptr, indices = data['indptr'], data['indices']
    chroms, starts = data['interval_chroms'], data['interval_starts']
    n_loci, n_groups = len(indptr) - 1, len(data['group_weight'])
    sizes = np.bincount(indices, minlength=n_groups)
    # locus -> community of the graph node on it (-1: no graph node)
    nodes = pd.read_csv(nodes_path, sep='\t', usecols=['chrom', 'start', 'community'])
    key = pd.Series(np.arange(n_loci), index=pd.MultiIndex.from_arrays([chroms, starts]))
    locus = key.reindex(pd.MultiIndex.from_arrays([nodes['chrom'], nodes['start']]))
    community = np.full(n_loci, -1, dtype=np.int64)
    found = locus.notna().to_numpy()
    community[locus[found].astype(int).to_numpy()] = nodes.loc[found, 'community'].to_numpy()
    print('graph nodes on SPRITE fragments: {:,} of {:,}'.format(int(found.sum()), len(nodes)), flush=True)
    chrom_codes, chrom_index = np.unique(chroms, return_inverse=True)
    chrom_start = np.array([np.argmax(chrom_index == c) for c in range(len(chrom_codes))])
    chrom_length = np.bincount(chrom_index)
    # entries of large clusters only, grouped by cluster
    rows = np.repeat(np.arange(n_loci, dtype=np.int32), np.diff(indptr))
    keep = sizes[indices] >= CLASSES[0][0]
    rows, groups = rows[keep], indices[keep]
    order = np.argsort(groups, kind='stable')
    rows, groups = rows[order], groups[order]
    del keep, order
    boundaries = np.flatnonzero(np.diff(groups)) + 1
    starts_g = np.concatenate([[0], boundaries])
    ends_g = np.concatenate([boundaries, [len(groups)]])
    cluster_ids = groups[starts_g]
    rng = np.random.default_rng(0)
    records = []
    for low, high in CLASSES:
        members = np.flatnonzero((sizes[cluster_ids] >= low) & (sizes[cluster_ids] <= high))
        if len(members) > per_class:
            members = rng.choice(members, per_class, replace=False)
        evaluated = 0
        for m in members:
            loci = rows[starts_g[m]:ends_g[m]]
            c = np.bincount(chrom_index[loci]).argmax()
            loci = loci[chrom_index[loci] == c]
            observed = measures(community[loci], starts[loci])
            if observed is None:
                continue
            evaluated += 1
            relative = loci - chrom_start[c]
            shifted = []
            for _ in range(N_SHIFTS):
                moved = chrom_start[c] + (relative + rng.integers(1, chrom_length[c])) % chrom_length[c]
                value = measures(community[moved], starts[moved])
                if value is not None:
                    shifted.append(value)
            if not shifted:
                continue
            shifted = np.array(shifted)
            records.append(dict(size_class='{}-{}'.format(low, high), fragments=int(sizes[cluster_ids[m]]),
                                fragments_on_chromosome=len(loci), graph_nodes=observed[3],
                                comembership=observed[0], largest_share=observed[1], communities=observed[2],
                                shifted_comembership=shifted[:, 0].mean(), shifted_largest_share=shifted[:, 1].mean(),
                                shifted_communities=shifted[:, 2].mean(), local_pairs=observed[5],
                                local_comembership=observed[4],
                                shifted_local_comembership=np.nanmean(shifted[:, 4]) if np.isfinite(shifted[:, 4]).any()
                                else np.nan,
                                span_mb=(starts[loci].max() - starts[loci].min()) / 1e6))
        print('{}-{}: {:,} sampled, {:,} with >= {} graph nodes'.format(low, high, len(members), evaluated,
                                                                       MIN_NODES), flush=True)
    table = pd.DataFrame(records)
    table.to_csv(prefix + '_clusters.tsv.gz', sep='\t', index=False)
    summary = []
    for size_class, sub in table.groupby('size_class', sort=False):
        local = sub.dropna(subset=['local_comembership', 'shifted_local_comembership'])
        summary.append(dict(
            clusters_with_local_pairs=len(local), local_comembership=local['local_comembership'].mean(),
            shifted_local_comembership=local['shifted_local_comembership'].mean(),
            local_concordance=float((np.sign(local['local_comembership'] - local['shifted_local_comembership'])
                                     + 1).mean() / 2),
            local_wilcoxon_p=stats.wilcoxon(local['local_comembership'],
                                            local['shifted_local_comembership']).pvalue if len(local) > 10 else np.nan,
            size_class=size_class, clusters=len(sub), median_graph_nodes=sub['graph_nodes'].median(),
            median_span_mb=sub['span_mb'].median(),
            comembership=sub['comembership'].mean(), shifted_comembership=sub['shifted_comembership'].mean(),
            concordance=float((np.sign(sub['comembership'] - sub['shifted_comembership']) + 1).mean() / 2),
            wilcoxon_p=stats.wilcoxon(sub['comembership'], sub['shifted_comembership']).pvalue,
            largest_share=sub['largest_share'].mean(), shifted_largest_share=sub['shifted_largest_share'].mean(),
            fraction_one_community_half=float((sub['largest_share'] >= 0.5).mean()),
            shifted_fraction_one_community_half=float((sub['shifted_largest_share'] >= 0.5).mean()),
            communities=sub['communities'].median(), shifted_communities=sub['shifted_communities'].median()))
    summary = pd.DataFrame(summary)
    summary.to_csv(prefix + '_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
