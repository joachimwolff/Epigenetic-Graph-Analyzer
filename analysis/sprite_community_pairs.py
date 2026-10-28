"""Do several graph communities come together in one large SPRITE cluster?

A community is the set of SPRITE fragments (loci) of its nodes. For SPRITE clusters with at least 6 fragments, every
pair of communities is scored by the number of clusters that contain loci of both. The reference is the same count
after shifting the complete community layout of each chromosome circularly by one random offset: every community
keeps its number of loci and their spacing, and cis pairs keep their distance, but the loci move to other fragments.
Communities sit in captured, active chromatin that is touched by more clusters than average chromatin, so the counts
are also compared as observed / expected ratios (n_ab N / (n_a n_b), with the marginal counts n_a, n_b of the same
layout), which measure whether two communities meet more often than their own frequencies predict. Reported per size class: the fraction of clusters touching two or more communities, and for cis and
trans pairs by distance the fraction of pairs above their shifted reference (concordance, 0.5 without relation) and
the median ratio. Groups of communities ("super-communities") are the connected components of pairs with at least
twice the shifted observed / expected ratio in at least MIN_CLUSTERS clusters; their sizes, number of chromosomes and genomic
span are written with their genes.

Usage: python sprite_community_pairs.py <sprite fragments .npz> <graph>_nodes.tsv <output prefix> [shifts, 5]
           [threads, default: available CPUs, at most 24; the shifted layouts are counted in parallel]
"""

import sys

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse

from epigeneticgraphanalyzer import parallel

CLASSES = ((6, 100), (101, 1000))
DISTANCE_BINS = np.array([0, 2, 5, 10, 20, 50, 300]) * 1e6
MIN_CLUSTERS = 5


def cluster_by_community(cluster_locus, members, n_comm):
    """Binary cluster x community incidence for communities given as (locus, community) pairs."""
    loci, comms = members
    assign = sparse.csr_matrix((np.ones(len(loci), dtype=np.float32), (loci, comms)),
                               shape=(cluster_locus.shape[1], n_comm))
    incidence = (cluster_locus @ assign).tocsr()
    incidence.data[:] = 1.0
    return incidence


def pair_counts(incidence):
    """Clusters with both communities, and the observed / expected ratio given how many clusters touch each."""
    touched = np.asarray(incidence.sum(axis=0)).ravel()
    together = (incidence.T @ incidence).tocoo()
    keep = together.row < together.col
    i, j, n = together.row[keep], together.col[keep], together.data[keep]
    index = pd.MultiIndex.from_arrays([i, j])
    return (pd.Series(n, index=index),
            pd.Series(n * incidence.shape[0] / (touched[i] * touched[j]), index=index))


def _shifted_counts(members):
    """pair_counts of one shifted community layout (worker function)."""
    return pair_counts(cluster_by_community(parallel.shared('cluster_locus'), members, parallel.shared('n_comm')))


def main():
    npz, nodes_path, prefix = sys.argv[1:4]
    n_shifts = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    threads = parallel.resolve_threads(int(sys.argv[5]) if len(sys.argv) > 5 else 0)
    data = np.load(npz)
    indptr, indices = data['indptr'], data['indices']
    chroms, starts = data['interval_chroms'], data['interval_starts']
    n_loci, n_groups = len(indptr) - 1, len(data['group_weight'])
    sizes = np.bincount(indices, minlength=n_groups)
    nodes = pd.read_csv(nodes_path, sep='\t', usecols=['node', 'chrom', 'start', 'end', 'community', 'genes'])
    key = pd.Series(np.arange(n_loci), index=pd.MultiIndex.from_arrays([chroms, starts]))
    locus = key.reindex(pd.MultiIndex.from_arrays([nodes['chrom'], nodes['start']])).to_numpy()
    found = ~np.isnan(locus)
    node_loci = locus[found].astype(np.int64)
    node_comms = nodes.loc[found, 'community'].to_numpy().astype(np.int64)
    n_comm = int(nodes['community'].max()) + 1
    info = nodes.groupby('community').agg(chrom=('chrom', 'first'), start=('start', 'min'), end=('end', 'max'),
                                          n_nodes=('node', 'size'),
                                          genes=('genes', lambda g: ','.join(sorted({x for v in g.dropna()
                                                                                     for x in str(v).split(',')}))))
    info = info.reindex(range(n_comm))
    centre = ((info['start'] + info['end']) / 2).to_numpy()
    comm_chrom = info['chrom'].to_numpy()
    chrom_codes, chrom_index = np.unique(chroms, return_inverse=True)
    chrom_start = np.array([np.argmax(chrom_index == c) for c in range(len(chrom_codes))])
    chrom_length = np.bincount(chrom_index)
    rng = np.random.default_rng(0)
    shifted_members = []
    for _ in range(n_shifts):
        offsets = rng.integers(1, chrom_length)
        c = chrom_index[node_loci]
        moved = chrom_start[c] + (node_loci - chrom_start[c] + offsets[c]) % chrom_length[c]
        shifted_members.append((moved, node_comms))
    # cluster x locus matrix (transpose of the stored locus x cluster matrix)
    locus_cluster = sparse.csr_matrix((np.ones(len(indices), dtype=np.float32), indices, indptr),
                                      shape=(n_loci, n_groups))
    summaries, tables = [], []
    for low, high in CLASSES:
        name = '{}-{}'.format(low, high)
        selected = np.flatnonzero((sizes >= low) & (sizes <= high))
        cluster_locus = locus_cluster[:, selected].T.tocsr()
        observed_inc = cluster_by_community(cluster_locus, (node_loci, node_comms), n_comm)
        touching = np.asarray(observed_inc.sum(axis=1)).ravel()
        observed, observed_oe = pair_counts(observed_inc)
        parallel.share(cluster_locus=cluster_locus, n_comm=n_comm)
        shifted_runs = parallel.map_chunks(_shifted_counts, shifted_members, threads)
        shifted = pd.concat([r[0] for r in shifted_runs], axis=1).reindex(observed.index).fillna(0).mean(axis=1)
        # observed / expected of the shifted layouts; pairs that never meet after a shift count with ratio 0
        shifted_oe = pd.concat([r[1] for r in shifted_runs], axis=1).reindex(observed.index).fillna(0).mean(axis=1)
        del cluster_locus
        a = observed.index.get_level_values(0).to_numpy()
        b = observed.index.get_level_values(1).to_numpy()
        table = pd.DataFrame({'community_a': a, 'community_b': b, 'clusters_both': observed.to_numpy().astype(int),
                              'shifted': shifted.to_numpy(), 'obs_exp': observed_oe.to_numpy(),
                              'shifted_obs_exp': shifted_oe.to_numpy()})
        table['ratio'] = (table['clusters_both'] + 1) / (table['shifted'] + 1)
        table['oe_ratio'] = table['obs_exp'] / table['shifted_obs_exp'].replace(0, np.nan)
        table['cis'] = comm_chrom[a] == comm_chrom[b]
        table['distance'] = np.where(table['cis'], np.abs(centre[a] - centre[b]), np.nan)
        table = table[(table['clusters_both'] >= MIN_CLUSTERS) | (table['shifted'] >= MIN_CLUSTERS)].copy()
        table['size_class'] = name
        tables.append(table)
        rows = []
        for label, sub in [('trans', table[~table['cis']])] + [
                ('cis {:.0f}-{:.0f} Mb'.format(DISTANCE_BINS[k] / 1e6, DISTANCE_BINS[k + 1] / 1e6),
                 table[table['cis'] & (table['distance'] >= DISTANCE_BINS[k])
                       & (table['distance'] < DISTANCE_BINS[k + 1])]) for k in range(len(DISTANCE_BINS) - 1)]:
            if len(sub):
                above = (np.sign(sub['clusters_both'] - sub['shifted']) + 1).mean() / 2
                above_oe = (np.sign(sub['obs_exp'] - sub['shifted_obs_exp']) + 1).mean() / 2
                rows.append('{}: {:,} pairs, counts {:.3f} (x{:.2f}), obs/exp {:.3f} (x{:.2f})'.format(
                    label, len(sub), above, sub['ratio'].median(), above_oe, sub['oe_ratio'].median()))
        strong = table[(table['oe_ratio'] >= 2) & (table['clusters_both'] >= MIN_CLUSTERS)]
        groups = nx.Graph()
        groups.add_edges_from(zip(strong['community_a'], strong['community_b']))
        components = sorted(nx.connected_components(groups), key=len, reverse=True)
        group_rows = []
        for k, members in enumerate(components):
            members = sorted(members)
            sub = info.loc[members]
            group_rows.append(dict(size_class=name, group=k, communities=len(members),
                                   chromosomes=sub['chrom'].nunique(), nodes=int(sub['n_nodes'].sum()),
                                   span_mb=float((sub['end'].max() - sub['start'].min()) / 1e6)
                                   if sub['chrom'].nunique() == 1 else np.nan,
                                   chrom_list=','.join(sorted(sub['chrom'].unique())),
                                   community_ids=','.join(str(m) for m in members),
                                   genes=','.join(g for g in ','.join(sub['genes'].fillna('')).split(',') if g)[:2000]))
        pd.DataFrame(group_rows).to_csv('{}_groups_{}.tsv'.format(prefix, name), sep='\t', index=False)
        summary = dict(size_class=name, clusters=len(selected),
                       fraction_touching_2plus=float((touching >= 2).mean()),
                       median_communities_touched=float(np.median(touching[touching >= 1])),
                       pairs=len(table), cis_pairs=int(table['cis'].sum()),
                       concordance_all=float((np.sign(table['clusters_both'] - table['shifted']) + 1).mean() / 2),
                       concordance_obs_exp=float((np.sign(table['obs_exp'] - table['shifted_obs_exp']) + 1).mean() / 2),
                       pairs_ratio_2plus=len(strong), trans_pairs_ratio_2plus=int((~strong['cis']).sum()),
                       groups=len(components), group_sizes=','.join(str(len(c)) for c in components[:8]),
                       by_distance=' | '.join(rows))
        summaries.append(summary)
        for k, v in summary.items():
            print('{}: {}'.format(k, v))
        print(flush=True)
    pd.concat(tables).to_csv(prefix + '_community_pairs.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(prefix + '_community_pairs_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
