"""Functional validation of graph paths with expression quantitative trait loci (eQTLs).

An eQTL pair is a variant-bearing fragment and the baited promoter fragment(s) of the gene it regulates. Every pair
is placed on a fragment map (restriction fragments or bins of the assembly of the graph), and its graph distance is
measured: a direct edge (1 step), a path over one mediator (2 steps), a path over two mediators (3 steps), a longer
path, or no connection. Controls are drawn per eQTL pair at a matched genomic distance from two pools:

  control_gene      the same variant fragment with the promoter of another gene that is not an eGene of any
                    significant variant on that fragment;
  control_fragment  the same gene with another fragment that carries a significant variant of the data set (for
                    any gene) but none for this gene.

Positions are 1-based throughout (eQTL tables, rmap fragments); TSSs of `io.read_tss_from_gtf` (0-based) are shifted
by one when promoters are built.
"""

import gzip
import io as _io
import logging
import tarfile
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from scipy.sparse.csgraph import connected_components

from .genome import IntervalIndex, normalize_chrom
from .validation import block_resample_items, block_resamples

log = logging.getLogger(__name__)

DISTANCE_EDGES = (0, 20000, 50000, 100000, 200000, 500000, 1000001)
AUTOSOMES = tuple('chr{}'.format(i) for i in range(1, 23))
GROUPS = ('eqtl', 'control', 'control_gene', 'control_fragment')
STEP_LABELS = ('1', '2', '3', '>3', 'none')


def distance_bin(distance):
    labels = np.array(['{}-{}kb'.format(DISTANCE_EDGES[i] // 1000, DISTANCE_EDGES[i + 1] // 1000)
                       for i in range(len(DISTANCE_EDGES) - 1)])
    index = np.clip(np.searchsorted(DISTANCE_EDGES, np.asarray(distance), side='right') - 1, 0, len(labels) - 1)
    return labels[index]


def distance_bin_labels():
    return list(distance_bin(np.array(DISTANCE_EDGES[:-1])))


# ---------------------------------------------------------------------------------------------------------------
# eQTL readers: every reader returns chrom, pos (1-based), gene_id (without version) and a ranking column


def _gene(values):
    return pd.Series(values).astype(str).str.split('.').str[0].to_numpy()


def read_eqtlgen(path):
    """eQTLGen cis-eQTL file (hg19): chrom, pos, gene_id, pvalue, abs_z."""
    df = pd.read_csv(path, sep='\t', usecols=['Pvalue', 'SNPChr', 'SNPPos', 'Gene', 'Zscore'],
                     dtype={'SNPChr': str})
    return pd.DataFrame({'chrom': df['SNPChr'].map(normalize_chrom).to_numpy(), 'pos': df['SNPPos'].astype(np.int64),
                         'gene_id': _gene(df['Gene']), 'pvalue': df['Pvalue'].astype(float),
                         'abs_z': df['Zscore'].abs().astype(float)})


def _split_variant(variant_ids):
    parts = pd.Series(variant_ids).astype(str).str.split('_', expand=True)
    return parts[0].map(normalize_chrom).to_numpy(), parts[1].astype(np.int64).to_numpy()


def read_gtex_pairs(path, tissue=None):
    """GTEx v8 signif_variant_gene_pairs (GRCh38): chrom, pos, gene_id, pvalue, abs_z (slope/se).

    `path` is the per-tissue .txt.gz or the GTEx_Analysis_v8_eQTL.tar holding it (read without extraction)."""
    if str(path).endswith('.tar'):
        with tarfile.open(path) as archive:
            member = next(m for m in archive.getmembers()
                          if m.name.endswith('/{}.v8.signif_variant_gene_pairs.txt.gz'.format(tissue)))
            raw = archive.extractfile(member).read()
        df = pd.read_csv(_io.BytesIO(gzip.decompress(raw)), sep='\t')
    else:
        df = pd.read_csv(path, sep='\t')
    chrom, pos = _split_variant(df['variant_id'])
    return pd.DataFrame({'chrom': chrom, 'pos': pos, 'gene_id': _gene(df['gene_id']),
                         'pvalue': df['pval_nominal'].astype(float),
                         'abs_z': (df['slope'] / df['slope_se']).abs().astype(float)})


def read_susie(path):
    """eQTL Catalogue SuSiE credible sets (GRCh38): chrom, pos, gene_id, pip, cs_id, pvalue."""
    df = pd.read_csv(path, sep='\t', usecols=['gene_id', 'cs_id', 'variant', 'pip', 'pvalue'])
    chrom, pos = _split_variant(df['variant'])
    return pd.DataFrame({'chrom': chrom, 'pos': pos, 'gene_id': _gene(df['gene_id']), 'pip': df['pip'].astype(float),
                         'cs_id': df['cs_id'].astype(str), 'pvalue': df['pvalue'].astype(float)})


def lead_per_gene(df):
    """The variant with the smallest P value per gene (ties: largest |z|, then position)."""
    order = df.sort_values(['gene_id', 'pvalue', 'abs_z', 'pos'], ascending=[True, True, False, True])
    return order.drop_duplicates('gene_id').reset_index(drop=True)


def fine_mapped(df, min_pip=0.5):
    return df[df['pip'] >= min_pip].reset_index(drop=True)


def lift_positions(df, chain_file):
    """Lift chrom/pos (1-based) with a chain file; unliftable rows and rows that change chromosome are dropped."""
    from pyliftover import LiftOver
    lift = LiftOver(chain_file)
    keys = df[['chrom', 'pos']].drop_duplicates()
    new_chrom, new_pos = [], []
    for chrom, pos in zip(keys['chrom'].to_numpy(), keys['pos'].to_numpy()):
        hits = lift.convert_coordinate(chrom, int(pos) - 1)
        if hits and hits[0][0] == chrom:
            new_chrom.append(hits[0][0])
            new_pos.append(int(hits[0][1]) + 1)
        else:
            new_chrom.append(None)
            new_pos.append(-1)
    keys = keys.assign(new_chrom=new_chrom, new_pos=new_pos)
    out = df.merge(keys, on=['chrom', 'pos'], how='left')
    out = out[out['new_pos'] > 0].copy()
    out['chrom'] = out.pop('new_chrom')
    out['pos'] = out.pop('new_pos').astype(np.int64)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------------------------------------------
# fragment maps and graphs


def read_rmap(path):
    """rmap / BED-like fragment file: chrom, start, end (1-based closed as in CHiCAGO designs)."""
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], comment='#', dtype={0: str})
    return pd.DataFrame({'chrom': df[0].map(normalize_chrom), 'start': df[1].astype(np.int64),
                         'end': df[2].astype(np.int64)})


def tile_bins(chrom_sizes, bin_size):
    """Fixed bins as fragments, 1-based closed (start k*size+1, end (k+1)*size), as the virtual capture rmap."""
    rows = []
    for chrom, size in chrom_sizes.items():
        starts = np.arange(0, size, bin_size, dtype=np.int64)
        rows.append(pd.DataFrame({'chrom': chrom, 'start': starts + 1, 'end': np.minimum(starts + bin_size, size)}))
    return pd.concat(rows, ignore_index=True)


def lift_fragments(fragments, chain_file, max_length_change=0.1):
    """Lift fragment ends with a chain file; fragments whose ends land on another chromosome, invert, or change
    length by more than `max_length_change` (relative) are dropped."""
    from pyliftover import LiftOver
    lift = LiftOver(chain_file)
    keep = []
    for chrom, start, end in zip(fragments['chrom'].to_numpy(), fragments['start'].to_numpy(),
                                 fragments['end'].to_numpy()):
        a = lift.convert_coordinate(chrom, int(start) - 1)
        b = lift.convert_coordinate(chrom, int(end) - 1)
        if not a or not b or a[0][0] != chrom or b[0][0] != chrom:
            continue
        new_start, new_end = int(a[0][1]) + 1, int(b[0][1]) + 1
        length = end - start
        if new_end <= new_start or abs((new_end - new_start) - length) > max_length_change * length + 10:
            continue
        keep.append((chrom, new_start, new_end))
    return pd.DataFrame(keep, columns=['chrom', 'start', 'end'])


def _drop_overlaps(fragments):
    fragments = fragments.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    keep = np.ones(len(fragments), dtype=bool)
    for _, index in fragments.groupby('chrom').indices.items():
        last_end = -1
        for i in index:
            if fragments.at[i, 'start'] <= last_end:
                keep[i] = False
            else:
                last_end = fragments.at[i, 'end']
    return fragments[keep].reset_index(drop=True)


def overlay_nodes(fragments, nodes):
    """Fragment map with the graph nodes in place of the fragments they overlap (nodes win)."""
    nodes = nodes[['chrom', 'start', 'end']].drop_duplicates()
    nodes = _drop_overlaps(nodes)
    index = FragmentMap(nodes)
    rest = []
    for chrom, group in fragments.groupby('chrom'):
        if chrom not in index.index.offsets:
            rest.append(group)
            continue
        first = index.index.offsets[chrom]
        count = index.index.n_chrom_bins(chrom)
        starts = index.start[first:first + count]
        ends = index.end[first:first + count]
        # a fragment overlaps a node if the first node ending at/after its start begins at/before its end
        j = np.searchsorted(ends, group['start'].to_numpy(), side='left')
        overlaps = (j < count) & (starts[np.clip(j, 0, count - 1)] <= group['end'].to_numpy())
        rest.append(group[~overlaps])
    return _drop_overlaps(pd.concat([nodes] + rest, ignore_index=True))


class FragmentMap:
    """Non-overlapping 1-based closed fragments with a global index ordered by chromosome and start."""

    def __init__(self, fragments):
        fragments = fragments.sort_values(['chrom', 'start']).reset_index(drop=True)
        self.index = IntervalIndex(fragments['chrom'].to_numpy(), fragments['start'].to_numpy(),
                                   fragments['end'].to_numpy() + 1)
        self.start = self.index.starts
        self.end = self.index.ends - 1
        self.chrom = np.array(self.index.chrom_names)[self.index.chrom_index_of_bins(np.arange(self.index.n_bins))]
        self.n = self.index.n_bins

    def locate(self, chroms, positions):
        """Fragment index of every position (-1 outside all fragments)."""
        return self.index.bins_of(np.asarray(chroms).astype(str), np.asarray(positions, dtype=np.int64))

    def midpoint(self, fragments):
        fragments = np.asarray(fragments, dtype=np.int64)
        return (self.start[fragments] + self.end[fragments]) // 2

    def length(self, fragments):
        fragments = np.asarray(fragments, dtype=np.int64)
        return self.end[fragments] - self.start[fragments] + 1


def locate_nodes(fmap, chroms, starts, ends):
    """Fragment index of graph nodes: the fragment holding the node midpoint; also returns how many nodes match a
    fragment exactly."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    index = fmap.locate(chroms, (starts + ends) // 2)
    exact = (index >= 0) & (fmap.start[np.clip(index, 0, None)] == starts) & (fmap.end[np.clip(index, 0, None)] == ends)
    return index, int(exact.sum())


class FragmentGraph:
    """Undirected simple graph on fragment indices (CSR adjacency)."""

    def __init__(self, n, a, b):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        keep = (a >= 0) & (b >= 0) & (a != b)
        a, b = a[keep], b[keep]
        matrix = sp.coo_matrix((np.ones(2 * len(a), dtype=np.int8), (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()
        matrix.sum_duplicates()
        matrix.data[:] = 1
        self.adjacency = matrix
        self.n = n
        self.degree = np.diff(matrix.indptr)
        self.n_edges = int(matrix.nnz // 2)
        _, self.component = connected_components(matrix, directed=False)
        self._sets = {}

    def neighbours(self, i):
        return self.adjacency.indices[self.adjacency.indptr[i]:self.adjacency.indptr[i + 1]]

    def neighbour_set(self, i):
        result = self._sets.get(i)
        if result is None:
            result = frozenset(self.neighbours(i).tolist())
            self._sets[i] = result
        return result

    def edges(self):
        upper = sp.triu(self.adjacency, k=1).tocoo()
        return upper.row.astype(np.int64), upper.col.astype(np.int64)

    def mirrored(self, fmap, bait):
        """Control graph with every edge's partner end reflected around its anchor (the bait end; the lower index
        if both or neither end is a bait), so that baits keep their degree and the edges their genomic length
        while the partners move to loci the data do not connect."""
        a, b = self.edges()
        swap = ~bait[a] & bait[b]
        anchor = np.where(swap, b, a)
        partner = np.where(swap, a, b)
        position = 2 * fmap.midpoint(anchor) - fmap.midpoint(partner)
        target = np.full(len(anchor), -1, dtype=np.int64)
        valid = position > 0
        target[valid] = fmap.locate(fmap.chrom[anchor[valid]], position[valid])
        return FragmentGraph(self.n, anchor, target)


def read_graph_tables(nodes_path, edges_path):
    nodes = pd.read_csv(nodes_path, sep='\t', usecols=['node', 'chrom', 'start', 'end', 'bait'])
    nodes['chrom'] = nodes['chrom'].astype(str).map(normalize_chrom)
    nodes['bait'] = nodes['bait'].astype(str).str.lower().isin(['true', '1'])
    edges = pd.read_csv(edges_path, sep='\t', usecols=['node_a', 'node_b'])
    return nodes, edges


def graph_on_fragments(fmap, nodes, edges):
    """FragmentGraph of a node/edge table pair, plus mapping statistics."""
    node_index, exact = locate_nodes(fmap, nodes['chrom'], nodes['start'], nodes['end'])
    lookup = dict(zip(nodes['node'], node_index))
    a = edges['node_a'].map(lookup).fillna(-1).astype(np.int64).to_numpy()
    b = edges['node_b'].map(lookup).fillna(-1).astype(np.int64).to_numpy()
    graph = FragmentGraph(fmap.n, a, b)
    info = {'graph_nodes': len(nodes), 'graph_nodes_exact_fragment': exact,
            'graph_nodes_unmapped': int((node_index < 0).sum()), 'graph_edges': len(edges),
            'graph_edges_on_fragments': graph.n_edges}
    return graph, node_index, info


# ---------------------------------------------------------------------------------------------------------------
# promoters and pairs


def gene_promoters(tss, fmap, bait):
    """Baited promoter fragments per gene: {gene_id: (chrom, sorted baited TSS positions, baited fragments)}.

    `tss` as from io.read_tss_from_gtf (0-based TSS); genes without a TSS on a baited fragment are left out."""
    chroms = tss['chrom'].to_numpy().astype(str)
    positions = tss['tss'].to_numpy().astype(np.int64) + 1
    fragments = fmap.locate(chroms, positions)
    baited = fragments >= 0
    baited[baited] = bait[fragments[baited]]
    table = pd.DataFrame({'gene_id': tss['gene_id'].to_numpy(), 'chrom': chroms, 'pos': positions,
                          'fragment': fragments})[baited]
    result = {}
    for gene, group in table.groupby('gene_id'):
        if group['chrom'].nunique() > 1:
            group = group[group['chrom'] == group['chrom'].mode().iloc[0]]
        result[gene] = (group['chrom'].iloc[0], np.sort(group['pos'].to_numpy()),
                        np.unique(group['fragment'].to_numpy()))
    return result


def _nearest(positions, sorted_targets):
    """Distance from every position to the nearest of the sorted targets."""
    positions = np.asarray(positions, dtype=np.int64)
    j = np.searchsorted(sorted_targets, positions)
    left = np.abs(positions - sorted_targets[np.clip(j - 1, 0, len(sorted_targets) - 1)])
    right = np.abs(sorted_targets[np.clip(j, 0, len(sorted_targets) - 1)] - positions)
    return np.minimum(left, right)


def map_pairs(eqtls, fmap, promoters, all_genes=None, max_distance=1000000, chroms=AUTOSOMES):
    """Map variant-gene pairs to (variant fragment, gene) pairs.

    Returns the pairs (chrom, var_frag, var_pos, gene_id, distance; one row per variant fragment and gene, the
    variant closest to the promoter kept) and a table of how many input rows were dropped at each step.
    `all_genes` (gene ids with any annotated TSS) separates unannotated genes from genes without a baited promoter."""
    counts = [('input_pairs', len(eqtls))]
    df = eqtls[['chrom', 'pos', 'gene_id']].copy()
    keep = df['chrom'].isin(chroms)
    counts.append(('not_autosomal', int((~keep).sum())))
    df = df[keep]
    df['var_frag'] = fmap.locate(df['chrom'], df['pos'])
    keep = df['var_frag'] >= 0
    counts.append(('variant_not_on_fragment', int((~keep).sum())))
    df = df[keep]
    has_promoter = df['gene_id'].isin(promoters.keys())
    if all_genes is not None:
        annotated = df['gene_id'].isin(all_genes)
        counts.append(('gene_not_annotated', int((~annotated).sum())))
        counts.append(('gene_promoter_not_baited', int((annotated & ~has_promoter).sum())))
    else:
        counts.append(('gene_without_baited_promoter', int((~has_promoter).sum())))
    df = df[has_promoter]
    rows = []
    status = {'other_chromosome': 0, 'variant_in_promoter_fragment': 0, 'beyond_max_distance': 0}
    for chrom, pos, gene, frag in zip(df['chrom'].to_numpy(), df['pos'].to_numpy(), df['gene_id'].to_numpy(),
                                      df['var_frag'].to_numpy()):
        p_chrom, p_pos, p_frag = promoters[gene]
        if p_chrom != chrom:
            status['other_chromosome'] += 1
            continue
        if frag in p_frag:
            status['variant_in_promoter_fragment'] += 1
            continue
        distance = int(_nearest([pos], p_pos)[0])
        if distance > max_distance:
            status['beyond_max_distance'] += 1
            continue
        rows.append((chrom, int(frag), int(pos), gene, distance))
    counts += list(status.items())
    pairs = pd.DataFrame(rows, columns=['chrom', 'var_frag', 'var_pos', 'gene_id', 'distance'])
    before = len(pairs)
    pairs = pairs.sort_values('distance').drop_duplicates(['var_frag', 'gene_id']).sort_values(
        ['chrom', 'var_pos', 'gene_id']).reset_index(drop=True)
    counts.append(('duplicate_fragment_gene', before - len(pairs)))
    counts.append(('pairs', len(pairs)))
    return pairs, pd.DataFrame(counts, columns=['step', 'n'])


class PairSet:
    """Membership test for (fragment, gene) pairs, e.g. all significant pairs of a data set."""

    def __init__(self, fragments, genes):
        self.genes = {g: i for i, g in enumerate(pd.unique(np.asarray(genes)))}
        self.codes = np.unique(self._code(fragments, genes))

    def _code(self, fragments, genes):
        gene_index = np.array([self.genes.get(g, -1) for g in genes], dtype=np.int64)
        fragments = np.asarray(fragments, dtype=np.int64)
        return np.where(gene_index >= 0, fragments * (len(self.genes) + 1) + gene_index, -1)

    def contains(self, fragments, genes):
        codes = self._code(np.atleast_1d(fragments), np.atleast_1d(genes))
        j = np.searchsorted(self.codes, codes)
        found = (j < len(self.codes)) & (self.codes[np.clip(j, 0, len(self.codes) - 1)] == codes)
        return found & (codes >= 0)


def significant_pairs(significant, fmap):
    """PairSet of the fragments of all significant variants of a data set and their genes."""
    frag = fmap.locate(significant['chrom'], significant['pos'])
    keep = frag >= 0
    return PairSet(frag[keep], significant['gene_id'].to_numpy()[keep])


def variant_fragments(significant, fmap, chroms=AUTOSOMES):
    """Fragments carrying a significant variant (any gene), with the leftmost such variant as position."""
    df = significant[significant['chrom'].isin(chroms)][['chrom', 'pos']].drop_duplicates()
    df = df.assign(frag=fmap.locate(df['chrom'], df['pos']))
    df = df[df['frag'] >= 0].sort_values('pos').drop_duplicates('frag')
    return df.sort_values(['chrom', 'pos']).reset_index(drop=True)


def sample_controls(pairs, promoters, universe, candidates, significant, n_controls=5, window=0.1, seed=0):
    """Distance-matched non-eQTL pairs for every eQTL pair (row index of `pairs` = set_id).

    control_gene: the variant fragment with genes of `universe` (gene ids with a baited promoter) whose nearest
    baited TSS lies within a factor 10**window of the eQTL distance; control_fragment: fragments of `candidates`
    (chrom, pos, frag) within that distance window of the gene's baited TSSs. Pairs in `significant` (PairSet), and
    fragments that are a promoter fragment of the gene, are never controls."""
    rng = np.random.default_rng(seed)
    factor = 10.0 ** window
    tss_rows = [(promoters[g][0], p, g) for g in universe if g in promoters for p in promoters[g][1]]
    tss = pd.DataFrame(tss_rows, columns=['chrom', 'pos', 'gene_id']).sort_values(['chrom', 'pos'])
    tss_by_chrom = {c: (g['pos'].to_numpy(), g['gene_id'].to_numpy()) for c, g in tss.groupby('chrom')}
    cand_by_chrom = {c: (g['pos'].to_numpy(), g['frag'].to_numpy()) for c, g in candidates.groupby('chrom')}
    rows = []
    for set_id, (chrom, frag, pos, gene, distance) in enumerate(zip(
            pairs['chrom'].to_numpy(), pairs['var_frag'].to_numpy(), pairs['var_pos'].to_numpy(),
            pairs['gene_id'].to_numpy(), pairs['distance'].to_numpy())):
        low, high = distance / factor, distance * factor
        # same variant fragment, other genes
        if chrom in tss_by_chrom:
            t_pos, t_gene = tss_by_chrom[chrom]
            a, b = np.searchsorted(t_pos, [pos - high, pos + high + 1])
            if b > a:
                window_genes = pd.DataFrame({'gene_id': t_gene[a:b], 'd': np.abs(t_pos[a:b] - pos)})
                nearest = window_genes.groupby('gene_id', sort=True)['d'].min()
                nearest = nearest[(nearest >= low) & (nearest <= high) & (nearest.index != gene)]
                if len(nearest):
                    names = nearest.index.to_numpy()
                    ok = ~significant.contains(np.full(len(names), frag), names)
                    ok &= np.array([frag not in promoters[g][2] for g in names])
                    names, dists = names[ok], nearest.to_numpy()[ok]
                    if len(names) > n_controls:
                        pick = np.sort(rng.choice(len(names), n_controls, replace=False))
                        names, dists = names[pick], dists[pick]
                    rows += [(set_id, 'control_gene', chrom, frag, pos, g, int(d)) for g, d in zip(names, dists)]
        # same gene, other variant-bearing fragments
        if chrom in cand_by_chrom:
            c_pos, c_frag = cand_by_chrom[chrom]
            g_pos, g_frags = promoters[gene][1], promoters[gene][2]
            a, b = np.searchsorted(c_pos, [g_pos[0] - high, g_pos[-1] + high + 1])
            if b > a:
                d = _nearest(c_pos[a:b], g_pos)
                ok = (d >= low) & (d <= high) & (c_frag[a:b] != frag) & ~np.isin(c_frag[a:b], g_frags)
                idx = np.nonzero(ok)[0]
                if len(idx):
                    idx = idx[~significant.contains(c_frag[a:b][idx], np.full(len(idx), gene, dtype=object))]
                if len(idx) > n_controls:
                    idx = np.sort(rng.choice(idx, n_controls, replace=False))
                rows += [(set_id, 'control_fragment', chrom, int(c_frag[a + i]), int(c_pos[a + i]), gene, int(d[i]))
                         for i in idx]
    controls = pd.DataFrame(rows, columns=['set_id', 'group', 'chrom', 'var_frag', 'var_pos', 'gene_id', 'distance'])
    eqtl = pairs.assign(set_id=np.arange(len(pairs)), group='eqtl')[controls.columns]
    table = pd.concat([eqtl, controls], ignore_index=True)
    return table.astype({'set_id': np.int64, 'var_frag': np.int64, 'var_pos': np.int64, 'distance': np.int64})


# ---------------------------------------------------------------------------------------------------------------
# graph features


def path_features(graph, source, targets):
    """(direct, two_step, three_step, shortest) between fragment `source` and the fragment set `targets`.

    direct: an edge source-target; two_step: a path source-m-target over a mediator m; three_step: a simple path
    source-a-b-target; the three are not exclusive. shortest: '1', '2', '3', '>3' (connected) or 'none'."""
    targets = [int(t) for t in targets]
    n_source = graph.neighbour_set(source)
    direct = any(t in n_source for t in targets)
    two = three = False
    for t in targets:
        n_target = graph.neighbour_set(t)
        if not two and len(n_source) and len(n_target):
            two = bool((n_source & n_target) - {source, t})
        if not three and len(n_source) and len(n_target):
            reach = n_target - {source}
            for a in n_source:
                if a == t:
                    continue
                if not reach.isdisjoint(graph.neighbour_set(a)):
                    three = True
                    break
        if two and three:
            break
    if direct:
        shortest = '1'
    elif two:
        shortest = '2'
    elif three:
        shortest = '3'
    elif any(graph.component[source] == graph.component[t] for t in targets):
        shortest = '>3'
    else:
        shortest = 'none'
    return direct, two, three, shortest


def add_features(table, graph, promoters, fmap, bait):
    """Graph features of every pair of a table with var_frag and gene_id."""
    out = table.copy()
    values = [path_features(graph, int(f), promoters[g][2]) for f, g in zip(out['var_frag'], out['gene_id'])]
    out['direct'] = np.array([v[0] for v in values], dtype=int)
    out['two_step'] = np.array([v[1] for v in values], dtype=int)
    out['three_step'] = np.array([v[2] for v in values], dtype=int)
    out['shortest'] = [v[3] for v in values]
    out['bait_degree'] = [int(graph.degree[promoters[g][2]].max()) for g in out['gene_id']]
    out['variant_degree'] = graph.degree[out['var_frag'].to_numpy()]
    out['variant_on_bait'] = bait[out['var_frag'].to_numpy()].astype(int)
    out['fragment_length'] = fmap.length(out['var_frag'].to_numpy())
    out['distance_bin'] = distance_bin(out['distance'].to_numpy())
    return out


# ---------------------------------------------------------------------------------------------------------------
# statistics


INDICATORS = {
    'step1': lambda t: t['shortest'] == '1',
    'within2': lambda t: t['shortest'].isin(['1', '2']),
    'within3': lambda t: t['shortest'].isin(['1', '2', '3']),
    'connected': lambda t: t['shortest'] != 'none',
    'exact2': lambda t: t['shortest'] == '2',
    'exact3': lambda t: t['shortest'] == '3',
}

# odds ratio features: (indicator, subset of pairs it is evaluated on)
CONTRASTS = {
    'direct': ('direct', None),
    'within2': ('within2', None),
    'within3': ('within3', None),
    'two_step_any': ('two_step', None),
    'three_step_any': ('three_step', None),
    'two_step_without_direct': ('two_step', lambda t: t['direct'] == 0),
    'three_step_without_1_2': ('three_step', lambda t: ~t['shortest'].isin(['1', '2'])),
}


def _group_mask(table, group):
    if group == 'control':
        return table['group'] != 'eqtl'
    return table['group'] == group


def _block_weights(blocks, n_boot, seed):
    """(labels, n_boot x n_blocks matrix of how often each block is drawn) of a chromosome block bootstrap."""
    labels, resamples = block_resamples(blocks, n_boot, seed)
    weights = np.array([np.bincount(chosen, minlength=len(labels)) for chosen in resamples], dtype=float)
    return labels, weights


def _per_block_sums(table, values, labels):
    index = np.searchsorted(labels, table['chrom'].to_numpy())
    return np.bincount(index, weights=np.asarray(values, dtype=float), minlength=len(labels))


def fraction_table(table, n_boot=1000, seed=0):
    """Fraction of pairs at each graph distance per group and distance bin, with block-bootstrap intervals for the
    cumulative fractions."""
    labels, weights = _block_weights(table['chrom'].to_numpy(), n_boot, seed)
    rows = []
    for group in GROUPS:
        g = table[_group_mask(table, group)]
        for bin_label in distance_bin_labels() + ['all']:
            sub = g if bin_label == 'all' else g[g['distance_bin'] == bin_label]
            row = {'group': group, 'distance_bin': bin_label, 'pairs': len(sub)}
            n_block = _per_block_sums(sub, np.ones(len(sub)), labels)
            n_boot_total = weights @ n_block
            for step in STEP_LABELS:
                row['fraction_' + step] = float(np.mean(sub['shortest'] == step)) if len(sub) else np.nan
            for name in ('step1', 'within2', 'within3', 'connected'):
                hit = INDICATORS[name](sub).to_numpy()
                row[name] = float(hit.mean()) if len(sub) else np.nan
                with np.errstate(invalid='ignore', divide='ignore'):
                    boot = (weights @ _per_block_sums(sub, hit, labels)) / n_boot_total
                row[name + '_low'] = _percentile(boot, 2.5)
                row[name + '_high'] = _percentile(boot, 97.5)
            for name in ('two_step', 'three_step'):
                row[name + '_any'] = float(sub[name].mean()) if len(sub) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def _percentile(values, q):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if len(values) else np.nan


def _odds(a, b, c, d):
    """Odds ratio with 0.5 added to every cell if any cell is empty."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a * d) / (b * c) if b * c > 0 else np.nan


def _mantel_haenszel(cells):
    """Mantel-Haenszel odds ratio over strata; cells: array strata x 4 (a, b, c, d)."""
    n = cells.sum(axis=1)
    use = n > 0
    a, b, c, d = cells[use].T
    num = np.sum(a * d / n[use])
    den = np.sum(b * c / n[use])
    return num / den if den > 0 else np.nan


def odds_table(table, n_boot=1000, seed=0, controls=('control', 'control_gene', 'control_fragment')):
    """eQTL pairs against controls: odds ratio of each contrast overall (crude and Mantel-Haenszel over distance
    bins) and per distance bin, with chromosome block-bootstrap intervals; Fisher's exact P for the crude table."""
    labels, weights = _block_weights(table['chrom'].to_numpy(), n_boot, seed)
    bins = distance_bin_labels()
    rows = []
    for control in controls:
        sub_all = table[(table['group'] == 'eqtl') | _group_mask(table, control)]
        for name, (feature, subset) in CONTRASTS.items():
            sub = sub_all if subset is None else sub_all[subset(sub_all)]
            hit = (INDICATORS[feature](sub) if feature in INDICATORS else sub[feature] == 1).to_numpy()
            case = (sub['group'] == 'eqtl').to_numpy()
            bin_index = np.array([bins.index(b) for b in sub['distance_bin']], dtype=int)
            block_index = np.searchsorted(labels, sub['chrom'].to_numpy())
            # per block and distance bin 2x2 cells
            cells = np.zeros((len(labels), len(bins), 4))
            for k, mask in enumerate((case & hit, case & ~hit, ~case & hit, ~case & ~hit)):
                np.add.at(cells[:, :, k], (block_index[mask], bin_index[mask]), 1)
            total = cells.sum(axis=0)
            boot_cells = np.einsum('rb,bsk->rsk', weights, cells)
            for bin_label in ['all', 'all_mh'] + bins:
                if bin_label == 'all':
                    a, b, c, d = total.sum(axis=0)
                    estimate = _odds(a, b, c, d)
                    boot = np.array([_odds(*x.sum(axis=0)) for x in boot_cells])
                elif bin_label == 'all_mh':
                    a, b, c, d = total.sum(axis=0)
                    estimate = _mantel_haenszel(total)
                    boot = np.array([_mantel_haenszel(x) for x in boot_cells])
                else:
                    s = bins.index(bin_label)
                    a, b, c, d = total[s]
                    estimate = _odds(a, b, c, d)
                    boot = np.array([_odds(*x[s]) for x in boot_cells])
                p = stats.fisher_exact([[a, b], [c, d]])[1] if a + b > 0 and c + d > 0 else np.nan
                rows.append({'control': control, 'contrast': name, 'distance_bin': bin_label,
                             'eqtl_pairs': int(a + b), 'eqtl_hits': int(a), 'control_pairs': int(c + d),
                             'control_hits': int(c),
                             'eqtl_fraction': a / (a + b) if a + b else np.nan,
                             'control_fraction': c / (c + d) if c + d else np.nan,
                             'odds_ratio': estimate, 'or_low': _percentile(boot, 2.5),
                             'or_high': _percentile(boot, 97.5), 'fisher_p': p})
    return pd.DataFrame(rows)


def _design(table, terms):
    columns = {'log10_distance': np.log10(table['distance'].clip(lower=1)),
               'log1p_bait_degree': np.log1p(table['bait_degree']),
               'log1p_variant_degree': np.log1p(table['variant_degree']),
               'log10_fragment_length': np.log10(table['fragment_length'].clip(lower=1))}
    X = pd.DataFrame({t: columns[t] if t in columns else table[t] for t in terms}, index=table.index)
    return X.astype(float)


def _fit(y, X, start=None):
    import statsmodels.api as sm
    model = sm.Logit(y, sm.add_constant(X, has_constant='add'))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            return model.fit(disp=0, method='newton', maxiter=100, start_params=start)
        except Exception:  # perfect separation or singular design in a resample
            return model.fit(disp=0, method='lbfgs', maxiter=500, start_params=start)


def _safe_pvalue(fit, k):
    try:
        return float(np.asarray(fit.pvalues)[k])
    except Exception:  # no covariance (separation); the likelihood-ratio test still applies
        return np.nan


MODEL_SETS = {
    'all_pairs': (None, ['log10_distance', 'direct', 'two_step', 'three_step']),
    'all_pairs_adjusted': (None, ['log10_distance', 'direct', 'two_step', 'three_step', 'log1p_bait_degree',
                                  'log1p_variant_degree', 'log10_fragment_length']),
    'without_direct': (lambda t: t['direct'] == 0, ['log10_distance', 'two_step', 'three_step']),
    'without_direct_adjusted': (lambda t: t['direct'] == 0, ['log10_distance', 'two_step', 'three_step',
                                                              'log1p_bait_degree', 'log1p_variant_degree',
                                                              'log10_fragment_length']),
}
TESTED_TERMS = ('direct', 'two_step', 'three_step')


def logistic_models(table, n_boot=200, seed=0, model_sets=None):
    """eQTL status ~ log distance + path terms (+ degrees, fragment length), on all pairs and on pairs without a
    direct edge. Per term: odds ratio, chromosome block-bootstrap interval (model refitted per resample) and the
    likelihood-ratio test of dropping the term; terms without variation in a set are left out of its model."""
    model_sets = model_sets or MODEL_SETS
    rows = []
    for name, (subset, terms) in model_sets.items():
        sub = table if subset is None else table[subset(table)]
        y = (sub['group'] == 'eqtl').astype(float)
        if not len(sub) or y.nunique() < 2:
            continue
        X = _design(sub, terms)
        terms = [t for t in terms if X[t].nunique() > 1]
        if not terms:
            continue
        X = X[terms]
        full = _fit(y, X)
        boot = []
        for items in block_resample_items(sub['chrom'].to_numpy(), n_boot, seed):
            try:
                fit = _fit(y.iloc[items], X.iloc[items], start=full.params.to_numpy())
                boot.append(fit.params[terms].to_numpy())
            except Exception:
                continue
        boot = np.array(boot)
        for k, term in enumerate(terms):
            row = {'model': name, 'term': term, 'pairs': len(sub), 'eqtl_pairs': int(y.sum()),
                   'term_positive': int((sub[term] > 0).sum()) if term in TESTED_TERMS else np.nan,
                   'coefficient': float(full.params[term]), 'odds_ratio': float(np.exp(full.params[term])),
                   'or_low': float(np.exp(np.nanpercentile(boot[:, k], 2.5))) if len(boot) else np.nan,
                   'or_high': float(np.exp(np.nanpercentile(boot[:, k], 97.5))) if len(boot) else np.nan,
                   'wald_p': _safe_pvalue(full, list(full.params.index).index(term)),
                   'n_boot': len(boot)}
            if term in TESTED_TERMS:
                reduced = _fit(y, X.drop(columns=[term]))
                statistic = max(0.0, 2.0 * (full.llf - reduced.llf))
                row['lr_statistic'] = statistic
                row['lr_p'] = float(stats.chi2.sf(statistic, 1))
            rows.append(row)
    return pd.DataFrame(rows)


def conditional_models(table):
    """Conditional logistic regression within matched sets (eQTL pair + its controls): likelihood-ratio tests of
    the path terms given log distance, on all sets and on pairs without a direct edge."""
    from statsmodels.discrete.conditional_models import ConditionalLogit
    rows = []
    for name, (subset, terms) in (('conditional_all_pairs', MODEL_SETS['all_pairs']),
                                  ('conditional_without_direct', MODEL_SETS['without_direct'])):
        sub = table if subset is None else table[subset(table)]
        y = (sub['group'] == 'eqtl').astype(float)
        informative = sub.groupby('set_id')['group'].transform(lambda g: (g == 'eqtl').any() and (g != 'eqtl').any())
        sub, y = sub[informative], y[informative]
        X = _design(sub, terms)
        terms = [t for t in terms if X[t].nunique() > 1]
        if not terms:
            continue
        X = X[terms]

        def fit(columns):
            """(params, log-likelihood, Wald P values or None) of the conditional model on `columns`; without
            columns the log-likelihood of the null model (all parameters zero)."""
            model = ConditionalLogit(y.to_numpy(), X[columns or terms[:1]].to_numpy(), groups=sub['set_id'].to_numpy())
            if not columns:
                return np.zeros(0), float(model.loglike(np.zeros(1))), None
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                try:
                    result = model.fit(disp=0)
                    return np.asarray(result.params), float(result.llf), np.asarray(result.pvalues)
                except Exception:  # no covariance under separation: maximise the likelihood directly
                    from scipy.optimize import minimize
                    solution = minimize(lambda p: -model.loglike(p), np.zeros(len(columns)), method='BFGS')
                    return solution.x, float(-solution.fun), None
        params, llf, pvalues = fit(terms)
        for k, term in enumerate(terms):
            row = {'model': name, 'term': term, 'pairs': len(sub), 'eqtl_pairs': int(y.sum()),
                   'coefficient': float(params[k]), 'odds_ratio': float(np.exp(params[k])),
                   'wald_p': float(pvalues[k]) if pvalues is not None else np.nan}
            if term in TESTED_TERMS:
                reduced_llf = fit([t for t in terms if t != term])[1]
                statistic = max(0.0, 2.0 * (llf - reduced_llf))
                row['lr_statistic'] = statistic
                row['lr_p'] = float(stats.chi2.sf(statistic, 1))
            rows.append(row)
    return pd.DataFrame(rows)
