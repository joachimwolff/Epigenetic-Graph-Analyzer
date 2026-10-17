"""Graph domains (communities, promoter hubs, edges, cascade explanations) against Hi-C TADs.

For every graph and TAD set:
  boundaries   community boundaries = midpoints between consecutive nodes (along the chromosome) of different
               communities; fraction within `window` of a TAD boundary, against 1,000 circular shifts of all boundaries
               of a chromosome (within the node span) and against 1,000 node-matched draws (the same number of
               boundaries drawn from all gaps between consecutive nodes, which keeps node and gene density)
  communities  communities with >= 3 nodes: span (first node start - last node end) inside one TAD, against shifts
  hubs         baited promoter + its direct partners (>= 3 elements): span inside one TAD, against shifts
  edges        fraction of edges whose two ends lie in different TADs (or outside), against the same edges shifted
               along the chromosome (distance matched); P of a deficit of crossings
Cascade explanations (hops_n *_explanations_all_levels.tsv): for every differential gene and graph level (step 1-3)
the strongest node is compared with the gene promoter node: same TAD or not (pairs with both ends in a TAD),
separately for explaining and non-explaining levels, against the same pairs placed at random along the chromosome
(distance matched), and explaining against non-explaining pairs (Mann-Whitney on observed - expected per pair).

TAD sets (GRCh38): hicFindTADs 50 kb of the 4DN mcool files (18_hic_reference.sh), 4DN H1 Micro-C insulation domains;
hg19: Rao et al. 2014 GM12878 Arrowhead domains. Nested domains are flattened into segments between consecutive
boundaries (flatten_domains), so "inside one TAD" means "no TAD boundary in between" for every set. hg19 graphs are lifted to GRCh38 by node midpoint for the
GRCh38 TAD sets.

Usage: python domains_graph.py [n_shifts] [window_bp]
Output: tables/domains_boundaries.tsv, domains_spans.tsv, domains_edges.tsv, domains_explanations.tsv,
        domains_explanations_genes.tsv.gz
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402
from epigeneticgraphanalyzer.compartments import (boundary_enrichment, community_boundaries,  # noqa: E402
                                                  domain_boundaries, domain_id, edges_crossing, flatten_domains,
                                                  read_domains,
                                                  spans_within_domain)
from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402

GRAPHS = [  # name, node table, edge table, assembly, cell
    ('GM_chicago', 'graphs/GM_chicago_nodes.tsv', 'graphs/GM_chicago_edges.tsv', 'hg19', 'GM12878'),
    ('K562', 'graphs/K562_nodes.tsv', 'graphs/K562_edges.tsv', 'hg38', 'K562'),
    ('H1_jung', 'graphs/H1_nodes.tsv', 'graphs/H1_edges.tsv', 'hg19', 'H1'),
    ('hESC_H9', 'graphs/hESC_nodes.tsv', 'graphs/hESC_edges.tsv', 'hg19', 'H1'),
    ('hNPC_H9', 'graphs/hNPC_nodes.tsv', 'graphs/hNPC_edges.tsv', 'hg19', 'H1'),
]
TADS = {  # cell -> list of (name, path, assembly)
    'GM12878': [('hicFindTADs_GM12878_50kb', os.path.join(cc.REF, 'GM12878_50kb_domains.bed'), 'hg38'),
                ('rao2014_arrowhead', os.path.join(cc.RAW, 'domains', 'GM12878_rao2014_hg19',
                                                   'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz'),
                 'hg19')],
    'K562': [('hicFindTADs_K562_50kb', os.path.join(cc.REF, 'K562_50kb_domains.bed'), 'hg38')],
    'H1': [('hicFindTADs_H1_50kb', os.path.join(cc.REF, 'H1_50kb_domains.bed'), 'hg38'),
           ('4dn_H1_insulation', os.path.join(cc.RAW, 'domains', 'H1_4dn_hg38', 'H1_domains_hg38.bed'), 'hg38')],
}
CASCADES = [  # comparison, assembly, TAD cell (closest reference)
    ('hESC_vs_hNEC', 'hg19', 'H1'),
    ('siMM_vs_siRE', 'hg38', 'K562'),
    ('t821_vs_FLT3ITD', 'hg38', 'K562'),
]


class Lifter:
    def __init__(self):
        self.lift = LiftOver(os.path.join(cc.RAW, 'annotation', 'hg19ToHg38.over.chain.gz'))
        self.cache = {}

    def __call__(self, chrom, pos):
        key = (chrom, int(pos))
        if key not in self.cache:
            hit = self.lift.position(chrom, int(pos))
            self.cache[key] = hit[1] if hit is not None and hit[0] == chrom else -1
        return self.cache[key]


def parse_nodes(ids):
    p = pd.Series(ids).str.extract(r'^(chr[^:]+):(\d+)-(\d+)$')
    return p[0].values, p[1].astype(np.int64).values, p[2].astype(np.int64).values


def to_assembly(chrom, start, end, source, target, lifter):
    """Node coordinates in the TAD assembly (hg19 -> hg38 by midpoint, length kept); -1 start when unmapped."""
    if source == target:
        return start, end
    mid = (start + end) // 2
    half = (end - start) // 2
    new = np.array([lifter(c, m) for c, m in zip(chrom, mid)], dtype=np.int64)
    ok = new >= 0
    return np.where(ok, new - half, -1), np.where(ok, new + half, -1)


def node_spans(chrom, mid):
    out = {}
    for c in pd.unique(chrom):
        m = mid[chrom == c]
        out[c] = (int(m.min()), int(m.max()) + 1)
    return out


def analyse_graph(name, nodes_path, edges_path, assembly, cell, lifter, n_shifts, window):
    nodes = pd.read_csv(os.path.join(cc.RESULTS, nodes_path), sep='\t')
    edges = pd.read_csv(os.path.join(cc.RESULTS, edges_path), sep='\t', usecols=['node_a', 'node_b', 'weight'])
    rows_b, rows_s, rows_e = [], [], []
    for tad_name, tad_path, tad_asm in TADS[cell]:
        if not os.path.exists(tad_path):
            print('missing', tad_path)
            continue
        domains = flatten_domains(read_domains(tad_path))
        domains = domains[domains.chrom.isin(cc.CHROMS)].reset_index(drop=True)
        bounds = domain_boundaries(domains)
        start, end = to_assembly(nodes.chrom.values, nodes.start.values, nodes.end.values, assembly, tad_asm, lifter)
        nd = nodes.assign(start=start, end=end)
        nd = nd[(nd.start >= 0) & nd.chrom.isin(cc.CHROMS)].copy()
        nd['mid'] = (nd.start + nd.end) // 2
        spans = node_spans(nd.chrom.values, nd.mid.values)
        # community boundaries
        cb = community_boundaries(nd[['chrom', 'start', 'end', 'community']])
        query = {c: g.position.values for c, g in cb.groupby('chrom')}
        gaps = {}
        for c, g in nd.sort_values(['chrom', 'mid']).groupby('chrom'):
            gaps[c] = ((g.end.values[:-1] + g.start.values[1:]) // 2)
        res = boundary_enrichment(query, bounds, spans, window=window, n_shifts=n_shifts, null_positions=gaps)
        res.update(graph=name, tads=tad_name, window=window, n_domains=len(domains),
                   median_domain_kb=float(np.median(domains.end - domains.start) / 1000))
        rows_b.append(res)
        # community spans
        com = nd.groupby('community').agg(chrom=('chrom', 'first'), n_chrom=('chrom', 'nunique'),
                                          start=('start', 'min'), end=('end', 'max'), n=('node', 'size'))
        com = com[(com.n >= 3) & (com.n_chrom == 1)]
        r = spans_within_domain(com[['chrom', 'start', 'end']].reset_index(drop=True), domains, spans,
                                n_shifts=n_shifts)
        r.update(graph=name, tads=tad_name, unit='community', median_span_kb=float(np.median(com.end - com.start) / 1000))
        rows_s.append(r)
        # promoter hubs
        pos = nd.set_index('node')[['chrom', 'start', 'end', 'bait']]
        e = edges[edges.node_a.isin(pos.index) & edges.node_b.isin(pos.index)]
        both = pd.concat([e.rename(columns={'node_a': 'hub', 'node_b': 'member'}),
                          e.rename(columns={'node_b': 'hub', 'node_a': 'member'})])
        both = both[pos.loc[both.hub.values, 'bait'].values.astype(bool)]
        hub_rows = []
        for hub, g in both.groupby('hub'):
            members = np.concatenate([[hub], g.member.values])
            if len(members) < 3:
                continue
            p = pos.loc[members]
            if p.chrom.nunique() != 1:
                continue
            hub_rows.append((p.chrom.iat[0], int(p.start.min()), int(p.end.max())))
        hubs = pd.DataFrame(hub_rows, columns=['chrom', 'start', 'end'])
        r = spans_within_domain(hubs, domains, spans, n_shifts=n_shifts)
        r.update(graph=name, tads=tad_name, unit='promoter_hub', median_span_kb=float(np.median(hubs.end - hubs.start) / 1000))
        rows_s.append(r)
        # edges crossing TAD boundaries
        ea = pos.loc[e.node_a.values]
        eb = pos.loc[e.node_b.values]
        ed = pd.DataFrame({'chrom': ea.chrom.values, 'mid_a': ((ea.start + ea.end) // 2).values,
                           'mid_b': ((eb.start + eb.end) // 2).values})
        ed = ed[ea.chrom.values == eb.chrom.values]
        r = edges_crossing(ed, domains, spans, n_shifts=n_shifts)
        dist = np.abs(ed.mid_b - ed.mid_a)
        r.update(graph=name, tads=tad_name, median_distance_kb=float(np.median(dist) / 1000))
        rows_e.append(r)
        print(name, tad_name, rows_b[-1], rows_s[-2], rows_s[-1], rows_e[-1], sep='\n  ', flush=True)
    return rows_b, rows_s, rows_e


def analyse_cascade(comparison, assembly, cell, lifter, n_shifts, seed=0):
    path = os.path.join(cc.HOPS, '{}_explanations_all_levels.tsv'.format(comparison))
    table = pd.read_csv(path, sep='\t')
    table = table[table.hops > 0]
    rows, genes = [], []
    rng = np.random.default_rng(seed)
    for tad_name, tad_path, tad_asm in TADS[cell]:
        if tad_asm != 'hg38' or not os.path.exists(tad_path):
            continue
        domains = flatten_domains(read_domains(tad_path))
        for roles in ('declared', 'learned'):
            sub = table.dropna(subset=['strongest_node_' + roles]).copy()
            ca, sa, ea = parse_nodes(sub.node.values)
            cb, sb, eb = parse_nodes(sub['strongest_node_' + roles].values)
            sa, ea = to_assembly(ca, sa, ea, assembly, tad_asm, lifter)
            sb, eb = to_assembly(cb, sb, eb, assembly, tad_asm, lifter)
            ok = (sa >= 0) & (sb >= 0) & (ca == cb)
            sub = sub[ok]
            ma, mb, ch = ((sa + ea) // 2)[ok], ((sb + eb) // 2)[ok], ca[ok]
            same = np.zeros(len(sub), dtype=bool)
            inside = np.zeros(len(sub), dtype=bool)
            for c in pd.unique(ch):
                k = ch == c
                ia, ib = domain_id(domains, c, ma[k]), domain_id(domains, c, mb[k])
                same[k] = (ia == ib) & (ia >= 0)
                inside[k] = (ia >= 0) & (ib >= 0)
            # distance-matched expectation: each pair placed at random along its chromosome, counted where both ends
            # fall in a domain (as for the observed pairs)
            dom_span = {c: (int(g.start.min()), int(g.end.max())) for c, g in domains.groupby('chrom')}
            hits = np.zeros(len(sub))
            valid = np.zeros(len(sub))
            for s in range(n_shifts):
                for c in pd.unique(ch):
                    if c not in dom_span:
                        continue
                    k = np.flatnonzero(ch == c)
                    lo, hi = dom_span[c]
                    a = rng.integers(lo, hi, size=len(k))
                    b = a + (mb[k] - ma[k])
                    ia, ib = domain_id(domains, c, a), domain_id(domains, c, b)
                    hits[k] += (ia == ib) & (ia >= 0)
                    valid[k] += (ia >= 0) & (ib >= 0)
            with np.errstate(invalid='ignore', divide='ignore'):
                expected = hits / valid
            keep = inside & np.isfinite(expected)
            sub, same, expected, ma, mb = sub[keep], same[keep], expected[keep], ma[keep], mb[keep]
            explains = sub['explains_{}_mean'.format(roles)].values.astype(bool)
            dist = np.abs(mb - ma)
            for level, g in pd.DataFrame({'level': sub.level.values, 'explains': explains, 'same': same,
                                          'expected': expected, 'dist': dist}).groupby('level'):
                # explaining vs not explaining: excess over the distance-matched expectation, Mann-Whitney on
                # (same - expected) per pair
                ex, ne = g[g.explains], g[~g.explains]
                p_vs = mannwhitneyu(ex.same - ex.expected, ne.same - ne.expected).pvalue \
                    if len(ex) and len(ne) else np.nan
                for label, h in (('explaining', ex), ('not_explaining', ne)):
                    if len(h) == 0:
                        continue
                    # P: probability that shifted pairs are in one TAD at least as often (Poisson-binomial, normal approx.)
                    mu = h.expected.sum()
                    var = np.sum(h.expected * (1 - h.expected))
                    z = (h.same.sum() - mu) / np.sqrt(var) if var > 0 else np.nan
                    rows.append({'comparison': comparison, 'tads': tad_name, 'roles': roles, 'level': level,
                                 'class': label, 'n': len(h), 'same_tad_fraction': h.same.mean(),
                                 'expected_fraction': h.expected.mean(),
                                 'ratio': h.same.mean() / h.expected.mean() if h.expected.mean() > 0 else np.nan,
                                 'z': z, 'p_enrichment': norm.sf(z) if np.isfinite(z) else np.nan,
                                 'median_distance_kb': h.dist.median() / 1000,
                                 'p_explaining_vs_not': p_vs})
            genes.append(pd.DataFrame({'comparison': comparison, 'tads': tad_name, 'roles': roles,
                                       'node': sub.node.values, 'genes': sub.genes.values, 'level': sub.level.values,
                                       'strongest_node': sub['strongest_node_' + roles].values,
                                       'explains': explains, 'same_tad': same, 'expected_same_tad': expected,
                                       'distance': dist}))
    return rows, genes


def main():
    n_shifts = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    window = int(sys.argv[2]) if len(sys.argv) > 2 else 50000
    os.makedirs(cc.TABLES, exist_ok=True)
    lifter = Lifter()
    all_b, all_s, all_e = [], [], []
    for args in GRAPHS:
        b, s, e = analyse_graph(*args, lifter=lifter, n_shifts=n_shifts, window=window)
        all_b += b
        all_s += s
        all_e += e
    pd.DataFrame(all_b).to_csv(os.path.join(cc.TABLES, 'domains_boundaries.tsv'), sep='\t', index=False,
                               float_format='%.4g')
    pd.DataFrame(all_s).to_csv(os.path.join(cc.TABLES, 'domains_spans.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(all_e).to_csv(os.path.join(cc.TABLES, 'domains_edges.tsv'), sep='\t', index=False, float_format='%.4g')
    rows, genes = [], []
    for comparison, assembly, cell in CASCADES:
        r, g = analyse_cascade(comparison, assembly, cell, lifter, n_shifts=min(n_shifts, 200))
        rows += r
        genes += g
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(cc.TABLES, 'domains_explanations.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.concat(genes).to_csv(os.path.join(cc.TABLES, 'domains_explanations_genes.tsv.gz'), sep='\t', index=False,
                            float_format='%.4g')
    print(out.to_string())


if __name__ == '__main__':
    main()
