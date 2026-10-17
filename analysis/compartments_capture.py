"""A/B compartments from promoter capture Hi-C interactions (the edges of the EGA graphs) against Hi-C.

Per chromosome the graph edges are binned (250 kb, 500 kb, 1 Mb) into a symmetric matrix of summed interaction
scores (edge weight) or of edge counts; observed / expected, Pearson correlation and the first eigenvector
(contact_eigenvector; bins with fewer than 3 non-empty entries are left out); as for the Hi-C reference
(hic_eigenvectors.py) the one of the 3 leading eigenvectors that correlates best with GC content, oriented by GC
(capture_eig_score / _count); capture_eig1_score is the first eigenvector only. A simpler
estimate is the contact density (log2 of the summed edge weight per bin). Combined estimate: mean of the
per-chromosome standardised chromatin PC1 (ENCODE 4 marks, hg38, compartments_chromatin.py) and the standardised
contact eigenvector where it exists (chromatin alone elsewhere).
Graphs in hg19 (GM12878 CHiCAGO, GM12878 Mifsud read counts, H1 Jung 2019, hESC H9) are lifted to hg38 by node
midpoint; all evaluation is on GRCh38 against the Hi-C eigenvector of the matching 4DN mcool (hic_eigenvectors.py, 100 kb, averaged
within the coarser bins) and, for H1, the published 4DN eigenvector.
Coverage = fraction of bins with a reference value that also have an estimate.

Usage: python compartments_capture.py   (after compartments_chromatin.py)
Output: tables/capture_eval.tsv, capture_eval_summary.tsv, capture_bins_<graph>_<res>.tsv.gz
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402
from epigeneticgraphanalyzer.compartments import (compare_to_reference, contact_eigenvector,  # noqa: E402
                                                  contacts_to_matrix, standardise_per_chrom)
from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402

GRAPHS = [
    # name, cell, edge table, assembly
    ('GM_chicago', 'GM12878', 'graphs/GM_chicago_edges.tsv', 'hg19'),
    ('GM_mifsud_counts', 'GM12878', 'graphs/GM_mifsud_c20_r1_edges.tsv', 'hg19'),
    ('K562', 'K562', 'graphs/K562_edges.tsv', 'hg38'),
    ('H1_jung', 'H1', 'graphs/H1_edges.tsv', 'hg19'),
    ('hESC_H9', 'H1', 'graphs/hESC_edges.tsv', 'hg19'),
]
CHROMATIN_SET = {'GM12878': 'GM12878_4marks_hg38', 'K562': 'K562_4marks_hg38', 'H1': 'H1_4marks_hg38'}
RESOLUTIONS = (250000, 500000, 1000000)


def edge_midpoints(path, assembly, lift):
    edges = pd.read_csv(path, sep='\t', usecols=['node_a', 'node_b', 'weight'])
    nodes = pd.unique(np.concatenate([edges.node_a.values, edges.node_b.values]))
    parsed = pd.Series(nodes).str.extract(r'^(chr[^:]+):(\d+)-(\d+)$')
    mid = (parsed[1].astype(np.int64) + parsed[2].astype(np.int64)) // 2
    chrom = parsed[0].values
    if assembly == 'hg19':
        lifted_c = np.empty(len(nodes), dtype=object)
        lifted_p = np.full(len(nodes), -1, dtype=np.int64)
        for i, (c, m) in enumerate(zip(chrom, mid.values)):
            hit = lift.position(c, int(m))
            if hit is not None and hit[0] == c:
                lifted_c[i], lifted_p[i] = hit
        chrom, pos = lifted_c, lifted_p
    else:
        pos = mid.values
    table = pd.DataFrame({'node': nodes, 'chrom': chrom, 'pos': pos}).set_index('node')
    a = table.loc[edges.node_a.values]
    b = table.loc[edges.node_b.values]
    out = pd.DataFrame({'chrom': a.chrom.values, 'pos_a': a.pos.values, 'chrom_b': b.chrom.values,
                        'pos_b': b.pos.values, 'weight': edges.weight.values})
    out = out[(out.pos_a >= 0) & (out.pos_b >= 0) & (out.chrom == out.chrom_b)]
    return out


def reference_coarse(cell, bins, fine_bins, fine_values):
    """Mean of the 100 kb reference (zeros masked) within each coarse bin."""
    size = int((bins.end - bins.start).max())
    fine = fine_bins.assign(v=np.where(fine_values == 0, np.nan, fine_values))
    fine['coarse'] = fine.start // size
    means = fine.groupby(['chrom', 'coarse'])['v'].mean()
    key = list(zip(bins.chrom.values, bins.start.values // size))
    return np.array([means.get(k, np.nan) for k in key])


def main():
    lift = LiftOver(os.path.join(cc.RAW, 'annotation', 'hg19ToHg38.over.chain.gz'))
    fine_bins = cc.bins_for('hg38', 100000)
    rows = []
    for name, cell, rel, assembly in GRAPHS:
        edges = edge_midpoints(os.path.join(cc.RESULTS, rel), assembly, lift)
        print(name, len(edges), 'edges', flush=True)
        chrom_frame = pd.read_csv(os.path.join(cc.TABLES, 'chromatin_bins_{}.tsv.gz'.format(CHROMATIN_SET[cell])),
                                  sep='\t', usecols=['chrom', 'start', 'marks_pc1'])
        refs_fine = {'hic_{}'.format(cell): cc.read_bedgraph_on_bins(
            os.path.join(cc.REF, '{}_eig_100kb.bedgraph'.format(cell)), fine_bins)}
        if cell == 'H1':
            refs_fine['4dn_H1'] = cc.read_bedgraph_on_bins(
                os.path.join(cc.RAW, 'domains', 'H1_4dn_hg38', '4DNFI475YIT8.bw'), fine_bins)
        for res in RESOLUTIONS:
            bins = cc.bins_for('hg38', res)
            genes = np.log2(cc.gene_density(bins, 'hg38') + 1)
            gc = cc.gc_content(bins, 'hg38')
            est = {k: np.full(len(bins), np.nan) for k in ('capture_eig_score', 'capture_eig_count', 'capture_eig1_score')}
            density = np.zeros(len(bins))
            for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
                sub = edges[edges.chrom == chrom]
                if sub.empty:
                    continue
                ba = (sub.pos_a.values // res).astype(np.int64)
                bb = (sub.pos_b.values // res).astype(np.int64)
                n = len(idx)
                keep = (ba < n) & (bb < n)
                ba, bb, w = ba[keep], bb[keep], sub.weight.values[keep]
                for key, values, n_cand in (('capture_eig_score', w, 3), ('capture_eig_count', np.ones(len(w)), 3),
                                             ('capture_eig1_score', w, 1)):
                    m = contacts_to_matrix(n, ba, bb, values)
                    vec, _ = contact_eigenvector(m, orient_by=gc[idx], min_contacts=3, n_candidates=n_cand)
                    est[key][idx] = vec
                np.add.at(density, idx[0] + ba, w)
                np.add.at(density, idx[0] + bb, w)
            est['capture_density'] = np.where(density > 0, np.log2(density + 1), np.nan)
            chrom_coarse = chrom_frame.assign(c=chrom_frame.start // res).groupby(['chrom', 'c'])['marks_pc1'].mean()
            chromatin = np.array([chrom_coarse.get(k, np.nan) for k in zip(bins.chrom.values, bins.start.values // res)])
            z_chrom = standardise_per_chrom(chromatin, bins.chrom.values)
            z_cap = standardise_per_chrom(est['capture_eig_score'], bins.chrom.values)
            est['chromatin_pc1'] = chromatin
            est['combined'] = np.where(np.isfinite(z_cap), (z_chrom + z_cap) / 2, z_chrom)
            est['gene_density'] = genes
            frame = bins.copy()
            for k, v in est.items():
                frame[k] = v
            for ref_name, fine_values in refs_fine.items():
                ref = reference_coarse(cell, bins, fine_bins, fine_values)
                frame['ref_' + ref_name] = ref
                has_ref = np.isfinite(ref) & (ref != 0)
                for method, values in est.items():
                    threshold = 'zero' if method.startswith('capture_eig') else 'median'
                    covered = np.isfinite(values) & has_ref
                    for scope, mask in (('covered_bins', np.isfinite(est['capture_eig_score'])),
                                        ('all_bins', np.ones(len(bins), dtype=bool))):
                        table = compare_to_reference(np.where(mask, values, np.nan), ref, bins.chrom.values,
                                                     threshold=threshold, min_bins=5)
                        table.insert(0, 'scope', scope)
                        table.insert(0, 'coverage', covered.sum() / max(has_ref.sum(), 1))
                        table.insert(0, 'reference', ref_name)
                        table.insert(0, 'method', method)
                        table.insert(0, 'resolution', res)
                        table.insert(0, 'cell', cell)
                        table.insert(0, 'graph', name)
                        rows.append(table)
            frame.to_csv(os.path.join(cc.TABLES, 'capture_bins_{}_{}kb.tsv.gz'.format(name, res // 1000)), sep='\t',
                         index=False, float_format='%.5g')
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(os.path.join(cc.TABLES, 'capture_eval.tsv'), sep='\t', index=False, float_format='%.4f')
    summary = result[result.chrom.isin(['genome', 'median_of_chroms'])]
    summary.to_csv(os.path.join(cc.TABLES, 'capture_eval_summary.tsv'), sep='\t', index=False, float_format='%.4f')
    show = summary[(summary.chrom == 'genome') & (summary.scope == 'covered_bins') & summary.reference.str.startswith('hic')]
    print(show.pivot_table(index=['graph', 'resolution'], columns='method', values='pearson').round(3).to_string())
    print(show.pivot_table(index=['graph', 'resolution'], columns='method', values='coverage').round(3).to_string())


if __name__ == '__main__':
    main()
