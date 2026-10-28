"""A/B compartments from chromatin tracks (no Hi-C) evaluated against Hi-C eigenvectors.

Per 100 kb bin the tracks are binned (mean signal) and turned into log2 enrichment over the genome-wide mean. Estimates:
  marks_pc1            first principal component of all tracks, oriented by the active-role tracks (unsupervised)
  marks_role           mean z of active tracks minus mean z of repressive tracks (declared roles)
  marks_active         mean z of the active tracks only
  marks_profile_eig    Fortin-Hansen style: per chromosome bin-bin correlation of the track profiles, eigenvector 1
  marks_pc1_smooth1    marks_pc1 smoothed over 3 bins
  single_<mark>        one mark alone (H3K27ac where present, else H3K4me1)
  gene_density, gc     sequence / annotation baselines (no chromatin)
  access_level         binned accessibility (DNase / ATAC) of the cell type, log enrichment (if the track exists)
  access_fh_eig        Fortin and Hansen 2015 style: per chromosome bin-bin correlation over all accessibility tracks
                       of the assembly (>= 3; here the DNase of GM12878, K562 and H1 in GRCh38), first eigenvector
                       oriented by gene density
  marks_plus_access_pc1  PC1 of the marks together with the accessibility track
  loco_ridge           leave-one-chromosome-out ridge regression on the reference (supervised, needs Hi-C; upper bound)
  loco_ridge_gc(_marks)  the same with GC content alone and with GC plus all tracks (added value of chromatin)
Cell-type differences (cell_differences.tsv): for each pair of GM12878, K562, H1 (hg38 ENCODE sets) the difference
of the Hi-C eigenvectors against the difference of the estimates (per-chromosome quantile-normalised pair, as in
compartment_switch.py): Pearson r, and recall / precision of the estimate's switches (sign change, |delta| >= 0.5 SD)
for the Hi-C switches (sign change of the Hi-C eigenvector, |delta| >= 0.5 SD).
References: Hi-C eigenvector (hic_eigenvectors.py: of the 3 leading eigenvectors the one phased by GC, hg38, 100 kb and
250 kb) of the matching 4DN mcool, the hicPCA first eigenvector (100 kb), the published 4DN H1 Micro-C
eigenvector (250 kb), the Rao et al. 2014 GM12878 subcompartments (hg19, A1/A2 = A). hg19 track sets (the inputs EGA
uses: ENCODE GM12878, Roadmap E003 / E008) are compared with the hg38 references after midpoint liftover of the bins.
Cross-cell rows (marks of one cell type against the Hi-C of another) show how much agreement is cell-type invariant.

Usage: python compartments_chromatin.py [bin_size]
Output: <results>/compartments_domains/tables/chromatin_eval.tsv (all metrics per set / method / reference / chromosome),
        chromatin_eval_summary.tsv (genome and median-of-chromosome rows), chromatin_loadings.tsv,
        chromatin_bins_<set>.tsv.gz (per-bin estimates and references)
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402
from epigeneticgraphanalyzer.compartments import (bin_bigwig, chromatin_pc1, compare_to_reference,  # noqa: E402
                                                  compartment_switches, log_enrichment, loco_regression, per_chrom,
                                                  profile_eigenvector, quantile_normalise_pair, role_score,
                                                  smooth_per_chrom, standardise_per_chrom)

ROADMAP = os.path.join(cc.RAW, 'roadmap')
ENC19 = os.path.join(cc.RAW, 'encode_gm12878')


def hg38_set(cell, marks=None):
    folder = cc.ENCODE_HG38[cell]
    names = sorted(f[:-len('.bigwig')] for f in os.listdir(folder) if f.endswith('.bigwig'))
    if marks is not None:
        names = [m for m in marks if m in names]
    return {m: os.path.join(folder, m + '.bigwig') for m in names}


def roadmap_set(eid, marks=('H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K27me3')):
    out = {}
    for m in marks:
        path = os.path.join(ROADMAP, '{}-{}.fc.signal.bigwig'.format(eid, m))
        if os.path.exists(path):
            out[m] = path
    return out


SETS = [
    # name, assembly, cell (for references), tracks
    ('GM12878_encode_hg38', 'hg38', 'GM12878', lambda: hg38_set('GM12878')),
    ('K562_encode_hg38', 'hg38', 'K562', lambda: hg38_set('K562')),
    ('H1_encode_hg38', 'hg38', 'H1', lambda: hg38_set('H1')),
    ('GM12878_4marks_hg38', 'hg38', 'GM12878',
     lambda: hg38_set('GM12878', ['H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K27me3'])),
    ('K562_4marks_hg38', 'hg38', 'K562', lambda: hg38_set('K562', ['H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K27me3'])),
    ('H1_4marks_hg38', 'hg38', 'H1', lambda: hg38_set('H1', ['H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K27me3'])),
    ('H1_3marks_hg38', 'hg38', 'H1', lambda: hg38_set('H1', ['H3K4me1', 'H3K4me3', 'H3K27me3'])),
    ('GM12878_ega_hg19', 'hg19', 'GM12878',
     lambda: {m: os.path.join(ENC19, 'GM12878_{}.bigWig'.format(f)) for m, f in
              [('H3K27ac', 'H3k27ac'), ('H3K4me1', 'H3k4me1'), ('H3K4me3', 'H3k4me3'), ('H3K27me3', 'H3k27me3'),
               ('CTCF', 'Ctcf')]}),
    ('E003_H1_hg19', 'hg19', 'H1', lambda: roadmap_set('E003')),
    ('E008_H9_hg19', 'hg19', 'H1', lambda: roadmap_set('E008')),
    ('E008_H9_3marks_hg19', 'hg19', 'H1', lambda: roadmap_set('E008', ('H3K4me1', 'H3K4me3', 'H3K27me3'))),
]
CELLS = ['GM12878', 'K562', 'H1']
# accessibility (DNase / ATAC) signal: <raw>/accessibility/<cell>_<assay>.<assembly>[.<accession>].bigWig, assembly
# GRCh38 / hg38 / hg19 / mm10; the key is the cell (H1 = E003, H9 = E008 for the Roadmap tracks in hg19)
ACCESS_DIR = os.path.join(cc.RAW, 'accessibility')
ACCESS_KEY = {'GM12878_ega_hg19': 'GM12878', 'E003_H1_hg19': 'H1', 'E008_H9_hg19': 'H9', 'E008_H9_3marks_hg19': 'H9'}
ASSEMBLY_ALIAS = {'GRCh38': 'hg38', 'hg38': 'hg38', 'GRCh37': 'hg19', 'hg19': 'hg19', 'mm10': 'mm10'}


def accessibility_files(assembly, directory=ACCESS_DIR):
    """dict cell -> bigWig of the accessibility tracks of one assembly."""
    out = {}
    if os.path.isdir(directory):
        for f in sorted(os.listdir(directory)):
            parts = f.split('.')
            if len(parts) < 3 or parts[-1].lower() not in ('bigwig', 'bw'):
                continue
            if ASSEMBLY_ALIAS.get(parts[1]) != assembly:
                continue
            out.setdefault(parts[0].split('_')[0], os.path.join(directory, f))
    return out


def references(cell, bins, assembly, chain):
    """dict reference name -> per-bin values (hg38 references looked up at lifted midpoints for hg19 bins)."""
    if assembly == 'hg38':
        chroms, pos = bins.chrom.values, ((bins.start.values + bins.end.values) // 2)
    else:
        chroms, pos = cc.lift_bin_midpoints(bins, chain)
    ok = pos >= 0
    out = {}

    def lookup(path):
        values = np.full(len(bins), np.nan)
        values[ok] = cc.reference_on_positions(path, chroms[ok], pos[ok])
        return values

    for other in CELLS:
        for res in ('100kb', '250kb'):
            if other != cell and res != '100kb':
                continue
            path = os.path.join(cc.REF, '{}_eig_{}.bedgraph'.format(other, res))   # hic_eigenvectors.py
            if os.path.exists(path):
                out['hic_{}_{}'.format(other, res)] = lookup(path)
    path = os.path.join(cc.REF, '{}_pc1_100kb.bedgraph'.format(cell))            # hicPCA, first eigenvector
    if os.path.exists(path):
        out['hicpca_{}_100kb'.format(cell)] = lookup(path)
    if cell == 'H1':
        out['4dn_H1_250kb'] = lookup(os.path.join(cc.RAW, 'domains', 'H1_4dn_hg38', '4DNFI475YIT8.bw'))
    if cell == 'GM12878' and assembly == 'hg19':
        sub = pd.read_csv(os.path.join(cc.RAW, 'domains', 'GM12878_rao2014_hg19', 'GSE63525_GM12878_subcompartments.bed.gz'),
                          sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'label'])
        sub['value'] = np.where(sub.label.str.startswith('A'), 1.0, np.where(sub.label.str.startswith('B'), -1.0, np.nan))
        tmp = os.path.join(cc.REF, 'GM12878_rao2014_subcompartments_hg19.bedgraph')
        sub.dropna()[['chrom', 'start', 'end', 'value']].to_csv(tmp, sep='\t', header=False, index=False)
        out['rao2014_subcompartments'] = cc.read_bedgraph_on_bins(tmp, bins)
    return out


def accessibility_estimates(est, feats, active, bins, own, panel, gene_density):
    """Accessibility level of the own track, Fortin-Hansen eigenvector of the panel of all accessibility tracks of
    the assembly (>= 3 tracks), and the PC1 of marks plus accessibility."""
    if own is not None:
        level = log_enrichment(own)
        est['access_level'] = level
        est['access_level_smooth1'] = smooth_per_chrom(level, bins.chrom.values, 1)
        both = np.column_stack([feats, level])
        est['marks_plus_access_pc1'] = chromatin_pc1(both, active=np.append(active, True))[0]
    if panel is not None and panel.shape[1] >= 3:
        est['access_fh_eig'] = per_chrom(profile_eigenvector, panel, bins, orient_by=np.log2(gene_density + 1))


def estimates(bins, tracks, assembly, gene_density, gc):
    names = list(tracks)
    raw = np.column_stack([bin_bigwig(tracks[m], bins) for m in names])
    feats = np.column_stack([log_enrichment(raw[:, j]) for j in range(raw.shape[1])])
    roles = [cc.ROLES.get(m, 'other') for m in names]
    active = np.array([r in ('active', 'enhancer', 'transcription') for r in roles])
    pc1, loadings, explained = chromatin_pc1(feats, active=active)
    est = {'marks_pc1': pc1, 'marks_role': role_score(feats, roles)}
    act_only = [r if r in ('active', 'enhancer', 'transcription') else 'other' for r in roles]
    est['marks_active'] = role_score(feats, act_only)
    est['marks_profile_eig'] = per_chrom(profile_eigenvector, feats, bins, orient_by=pc1) if feats.shape[1] >= 3 \
        else np.full(len(bins), np.nan)
    est['marks_pc1_smooth1'] = smooth_per_chrom(pc1, bins.chrom.values, 1)
    single = 'H3K27ac' if 'H3K27ac' in names else 'H3K4me1'
    est['single_' + single] = feats[:, names.index(single)]
    est['gene_density'] = np.log2(gene_density + 1)
    est['gc'] = gc
    load = pd.DataFrame({'track': names, 'role': roles, 'loading': loadings, 'explained_pc1': explained})
    return est, feats, load, active


def scale_per_chrom(values, chroms):
    out = np.full(len(values), np.nan)
    for chrom in pd.unique(chroms):
        k = chroms == chrom
        sd = np.nanstd(values[k])
        if sd > 0:
            out[k] = values[k] / sd
    return out


def cell_differences(threshold=0.5):
    """Do the estimates follow the compartment differences between cell types?"""
    frames = {}
    for cell in CELLS:
        path = os.path.join(cc.TABLES, 'chromatin_bins_{}_encode_hg38.tsv.gz'.format(cell))
        if os.path.exists(path):
            frames[cell] = pd.read_csv(path, sep='\t')
    rows = []
    for i, c1 in enumerate(CELLS):
        for c2 in CELLS[i + 1:]:
            if c1 not in frames or c2 not in frames:
                continue
            f1, f2 = frames[c1], frames[c2]
            chroms = f1.chrom.values
            h1 = np.where(f1['ref_hic_{}_100kb'.format(c1)] == 0, np.nan, f1['ref_hic_{}_100kb'.format(c1)])
            h2 = np.where(f2['ref_hic_{}_100kb'.format(c2)] == 0, np.nan, f2['ref_hic_{}_100kb'.format(c2)])
            # Hi-C: the sign of the eigenvector is meaningful; scale each chromosome by its standard deviation
            s1, s2 = scale_per_chrom(h1, chroms), scale_per_chrom(h2, chroms)
            hic_delta = s2 - s1
            hic_switch = np.isin(compartment_switches(s1, s2, threshold), ['B->A', 'A->B'])
            for method in ('marks_pc1', 'marks_pc1_smooth1', 'marks_active', 'single_H3K27ac', 'access_level',
                           'access_level_smooth1', 'marks_plus_access_pc1'):
                if method not in f1 or method not in f2:
                    continue
                e1, e2 = quantile_normalise_pair(f1[method].values, f2[method].values, chroms)
                e1, e2 = standardise_per_chrom(e1, chroms), standardise_per_chrom(e2, chroms)
                delta = e2 - e1
                ok = np.isfinite(delta) & np.isfinite(hic_delta)
                est_switch = np.isin(compartment_switches(e1, e2, threshold), ['B->A', 'A->B'])
                both = (est_switch & hic_switch & ok)
                same_dir = both & (np.sign(delta) == np.sign(hic_delta))
                rows.append({'cell_A': c1, 'cell_B': c2, 'method': method, 'n_bins': int(ok.sum()),
                             'pearson_delta': float(np.corrcoef(delta[ok], hic_delta[ok])[0, 1]),
                             'hic_switch_bins': int((hic_switch & ok).sum()),
                             'estimate_switch_bins': int((est_switch & ok).sum()),
                             'recall': float(same_dir.sum() / max((hic_switch & ok).sum(), 1)),
                             'precision': float(same_dir.sum() / max((est_switch & ok).sum(), 1)),
                             'expected_precision': float((hic_switch & ok).sum() / ok.sum() / 2)})
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(cc.TABLES, 'cell_differences.tsv'), sep='\t', index=False, float_format='%.4f')
    print(table.round(3).to_string())


def main():
    bin_size = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    os.makedirs(cc.TABLES, exist_ok=True)
    chain = os.path.join(cc.RAW, 'annotation', 'hg19ToHg38.over.chain.gz')
    long_rows, loads = [], []
    base = {}
    access = {}
    for assembly in ('hg38', 'hg19'):
        bins = cc.bins_for(assembly, bin_size)
        base[assembly] = (bins, cc.gene_density(bins, assembly), cc.gc_content(bins, assembly))
        files = accessibility_files(assembly)
        binned = {k: bin_bigwig(v, bins) for k, v in files.items()}
        panel = np.column_stack([log_enrichment(v) for v in binned.values()]) if len(binned) >= 3 else None
        access[assembly] = (binned, panel)
        print(assembly, 'accessibility tracks:', list(files), flush=True)
    for name, assembly, cell, make in SETS:
        tracks = make()
        if len(tracks) < 2:
            print('skip', name, 'tracks', tracks)
            continue
        print(name, list(tracks), flush=True)
        bins, genes, gc = base[assembly]
        est, feats, load, active = estimates(bins, tracks, assembly, genes, gc)
        binned, panel = access[assembly]
        key = ACCESS_KEY.get(name, cell)
        accessibility_estimates(est, feats, active, bins, binned.get(key), panel, genes)
        load['set'] = name
        loads.append(load)
        refs = references(cell, bins, assembly, chain)
        own = 'hic_{}_100kb'.format(cell)
        if own in refs:
            est['loco_ridge'] = loco_regression(feats, refs[own], bins.chrom.values)
            # does chromatin add to sequence? GC alone and GC + tracks, same supervised model
            est['loco_ridge_gc'] = loco_regression(gc[:, None], refs[own], bins.chrom.values)
            est['loco_ridge_gc_marks'] = loco_regression(np.column_stack([feats, gc]), refs[own], bins.chrom.values)
        frame = bins.copy()
        for k, v in est.items():
            frame[k] = v
        for k, v in refs.items():
            frame['ref_' + k] = v
        frame.to_csv(os.path.join(cc.TABLES, 'chromatin_bins_{}.tsv.gz'.format(name)), sep='\t', index=False,
                     float_format='%.5g')
        for ref_name, ref in refs.items():
            for method, values in est.items():
                if method.startswith('loco_ridge') and ref_name != own:
                    continue
                threshold = 'zero' if method in ('marks_profile_eig', 'access_fh_eig') or method.startswith('loco_ridge') \
                    else 'median'
                table = compare_to_reference(values, ref, bins.chrom.values, threshold=threshold)
                table.insert(0, 'reference', ref_name)
                table.insert(0, 'method', method)
                table.insert(0, 'cell', cell)
                table.insert(0, 'set', name)
                table['cross_cell'] = ref_name.startswith('hic_') and not ref_name.startswith('hic_' + cell + '_')
                long_rows.append(table)
    result = pd.concat(long_rows, ignore_index=True)
    result.to_csv(os.path.join(cc.TABLES, 'chromatin_eval.tsv'), sep='\t', index=False, float_format='%.4f')
    summary = result[result.chrom.isin(['genome', 'median_of_chroms'])]
    summary.to_csv(os.path.join(cc.TABLES, 'chromatin_eval_summary.tsv'), sep='\t', index=False, float_format='%.4f')
    pd.concat(loads).to_csv(os.path.join(cc.TABLES, 'chromatin_loadings.tsv'), sep='\t', index=False,
                            float_format='%.4f')
    cell_differences()
    show = summary[(summary.chrom == 'genome') & ~summary.cross_cell]
    print(show.pivot_table(index=['set', 'method'], columns='reference', values='pearson').round(3).to_string())


if __name__ == '__main__':
    main()
