"""Boundaries over keratinocyte differentiation (Rubin et al. 2017, hg19, 50 kb; 27_boundaries.sh after
boundaries_evaluate.py).

Per day (D0, D3, D6): capture insulation of the day graph, CTCF (one sample) and SMC1A (two donors) sites, Hi-C
insulation of the merged day matrix (boundaries_hic.py; local minima with prominence 0.3 = Hi-C boundaries).
Evaluation is restricted to chr2-chr22 (the day 0 matrices lack chr1, no chrX matrices).
Boundary calls per day:
  capture_minima   local minima of the capture insulation (prominence 0.3)
  <method>_encode  the three ENCODE cross-cell models (GM12878, K562, H1, 50 kb, boundaries_evaluate.py) applied to
                   the day features (z-scored within the day); mean probability, mean threshold
  <method>_split   within-day chromosome split (odd / even) trained on the day's own Hi-C boundaries
Changes between days (D0-D3, D3-D6, D0-D6): per bin, Spearman correlation of the change of the Hi-C insulation with
the change of the capture insulation (bins covered on both days), of the CTCF and of the SMC1A site strength and of the
combined score; boundary gains / losses: a boundary of one day without a boundary of the other day within 2 bins and
an insulation change of at least 0.2 (Hi-C: Hi-C insulation; calls: the score of that caller), overlap of the
Hi-C changes with the changes of each caller within 1 bin against 200 circular shifts of the caller's changes per
chromosome. Expression: protein-coding genes (TSS) in the two domains flanking a changed Hi-C boundary (domains of the
later day for gains, of the earlier day for losses) against genes flanking only unchanged boundaries: |log2 fold
change| (PyDESeq2, later vs earlier day) by Mann-Whitney and the fraction with padj < 0.05.

Usage: python boundaries_timecourse.py <results dir> [n_shifts]
Output: <results>/boundaries/tables/timecourse_*.tsv, cache/timecourse_bins_50kb.tsv.gz
"""

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402

DAYS = ['D0', 'D3', 'D6']
PAIRS = [('D0', 'D3', 'D3_vs_D0'), ('D3', 'D6', 'D6_vs_D3'), ('D0', 'D6', 'D6_vs_D0')]
CHROMS = ['chr{}'.format(c) for c in range(2, 23)]
RES = 50000
PROM = 0.3
METHODS = ['density', 'capture', 'ctcf_cohesin', 'ctcf_orientation', 'ctcf_orientation_density', 'combined']
MIN_CHANGE = 0.2


def main():
    results = sys.argv[1]
    n_shifts = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    base = os.path.join(results, 'boundaries')
    tables = os.path.join(base, 'tables')
    rng = np.random.default_rng(0)
    models = {}
    for method in METHODS:
        models[method] = []
        for cell in ('GM12878', 'K562', 'H1'):
            with open(os.path.join(base, 'models', '{}_50kb_{}.json'.format(cell, method))) as h:
                spec = json.load(h)
            models[method].append((bd.LogisticModel.from_dict(spec), spec['threshold']))
    radius = 1
    data, eval_rows = {}, []
    for day in DAYS:
        t = pd.read_csv(os.path.join(base, 'cache', 'features_rubin_{}_50kb.tsv.gz'.format(day)), sep='\t')
        t = t[t.chrom.isin(CHROMS)].reset_index(drop=True)
        codes = bd.chrom_codes(t)
        covered = t.covered.values > 0
        hic = bd.local_minima(t.hic_insulation.values, codes, PROM)
        calls, scores = {}, {}
        ins = t.capture_insulation.values
        calls['capture_minima'] = bd.local_minima(ins, codes, PROM)
        scores['capture_minima'] = np.where(np.isfinite(ins), -ins, np.nanmin(-ins) - 1)
        for method in METHODS:
            s = np.mean([m.predict(t) for m, _ in models[method]], axis=0)
            thr = float(np.mean([th for _, th in models[method]]))
            scores[method + '_encode'] = s
            calls[method + '_encode'] = bd.peak_calls(s, codes, thr, radius)
            label = bd.bin_labels(len(t), hic)
            num = t.chrom.str.replace('chr', '').astype(int).values
            fold = num % 2
            sc = np.zeros(len(t))
            cl = []
            for f in (0, 1):
                tr, te = fold != f, fold == f
                m = bd.LogisticModel(bd.METHODS[method]).fit(t[tr], label[tr])
                p = m.predict(t)
                th, _ = bd.best_threshold(p[tr], np.flatnonzero(label[tr]), codes[tr], radius=radius)
                sc[te] = p[te]
                cl.append(bd.peak_calls(np.where(te, p, -np.inf), codes, th, radius))
            scores[method + '_split'] = sc
            calls[method + '_split'] = np.sort(np.concatenate(cl))
        for name in calls:
            for region, mask in (('all', None), ('covered', covered)):
                sel = np.ones(len(t), dtype=bool) if mask is None else mask
                row = {'day': day, 'method': name, 'region': region}
                for tol in (1, 2):
                    m = bd.match_boundaries(calls[name], hic, codes, tol, mask)
                    row.update({'{}_tol{}'.format(k, tol): v for k, v in m.items()})
                lab = bd.bin_labels(len(t), hic)
                row['auprc'] = bd.average_precision(scores[name][sel], lab[sel])
                row['positive_rate'] = float(lab[sel].mean())
                eval_rows.append(row)
        data[day] = dict(t=t, codes=codes, covered=covered, hic=hic, calls=calls, scores=scores)
        print(day, 'covered', round(float(covered.mean()), 3), 'hic boundaries', len(hic), flush=True)
    ev = pd.DataFrame(eval_rows)
    ev.to_csv(os.path.join(tables, 'timecourse_eval.tsv'), sep='\t', index=False, float_format='%.4g')

    # per-bin table
    t0 = data['D0']['t'][['chrom', 'start', 'end']].copy()
    for day in DAYS:
        d = data[day]
        t0['hic_insulation_' + day] = d['t'].hic_insulation.values
        t0['capture_insulation_' + day] = d['t'].capture_insulation.values
        t0['ctcf_' + day] = d['t'].ctcf_max.values
        t0['cohesin_' + day] = d['t'].cohesin_max.values
        t0['combined_' + day] = d['scores']['combined_encode']
        t0['hic_boundary_' + day] = np.isin(np.arange(len(t0)), d['hic'])
        t0['combined_call_' + day] = np.isin(np.arange(len(t0)), d['calls']['combined_encode'])
        t0['capture_call_' + day] = np.isin(np.arange(len(t0)), d['calls']['capture_minima'])
    t0.to_csv(os.path.join(base, 'cache', 'timecourse_bins_50kb.tsv.gz'), sep='\t', index=False, float_format='%.4g')
    codes = data['D0']['codes']

    def changes(b_a, b_b, v_a, v_b, sign):
        """Boundaries of day b absent (+-2 bins) on day a with a deeper value on b by >= MIN_CHANGE (sign=-1 for
        insulation, where lower = stronger boundary; +1 for scores)."""
        lab_a = bd.bin_labels(len(codes), b_a, 2, codes)
        new = b_b[~lab_a[b_b]]
        delta = (v_b - v_a)[new] * sign
        return new[np.nan_to_num(delta, nan=0) >= MIN_CHANGE]

    def shifted_overlap(query, target, n):
        obs = bd.match_boundaries(target, query, codes, 1)['recall'] if len(query) and len(target) else np.nan
        null = []
        for _ in range(n):
            shifted = []
            for c in np.unique(codes):
                idx = np.flatnonzero(codes == c)
                k = target[codes[target] == c] - idx[0]
                if len(k) == 0:
                    continue
                shifted.append(idx[0] + (k + rng.integers(len(idx))) % len(idx))
            sh = np.concatenate(shifted) if shifted else np.array([], dtype=np.int64)
            null.append(bd.match_boundaries(sh, query, codes, 1)['recall'] if len(query) and len(sh) else np.nan)
        null = np.array(null)
        return obs, float(np.nanmean(null)), float(np.mean(null >= obs)) if np.isfinite(obs) else np.nan

    genes = pd.read_csv(os.path.join(results, 'compartments_domains', 'hic_reference', 'genes_hg19_protein_coding.bed'),
                        sep='\t', header=None, names=['chrom', 'start', 'end', 'gene', 'score', 'strand'])
    genes['tss'] = np.where(genes.strand == '+', genes.start, genes.end - 1)
    genes = genes[genes.chrom.isin(CHROMS)].drop_duplicates('gene')
    first = {c: int(np.flatnonzero(t0.chrom.values == c)[0]) for c in CHROMS}
    genes['bin'] = [first[c] + int(p // RES) for c, p in zip(genes.chrom, genes.tss)]
    change_rows, expr_rows, corr_rows = [], [], []
    for a, b, contrast in PAIRS:
        da, db = data[a], data[b]
        ok = da['covered'] & db['covered']
        for name, va, vb in (('capture_insulation', da['t'].capture_insulation.values,
                              db['t'].capture_insulation.values),
                             ('ctcf', da['t'].ctcf_max.values, db['t'].ctcf_max.values),
                             ('cohesin', da['t'].cohesin_max.values, db['t'].cohesin_max.values),
                             ('combined_encode', da['scores']['combined_encode'], db['scores']['combined_encode'])):
            dh = db['t'].hic_insulation.values - da['t'].hic_insulation.values
            dv = vb - va
            for region, sel in (('all', np.ones(len(dh), dtype=bool)), ('covered_both', ok)):
                m = sel & np.isfinite(dh) & np.isfinite(dv)
                rho, p = spearmanr(dh[m], dv[m])
                corr_rows.append({'pair': '{}_{}'.format(a, b), 'track': name, 'region': region, 'n_bins': int(m.sum()),
                                  'spearman_vs_hic_change': rho, 'p': p})
        ha, hb = da['t'].hic_insulation.values, db['t'].hic_insulation.values
        hic_gain = changes(da['hic'], db['hic'], ha, hb, -1)
        hic_loss = changes(db['hic'], da['hic'], hb, ha, -1)
        for name in ('capture_minima', 'combined_encode', 'ctcf_orientation_encode', 'ctcf_cohesin_encode'):
            if name == 'capture_minima':
                va, vb, sign = da['t'].capture_insulation.values, db['t'].capture_insulation.values, -1
            else:
                va, vb, sign = da['scores'][name], db['scores'][name], 1
                sd = np.nanstd(np.concatenate([va, vb]))
                va, vb = va / sd * 0.5, vb / sd * 0.5      # score changes in units of 2 standard deviations
            gain = changes(da['calls'][name], db['calls'][name], va, vb, sign)
            loss = changes(db['calls'][name], da['calls'][name], vb, va, sign)
            for kind, hq, cq in (('gain', hic_gain, gain), ('loss', hic_loss, loss)):
                obs, null, p = shifted_overlap(hq, cq, n_shifts)
                change_rows.append({'pair': '{}_{}'.format(a, b), 'caller': name, 'change': kind,
                                    'n_hic_changes': len(hq), 'n_caller_changes': len(cq),
                                    'fraction_hic_changes_matched': obs, 'expected_shifted': null, 'p_shift': p})
        # expression in domains flanking changed Hi-C boundaries
        lfc = pd.read_csv(os.path.join(results, 'timecourse_rubin2017', 'expression', 'deseq2_{}.tsv'.format(contrast)),
                          sep='\t').drop_duplicates('gene_name').set_index('gene_name')
        g = genes[genes.gene.isin(lfc.index)].copy()
        g['abs_lfc'] = lfc.loc[g.gene, 'log2FoldChange'].abs().values
        g['de'] = (lfc.loc[g.gene, 'padj'] < 0.05).values
        for source, gain, loss in (('hic', hic_gain, hic_loss),):
            for kind, changed, part_bounds in (('gain', gain, db['hic']), ('loss', loss, da['hic'])):
                lab = bd.segment_labels(codes, part_bounds)
                changed_dom = set()
                for k in changed:
                    changed_dom.update([lab[k] - 1, lab[k]])   # the domains left and right of the boundary
                stable_dom = set()
                stable = np.setdiff1d(part_bounds, changed)
                for k in stable:
                    stable_dom.update([lab[k] - 1, lab[k]])
                stable_dom -= changed_dom
                gl = lab[g.bin.values]
                in_c = np.isin(gl, list(changed_dom))
                in_s = np.isin(gl, list(stable_dom))
                x, y = g.abs_lfc.values[in_c], g.abs_lfc.values[in_s]
                p = mannwhitneyu(x, y, alternative='greater').pvalue if len(x) and len(y) else np.nan
                expr_rows.append({'pair': '{}_{}'.format(a, b), 'contrast': contrast, 'boundaries': source,
                                  'change': kind, 'n_changed_boundaries': len(changed), 'n_genes_changed': len(x),
                                  'n_genes_stable': len(y),
                                  'median_abs_lfc_changed': np.median(x) if len(x) else np.nan,
                                  'median_abs_lfc_stable': np.median(y) if len(y) else np.nan,
                                  'de_fraction_changed': g.de.values[in_c].mean() if len(x) else np.nan,
                                  'de_fraction_stable': g.de.values[in_s].mean() if len(y) else np.nan,
                                  'p_mannwhitney_greater': p})
    pd.DataFrame(corr_rows).to_csv(os.path.join(tables, 'timecourse_change_correlation.tsv'), sep='\t', index=False,
                                   float_format='%.4g')
    pd.DataFrame(change_rows).to_csv(os.path.join(tables, 'timecourse_changes.tsv'), sep='\t', index=False,
                                     float_format='%.4g')
    pd.DataFrame(expr_rows).to_csv(os.path.join(tables, 'timecourse_expression.tsv'), sep='\t', index=False,
                                   float_format='%.4g')
    pd.set_option('display.width', 250)
    print(ev[ev.region == 'all'][['day', 'method', 'precision_tol1', 'recall_tol1', 'f1_tol1', 'auprc',
                                  'positive_rate']].to_string())
    print(pd.DataFrame(corr_rows).to_string())
    print(pd.DataFrame(change_rows).to_string())
    print(pd.DataFrame(expr_rows).to_string())


if __name__ == '__main__':
    main()
