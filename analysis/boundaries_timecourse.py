"""Boundaries over keratinocyte differentiation (Rubin et al. 2017, hg19, 50 kb; 27_boundaries.sh after
boundaries_features.py).

Per day (D0, D3, D6): capture insulation of the day graph, CTCF (one sample) and SMC1A (two donors) sites, convergent
loops, Hi-C insulation of the merged day matrix (boundaries_hic.py; local minima with prominence 0.3 = Hi-C
boundaries). Evaluation is restricted to chr2-chr22 (the day 0 matrices lack chr1, no chrX matrices).
Boundary calls per day: the rules of boundaries_evaluate.py with the a priori thresholds (dip prominence 0.3, dip
value <= 0, cohesin-bound CTCF enrichment >= 10, window +-1 bin); no fitted model.
Changes between days (D0-D3, D3-D6, D0-D6): per bin, Spearman correlation of the change of the Hi-C insulation with
the change of the capture insulation, of the CTCF and of the SMC1A site strength (bins covered on both days and all
bins); boundary gains / losses: Hi-C: a boundary of one day without a boundary of the other day within 2 bins and an
insulation change of at least 0.2; rules: a call of one day without a call of the other day within 2 bins. Overlap of
the Hi-C changes with the changes of each rule within 1 bin against 200 circular shifts of the rule's changes per
chromosome. Expression: protein-coding genes (TSS) in the two domains flanking a changed boundary (domains of the
later day for gains, of the earlier day for losses; Hi-C boundaries and the calls of the rules dip and
dip_and_ctcf_pattern) against genes flanking only unchanged boundaries: |log2 fold change| (PyDESeq2, later vs
earlier day) by Mann-Whitney and the fraction with padj < 0.05.

Usage: python boundaries_timecourse.py <results dir> [n_shifts]
Output: <results>/boundaries/tables/timecourse_*.tsv, cache/timecourse_bins_50kb.tsv.gz
"""

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
RULES = list(bd.RULES)
CHANGE_RULES = ['dip', 'ctcf', 'ctcf_pattern', 'dip_and_ctcf', 'dip_and_ctcf_pattern', 'dip_or_ctcf']
EXPRESSION_SOURCES = ['hic', 'dip', 'dip_and_ctcf_pattern']
MIN_CHANGE = 0.2


def main():
    results = sys.argv[1]
    n_shifts = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    base = os.path.join(results, 'boundaries')
    tables = os.path.join(base, 'tables')
    rng = np.random.default_rng(0)
    data, eval_rows = {}, []
    for day in DAYS:
        t = pd.read_csv(os.path.join(base, 'cache', 'features_rubin_{}_50kb.tsv.gz'.format(day)), sep='\t')
        keep = t.chrom.isin(CHROMS).values
        sites = pd.read_csv(os.path.join(base, 'cache', 'sites_rubin_{}.tsv.gz'.format(day)), sep='\t',
                            keep_default_na=False, na_values=['nan', 'NaN'])
        loops = pd.read_csv(os.path.join(base, 'cache', 'loops_rubin_{}.tsv.gz'.format(day)), sep='\t')
        t = t[keep].reset_index(drop=True)
        codes = bd.chrom_codes(t)
        covered = t.covered.values.astype(bool)
        ev = bd.bin_evidence(t[['chrom', 'start', 'end']], sites, t[['insulation', 'covered']].astype(
            {'covered': bool}), loops)
        hic = bd.local_minima(t.hic_insulation.values, codes, PROM)
        lab = bd.bin_labels(len(t), hic)
        calls, scores = {}, {}
        for rule in RULES:
            flags, score = bd.rule_flags(ev, rule)
            calls[rule] = bd.calls_from_flags(flags, score, codes)
            scores[rule] = score
            for region, mask in (('all', None), ('covered', covered)):
                sel = np.ones(len(t), dtype=bool) if mask is None else mask
                row = {'day': day, 'rule': rule, 'region': region}
                for tol in (1, 2):
                    m = bd.match_boundaries(calls[rule], hic, codes, tol, mask)
                    row.update({'{}_tol{}'.format(k, tol): v for k, v in m.items()})
                row['auprc'] = bd.average_precision(score[sel], lab[sel])
                row['positive_rate'] = float(lab[sel].mean())
                eval_rows.append(row)
        data[day] = dict(t=t, ev=ev, codes=codes, covered=covered, hic=hic, calls=calls, scores=scores)
        print(day, 'covered', round(float(covered.mean()), 3), 'hic boundaries', len(hic), flush=True)
    ev_table = pd.DataFrame(eval_rows)
    ev_table.to_csv(os.path.join(tables, 'timecourse_eval.tsv'), sep='\t', index=False, float_format='%.4g')

    # per-bin table
    t0 = data['D0']['t'][['chrom', 'start', 'end']].copy()
    for day in DAYS:
        d = data[day]
        t0['hic_insulation_' + day] = d['t'].hic_insulation.values
        t0['capture_insulation_' + day] = d['ev'].capture_insulation.values
        t0['ctcf_strength_' + day] = d['ev'].ctcf_strength.values
        t0['ctcf_pattern_' + day] = d['ev'].ctcf_pattern.values
        t0['hic_boundary_' + day] = np.isin(np.arange(len(t0)), d['hic'])
        for rule in ('dip', 'ctcf_pattern', 'dip_and_ctcf_pattern'):
            t0['call_{}_{}'.format(rule, day)] = np.isin(np.arange(len(t0)), d['calls'][rule])
    t0.to_csv(os.path.join(base, 'cache', 'timecourse_bins_50kb.tsv.gz'), sep='\t', index=False, float_format='%.4g')
    codes = data['D0']['codes']

    def changes(b_a, b_b, v_a=None, v_b=None):
        """Boundaries of day b absent (+-2 bins) on day a; with values (Hi-C insulation, lower = stronger boundary)
        also a decrease by >= MIN_CHANGE from day a to day b."""
        lab_a = bd.bin_labels(len(codes), b_a, 2, codes)
        new = b_b[~lab_a[b_b]]
        if v_a is None:
            return new
        delta = (v_a - v_b)[new]
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
        dh = db['t'].hic_insulation.values - da['t'].hic_insulation.values
        for name in ('capture_insulation', 'ctcf_strength', 'pattern_score'):
            dv = db['ev'][name].values - da['ev'][name].values
            for region, sel in (('all', np.ones(len(dh), dtype=bool)), ('covered_both', ok)):
                m = sel & np.isfinite(dh) & np.isfinite(dv)
                rho, p = spearmanr(dh[m], dv[m])
                corr_rows.append({'pair': '{}_{}'.format(a, b), 'track': name, 'region': region,
                                  'n_bins': int(m.sum()), 'spearman_vs_hic_change': rho, 'p': p})
        ha, hb = da['t'].hic_insulation.values, db['t'].hic_insulation.values
        hic_gain = changes(da['hic'], db['hic'], ha, hb)
        hic_loss = changes(db['hic'], da['hic'], hb, ha)
        rule_changes = {}
        for name in CHANGE_RULES:
            gain = changes(da['calls'][name], db['calls'][name])
            loss = changes(db['calls'][name], da['calls'][name])
            rule_changes[name] = (gain, loss)
            for kind, hq, cq in (('gain', hic_gain, gain), ('loss', hic_loss, loss)):
                obs, null, p = shifted_overlap(hq, cq, n_shifts)
                change_rows.append({'pair': '{}_{}'.format(a, b), 'rule': name, 'change': kind,
                                    'n_hic_changes': len(hq), 'n_rule_changes': len(cq),
                                    'fraction_hic_changes_matched': obs, 'expected_shifted': null, 'p_shift': p})
        # expression in domains flanking changed Hi-C boundaries
        lfc = pd.read_csv(os.path.join(results, 'timecourse_rubin2017', 'expression', 'deseq2_{}.tsv'.format(contrast)),
                          sep='\t').drop_duplicates('gene_name').set_index('gene_name')
        g = genes[genes.gene.isin(lfc.index)].copy()
        g['abs_lfc'] = lfc.loc[g.gene, 'log2FoldChange'].abs().values
        g['de'] = (lfc.loc[g.gene, 'padj'] < 0.05).values
        sources = [('hic', hic_gain, hic_loss, da['hic'], db['hic'])] + \
            [(r, rule_changes[r][0], rule_changes[r][1], da['calls'][r], db['calls'][r])
             for r in EXPRESSION_SOURCES if r != 'hic']
        for source, gain, loss, bounds_a, bounds_b in sources:
            for kind, changed, part_bounds in (('gain', gain, bounds_b), ('loss', loss, bounds_a)):
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
    print(ev_table[ev_table.region == 'all'][['day', 'rule', 'precision_tol1', 'recall_tol1', 'f1_tol1', 'auprc',
                                              'positive_rate']].to_string())
    print(pd.DataFrame(corr_rows).to_string())
    print(pd.DataFrame(change_rows).to_string())
    print(pd.DataFrame(expr_rows).to_string())


if __name__ == '__main__':
    main()
