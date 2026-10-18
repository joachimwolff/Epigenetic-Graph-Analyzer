"""Depth-aware differential test of the Hoxa2 contacts into HIRE1, Md versus PA2 (Kessler et al. 2023, mm10).

usage: hoxa2_depth_test.py <Md 25 kb cool> <PA2 25 kb cool> <mm10 RefSeq GTF> <capture washU dir> <output dir>

Hi-C (GSE211901, one biological replicate per tissue; raw integer counts of the .cool converted from the .hic):
  virtual 4C from the 25 kb bin holding the Hoxa2 promoter (chr6:52,150,000-52,175,000). Counts into HIRE1
  (chr6:50,925,000-51,175,000, 10 bins; and the Kessler HIRE1 coordinates 50,913,170-51,087,888, 8 bins) and per
  25 kb bin across 50.9-51.2 Mb. Conditional on n = c_Md + c_PA2, c_PA2 is tested against the PA2 share expected from
  the viewpoint cis totals (cis counts at >= 2 bins distance) and from a distance-matched background: windows of the
  same bin offsets from every other anchor bin (both directions; all valid anchors and RefSeq NM_ promoter anchors).
  The background fixes the residual PA2 / Md scaling at that distance (beta0, all windows) and the overdispersion
  (rho, the 500 windows closest in n) of a beta-binomial model; P values: binomial (no overdispersion,
  anti-conservative), beta-binomial, a robust z against the median / MAD of the log2 enrichment of the 500
  count-matched windows, and the empirical rank among the background windows (all, and n within 0.5-2x).
  log2 ratio = log2[(c_PA2 / c_Md) / (expected PA2 / Md odds)], CI from the logit Wald interval inflated by the
  beta-binomial design effect 1 + (n - 1) rho (binomial CI: Clopper-Pearson; robust CI: +- 1.96 MAD-SD).
Capture Hi-C (GSE211902 CHiCAGO washU text, 2 replicates per tissue): only interactions with score >= 5 are
  published, read counts are not; per replicate the number of Hoxa2-bait calls into HIRE1 and genome-wide calls.
Outputs: hic_windows.tsv, hic_bins.tsv, hic_virtual4c.tsv, hic_depth.tsv, capture_replicates.tsv,
  virtual4c_Md_PA2.pdf, summary.json.
"""

import glob
import gzip
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import optimize, special, stats

BIN = 25000
CHROM = '6'
VIEW = 52150000                      # bin holding the Hoxa2 TSS (chr6:52,164,831, minus strand)
BAIT = 'chr6,52162863,52166906'      # Hoxa2 capture fragment
WINDOWS = {'HIRE1_partners': (50925000, 51175000), 'HIRE1_Kessler': (50900000, 51100000)}
SCAN = (50900000, 51200000)
MIN_SEP = 2                          # viewpoint cis totals from distances >= 2 bins
AUTOSOMES_X = [str(i) for i in range(1, 20)] + ['X']


def load(path):
    with h5py.File(path, 'r') as f:
        names = [x.decode() for x in f['chroms/name'][:]]
        chrom = np.asarray(f['bins/chrom'][:])
        start = np.asarray(f['bins/start'][:])
        weight = np.asarray(f['bins/weight'][:])
        b1 = np.asarray(f['pixels/bin1_id'][:])
        b2 = np.asarray(f['pixels/bin2_id'][:])
        cnt = np.asarray(f['pixels/count'][:]).astype(np.int64)
    cis = chrom[b1] == chrom[b2]
    d = b2 - b1
    depth = {'total': int(cnt.sum()), 'cis': int(cnt[cis & (d > 0)].sum()),
             'cis_ge2bins': int(cnt[cis & (d >= MIN_SEP)].sum())}
    keep = cis & (d >= MIN_SEP)
    rowcis = np.bincount(b1[keep], cnt[keep], minlength=len(chrom)) + \
        np.bincount(b2[keep], cnt[keep], minlength=len(chrom))
    return dict(names=names, chrom=chrom, start=start, weight=weight, b1=b1[cis], d=d[cis], cnt=cnt[cis],
                rowcis=rowcis, depth=depth)


def band(m, dmax):
    """band[i, k] = count(bin i, bin i + k), k <= dmax (cis, upper triangle)."""
    nb = len(m['chrom'])
    sel = m['d'] <= dmax
    flat = np.bincount(m['b1'][sel] * (dmax + 1) + m['d'][sel], m['cnt'][sel], minlength=nb * (dmax + 1))
    return flat.reshape(nb, dmax + 1).astype(np.int64)


def window_counts(bands, valid, chrom, anchors, offsets):
    """Counts of the windows of bin offsets `offsets` upstream (-) and downstream (+) of each anchor.
    Returns arrays anchor, direction, counts per tissue; windows with invalid or off-chromosome bins dropped."""
    nb = len(chrom)
    rows = []
    for sign in (-1, 1):
        tgt = anchors[:, None] + sign * offsets[None, :]
        inside = (tgt >= 0) & (tgt < nb)
        tgt_c = np.clip(tgt, 0, nb - 1)
        ok = inside.all(1) & (chrom[tgt_c] == chrom[anchors][:, None]).all(1) & valid[tgt_c].all(1)
        a, t = anchors[ok], tgt_c[ok]
        lo = np.minimum(a[:, None], t)
        cs = [b[lo, np.abs(t - a[:, None])].sum(1) for b in bands]
        rows.append((a, np.full(len(a), sign), cs[0], cs[1]))
    return [np.concatenate([r[k] for r in rows]) for k in range(4)]


def bb_negll(params, y, n, base_logit):
    beta0, log_s = params
    p = special.expit(base_logit + beta0)
    s = np.exp(log_s)
    a, b = p * s, (1 - p) * s
    return -np.sum(special.betaln(y + a, n - y + b) - special.betaln(a, b))


def fit_bb(y, n, base_logit):
    m = n > 0
    y, n, base_logit = y[m], n[m], base_logit[m]
    b0 = special.logit(np.clip(y.sum() / n.sum(), 1e-6, 1 - 1e-6)) - np.median(base_logit)
    res = optimize.minimize(bb_negll, [b0, np.log(50.0)], args=(y, n, base_logit), method='Nelder-Mead',
                            options={'xatol': 1e-6, 'fatol': 1e-6, 'maxiter': 4000})
    beta0, log_s = res.x
    return beta0, 1.0 / (np.exp(log_s) + 1.0), len(y)


def bb_sf(y, n, p, rho):
    """P(Y >= y) under beta-binomial(n, p, rho), summed in log space."""
    if y <= 0:
        return 1.0
    if y > n:
        return 0.0
    s = 1.0 / max(rho, 1e-12) - 1.0
    a, b = p * s, (1 - p) * s
    k = np.arange(y, n + 1)
    lp = (special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1)
          + special.betaln(k + a, n - k + b) - special.betaln(a, b))
    return float(min(1.0, np.exp(special.logsumexp(lp))))


def fit_rho_local(y, n, base_logit, beta0, n_obs, k=500):
    """Overdispersion from the k background windows closest in total count to the tested window (beta0 fixed)."""
    m = n > 0
    y, n, base_logit = y[m], n[m], base_logit[m]
    idx = np.argsort(np.abs(np.log(n) - np.log(max(n_obs, 1))))[:k]
    f = lambda ls: bb_negll([beta0, ls], y[idx], n[idx], base_logit[idx])  # noqa: E731
    res = optimize.minimize_scalar(f, bounds=(np.log(0.5), np.log(1e7)), method='bounded')
    return 1.0 / (np.exp(res.x) + 1.0), int(n[idx].min()), int(n[idx].max())


def test_window(y, n, p_exp, rho, null_l2, null_n, k=500):
    """y = PA2 count, n = Md + PA2 count, p_exp = expected PA2 share, rho = count-matched overdispersion."""
    odds_exp = p_exp / (1 - p_exp)
    l2 = np.log2(((y + 0.5) / (n - y + 0.5)) / odds_exp)
    ph = (y + 0.5) / (n + 1.0)
    se = np.sqrt(1.0 / ((n + 1.0) * ph * (1 - ph)))
    deff = 1 + (n - 1) * rho
    z = stats.norm.ppf(0.975)
    lg = special.logit(ph) - np.log(odds_exp)
    ci_bb = (lg - z * se * np.sqrt(deff)) / np.log(2), (lg + z * se * np.sqrt(deff)) / np.log(2)
    cp = stats.beta.ppf([0.025, 0.975], [y, y + 1], [n - y + 1, n - y])
    cp = np.nan_to_num(cp, nan=0.0)
    ci_bin = [np.log2(max(c, 1e-12) / max(1 - c, 1e-12) / odds_exp) for c in cp]
    p_bin_up = float(stats.binom.sf(y - 1, n, p_exp))
    p_bin_two = float(stats.binomtest(int(y), int(n), p_exp).pvalue)
    p_bb_up = bb_sf(y, n, p_exp, rho)
    p_bb_two = min(1.0, 2 * min(p_bb_up, 1 - bb_sf(y + 1, n, p_exp, rho)))
    cm = (null_n >= 0.5 * n) & (null_n <= 2 * n)
    emp_all = (1 + np.sum(null_l2 >= l2)) / (1 + len(null_l2))
    emp_cm = (1 + np.sum(null_l2[cm] >= l2)) / (1 + cm.sum())
    near = np.argsort(np.abs(np.log(null_n) - np.log(max(n, 1))))[:k]
    med = np.median(null_l2[near])
    sd = 1.4826 * np.median(np.abs(null_l2[near] - med))
    z_rob = (l2 - med) / sd
    p_rob_two = float(2 * stats.norm.sf(abs(z_rob)))
    return dict(log2_ratio=l2, ci_bb_lo=ci_bb[0], ci_bb_hi=ci_bb[1], ci_binom_lo=ci_bin[0], ci_binom_hi=ci_bin[1],
                p_binom_upper=p_bin_up, p_binom_two=p_bin_two, p_bb_upper=p_bb_up, p_bb_two=p_bb_two,
                p_empirical_all=emp_all, p_empirical_countmatched=emp_cm, n_null=len(null_l2),
                n_null_countmatched=int(cm.sum()), robust_null_median=med, robust_null_sd=sd, robust_z=z_rob,
                p_robust_z_two=p_rob_two, robust_ci_lo=l2 - 1.96 * sd, robust_ci_hi=l2 + 1.96 * sd)


def promoter_bins(gtf, m):
    idx = {(c, s): i for i, (c, s) in enumerate(zip(np.asarray(m['names'])[m['chrom']], m['start']))}
    out = set()
    with open(gtf) as f:
        for line in f:
            p = line.split('\t')
            if len(p) < 9 or p[2] != 'transcript' or 'transcript_id "NM_' not in p[8]:
                continue
            c = p[0][3:] if p[0].startswith('chr') else p[0]
            tss = int(p[3]) - 1 if p[6] == '+' else int(p[4]) - 1
            i = idx.get((c, tss // BIN * BIN))
            if i is not None:
                out.add(i)
    return np.array(sorted(out))


def capture(washu_dir, lo, hi):
    rows = []
    for f in sorted(glob.glob(os.path.join(washu_dir, 'GSM65052*_Chicago_output_washU_text.txt.gz'))):
        name = os.path.basename(f).split('_Chicago')[0]
        total = bait = hire = 0
        scores = []
        with gzip.open(f, 'rt') as fh:
            for line in fh:
                a, b, s = line.rstrip('\n').split('\t')
                total += 1
                if BAIT in (a, b):
                    bait += 1
                    o = b if a == BAIT else a
                    c, st, _ = o.split(',')
                    if c == 'chr6' and lo <= int(st) < hi:
                        hire += 1
                        scores.append(float(s))
        rows.append(dict(sample=name, tissue='Md' if '_Md_' in name else 'PA2', calls_genome=total,
                         calls_hoxa2_bait=bait, calls_hoxa2_HIRE1=hire,
                         frac_bait_calls_in_HIRE1=hire / bait if bait else np.nan,
                         HIRE1_calls_per_100k_genome_calls=1e5 * hire / total,
                         median_score_HIRE1=float(np.median(scores)) if scores else np.nan))
    return pd.DataFrame(rows)


def main():
    md_path, pa2_path, gtf, washu_dir, out = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    mats = {'Md': load(md_path), 'PA2': load(pa2_path)}
    md, pa = mats['Md'], mats['PA2']
    assert np.array_equal(md['start'], pa['start']) and np.array_equal(md['chrom'], pa['chrom'])
    names = np.asarray(md['names'])
    chrom, start = md['chrom'], md['start']
    dmax = 60
    bands = [band(md, dmax), band(pa, dmax)]
    valid = np.isfinite(md['weight']) & np.isfinite(pa['weight']) & (md['rowcis'] > 0) & (pa['rowcis'] > 0)
    valid &= np.isin(names[chrom], AUTOSOMES_X)
    view = int(np.where((names[chrom] == CHROM) & (start == VIEW))[0][0])
    assert valid[view]
    depth = pd.DataFrame([dict(tissue=t, **mats[t]['depth'], viewpoint_cis_ge2bins=int(mats[t]['rowcis'][view]))
                          for t in ('Md', 'PA2')])
    depth.to_csv(os.path.join(out, 'hic_depth.tsv'), sep='\t', index=False)

    anchors_all = np.where(valid)[0]
    anchors_all = anchors_all[anchors_all != view]
    prom = promoter_bins(gtf, md)
    anchors_prom = prom[valid[prom] & (prom != view)]
    base_view = np.log(pa['rowcis'][view] / md['rowcis'][view])
    lib_ratio = np.log(pa['depth']['cis_ge2bins'] / md['depth']['cis_ge2bins'])

    def analyse(offsets, label, region):
        a_md = sum(int(bands[0][view - k, k]) for k in offsets)
        a_pa = sum(int(bands[1][view - k, k]) for k in offsets)
        mir_md = sum(int(bands[0][view, k]) for k in offsets)
        mir_pa = sum(int(bands[1][view, k]) for k in offsets)
        res = []
        for nullname, anchors in (('all_anchors', anchors_all), ('promoter_anchors', anchors_prom)):
            a, sign, c_md, c_pa = window_counts(bands, valid, chrom, anchors, offsets)
            not_hox = ~((names[chrom[a]] == CHROM) & (np.abs(a - view) <= 2))
            a, sign, c_md, c_pa = a[not_hox], sign[not_hox], c_md[not_hox], c_pa[not_hox]
            for norm in ('viewpoint_cis', 'library_cis'):
                if norm == 'viewpoint_cis':
                    base = np.log(pa['rowcis'][a] / md['rowcis'][a])
                    bv = base_view
                else:
                    base = np.full(len(a), lib_ratio)
                    bv = lib_ratio
                n = c_md + c_pa
                beta0, rho_global, nfit = fit_bb(c_pa, n, base)
                rho, nlo, nhi = fit_rho_local(c_pa, n, base, beta0, a_md + a_pa)
                p_exp = float(special.expit(bv + beta0))
                m = n > 0
                null_l2 = (np.log2((c_pa[m] + 0.5) / (c_md[m] + 0.5)) - (base[m] + beta0) / np.log(2))
                r = test_window(a_pa, a_md + a_pa, p_exp, rho, null_l2, n[m])
                mir_l2 = np.log2(((mir_pa + 0.5) / (mir_md + 0.5)) / (p_exp / (1 - p_exp)))
                res.append(dict(window=label, region=region, null=nullname, normalisation=norm,
                                count_Md=a_md, count_PA2=a_pa,
                                per1e4_viewpoint_Md=1e4 * a_md / md['rowcis'][view],
                                per1e4_viewpoint_PA2=1e4 * a_pa / pa['rowcis'][view],
                                per1e6_library_Md=1e6 * a_md / md['depth']['cis_ge2bins'],
                                per1e6_library_PA2=1e6 * a_pa / pa['depth']['cis_ge2bins'],
                                expected_PA2_share=p_exp, null_beta0_log2=beta0 / np.log(2), null_rho_global=rho_global,
                                null_rho_countmatched=rho, null_countmatched_n_range=f'{nlo}-{nhi}',
                                null_windows_fitted=nfit,
                                mirror_count_Md=mir_md, mirror_count_PA2=mir_pa, mirror_log2_ratio=mir_l2, **r))
        return res

    rows = []
    for label, (lo, hi) in WINDOWS.items():
        offsets = np.arange((VIEW - hi) // BIN + 1, (VIEW - lo) // BIN + 1)
        rows += analyse(offsets, label, f'chr6:{lo}-{hi}')
    pd.DataFrame(rows).to_csv(os.path.join(out, 'hic_windows.tsv'), sep='\t', index=False)

    brow = []
    for s in range(SCAN[0], SCAN[1], BIN):
        k = (VIEW - s) // BIN
        for r in analyse(np.array([k]), f'bin_{s}', f'chr6:{s}-{s + BIN}'):
            brow.append(r)
    bins = pd.DataFrame(brow)
    for (nl, nm), g in bins.groupby(['null', 'normalisation']):
        for col in ('p_bb_two', 'p_binom_two', 'p_robust_z_two'):
            q = stats.false_discovery_control(g[col].values)
            bins.loc[g.index, col.replace('p_', 'q_')] = q
    bins.to_csv(os.path.join(out, 'hic_bins.tsv'), sep='\t', index=False)

    lo_v, hi_v = 50500000, 53800000
    sel = np.where((names[chrom] == CHROM) & (start >= lo_v) & (start < hi_v))[0]
    v4c = []
    for i in sel:
        k = abs(i - view)
        lo_i = min(i, view)
        c = [int(b[lo_i, k]) if k <= dmax else np.nan for b in bands]
        v4c.append(dict(start=int(start[i]), distance=int(start[i] - VIEW), count_Md=c[0], count_PA2=c[1],
                        per1e4_viewpoint_Md=1e4 * c[0] / md['rowcis'][view],
                        per1e4_viewpoint_PA2=1e4 * c[1] / pa['rowcis'][view]))
    v4c = pd.DataFrame(v4c)
    v4c.to_csv(os.path.join(out, 'hic_virtual4c.tsv'), sep='\t', index=False)

    cap = capture(washu_dir, 50926155, 51171572)
    cap.to_csv(os.path.join(out, 'capture_replicates.tsv'), sep='\t', index=False)
    tab = [[int(cap.loc[cap.tissue == t, 'calls_hoxa2_HIRE1'].sum()),
            int((cap.loc[cap.tissue == t, 'calls_hoxa2_bait'] - cap.loc[cap.tissue == t, 'calls_hoxa2_HIRE1']).sum())]
           for t in ('PA2', 'Md')]
    fisher = stats.fisher_exact(tab)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    v = v4c.dropna()
    v = v[v.distance != 0]
    ax.step(v.start / 1e6, v.per1e4_viewpoint_Md, where='post', color='#4c72b0', label='Md', lw=1)
    ax.step(v.start / 1e6, v.per1e4_viewpoint_PA2, where='post', color='#c44e52', label='PA2', lw=1)
    ax.axvspan(WINDOWS['HIRE1_partners'][0] / 1e6, WINDOWS['HIRE1_partners'][1] / 1e6, color='0.85', zorder=0)
    ax.axvline(51175000 / 1e6, color='0.4', ls='--', lw=0.8)
    ax.axvline(VIEW / 1e6, color='k', lw=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('chr6 (Mb, mm10)')
    ax.set_ylabel('contacts per 10$^4$\nviewpoint cis contacts')
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'virtual4c_Md_PA2.pdf'))

    summary = dict(depth=depth.to_dict('records'), capture_fisher_pooled=dict(table_PA2_Md=tab,
                   odds_ratio=float(fisher.statistic), p=float(fisher.pvalue)),
                   n_promoter_anchors=int(len(anchors_prom)), n_all_anchors=int(len(anchors_all)))
    with open(os.path.join(out, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1, default=float)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    print(depth)
    print(pd.DataFrame(rows))
    print(bins[(bins.null == 'all_anchors') & (bins.normalisation == 'viewpoint_cis')][
        ['region', 'count_Md', 'count_PA2', 'log2_ratio', 'ci_bb_lo', 'ci_bb_hi', 'null_rho_countmatched',
         'p_binom_two', 'p_bb_two', 'q_bb_two', 'p_robust_z_two', 'q_robust_z_two', 'p_empirical_countmatched']])
    print(cap)
    print(summary['capture_fisher_pooled'])


if __name__ == '__main__':
    main()
