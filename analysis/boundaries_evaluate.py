"""Boundary callers against Hi-C (27_boundaries.sh, after boundaries_features.py).

Methods (epigeneticgraphanalyzer.boundaries.METHODS; L2 logistic regression on per-cell z-scored features):
  density           bait and edge density only (covered, log expected spanning edges, number of baits): a control
                    that knows where the baits are but not to which side their edges point
  capture           density + capture insulation features
  ctcf_cohesin      CTCF / cohesin site strength and counts
  ctcf_orientation  + CTCF motif orientation (boundary pattern: reverse sites left, forward sites right)
  ctcf_orientation_density  + density
  combined          ctcf_orientation + density + capture insulation + convergent cohesin-CTCF loops among the edges
  capture_minima    unsupervised: local minima of the capture insulation with prominence 0.3 (no training)
Training label: Hi-C insulation boundaries of the training cell at the same resolution (local minima of the diamond
insulation with prominence 0.3, boundaries_hic.py). Calls: local maxima of the predicted probability (+-50 kb) above
the threshold that maximises F1 (tolerance 1 bin) on the training data.
Designs: cross-cell (train on one cell type, all chromosomes; test on another) and within-cell chromosome split
(odd vs even chromosomes, chrX with the even ones; each half predicted by the model of the other half).
References: hic_insulation (same pipeline as the label), hicFindTADs 50 kb (4DN mcool, 23_hic_reference.sh),
Rao et al. 2014 Arrowhead (GM12878, hg19 lifted to GRCh38 by position), 4DN H1 Micro-C boundaries.
Metrics: precision, recall, F1 at tolerance 1 and 2 bins; AUPRC of the per-bin score against the bins holding a
reference boundary (tolerance 0) with the positive rate; all bins and covered bins (capture expectation >= 5 edges).
Contacts beyond CTCF: combined - ctcf_orientation, combined - ctcf_orientation_density (the contribution of the
direction of the contacts beyond bait density), ctcf_orientation_density - ctcf_orientation and capture - density;
delta AUPRC and delta F1 with chromosome block-bootstrap intervals (resampling chromosomes with replacement), and
likelihood-ratio statistics of the nested models on the training cell (bins are autocorrelated, so the P values are
optimistic). Domains: segments between consecutive calls against
segments between reference boundaries, normalised variation of information and the agreement on same-domain status
of bin pairs <= 1 Mb apart (local Rand), each against random boundaries of the same number per chromosome.

Usage: python boundaries_evaluate.py <results dir> <raw dir> [threads] [n_boot]
Output: <results>/boundaries/tables/eval_*.tsv, models/*.json, calls/*.bed, cache/predictions_*.tsv.gz
"""

import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402
from epigeneticgraphanalyzer.compartments import domain_boundaries, flatten_domains, read_domains  # noqa: E402
from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402

CELLS = ['GM12878', 'K562', 'H1']
RESOLUTIONS = (10000, 25000, 50000)
METHODS = ['density', 'capture', 'ctcf_cohesin', 'ctcf_orientation', 'ctcf_orientation_density', 'combined']
COMPARISONS = [('combined', 'ctcf_orientation'), ('combined', 'ctcf_orientation_density'),
               ('ctcf_orientation_density', 'ctcf_orientation'), ('capture', 'density')]
HIC_PROMINENCE = 0.3
CAPTURE_PROMINENCE = 0.3
CONTEXT_BP = 50000
RAW = None
RESULTS = None


def out_dir(*parts):
    return os.path.join(RESULTS, 'boundaries', *parts)


def load_features(cell, res):
    t = pd.read_csv(out_dir('cache', 'features_{}_{}kb.tsv.gz'.format(cell, res // 1000)), sep='\t')
    t['covered'] = t.covered.astype(float)
    return t


def references(cell, table):
    """Reference boundary bins of a cell: dict name -> bin indices."""
    codes = bd.chrom_codes(table)
    refs = {'hic_insulation': bd.local_minima(table.hic_insulation.values, codes, HIC_PROMINENCE)}
    ref_dir = os.path.join(RESULTS, 'compartments_domains', 'hic_reference')
    dom = flatten_domains(read_domains(os.path.join(ref_dir, '{}_50kb_domains.bed'.format(cell))))
    refs['hicFindTADs_50kb'] = bd.positions_to_bins(table, domain_boundaries(dom))
    if cell == 'GM12878':
        path = os.path.join(RAW, 'domains', 'GM12878_rao2014_hg19',
                            'GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz')
        bounds = domain_boundaries(flatten_domains(read_domains(path)))
        lift = LiftOver(os.path.join(RAW, 'annotation', 'hg19ToHg38.over.chain.gz'))
        lifted = {}
        for chrom, pos in bounds.items():
            hits = [lift.position(chrom, int(p)) for p in pos]
            lifted[chrom] = np.array([h[1] for h in hits if h is not None and h[0] == chrom])
        refs['rao2014_arrowhead'] = bd.positions_to_bins(table, lifted)
    if cell == 'H1':
        z = pd.read_csv(os.path.join(RAW, 'domains', 'H1_4dn_hg38', '4DNFIED5HLDC.bed.gz'), sep='\t', header=None)
        refs['4dn_H1_microc'] = bd.positions_to_bins(table, {c: ((g[1] + g[2]) // 2).values
                                                             for c, g in z.groupby(0)})
    return refs


def radius(res):
    return max(1, int(round(CONTEXT_BP / res)))


def fit_method(method, table, label, codes, res):
    model = bd.LogisticModel(bd.METHODS[method]).fit(table, label)
    score = model.predict(table)
    thr, f1 = bd.best_threshold(score, np.flatnonzero(label), codes, tolerance=1, radius=radius(res))
    return model, thr, f1


def evaluate(score, calls, refs, codes, covered, res):
    rows = []
    for ref_name, ref in refs.items():
        for region, mask in (('all', None), ('covered', covered)):
            sel = np.ones(len(codes), dtype=bool) if mask is None else mask
            row = {'reference': ref_name, 'region': region}
            for tol in (1, 2):
                m = bd.match_boundaries(calls, ref, codes, tolerance=tol, mask=mask)
                row.update({'{}_tol{}'.format(k, tol): v for k, v in m.items()})
            label = bd.bin_labels(len(codes), ref)
            row['auprc'] = bd.average_precision(score[sel], label[sel])
            row['positive_rate'] = float(label[sel].mean())
            row['auprc_tol1'] = bd.average_precision(score[sel], bd.bin_labels(len(codes), ref, 1, codes)[sel])
            row['positive_rate_tol1'] = float(bd.bin_labels(len(codes), ref, 1, codes)[sel].mean())
            rows.append(row)
    return rows


def per_chrom_counts(calls, ref, codes, tol, mask=None):
    """Per chromosome: called, called matched, reference, reference matched."""
    out = np.zeros((codes.max() + 1, 4))
    for c in range(codes.max() + 1):
        sel = codes == c
        k_calls = calls[sel[calls]]
        k_ref = ref[sel[ref]]
        m = bd.match_boundaries(k_calls, k_ref, codes, tolerance=tol, mask=mask)
        out[c] = [m['n_called'], m['n_called'] * np.nan_to_num(m['precision']), m['n_reference'],
                  m['n_reference'] * np.nan_to_num(m['recall'])]
    return out


def f1_from_counts(cnt):
    p = cnt[1] / cnt[0] if cnt[0] else 0.0
    r = cnt[3] / cnt[2] if cnt[2] else 0.0
    return 2 * p * r / (p + r) if p + r > 0 else 0.0


def bootstrap_delta(args):
    key, s_full, s_red, calls_full, calls_red, ref, codes, mask, n_boot, seed = args
    rng = np.random.default_rng(seed)
    sel = np.ones(len(codes), dtype=bool) if mask is None else mask
    label = bd.bin_labels(len(codes), ref)
    n_chrom = codes.max() + 1
    by_chrom = [np.flatnonzero((codes == c) & sel) for c in range(n_chrom)]
    cf = per_chrom_counts(calls_full, ref, codes, 1, mask)
    cr = per_chrom_counts(calls_red, ref, codes, 1, mask)
    obs_ap = bd.average_precision(s_full[sel], label[sel]) - bd.average_precision(s_red[sel], label[sel])
    obs_f1 = f1_from_counts(cf.sum(0)) - f1_from_counts(cr.sum(0))
    d_ap, d_f1 = [], []
    for _ in range(n_boot):
        pick = rng.integers(0, n_chrom, n_chrom)
        idx = np.concatenate([by_chrom[c] for c in pick])
        d_ap.append(bd.average_precision(s_full[idx], label[idx]) - bd.average_precision(s_red[idx], label[idx]))
        d_f1.append(f1_from_counts(cf[pick].sum(0)) - f1_from_counts(cr[pick].sum(0)))
    d_ap, d_f1 = np.array(d_ap), np.array(d_f1)
    row = dict(key)
    row.update({'delta_auprc': obs_ap, 'delta_auprc_lo': np.nanpercentile(d_ap, 2.5),
                'delta_auprc_hi': np.nanpercentile(d_ap, 97.5), 'p_delta_auprc_le0': float(np.mean(d_ap <= 0)),
                'delta_f1_tol1': obs_f1, 'delta_f1_lo': np.percentile(d_f1, 2.5),
                'delta_f1_hi': np.percentile(d_f1, 97.5), 'p_delta_f1_le0': float(np.mean(d_f1 <= 0)),
                'n_boot': n_boot})
    return row


def domain_agreement(calls, ref, codes, res, rng, n_null=20):
    lab_c = bd.segment_labels(codes, calls)
    lab_r = bd.segment_labels(codes, ref)
    max_off = int(1_000_000 // res)
    vi = bd.variation_of_information(lab_c, lab_r)
    lr, same_c, same_r = bd.local_rand(lab_c, lab_r, codes, max_off)
    counts = pd.Series(codes[calls]).value_counts().to_dict()
    vi0, lr0 = [], []
    for _ in range(n_null):
        rb = bd.random_boundaries(codes, counts, rng)
        lab0 = bd.segment_labels(codes, rb)
        vi0.append(bd.variation_of_information(lab0, lab_r))
        lr0.append(bd.local_rand(lab0, lab_r, codes, max_off)[0])
    return {'vi_norm': vi, 'vi_norm_random': float(np.mean(vi0)), 'local_rand_1mb': lr,
            'local_rand_1mb_random': float(np.mean(lr0)), 'same_domain_pairs_called': same_c,
            'same_domain_pairs_reference': same_r,
            'median_domain_kb_called': float(np.median(np.bincount(lab_c)[np.unique(lab_c)]) * res / 1000),
            'median_domain_kb_reference': float(np.median(np.bincount(lab_r)[np.unique(lab_r)]) * res / 1000)}


def main():
    global RAW, RESULTS
    RESULTS, RAW = sys.argv[1], sys.argv[2]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    n_boot = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    for d in ('tables', 'models', 'calls', 'cache'):
        os.makedirs(out_dir(d), exist_ok=True)
    rows, coefs, lrt, cover, dom_rows, boot_jobs = [], [], [], [], [], []
    rng = np.random.default_rng(0)
    for res in RESOLUTIONS:
        data = {}
        for cell in CELLS:
            t = load_features(cell, res)
            codes = bd.chrom_codes(t)
            refs = references(cell, t)
            covered = t.covered.values > 0
            data[cell] = (t, codes, refs, covered)
            label = bd.bin_labels(len(t), refs['hic_insulation'])
            cover.append({'cell': cell, 'resolution_kb': res // 1000, 'n_bins': len(t),
                          'covered_fraction': float(covered.mean()),
                          **{'ref_in_covered_' + k: float(covered[v].mean()) for k, v in refs.items()},
                          **{'n_ref_' + k: len(v) for k, v in refs.items()}})
            # nested likelihood ratio on the full training cell
            for big, small in COMPARISONS:
                full = bd.LogisticModel(bd.METHODS[big]).fit(t, label)
                red = bd.LogisticModel(bd.METHODS[small]).fit(t, label)
                stat = 2 * (full.loglik - red.loglik)
                df = len(bd.METHODS[big]) - len(bd.METHODS[small])
                lrt.append({'cell': cell, 'resolution_kb': res // 1000, 'full': big, 'reduced': small,
                            'lr_statistic': stat, 'df': df, 'p_naive': chi2.sf(stat, df), 'loglik_full': full.loglik,
                            'loglik_reduced': red.loglik, 'n_bins': len(t), 'n_boundaries': int(label.sum())})
        preds = {cell: pd.DataFrame(index=data[cell][0].index) for cell in CELLS}
        # within-cell chromosome split
        for cell in CELLS:
            t, codes, refs, covered = data[cell]
            label = bd.bin_labels(len(t), refs['hic_insulation'])
            num = t.chrom.str.replace('chr', '').replace({'X': '24'}).astype(int).values
            fold = (num % 2 == 0).astype(int)
            for method in METHODS:
                score = np.zeros(len(t))
                calls = []
                for f in (0, 1):
                    tr, te = fold != f, fold == f
                    model, thr, f1_train = fit_method(method, t[tr], label[tr], codes[tr], res)
                    s = model.predict(t)
                    score[te] = s[te]
                    calls.append(bd.peak_calls(np.where(te, s, -np.inf), codes, thr, radius(res)))
                calls = np.sort(np.concatenate(calls))
                preds[cell]['{}__{}'.format(method, 'chromsplit')] = score
                for r in evaluate(score, calls, refs, codes, covered, res):
                    rows.append({'design': 'chromosome_split', 'train': cell, 'test': cell, 'method': method,
                                 'resolution_kb': res // 1000, **r})
                if method == 'combined':
                    dom_rows.append({'design': 'chromosome_split', 'train': cell, 'test': cell, 'method': method,
                                     'resolution_kb': res // 1000, 'reference': 'hic_insulation',
                                     **domain_agreement(calls, refs['hic_insulation'], codes, res, rng)})
                preds[cell]['calls__{}__chromsplit'.format(method)] = np.isin(np.arange(len(t)), calls)
            for ref_name in ('hic_insulation', 'hicFindTADs_50kb'):
                for region in ('all', 'covered'):
                    boot_jobs.append(({'design': 'chromosome_split', 'train': cell, 'test': cell,
                                       'resolution_kb': res // 1000, 'reference': ref_name, 'region': region},
                                      cell, 'chromsplit', ref_name, region))
        # cross-cell
        for train in CELLS:
            t, codes, refs, covered = data[train]
            label = bd.bin_labels(len(t), refs['hic_insulation'])
            for method in METHODS:
                model, thr, f1_train = fit_method(method, t, label, codes, res)
                spec = model.to_dict()
                spec.update(threshold=thr, radius_bins=radius(res), resolution=res, train=train, method=method,
                            train_f1_tol1=f1_train, label='hic_insulation prominence {}'.format(HIC_PROMINENCE))
                with open(out_dir('models', '{}_{}kb_{}.json'.format(train, res // 1000, method)), 'w') as h:
                    json.dump(spec, h, indent=1)
                for f, c in zip(model.features, model.coef):
                    coefs.append({'train': train, 'resolution_kb': res // 1000, 'method': method, 'feature': f,
                                  'coef': c})
                for test in CELLS:
                    if test == train:
                        continue
                    tt, tcodes, trefs, tcov = data[test]
                    s = model.predict(tt)
                    calls = bd.peak_calls(s, tcodes, thr, radius(res))
                    preds[test]['{}__{}'.format(method, train)] = s
                    preds[test]['calls__{}__{}'.format(method, train)] = np.isin(np.arange(len(tt)), calls)
                    for r in evaluate(s, calls, trefs, tcodes, tcov, res):
                        rows.append({'design': 'cross_cell', 'train': train, 'test': test, 'method': method,
                                     'resolution_kb': res // 1000, 'threshold': thr, **r})
                    if method in ('combined', 'ctcf_orientation', 'capture', 'density'):
                        for ref_name in ('hic_insulation', 'hicFindTADs_50kb'):
                            dom_rows.append({'design': 'cross_cell', 'train': train, 'test': test, 'method': method,
                                             'resolution_kb': res // 1000, 'reference': ref_name,
                                             **domain_agreement(calls, trefs[ref_name], tcodes, res, rng)})
                    if method == 'combined':
                        cells = tt.loc[calls, ['chrom', 'start', 'end']].assign(name='boundary', score=s[calls])
                        name = '{}_{}kb_{}_from_{}.bed'.format(test, res // 1000, method, train)
                        cells.to_csv(out_dir('calls', name),
                                     sep='\t', header=False, index=False, float_format='%.4g')
            for test in CELLS:
                if test == train:
                    continue
                for ref_name in ('hic_insulation', 'hicFindTADs_50kb'):
                    for region in ('all', 'covered'):
                        boot_jobs.append(({'design': 'cross_cell', 'train': train, 'test': test,
                                           'resolution_kb': res // 1000, 'reference': ref_name, 'region': region},
                                          test, train, ref_name, region))
        # unsupervised capture minima and the Hi-C insulation itself (upper reference)
        for cell in CELLS:
            t, codes, refs, covered = data[cell]
            ins = t.capture_insulation.values
            calls = bd.local_minima(ins, codes, CAPTURE_PROMINENCE)
            score = np.where(np.isfinite(ins), -ins, np.nanmin(-ins) - 1)
            for r in evaluate(score, calls, refs, codes, covered, res):
                rows.append({'design': 'unsupervised', 'train': '-', 'test': cell, 'method': 'capture_minima',
                             'resolution_kb': res // 1000, **r})
            dom_rows.append({'design': 'unsupervised', 'train': '-', 'test': cell, 'method': 'capture_minima',
                             'resolution_kb': res // 1000, 'reference': 'hic_insulation',
                             **domain_agreement(calls, refs['hic_insulation'], codes, res, rng)})
            preds[cell]['capture_minima'] = np.isin(np.arange(len(t)), calls)
            for name, ref in refs.items():
                preds[cell]['ref__' + name] = np.isin(np.arange(len(t)), ref)
            pd.concat([t[['chrom', 'start', 'end', 'capture_insulation', 'hic_insulation', 'covered']], preds[cell]],
                      axis=1).to_csv(out_dir('cache', 'predictions_{}_{}kb.tsv.gz'.format(cell, res // 1000)),
                                     sep='\t', index=False, float_format='%.4g')
        # bootstrap jobs of this resolution
        jobs = []
        for i, (key, cell, src, ref_name, region) in enumerate(boot_jobs):
            t, codes, refs, covered = data[cell]
            for j, (big, small) in enumerate(COMPARISONS):
                pf = preds[cell]['{}__{}'.format(big, src)].values
                pr = preds[cell]['{}__{}'.format(small, src)].values
                cf = np.flatnonzero(preds[cell]['calls__{}__{}'.format(big, src)].values)
                cr = np.flatnonzero(preds[cell]['calls__{}__{}'.format(small, src)].values)
                jobs.append((dict(key, full=big, reduced=small), pf, pr, cf, cr, refs[ref_name], codes,
                             covered if region == 'covered' else None, n_boot, 100 * i + j))
        boot_jobs = []
        with Pool(threads) as pool:
            boot = pool.map(bootstrap_delta, jobs)
        pd.DataFrame(boot).to_csv(out_dir('tables', 'eval_delta_{}kb.tsv'.format(res // 1000)), sep='\t',
                                  index=False, float_format='%.4g')
        print('resolution', res, 'done', flush=True)
    pd.DataFrame(rows).to_csv(out_dir('tables', 'eval_boundaries.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(coefs).to_csv(out_dir('tables', 'model_coefficients.tsv'), sep='\t', index=False,
                               float_format='%.4g')
    pd.DataFrame(lrt).to_csv(out_dir('tables', 'eval_lrt.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(cover).to_csv(out_dir('tables', 'coverage.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(dom_rows).to_csv(out_dir('tables', 'eval_domains.tsv'), sep='\t', index=False, float_format='%.4g')
    delta = pd.concat([pd.read_csv(out_dir('tables', 'eval_delta_{}kb.tsv'.format(r // 1000)), sep='\t')
                       for r in RESOLUTIONS])
    delta.to_csv(out_dir('tables', 'eval_delta.tsv'), sep='\t', index=False, float_format='%.4g')
    ev = pd.DataFrame(rows)
    pd.set_option('display.width', 250)
    sub = ev[(ev.reference == 'hic_insulation')]
    print(sub.pivot_table(index=['design', 'train', 'test', 'resolution_kb'], columns=['method', 'region'],
                          values='auprc').round(3).to_string())
    print(delta.to_string())


if __name__ == '__main__':
    main()
