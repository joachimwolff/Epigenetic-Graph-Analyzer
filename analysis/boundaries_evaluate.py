"""Rule-based boundary calls against Hi-C (27_boundaries.sh, after boundaries_features.py).

No fitted model: every rule combines the per-bin evidence (epigeneticgraphanalyzer.boundaries.bin_evidence) with
fixed thresholds chosen before the evaluation.
Evidence (window = +-1 bin):
  dip           capture insulation (log2 observed / expected spanning edges, mirror null, covered bins) is a local
                minimum with prominence >= P and value <= 0
  ctcf          a cohesin-bound CTCF site with CTCF enrichment >= C in the bin (combined with a dip: within the
                window)
  ctcf_pattern  reverse site(s) in bins i-1..i and forward site(s) in bins i..i+1 (cohesin-bound, enrichment >= C)
  loop_anchor   an end of a convergent cohesin-CTCF loop (graph edge) within the window
Rules: dip, ctcf, ctcf_pattern, loop_anchor, dip_or_ctcf, dip_and_ctcf, dip_or_ctcf_pattern, dip_and_ctcf_pattern;
one call per run of consecutive flagged bins (the bin with the highest score).
Thresholds: a priori P = 0.3 (the prominence of the Hi-C insulation boundaries) and C = 10 (all called sites); a
grid P in {0.2, 0.3, 0.5} x C in {10, 20, 40} is reported in full and is not used to pick a setting.
Scores for AUPRC: dip = -insulation, ctcf = summed enrichment, ctcf_pattern = min(strongest reverse-left,
strongest forward-right site), loop_anchor = summed enrichment of the anchors; AND / OR = minimum / maximum of the
percentile ranks of the two component scores (unweighted).
References: hic_insulation (diamond insulation minima, prominence 0.3, boundaries_hic.py), hicFindTADs 50 kb
(23_hic_reference.sh), Rao et al. 2014 Arrowhead (GM12878, lifted to GRCh38), 4DN H1 Micro-C boundaries; the Hi-C
insulation boundaries are also evaluated against the other references (agreement between Hi-C callers).
Metrics: precision, recall, F1 at tolerance 1 and 2 bins, AUPRC (label: bin holds a reference boundary) and the
positive rate; all bins and covered bins (capture expectation >= 5 spanning edges; calls and reference boundaries
both restricted). Does the dip add to CTCF: dip_and_X - X and dip_or_X - X (X = ctcf, ctcf_pattern) in precision,
recall, F1 and AUPRC with 95 % chromosome block-bootstrap intervals (chromosomes resampled with replacement).
Domains: segments between consecutive calls against segments between reference boundaries: normalised variation of
information and the agreement on same-domain status of bin pairs <= 1 Mb apart (local Rand), each with random
boundaries of the same number per chromosome.

Usage: python boundaries_evaluate.py <results dir> <raw dir> [threads] [n_boot]
Output: <results>/boundaries/tables/eval_rules.tsv, eval_rules_grid.tsv, eval_rules_delta.tsv, eval_rules_domains.tsv,
        coverage.tsv, evidence_summary.tsv; cache/evidence_<cell>_<res>kb.tsv.gz; calls/<cell>_<res>kb_<rule>.bed
"""

import itertools
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402
from epigeneticgraphanalyzer.compartments import domain_boundaries, flatten_domains, read_domains  # noqa: E402
from epigeneticgraphanalyzer.genome import LiftOver  # noqa: E402

CELLS = ['GM12878', 'K562', 'H1']
RESOLUTIONS = (10000, 25000, 50000)
RULES = list(bd.RULES)
HIC_PROMINENCE = 0.3
DEFAULT = {'dip_prominence': 0.3, 'ctcf_min': 10.0}
GRID = list(itertools.product((0.2, 0.3, 0.5), (10.0, 20.0, 40.0)))
COMPARISONS = [('dip_and_ctcf', 'ctcf'), ('dip_or_ctcf', 'ctcf'), ('dip_and_ctcf_pattern', 'ctcf_pattern'),
               ('dip_or_ctcf_pattern', 'ctcf_pattern'), ('ctcf_pattern', 'ctcf'), ('dip_and_ctcf_pattern', 'dip')]
RAW = None
RESULTS = None


def out_dir(*parts):
    return os.path.join(RESULTS, 'boundaries', *parts)


def load(cell, res):
    t = pd.read_csv(out_dir('cache', 'features_{}_{}kb.tsv.gz'.format(cell, res // 1000)), sep='\t')
    sites = pd.read_csv(out_dir('cache', 'sites_{}.tsv.gz'.format(cell)), sep='\t', keep_default_na=False,
                        na_values=['nan', 'NaN'])
    loops = pd.read_csv(out_dir('cache', 'loops_{}.tsv.gz'.format(cell)), sep='\t')
    return t, sites, loops


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


def insulation_frame(t):
    return pd.DataFrame({'insulation': t.insulation.values, 'covered': t.covered.values.astype(bool)})


def evaluate(score, calls, refs, codes, covered):
    rows = []
    for ref_name, ref in refs.items():
        for region, mask in (('all', None), ('covered', covered)):
            sel = np.ones(len(codes), dtype=bool) if mask is None else mask
            row = {'reference': ref_name, 'region': region}
            for tol in (1, 2):
                m = bd.match_boundaries(calls, ref, codes, tolerance=tol, mask=mask)
                row.update({'{}_tol{}'.format(k, tol): v for k, v in m.items()})
            label = bd.bin_labels(len(codes), ref)
            row['auprc'] = bd.average_precision(score[sel], label[sel]) if score is not None else np.nan
            row['positive_rate'] = float(label[sel].mean())
            rows.append(row)
    return rows


def per_chrom_counts(calls, ref, codes, tol, mask=None):
    """Per chromosome: called, called matched, reference, reference matched."""
    out = np.zeros((codes.max() + 1, 4))
    for c in range(codes.max() + 1):
        sel = codes == c
        m = bd.match_boundaries(calls[sel[calls]], ref[sel[ref]], codes, tolerance=tol, mask=mask)
        out[c] = [m['n_called'], m['n_called'] * np.nan_to_num(m['precision']), m['n_reference'],
                  m['n_reference'] * np.nan_to_num(m['recall'])]
    return out


def prf(cnt):
    p = cnt[1] / cnt[0] if cnt[0] else 0.0
    r = cnt[3] / cnt[2] if cnt[2] else 0.0
    return p, r, (2 * p * r / (p + r) if p + r > 0 else 0.0)


def bootstrap_delta(args):
    key, s_a, s_b, calls_a, calls_b, ref, codes, mask, n_boot, seed = args
    rng = np.random.default_rng(seed)
    sel = np.ones(len(codes), dtype=bool) if mask is None else mask
    label = bd.bin_labels(len(codes), ref)
    n_chrom = codes.max() + 1
    by_chrom = [np.flatnonzero((codes == c) & sel) for c in range(n_chrom)]
    ca = per_chrom_counts(calls_a, ref, codes, 1, mask)
    cb = per_chrom_counts(calls_b, ref, codes, 1, mask)
    obs = np.array(prf(ca.sum(0))) - np.array(prf(cb.sum(0)))
    obs_ap = bd.average_precision(s_a[sel], label[sel]) - bd.average_precision(s_b[sel], label[sel])
    boot, boot_ap = [], []
    for _ in range(n_boot):
        pick = rng.integers(0, n_chrom, n_chrom)
        boot.append(np.array(prf(ca[pick].sum(0))) - np.array(prf(cb[pick].sum(0))))
        idx = np.concatenate([by_chrom[c] for c in pick])
        boot_ap.append(bd.average_precision(s_a[idx], label[idx]) - bd.average_precision(s_b[idx], label[idx]))
    boot, boot_ap = np.array(boot), np.array(boot_ap)
    row = dict(key)
    for k, name in enumerate(('precision', 'recall', 'f1')):
        row.update({'delta_' + name: obs[k], name + '_lo': np.percentile(boot[:, k], 2.5),
                    name + '_hi': np.percentile(boot[:, k], 97.5)})
    row.update({'delta_auprc': obs_ap, 'auprc_lo': np.nanpercentile(boot_ap, 2.5),
                'auprc_hi': np.nanpercentile(boot_ap, 97.5),
                'n_called_with': int(ca[:, 0].sum()), 'n_called_without': int(cb[:, 0].sum()), 'n_boot': n_boot})
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
        lab0 = bd.segment_labels(codes, bd.random_boundaries(codes, counts, rng))
        vi0.append(bd.variation_of_information(lab0, lab_r))
        lr0.append(bd.local_rand(lab0, lab_r, codes, max_off)[0])
    sizes_c = np.bincount(lab_c)[np.unique(lab_c)]
    sizes_r = np.bincount(lab_r)[np.unique(lab_r)]
    return {'n_calls': len(calls), 'vi_norm': vi, 'vi_norm_random': float(np.mean(vi0)), 'local_rand_1mb': lr,
            'local_rand_1mb_random': float(np.mean(lr0)), 'same_domain_pairs_called': same_c,
            'same_domain_pairs_reference': same_r, 'median_domain_kb_called': float(np.median(sizes_c) * res / 1000),
            'median_domain_kb_reference': float(np.median(sizes_r) * res / 1000)}


def evidence_summary(cell, res, ev, refs, codes):
    """How often each evidence occurs and how often it lies within 1 bin of a Hi-C insulation boundary."""
    near = bd.bin_labels(len(codes), refs['hic_insulation'], 1, codes)
    rows = []
    for name in ('dip', 'ctcf', 'ctcf_pattern', 'loop'):
        f = ev[name].values
        for region, sel in (('all', np.ones(len(f), dtype=bool)), ('covered', ev.covered.values)):
            rows.append({'cell': cell, 'resolution_kb': res // 1000, 'evidence': name, 'region': region,
                         'fraction_of_bins': float(f[sel].mean()),
                         'fraction_near_hic_boundary': float(near[sel & f].mean()) if (sel & f).any() else np.nan,
                         'fraction_near_hic_boundary_all_bins': float(near[sel].mean())})
    return rows


def main():
    global RAW, RESULTS
    RESULTS, RAW = sys.argv[1], sys.argv[2]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    n_boot = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    for d in ('tables', 'calls', 'cache'):
        os.makedirs(out_dir(d), exist_ok=True)
    rows, grid_rows, cover, dom_rows, summ, jobs = [], [], [], [], [], []
    rng = np.random.default_rng(0)
    for res in RESOLUTIONS:
        for cell in CELLS:
            t, sites, loops = load(cell, res)
            codes = bd.chrom_codes(t)
            refs = references(cell, t)
            covered = t.covered.values.astype(bool)
            ins = insulation_frame(t)
            cover.append({'cell': cell, 'resolution_kb': res // 1000, 'n_bins': len(t),
                          'covered_fraction': float(covered.mean()),
                          **{'ref_in_covered_' + k: float(covered[v].mean()) for k, v in refs.items()},
                          **{'n_ref_' + k: len(v) for k, v in refs.items()}})
            # agreement between the Hi-C references themselves
            others = {k: v for k, v in refs.items() if k != 'hic_insulation'}
            for r in evaluate(None, refs['hic_insulation'], others, codes, covered):
                rows.append({'cell': cell, 'resolution_kb': res // 1000, 'rule': 'hic_insulation_minima',
                             'dip_prominence': np.nan, 'ctcf_min': np.nan, **r})
            for prom, cmin in GRID:
                ev = bd.bin_evidence(t[['chrom', 'start', 'end']], sites, ins, loops, dip_prominence=prom,
                                     ctcf_min=cmin)
                default = prom == DEFAULT['dip_prominence'] and cmin == DEFAULT['ctcf_min']
                out_ev = t[['chrom', 'start', 'end', 'hic_insulation']].copy() if default else None
                for rule in RULES:
                    flags, score = bd.rule_flags(ev, rule)
                    calls = bd.calls_from_flags(flags, score, codes)
                    use = refs if default else {'hic_insulation': refs['hic_insulation']}
                    for r in evaluate(score, calls, use, codes, covered):
                        rec = {'cell': cell, 'resolution_kb': res // 1000, 'rule': rule, 'dip_prominence': prom,
                               'ctcf_min': cmin, **r}
                        grid_rows.append(rec)
                        if default:
                            rows.append(rec)
                    if default:
                        out_ev['call_' + rule] = np.isin(np.arange(len(t)), calls)
                        out_ev['score_' + rule] = score
                        bed = t.loc[calls, ['chrom', 'start', 'end']].assign(name=rule, score=score[calls])
                        bed.to_csv(out_dir('calls', '{}_{}kb_{}.bed'.format(cell, res // 1000, rule)), sep='\t',
                                   header=False, index=False, float_format='%.4g')
                        for ref_name in ('hic_insulation', 'hicFindTADs_50kb'):
                            dom_rows.append({'cell': cell, 'resolution_kb': res // 1000, 'rule': rule,
                                             'reference': ref_name,
                                             **domain_agreement(calls, refs[ref_name], codes, res, rng)})
                if default:
                    for c in ev.columns:
                        out_ev['ev_' + c] = ev[c].values
                    for name, ref in refs.items():
                        out_ev['ref_' + name] = np.isin(np.arange(len(t)), ref)
                    out_ev.to_csv(out_dir('cache', 'evidence_{}_{}kb.tsv.gz'.format(cell, res // 1000)), sep='\t',
                                  index=False, float_format='%.4g')
                    summ += evidence_summary(cell, res, ev, refs, codes)
                    for big, small in COMPARISONS:
                        fa, sa = bd.rule_flags(ev, big)
                        fb, sb = bd.rule_flags(ev, small)
                        ca, cb = bd.calls_from_flags(fa, sa, codes), bd.calls_from_flags(fb, sb, codes)
                        for ref_name in ('hic_insulation', 'hicFindTADs_50kb'):
                            for region in ('all', 'covered'):
                                key = {'cell': cell, 'resolution_kb': res // 1000, 'with': big, 'without': small,
                                       'reference': ref_name, 'region': region}
                                jobs.append((key, sa, sb, ca, cb, refs[ref_name], codes,
                                             covered if region == 'covered' else None, n_boot, len(jobs)))
            print(cell, res, 'done', flush=True)
    pd.DataFrame(rows).to_csv(out_dir('tables', 'eval_rules.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(grid_rows).to_csv(out_dir('tables', 'eval_rules_grid.tsv'), sep='\t', index=False,
                                   float_format='%.4g')
    pd.DataFrame(cover).to_csv(out_dir('tables', 'coverage.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(dom_rows).to_csv(out_dir('tables', 'eval_rules_domains.tsv'), sep='\t', index=False,
                                  float_format='%.4g')
    pd.DataFrame(summ).to_csv(out_dir('tables', 'evidence_summary.tsv'), sep='\t', index=False, float_format='%.4g')
    with Pool(threads) as pool:
        delta = pd.DataFrame(pool.map(bootstrap_delta, jobs))
    delta.to_csv(out_dir('tables', 'eval_rules_delta.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    ev = pd.DataFrame(rows)
    sub = ev[(ev.reference == 'hic_insulation')]
    print(sub.pivot_table(index=['cell', 'resolution_kb', 'region'], columns='rule',
                          values=['precision_tol1', 'recall_tol1', 'f1_tol1']).round(2).to_string())
    print(sub.pivot_table(index=['cell', 'resolution_kb', 'region'], columns='rule', values='auprc').round(3)
          .to_string())
    print(delta[delta.reference == 'hic_insulation'].round(3).to_string())


if __name__ == '__main__':
    main()
