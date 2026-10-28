"""Summary of the TAD-crossing runs (28_tad_crossing.sh).

usage: tad_crossing_summary.py <tad_crossing results dir>

Reads <dir>/<comparison>/<setting>_{summary,promoters,boundaries,genes}.tsv of egaTadCrossing and the egaBoundaries
runs in <dir>/egaboundaries/ and writes to <dir>/summary/:
  settings.tsv        per comparison and setting: edges crossing per condition and status, promoters tested and
                      significant (q < 0.05) for gained / lost, boundaries tested and significant, and the focus
                      promoter (Hoxa2 / EPHA7) and focus boundary (the strongest crossed by its gained contacts)
  top_promoters.tsv   per comparison and setting the 10 best promoters for gained and for lost contacts
  top_boundaries.tsv  per comparison and setting the 10 best boundaries for gained and for lost crossing contacts
  egaboundaries.tsv   per condition and egaBoundaries rule: calls, precision / recall / F1 against the hicFindTADs
                      boundaries of the same condition (tolerance 0, 1, 2 bins), and whether the rule calls a bin
                      within 1 bin of the focus boundary (chr6:51.175 Mb)
"""

import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402

FOCUS = {'Md_vs_PA2': 'Hoxa2', 'hESC_vs_hNEC': 'EPHA7'}
FOCUS_BOUNDARY = ('chr6', 51175000)
RULES = ['dip', 'ctcf', 'ctcf_pattern', 'loop_anchor', 'dip_or_ctcf', 'dip_and_ctcf', 'dip_or_ctcf_pattern',
         'dip_and_ctcf_pattern']


def read(path):
    return pd.read_csv(path, sep='\t', low_memory=False)


def has_gene(genes, gene):
    return genes.fillna('').astype(str).str.split(',').map(lambda g: gene in g)


def settings(root):
    rows, top_p, top_b = [], [], []
    for path in sorted(glob.glob(os.path.join(root, '*_vs_*', '*_promoters.tsv'))):
        comparison = os.path.basename(os.path.dirname(path))
        setting = os.path.basename(path)[:-len('_promoters.tsv')]
        stem = path[:-len('_promoters.tsv')]
        if not os.path.exists(stem + '_boundaries.tsv'):
            continue
        promoters, bounds, summary = read(path), read(stem + '_boundaries.tsv'), read(stem + '_summary.tsv')
        row = dict(comparison=comparison, setting=setting)
        for r in summary.itertuples():
            row['crossing_{}'.format(r.group)] = r.fraction_crossing
            row['edges_{}'.format(r.group)] = r.edges
        for s in ('gained', 'lost'):
            row['promoters_tested_{}'.format(s)] = int(promoters['{}_p'.format(s)].notna().sum())
            row['promoters_q05_{}'.format(s)] = int((promoters['{}_q'.format(s)] < 0.05).sum())
            row['boundaries_tested_{}'.format(s)] = int(bounds['{}_p'.format(s)].notna().sum())
            row['boundaries_q05_{}'.format(s)] = int((bounds['{}_q'.format(s)] < 0.05).sum())
        gene = FOCUS.get(comparison)
        if gene:
            hit = promoters[has_gene(promoters['genes'], gene)]
            if len(hit):
                best = hit.sort_values('gained_p').iloc[0]
                row.update(focus_gene=gene, focus_node=best['node'], focus_gained_crossing=best['gained_crossing'],
                           focus_gained_within=best['gained_within'], focus_lost_crossing=best['lost_crossing'],
                           focus_lost_within=best['lost_within'], focus_shared_crossing=best['shared_crossing'],
                           focus_shared_within=best['shared_within'], focus_gained_odds_ratio=best['gained_odds_ratio'],
                           focus_gained_p=best['gained_p'], focus_gained_q=best['gained_q'],
                           focus_gained_rank=best['gained_rank'])
                genes_path = stem + '_genes.tsv'
                if os.path.exists(genes_path):
                    genes = read(genes_path).set_index('node')
                    if best['node'] in genes.index:
                        g = genes.loc[best['node']]
                        row['focus_sentence'] = g['sentence']
                        row['focus_boundary'] = g.get('gained_boundary', '')
            fb = row.get('focus_boundary')
            if isinstance(fb, str) and fb in set(bounds['id']):
                b = bounds.set_index('id').loc[fb]
                row.update(focus_boundary_strength_rank=b['strength_rank'], focus_boundary_crossing=b['crossing'],
                           focus_boundary_gained=b.get('gained_crossing'),
                           focus_boundary_gained_expected=b.get('gained_expected'),
                           focus_boundary_lost=b.get('lost_crossing'),
                           focus_boundary_lost_expected=b.get('lost_expected'),
                           focus_boundary_gained_p=b['gained_p'], focus_boundary_gained_q=b['gained_q'],
                           focus_boundary_gained_rank=b['gained_rank'])
        rows.append(row)
        for s in ('gained', 'lost'):
            t = promoters.dropna(subset=['{}_p'.format(s)]).sort_values('{}_p'.format(s)).head(10)
            top_p.append(t.assign(comparison=comparison, setting=setting, direction=s))
            t = bounds.dropna(subset=['{}_p'.format(s)]).sort_values('{}_p'.format(s)).head(10)
            top_b.append(t.assign(comparison=comparison, setting=setting, direction=s))
    first = ['comparison', 'setting', 'direction']

    def order(frames):
        if not frames:
            return pd.DataFrame()
        df = pd.concat(frames, ignore_index=True)
        return df[first + [c for c in df.columns if c not in first]]
    return pd.DataFrame(rows), order(top_p), order(top_b)


def egaboundaries(root):
    rows = []
    for path in sorted(glob.glob(os.path.join(root, 'egaboundaries', '*_ega_bins.tsv'))):
        cond = os.path.basename(path).split('_ega_bins.tsv')[0]
        ref_path = os.path.join(root, 'egaboundaries', 'inputs', '{}_25kb_boundaries_chr.bed'.format(cond))
        if not os.path.exists(ref_path):
            continue
        bins = read(path)
        codes = bd.chrom_codes(bins)
        ref = pd.read_csv(ref_path, sep='\t', header=None, usecols=[0, 1, 2])
        ref = ref[ref[0].isin(set(bins.chrom))]
        positions = {c: ((g[1] + g[2]) // 2).values for c, g in ref.groupby(0)}
        rb = bd.positions_to_bins(bins, positions)
        focus = np.flatnonzero((bins.chrom.values == FOCUS_BOUNDARY[0]) &
                               (bins.start.values <= FOCUS_BOUNDARY[1]) & (bins.end.values > FOCUS_BOUNDARY[1]))
        covered = bins['covered'].astype(str).str.lower().eq('true').values
        for rule in RULES:
            col = 'call_' + rule
            if col not in bins:
                continue
            calls = np.flatnonzero(bins[col].astype(str).str.lower().eq('true').values)
            near = bool(len(focus) and np.any(np.abs(calls - focus[0]) <= 1))
            for tol in (0, 1, 2):
                m = bd.match_boundaries(calls, rb, codes, tolerance=tol)
                mc = bd.match_boundaries(calls, rb, codes, tolerance=tol, mask=covered)
                rows.append(dict(condition=cond, rule=rule, tolerance_bins=tol, calls=len(calls),
                                 reference=m['n_reference'], precision=m['precision'], recall=m['recall'],
                                 f1=m['f1'], precision_covered=mc['precision'], recall_covered=mc['recall'],
                                 f1_covered=mc['f1'], calls_focus_boundary=near))
        if len(focus):
            f = bins.iloc[max(focus[0] - 1, 0):focus[0] + 2]
            for r in f.itertuples():
                rows.append(dict(condition=cond, rule='focus_bin', tolerance_bins=np.nan,
                                 calls=np.nan, reference=np.nan,
                                 note='{}:{}-{} covered={} capture_insulation={} ctcf_n={} rev_left={} '
                                      'fwd_right={}'.format(r.chrom, r.start, r.end, r.covered,
                                                            getattr(r, 'capture_insulation', np.nan), r.ctcf_n,
                                                            r.rev_left, r.fwd_right)))
    return pd.DataFrame(rows)


def main(root):
    out = os.path.join(root, 'summary')
    os.makedirs(out, exist_ok=True)
    table, top_p, top_b = settings(root)
    table.to_csv(os.path.join(out, 'settings.tsv'), sep='\t', index=False)
    top_p.to_csv(os.path.join(out, 'top_promoters.tsv'), sep='\t', index=False)
    top_b.to_csv(os.path.join(out, 'top_boundaries.tsv'), sep='\t', index=False)
    ega = egaboundaries(root)
    ega.to_csv(os.path.join(out, 'egaboundaries.tsv'), sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40, 'display.max_colwidth', 400):
        cols = [c for c in ('comparison', 'setting', 'promoters_tested_gained', 'promoters_q05_gained',
                            'promoters_q05_lost', 'focus_gene', 'focus_gained_crossing', 'focus_gained_p',
                            'focus_gained_q', 'focus_gained_rank', 'focus_boundary', 'focus_boundary_gained',
                            'focus_boundary_gained_expected', 'focus_boundary_gained_q', 'focus_boundary_gained_rank')
                if c in table]
        print(table[cols].to_string(index=False))
        if len(ega):
            print(ega[ega.tolerance_bins == 1][['condition', 'rule', 'calls', 'precision', 'recall', 'f1',
                                                'calls_focus_boundary']].to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1])
