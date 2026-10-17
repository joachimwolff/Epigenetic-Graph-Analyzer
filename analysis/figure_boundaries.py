"""Figures of the boundary analysis (27_boundaries.sh); one PDF per panel.

  python figure_boundaries.py select <results dir>   prints the example region (GM12878, 25 kb): the 2 Mb window
                                                      (80 bins) with the highest capture coverage that holds >= 3
                                                      Hi-C insulation boundaries (deterministic, not chosen by result)
  python figure_boundaries.py plot <results dir> [region npz]
Panels (results/boundaries/figures): example_region.pdf (Hi-C triangle, capture insulation with dips, cohesin-bound
CTCF sites with motif orientation, calls of the rules and Hi-C boundaries), precision_recall_<res>kb_<region>.pdf
(every rule and cell type against the Hi-C insulation boundaries), auprc_<res>kb.pdf (rules with a score),
delta_dip.pdf (change of precision and recall when the dip is added to CTCF, bootstrap intervals), coverage.pdf,
profile_capture.pdf and profile_orientation.pdf (capture insulation and forward / reverse cohesin-bound CTCF sites
around Hi-C boundaries, 25 kb), loop_orientation.pdf, timecourse_insulation_change.pdf, timecourse_changes.pdf.
"""

import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402

CELL_COLOUR = {'GM12878': '#2a78d6', 'K562': '#eb6834', 'H1': '#1baf7a'}
CELL_MARKER = {'GM12878': 'o', 'K562': 's', 'H1': '^'}
RULE_COLOUR = {'dip': '#E69F00', 'ctcf': '#56B4E9', 'ctcf_pattern': '#0072B2', 'loop_anchor': '#009E73',
               'dip_or_ctcf': '#9ecae1', 'dip_and_ctcf': '#CC79A7', 'dip_or_ctcf_pattern': '#D55E00',
               'dip_and_ctcf_pattern': '#7a3b63'}
RULE_LABEL = {'dip': 'dip', 'ctcf': 'CTCF / cohesin', 'ctcf_pattern': 'CTCF pattern', 'loop_anchor': 'loop anchor',
              'dip_or_ctcf': 'dip OR CTCF', 'dip_and_ctcf': 'dip AND CTCF', 'dip_or_ctcf_pattern': 'dip OR pattern',
              'dip_and_ctcf_pattern': 'dip AND pattern'}
GREY = '#8a8984'
DARK = '#3d3c39'
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'axes.edgecolor': DARK, 'axes.labelcolor': DARK, 'xtick.color': DARK, 'ytick.color': DARK})
CELLS = ['GM12878', 'K562', 'H1']
RULES = list(bd.RULES)


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name), bbox_inches='tight')
    plt.close(fig)


def select_region(base):
    t = pd.read_csv(os.path.join(base, 'cache', 'evidence_GM12878_25kb.tsv.gz'), sep='\t')
    best = None
    for chrom, g in t.groupby('chrom', sort=False):
        cov = g.ev_covered.values.astype(float)
        ref = g.ref_hic_insulation.values
        for s in range(0, len(g) - 80, 10):
            c = cov[s:s + 80].mean()
            if ref[s:s + 80].sum() >= 3 and (best is None or c > best[0]):
                best = (c, chrom, int(g.start.values[s]), int(g.start.values[s]) + 80 * 25000)
    return best


def example_region(base, npz, out):
    z = np.load(npz)
    chrom, start, end, res = str(z['chrom']), int(z['start']), int(z['end']), int(z['res'])
    m = z['matrix']
    t = pd.read_csv(os.path.join(base, 'cache', 'evidence_GM12878_25kb.tsv.gz'), sep='\t')
    g = t[(t.chrom == chrom) & (t.start >= start) & (t.start < end)]
    sites = pd.read_csv(os.path.join(base, 'cache', 'sites_GM12878.tsv.gz'), sep='\t', keep_default_na=False,
                        na_values=['nan'])
    s = sites[(sites.chrom == chrom) & (sites.summit >= start) & (sites.summit < end) & sites.cohesin_bound &
              (sites.ctcf >= 10)]
    fig, axes = plt.subplots(4, 1, figsize=(5.0, 6.0), sharex=True,
                             gridspec_kw={'height_ratios': [2.2, 1, 0.8, 1.1], 'hspace': 0.12})
    ax = axes[0]
    n = m.shape[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        lm = np.log10(m)
    vmin, vmax = np.nanpercentile(lm[np.isfinite(lm)], [5, 99])
    x = np.arange(n + 1) * res + start
    xx, yy = np.meshgrid(x, x)
    ax.pcolormesh((xx + yy) / 2, (yy - xx) / 2, lm, cmap='Reds', vmin=vmin, vmax=vmax, shading='flat',
                  rasterized=True)
    ax.set_ylim(0, 1_000_000)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_ylabel('Hi-C\n(GM12878)')
    for b in g.start.values[g.ref_hic_insulation.values.astype(bool)] + res // 2:
        ax.axvline(b, color=DARK, lw=0.5, ls=':')
    ax = axes[1]
    mid = g.start.values + res // 2
    ax.axhline(0, color=GREY, lw=0.5)
    ax.plot(mid, g.hic_insulation.values, color=GREY, lw=0.8)
    ax.plot(mid, g.ev_capture_insulation.values, color=RULE_COLOUR['dip'], lw=1.2)
    d = g.ev_dip.values.astype(bool)
    ax.plot(mid[d], g.ev_capture_insulation.values[d], 'v', color=RULE_COLOUR['dip'], ms=4)
    ax.text(0.01, 0.97, 'capture insulation (dips marked); Hi-C insulation in grey', transform=ax.transAxes,
            va='top', fontsize=6.5)
    ax.set_ylabel('log2 obs / exp')
    ax = axes[2]
    for _, r in s.iterrows():
        h = np.log2(1 + r.ctcf)
        if r.orientation in ('+', '-'):
            dx = 25000 if r.orientation == '+' else -25000
            colour = '#0072B2' if r.orientation == '+' else '#D55E00'
            ax.annotate('', xy=(r.summit + dx, h), xytext=(r.summit, h),
                        arrowprops=dict(arrowstyle='-|>', color=colour, lw=0.8, mutation_scale=6))
        else:
            ax.plot([r.summit], [h], 'o', ms=2.5, color=GREY)
    ax.set_ylim(0, max(1, np.log2(1 + s.ctcf.max()) * 1.2) if len(s) else 1)
    ax.set_ylabel('CTCF log2\n(1 + enr.)')
    ax.text(0.01, 0.97, 'cohesin-bound CTCF: forward (blue, right), reverse (red, left), no / both motifs (dot)',
            transform=ax.transAxes, va='top', fontsize=6.5)
    ax = axes[3]
    tracks = [('ref_hic_insulation', 'Hi-C insulation', DARK), ('ref_hicFindTADs_50kb', 'hicFindTADs', GREY),
              ('call_dip', RULE_LABEL['dip'], RULE_COLOUR['dip']),
              ('call_ctcf_pattern', RULE_LABEL['ctcf_pattern'], RULE_COLOUR['ctcf_pattern']),
              ('call_dip_and_ctcf_pattern', RULE_LABEL['dip_and_ctcf_pattern'], RULE_COLOUR['dip_and_ctcf_pattern'])]
    for y, (col, lab, colour) in enumerate(tracks[::-1]):
        pos = g.start.values[g[col].values.astype(bool)] + res // 2
        ax.plot(pos, np.full(len(pos), y), '|', ms=8, mew=1.5, color=colour)
        ax.text(end + 20000, y, lab, va='center', fontsize=7)
    ax.set_ylim(-0.6, len(tracks) - 0.4)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(start, end)
    ax.set_xlabel('{} (Mb, GRCh38)'.format(chrom))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '{:.1f}'.format(v / 1e6)))
    save(fig, out, 'example_region.pdf')


def precision_recall(ev, res, region, out):
    sub = ev[(ev.reference == 'hic_insulation') & (ev.resolution_kb == res) & (ev.region == region)]
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    for f1 in (0.2, 0.4, 0.6):
        r = np.linspace(f1 / 2 + 1e-3, 1, 100)
        p = f1 * r / (2 * r - f1)
        ok = (p > 0) & (p <= 1)
        ax.plot(r[ok], p[ok], color=GREY, lw=0.4, ls=':')
        ax.text(1.0, p[ok][-1], 'F1 {}'.format(f1), fontsize=6, color=GREY, va='center')
    for rule in RULES:
        for cell in CELLS:
            g = sub[(sub.rule == rule) & (sub.cell == cell)]
            ax.scatter(g.recall_tol1, g.precision_tol1, color=RULE_COLOUR[rule], marker=CELL_MARKER[cell], s=18,
                       edgecolor='white', linewidth=0.5, zorder=3)
    hs = [plt.Line2D([], [], color=RULE_COLOUR[r], marker='o', ls='', label=RULE_LABEL[r]) for r in RULES]
    hs += [plt.Line2D([], [], color=DARK, marker=CELL_MARKER[c], ls='', label=c) for c in CELLS]
    ax.legend(handles=hs, frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(1.08, 1.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('recall (Hi-C boundaries within 1 bin of a call)')
    ax.set_ylabel('precision (calls within 1 bin of a Hi-C boundary)')
    ax.set_title('{} kb, {} bins'.format(res, region), fontsize=8)
    save(fig, out, 'precision_recall_{}kb_{}.pdf'.format(res, region))


def auprc_plot(ev, res, out):
    rules = ['dip', 'ctcf', 'ctcf_pattern', 'loop_anchor', 'dip_and_ctcf', 'dip_and_ctcf_pattern']
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    for k, cell in enumerate(CELLS):
        for region, alpha in (('all', 1.0), ('covered', 0.45)):
            g = ev[(ev.reference == 'hic_insulation') & (ev.resolution_kb == res) & (ev.region == region) &
                   (ev.cell == cell)].set_index('rule')
            xs = np.arange(len(rules)) + (k - 1) * 0.26 + (0.06 if region == 'covered' else -0.06)
            ax.scatter(xs, g.loc[rules, 'auprc'], color=CELL_COLOUR[cell], marker=CELL_MARKER[cell], s=14,
                       alpha=alpha, label='{} ({} bins)'.format(cell, region))
            ax.plot([xs.min() - 0.05, xs.max() + 0.05], [g.positive_rate.iloc[0]] * 2, color=CELL_COLOUR[cell],
                    lw=0.5, ls=':', alpha=alpha)
    ax.set_xticks(range(len(rules)))
    ax.set_xticklabels([RULE_LABEL[r] for r in rules], rotation=30, ha='right', fontsize=6.5)
    ax.set_ylabel('AUPRC (dotted: positive rate)')
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=5.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    ax.set_title('{} kb, reference: Hi-C insulation boundaries'.format(res), fontsize=8)
    save(fig, out, 'auprc_{}kb.pdf'.format(res))


def delta_plot(delta, out):
    d = delta[(delta.reference == 'hic_insulation') &
              delta['with'].isin(['dip_and_ctcf', 'dip_or_ctcf', 'dip_and_ctcf_pattern', 'dip_or_ctcf_pattern']) &
              delta.without.isin(['ctcf', 'ctcf_pattern'])].reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 4.2), sharey=True)
    labels = ['{} kb {} {} ({})'.format(int(r.resolution_kb), r.cell, RULE_LABEL[r['with']], r.region)
              for _, r in d.iterrows()]
    for ax, metric in zip(axes, ('precision', 'recall')):
        for i, r in d.iterrows():
            ax.errorbar(r['delta_' + metric], i, xerr=[[r['delta_' + metric] - r[metric + '_lo']],
                                                       [r[metric + '_hi'] - r['delta_' + metric]]],
                        fmt='o', ms=2.5, color=RULE_COLOUR[r['with']], elinewidth=0.7, capsize=0)
        ax.axvline(0, color=GREY, lw=0.6)
        ax.set_xlabel('change of {}\nwhen the dip is added'.format(metric))
    axes[0].set_yticks(range(len(d)))
    axes[0].set_yticklabels(labels, fontsize=4.5)
    axes[0].invert_yaxis()
    save(fig, out, 'delta_dip.pdf')


def coverage_plot(cov, out):
    fig, ax = plt.subplots(figsize=(3.0, 2.2))
    for k, cell in enumerate(CELLS):
        c = cov[cov.cell == cell].sort_values('resolution_kb')
        ax.bar(np.arange(len(c)) + (k - 1) * 0.27, c.covered_fraction, width=0.25, color=CELL_COLOUR[cell],
               label=cell)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['10 kb', '25 kb', '50 kb'])
    ax.set_ylabel('fraction of bins with\n>= 5 expected spanning edges')
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7)
    save(fig, out, 'coverage.pdf')


def profiles(base, out, res=25000, flank=40):
    fig1, ax1 = plt.subplots(figsize=(3.0, 2.2))
    fig2, ax2 = plt.subplots(figsize=(3.0, 2.2))
    offs = np.arange(-flank, flank + 1)
    for cell in CELLS:
        t = pd.read_csv(os.path.join(base, 'cache', 'evidence_{}_{}kb.tsv.gz'.format(cell, res // 1000)), sep='\t')
        sites = pd.read_csv(os.path.join(base, 'cache', 'sites_{}.tsv.gz'.format(cell)), sep='\t',
                            keep_default_na=False, na_values=['nan'])
        sites = sites[sites.cohesin_bound & (sites.ctcf >= 10)]
        codes = bd.chrom_codes(t)
        first = {c: int(np.flatnonzero(t.chrom.values == c)[0]) for c in pd.unique(t.chrom)}
        fwd = np.zeros(len(t))
        rev = np.zeros(len(t))
        for o, arr in (('+', fwd), ('-', rev)):
            ss = sites[(sites.orientation == o) & sites.chrom.isin(first)]
            k = np.array([first[c] for c in ss.chrom], dtype=np.int64) + (ss.summit.values // res)
            np.add.at(arr, k[k < len(t)], 1)
        ref = np.flatnonzero(t.ref_hic_insulation.values)
        ins = t.ev_capture_insulation.values
        prof, pf, pr = [], [], []
        for o in offs:
            k = ref + o
            ok = (k >= 0) & (k < len(t))
            k, b = k[ok], ref[ok]
            k = k[codes[k] == codes[b]]
            prof.append(np.nanmean(ins[k]))
            pf.append(fwd[k].mean())
            pr.append(rev[k].mean())
        ax1.plot(offs * res / 1000, prof, color=CELL_COLOUR[cell], lw=1.2, label=cell)
        ax2.plot(offs * res / 1000, pf, color=CELL_COLOUR[cell], lw=1.2, label=cell + ' forward')
        ax2.plot(offs * res / 1000, pr, color=CELL_COLOUR[cell], lw=1.2, ls='--', label=cell + ' reverse')
    for ax in (ax1, ax2):
        ax.axvline(0, color=GREY, lw=0.5)
        ax.set_xlabel('distance to Hi-C boundary (kb)')
    ax1.set_ylabel('mean capture insulation\n(covered bins, log2)')
    ax1.legend(frameon=False, fontsize=7)
    ax2.set_ylabel('cohesin-bound CTCF sites\nper {} kb bin'.format(res // 1000))
    ax2.legend(frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    save(fig1, out, 'profile_capture.pdf')
    save(fig2, out, 'profile_orientation.pdf')


def loop_plot(lo, out):
    sets = list(lo.set.unique())
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    classes = ['convergent', 'tandem_forward', 'tandem_reverse', 'divergent']
    colours = ['#7a3b63', '#56B4E9', '#0072B2', GREY]
    left = np.zeros(len(sets))
    for c, col in zip(classes, colours):
        v = np.array([lo[(lo.set == s) & (lo['class'] == c)].fraction.values[0] for s in sets])
        ax.barh(range(len(sets)), v, left=left, color=col, label=c.replace('_', ' '), height=0.7,
                edgecolor='white', linewidth=0.8)
        left += v
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels(sets)
    ax.set_xlabel('fraction of edges with oriented cohesin-bound\nCTCF at both ends (one orientation per end)')
    ax.axvline(0.25, color=DARK, lw=0.6, ls=':')
    ax.legend(frameon=False, fontsize=6.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    save(fig, out, 'loop_orientation.pdf')


def timecourse_plot(base, out):
    t = pd.read_csv(os.path.join(base, 'cache', 'timecourse_bins_50kb.tsv.gz'), sep='\t')
    dh = t.hic_insulation_D6 - t.hic_insulation_D0
    dc = t.capture_insulation_D6 - t.capture_insulation_D0
    ok = np.isfinite(dh) & np.isfinite(dc)
    fig, ax = plt.subplots(figsize=(2.6, 2.4))
    ax.hexbin(dh[ok], dc[ok], gridsize=40, cmap='Greys', mincnt=1, bins='log', extent=(-1.5, 1.5, -2, 2))
    ax.set_xlabel('Hi-C insulation change\nday 6 - day 0 (log2)')
    ax.set_ylabel('capture insulation change\nday 6 - day 0 (log2)')
    rho = spearmanr(dh[ok], dc[ok])[0]
    ax.text(0.03, 0.97, 'Spearman {:.2f}, n = {}'.format(rho, int(ok.sum())), transform=ax.transAxes, va='top',
            fontsize=7)
    save(fig, out, 'timecourse_insulation_change.pdf')
    ch = pd.read_csv(os.path.join(base, 'tables', 'timecourse_changes.tsv'), sep='\t')
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    labels = []
    for i, r in ch.iterrows():
        ax.barh(i, r.fraction_hic_changes_matched, color=RULE_COLOUR.get(r.rule, GREY), height=0.7)
        ax.plot([r.expected_shifted], [i], '|', color=DARK, ms=7)
        labels.append('{} {} {}'.format(r.pair.replace('_', '-'), r.change, RULE_LABEL.get(r.rule, r.rule)))
    ax.set_yticks(range(len(ch)))
    ax.set_yticklabels(labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel('fraction of Hi-C boundary changes matched\nby changes of the rule calls (| = shifted)')
    save(fig, out, 'timecourse_changes.pdf')


def main():
    mode, results = sys.argv[1], sys.argv[2]
    base = os.path.join(results, 'boundaries')
    if mode == 'select':
        c, chrom, start, end = select_region(base)
        print('{}:{}-{}'.format(chrom, start, end))
        return
    out = os.path.join(base, 'figures')
    npz = sys.argv[3] if len(sys.argv) > 3 else (sorted(glob.glob(os.path.join(base, 'hic', 'hic_region_*.npz')))
                                                 or [None])[0]
    if npz:
        example_region(base, npz, out)
    ev = pd.read_csv(os.path.join(base, 'tables', 'eval_rules.tsv'), sep='\t')
    for res in (10, 25, 50):
        for region in ('all', 'covered'):
            precision_recall(ev, res, region, out)
        auprc_plot(ev, res, out)
    delta_plot(pd.read_csv(os.path.join(base, 'tables', 'eval_rules_delta.tsv'), sep='\t'), out)
    coverage_plot(pd.read_csv(os.path.join(base, 'tables', 'coverage.tsv'), sep='\t'), out)
    profiles(base, out)
    loop_plot(pd.read_csv(os.path.join(base, 'tables', 'loop_orientation.tsv'), sep='\t'), out)
    if os.path.exists(os.path.join(base, 'cache', 'timecourse_bins_50kb.tsv.gz')):
        timecourse_plot(base, out)


if __name__ == '__main__':
    main()
