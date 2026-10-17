"""Figures of the boundary analysis (27_boundaries.sh); one PDF per panel.

  python figure_boundaries.py select <results dir>   prints the example region (GM12878, 25 kb): the 2 Mb window
                                                      (80 bins) with the highest capture coverage that holds >= 3
                                                      Hi-C insulation boundaries (deterministic, not chosen by result)
  python figure_boundaries.py plot <results dir> [region npz]
Panels (results/boundaries/figures): example_region.pdf (Hi-C triangle, capture insulation, cohesin-bound CTCF
sites with motif orientation, calls of the combined model trained on K562 and Hi-C boundaries),
auprc_<res>kb.pdf and f1_<res>kb.pdf (cross-cell and chromosome split per method, reference Hi-C insulation),
delta_auprc_combined_vs_ctcf[_density].pdf (bootstrap intervals), coverage.pdf, profile_capture.pdf and
profile_orientation.pdf (capture insulation and forward / reverse CTCF sites around Hi-C boundaries, 25 kb),
loop_orientation.pdf, timecourse_changes.pdf.
"""

import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402

CELL_COLOUR = {'GM12878': '#2a78d6', 'K562': '#eb6834', 'H1': '#1baf7a'}
METHOD_COLOUR = {'density': '#c9c8c3', 'capture': '#E69F00', 'ctcf_cohesin': '#56B4E9', 'ctcf_orientation': '#0072B2',
                 'ctcf_orientation_density': '#1c3f75', 'combined': '#CC79A7', 'capture_minima': '#8a8984'}
METHOD_LABEL = {'density': 'bait density (control)', 'capture': 'bait density + capture insulation',
                'ctcf_cohesin': 'CTCF + cohesin', 'ctcf_orientation': 'CTCF + cohesin + orientation',
                'ctcf_orientation_density': 'CTCF + cohesin + orientation + bait density',
                'combined': 'combined (all)', 'capture_minima': 'capture minima (unsupervised)'}
GREY = '#8a8984'
DARK = '#3d3c39'
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'axes.edgecolor': DARK, 'axes.labelcolor': DARK, 'xtick.color': DARK, 'ytick.color': DARK})
CELLS = ['GM12878', 'K562', 'H1']
METHODS = ['density', 'capture', 'ctcf_cohesin', 'ctcf_orientation', 'ctcf_orientation_density', 'combined']


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name), bbox_inches='tight')
    plt.close(fig)


def select_region(base):
    t = pd.read_csv(os.path.join(base, 'cache', 'predictions_GM12878_25kb.tsv.gz'), sep='\t')
    best = None
    for chrom, g in t.groupby('chrom', sort=False):
        cov = g.covered.values
        ref = g.ref__hic_insulation.values
        for s in range(0, len(g) - 80, 10):
            c = cov[s:s + 80].mean()
            if ref[s:s + 80].sum() >= 3 and (best is None or c > best[0]):
                best = (c, chrom, int(g.start.values[s]), int(g.start.values[s]) + 80 * 25000)
    return best


def example_region(base, npz, out):
    z = np.load(npz)
    chrom, start, end, res = str(z['chrom']), int(z['start']), int(z['end']), int(z['res'])
    m = z['matrix']
    t = pd.read_csv(os.path.join(base, 'cache', 'predictions_GM12878_25kb.tsv.gz'), sep='\t')
    g = t[(t.chrom == chrom) & (t.start >= start) & (t.start < end)]
    sites = pd.read_csv(os.path.join(base, 'cache', 'sites_GM12878.tsv.gz'), sep='\t', keep_default_na=False,
                        na_values=['nan'])
    s = sites[(sites.chrom == chrom) & (sites.summit >= start) & (sites.summit < end) & sites.cohesin_bound]
    fig, axes = plt.subplots(4, 1, figsize=(5.0, 5.4), sharex=True,
                             gridspec_kw={'height_ratios': [2.2, 1, 0.8, 0.6], 'hspace': 0.12})
    # Hi-C triangle (log balanced contacts, rotated)
    ax = axes[0]
    n = m.shape[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        lm = np.log10(m)
    vmin, vmax = np.nanpercentile(lm[np.isfinite(lm)], [5, 99])
    x = np.arange(n + 1) * res + start
    xx, yy = np.meshgrid(x, x)
    xr = (xx + yy) / 2
    yr = (yy - xx) / 2
    ax.pcolormesh(xr, yr, lm, cmap='Reds', vmin=vmin, vmax=vmax, shading='flat',
                  rasterized=True)
    ax.set_ylim(0, 1_000_000)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_ylabel('Hi-C\n(GM12878)')
    for b in g.start.values[g.ref__hic_insulation.values.astype(bool)] + res // 2:
        ax.axvline(b, color=DARK, lw=0.5, ls=':')
    # capture insulation
    ax = axes[1]
    mid = g.start.values + res // 2
    ax.axhline(0, color=GREY, lw=0.5)
    ax.plot(mid, g.capture_insulation.values, color=METHOD_COLOUR['capture'], lw=1.2)
    ax.plot(mid, g.hic_insulation.values, color=DARK, lw=0.8, alpha=0.6)
    ax.text(0.01, 0.95, 'capture insulation', color=DARK, transform=ax.transAxes, va='top', fontsize=7)
    ax.text(0.01, 0.05, 'Hi-C insulation (grey)', color=DARK, transform=ax.transAxes, va='bottom', fontsize=7)
    ax.set_ylabel('log2 obs / exp')
    # CTCF orientation arrows
    ax = axes[2]
    for _, r in s.iterrows():
        h = np.log2(1 + r.ctcf)
        if r.orientation == '+':
            ax.annotate('', xy=(r.summit + 25000, h), xytext=(r.summit, h),
                        arrowprops=dict(arrowstyle='-|>', color='#0072B2', lw=0.8, mutation_scale=6))
        elif r.orientation == '-':
            ax.annotate('', xy=(r.summit - 25000, h), xytext=(r.summit, h),
                        arrowprops=dict(arrowstyle='-|>', color='#D55E00', lw=0.8, mutation_scale=6))
        else:
            ax.plot([r.summit], [h], 'o', ms=2.5, color=GREY)
    ax.set_ylim(0, max(1, np.log2(1 + s.ctcf.max()) * 1.15) if len(s) else 1)
    ax.set_ylabel('CTCF log2\n(1 + enr.)')
    ax.text(0.01, 0.95, 'forward (blue, right) / reverse (red, left)', transform=ax.transAxes, va='top', fontsize=7)
    # calls
    ax = axes[3]
    for y, col, lab, colour in ((1, 'ref__hic_insulation', 'Hi-C insulation', DARK),
                                (0.5, 'ref__hicFindTADs_50kb', 'hicFindTADs', GREY),
                                (0, 'calls__combined__K562', 'combined (K562 model)', METHOD_COLOUR['combined'])):
        if col in g:
            pos = g.start.values[g[col].values.astype(bool)] + res // 2
            ax.plot(pos, np.full(len(pos), y), '|', ms=8, mew=1.5, color=colour)
            ax.text(end + 20000, y, lab, va='center', fontsize=7)
    ax.set_ylim(-0.4, 1.4)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlim(start, end)
    ax.set_xlabel('{} (Mb, GRCh38)'.format(chrom))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '{:.1f}'.format(v / 1e6)))
    save(fig, out, 'example_region.pdf')


def method_points(ev, metric, res, out, name, ylabel):
    sub = ev[(ev.reference == 'hic_insulation') & (ev.resolution_kb == res) & (ev.region == 'all')]
    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    designs = [('cross_cell', 'cross-cell'), ('chromosome_split', 'chromosome split')]
    xpos = 0
    ticks, labels = [], []
    for test in CELLS:
        for design, dl in designs:
            d = sub[(sub.test == test) & (sub.design == design)]
            for k, method in enumerate(METHODS):
                v = d[d.method == method][metric].values
                ax.scatter(np.full(len(v), xpos + (k - 2.5) * 0.13), v, s=12, color=METHOD_COLOUR[method],
                           edgecolor='white', linewidth=0.5, zorder=3, label=METHOD_LABEL[method]
                           if (test == 'GM12878' and design == 'cross_cell') else None)
            ticks.append(xpos)
            labels.append('{}\n{}'.format(test, dl))
            xpos += 1
        u = ev[(ev.reference == 'hic_insulation') & (ev.resolution_kb == res) & (ev.region == 'all') &
               (ev.test == test) & (ev.method == 'capture_minima')]
        if metric in u and len(u):
            ax.scatter([xpos - 1.5], u[metric].values, marker='x', s=14, color=METHOD_COLOUR['capture_minima'],
                       label=METHOD_LABEL['capture_minima'] if test == 'GM12878' else None, zorder=3)
        if metric == 'auprc':
            pr = sub[(sub.test == test)].positive_rate.mean()
            ax.plot([xpos - 2.4, xpos - 0.6], [pr, pr], color=GREY, lw=0.7, ls='--')
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=6.5, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    ax.set_title('{} kb, reference: Hi-C insulation boundaries'.format(res), fontsize=8)
    save(fig, out, name)


def delta_plot(delta, out, full, reduced, name):
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    d = delta[(delta.reference == 'hic_insulation') & (delta.full == full) &
              (delta.reduced == reduced)].reset_index(drop=True)
    labels = []
    for i, r in d.iterrows():
        colour = '#CC79A7' if r.region == 'all' else '#7a3b63'
        ax.errorbar(r.delta_auprc, i, xerr=[[r.delta_auprc - r.delta_auprc_lo], [r.delta_auprc_hi - r.delta_auprc]],
                    fmt='o', ms=3, color=colour, elinewidth=0.8, capsize=0)
        tr = 'split' if r.design == 'chromosome_split' else r.train
        labels.append('{} kb  {} -> {}  {}'.format(int(r.resolution_kb), tr, r.test, r.region))
    ax.axvline(0, color=GREY, lw=0.7)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel('AUPRC {} - AUPRC {}\n(95 % chromosome bootstrap interval)'.format(
        METHOD_LABEL[full], METHOD_LABEL[reduced]),
                  fontsize=7)
    save(fig, out, name)


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
        t = pd.read_csv(os.path.join(base, 'cache', 'features_{}_{}kb.tsv.gz'.format(cell, res // 1000)), sep='\t')
        codes = bd.chrom_codes(t)
        ref = bd.local_minima(t.hic_insulation.values, codes, 0.3)
        ins = t.capture_insulation.values
        fwd, rev = np.expm1(t.fwd_n.values), np.expm1(t.rev_n.values)
        prof, pf, pr = [], [], []
        for o in offs:
            k = ref + o
            ok = (k >= 0) & (k < len(t))
            k, b = k[ok], ref[ok]
            k, b = k[codes[k] == codes[b]], b[codes[k] == codes[b]]
            prof.append(np.nanmean(ins[k]))
            pf.append(fwd[k].mean())
            pr.append(rev[k].mean())
        ax1.plot(offs * res / 1000, prof, color=CELL_COLOUR[cell], lw=1.2, label=cell)
        ax2.plot(offs * res / 1000, pf, color=CELL_COLOUR[cell], lw=1.2, label=cell + ' forward')
        ax2.plot(offs * res / 1000, pr, color=CELL_COLOUR[cell], lw=1.2, ls='--', label=cell + ' reverse')
    ax1.axvline(0, color=GREY, lw=0.5)
    ax1.set_xlabel('distance to Hi-C boundary (kb)')
    ax1.set_ylabel('mean capture insulation\n(covered bins, log2)')
    ax1.legend(frameon=False, fontsize=7)
    ax2.axvline(0, color=GREY, lw=0.5)
    ax2.set_xlabel('distance to Hi-C boundary (kb)')
    ax2.set_ylabel('cohesin-bound CTCF sites\nper {} kb bin'.format(res // 1000))
    ax2.legend(frameon=False, fontsize=6, ncol=1, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    save(fig1, out, 'profile_capture.pdf')
    save(fig2, out, 'profile_orientation.pdf')


def loop_plot(lo, out):
    sets = [s for s in lo.set.unique()]
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    classes = ['convergent', 'tandem_forward', 'tandem_reverse', 'divergent']
    colours = ['#CC79A7', '#56B4E9', '#0072B2', GREY]
    left = np.zeros(len(sets))
    for c, col in zip(classes, colours):
        v = np.array([lo[(lo.set == s) & (lo['class'] == c)].fraction.values[0] for s in sets])
        ax.barh(range(len(sets)), v, left=left, color=col, label=c.replace('_', ' '), height=0.7,
                edgecolor='white', linewidth=0.8)
        left += v
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels(sets)
    ax.set_xlabel('fraction of edges with oriented cohesin-bound CTCF at both ends')
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
    from scipy.stats import spearmanr
    rho = spearmanr(dh[ok], dc[ok])[0]
    ax.text(0.03, 0.97, 'Spearman {:.2f}, n = {}'.format(rho, int(ok.sum())), transform=ax.transAxes, va='top',
            fontsize=7)
    save(fig, out, 'timecourse_insulation_change.pdf')
    ch = pd.read_csv(os.path.join(base, 'tables', 'timecourse_changes.tsv'), sep='\t')
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    sub = ch.reset_index(drop=True)
    labels = []
    for i, r in sub.iterrows():
        ax.barh(i, r.fraction_hic_changes_matched, color=METHOD_COLOUR.get(r.caller.replace('_encode', ''), GREY),
                height=0.7)
        ax.plot([r.expected_shifted], [i], '|', color=DARK, ms=7)
        labels.append('{} {} {}'.format(r.pair.replace('_', '-'), r.change, r.caller.replace('_encode', '')))
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel('fraction of Hi-C boundary changes\nmatched by caller changes (| = shifted)')
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
    ev = pd.read_csv(os.path.join(base, 'tables', 'eval_boundaries.tsv'), sep='\t')
    for res in (10, 25, 50):
        method_points(ev, 'auprc', res, out, 'auprc_{}kb.pdf'.format(res), 'AUPRC (boundary bins)')
        method_points(ev, 'f1_tol1', res, out, 'f1_{}kb.pdf'.format(res), 'F1 (tolerance 1 bin)')
    delta = pd.read_csv(os.path.join(base, 'tables', 'eval_delta.tsv'), sep='\t')
    delta_plot(delta, out, 'combined', 'ctcf_orientation', 'delta_auprc_combined_vs_ctcf.pdf')
    delta_plot(delta, out, 'combined', 'ctcf_orientation_density', 'delta_auprc_combined_vs_ctcf_density.pdf')
    coverage_plot(pd.read_csv(os.path.join(base, 'tables', 'coverage.tsv'), sep='\t'), out)
    profiles(base, out)
    loop_plot(pd.read_csv(os.path.join(base, 'tables', 'loop_orientation.tsv'), sep='\t'), out)
    if os.path.exists(os.path.join(base, 'cache', 'timecourse_bins_50kb.tsv.gz')):
        timecourse_plot(base, out)


if __name__ == '__main__':
    main()
