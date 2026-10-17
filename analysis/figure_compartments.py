"""Figures of the compartment and domain analyses, one PDF per panel.

Reads <results>/compartments_domains/{tables,switch} (compartments_chromatin.py, compartments_capture.py,
domains_graph.py, compartment_switch.py) and writes <results>/compartments_domains/figures/*.pdf.

Usage: python figure_compartments.py
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402

CELL_COLOUR = {'GM12878': '#2a78d6', 'K562': '#eb6834', 'H1': '#1baf7a'}
GREY = '#8a8984'
DARK = '#3d3c39'
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
                     'axes.edgecolor': DARK, 'axes.labelcolor': DARK, 'xtick.color': DARK, 'ytick.color': DARK})

METHOD_LABEL = {
    'marks_pc1': 'Marks PC1', 'marks_role': 'Active - repressive', 'marks_active': 'Active marks',
    'marks_profile_eig': 'Mark-profile eigenvector', 'marks_pc1_smooth1': 'Marks PC1, 300 kb smoothing',
    'single_H3K27ac': 'H3K27ac alone', 'single_H3K4me1': 'H3K4me1 alone', 'gene_density': 'Gene density',
    'gc': 'GC content', 'loco_ridge': 'Supervised (LOCO ridge)', 'access_level': 'Accessibility level',
    'access_fh_eig': 'Accessibility eigenvector (F&H)', 'marks_plus_access_pc1': 'Marks + accessibility PC1',
    'access_level_smooth1': 'Accessibility, 300 kb smoothing', 'loco_ridge_gc': 'Supervised, GC only',
    'loco_ridge_gc_marks': 'Supervised, GC + tracks',
}


def save(fig, name):
    os.makedirs(cc.FIGURES, exist_ok=True)
    fig.savefig(os.path.join(cc.FIGURES, name), bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


def chromatin_bars(summary, metric, label, name):
    sets = {'GM12878': 'GM12878_encode_hg38', 'K562': 'K562_encode_hg38', 'H1': 'H1_encode_hg38'}
    own = summary[summary.set.isin(sets.values()) & (summary.chrom == 'genome')]
    methods = [m for m in METHOD_LABEL if m in set(own.method) and own[own.method == m][metric].notna().any()]
    fig, ax = plt.subplots(figsize=(3.6, 0.32 * len(methods) + 0.6))
    height = 0.26
    for k, (cell, set_name) in enumerate(sets.items()):
        sub = summary[(summary.set == set_name) & (summary.reference == 'hic_{}_100kb'.format(cell)) &
                      (summary.chrom == 'genome')].set_index('method')
        y = np.arange(len(methods)) + (k - 1) * height
        values = [sub[metric].get(m, np.nan) for m in methods]
        ax.barh(y, values, height=height - 0.03, color=CELL_COLOUR[cell], label=cell)
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels([METHOD_LABEL[m] for m in methods])
    ax.invert_yaxis()
    ax.set_xlabel(label)
    if metric in ('sign_agreement', 'auc'):
        ax.axvline(0.5, color=GREY, lw=0.8, ls='--')
        ax.set_xlim(0.4, 1)
    else:
        ax.axvline(0, color=GREY, lw=0.8)
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.4, -0.12 if len(methods) > 12 else -0.2), ncol=3)
    save(fig, name)


def chromatin_crosscell(summary):
    cells = ['GM12878', 'K562', 'H1']
    rows = [('marks_pc1', c) for c in cells] + [('gene_density', 'GM12878')]
    mat = np.full((len(rows), 3), np.nan)
    for i, (method, cell) in enumerate(rows):
        for j, ref in enumerate(cells):
            sub = summary[(summary.set == '{}_encode_hg38'.format(cell)) & (summary.method == method) &
                          (summary.reference == 'hic_{}_100kb'.format(ref)) & (summary.chrom == 'genome')]
            if len(sub):
                mat[i, j] = sub.pearson.iloc[0]
    fig, ax = plt.subplots(figsize=(2.8, 2.3))
    im = ax.imshow(mat, cmap='Blues', vmin=0, vmax=1)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, '{:.2f}'.format(mat[i, j]), ha='center', va='center',
                        color='white' if mat[i, j] > 0.6 else DARK)
    ax.set_xticks(range(3))
    ax.set_xticklabels(cells)
    ax.set_xlabel('Hi-C eigenvector of')
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(['{} marks PC1'.format(c) for c in cells] + ['Gene density'])
    fig.colorbar(im, ax=ax, fraction=0.05, label='Pearson r (genome)')
    save(fig, 'chromatin_crosscell.pdf')


def chromatin_sets(summary):
    rows = [('GM12878_encode_hg38', 'hic_GM12878_100kb', 'GM12878 ENCODE, 12 tracks'),
            ('GM12878_4marks_hg38', 'hic_GM12878_100kb', 'GM12878 ENCODE, 4 marks'),
            ('GM12878_ega_hg19', 'hic_GM12878_100kb', 'GM12878 EGA input (hg19), 5 tracks'),
            ('GM12878_ega_hg19', 'rao2014_subcompartments', 'GM12878 EGA input vs Rao 2014'),
            ('K562_encode_hg38', 'hic_K562_100kb', 'K562 ENCODE, 12 tracks'),
            ('K562_4marks_hg38', 'hic_K562_100kb', 'K562 ENCODE, 4 marks'),
            ('H1_encode_hg38', 'hic_H1_100kb', 'H1 ENCODE, 12 tracks'),
            ('H1_4marks_hg38', 'hic_H1_100kb', 'H1 ENCODE, 4 marks'),
            ('H1_3marks_hg38', 'hic_H1_100kb', 'H1 ENCODE, 3 marks'),
            ('E003_H1_hg19', 'hic_H1_100kb', 'H1 Roadmap E003 (hg19)'),
            ('E008_H9_hg19', 'hic_H1_100kb', 'H9 Roadmap E008 vs H1'),
            ('E008_H9_3marks_hg19', 'hic_H1_100kb', 'H9 E008, 3 marks vs H1'),
            ('E003_H1_hg19', '4dn_H1_250kb', 'H1 E003 vs 4DN published')]
    values = []
    for set_name, ref, label in rows:
        sub = summary[(summary.set == set_name) & (summary.reference == ref) & (summary.method == 'marks_pc1') &
                      (summary.chrom == 'genome')]
        values.append((label, sub.auc.iloc[0] if len(sub) else np.nan, sub.pearson.iloc[0] if len(sub) else np.nan,
                       set_name.split('_')[0]))
    fig, ax = plt.subplots(figsize=(3.4, 0.25 * len(values) + 0.6))
    for i, (label, auc, r, cell) in enumerate(values):
        colour = CELL_COLOUR.get({'E003': 'H1', 'E008': 'H1'}.get(cell, cell), GREY)
        ax.barh(i, auc, height=0.7, color=colour)
        if np.isfinite(auc):
            ax.text(auc + 0.005, i, 'AUC {:.2f}, r {:.2f}'.format(auc, r), va='center', fontsize=6.5, color=DARK)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels([v[0] for v in values])
    ax.invert_yaxis()
    ax.set_xlim(0.5, 1.05)
    ax.set_xlabel('AUC of marks PC1 for Hi-C compartment A')
    save(fig, 'chromatin_sets_auc.pdf')


def cell_difference_panel():
    path = os.path.join(cc.TABLES, 'cell_differences.tsv')
    if not os.path.exists(path):
        return
    t = pd.read_csv(path, sep='\t')
    methods = [m for m in ('marks_pc1', 'marks_pc1_smooth1', 'access_level', 'marks_plus_access_pc1')
               if m in set(t.method)]
    pairs = list(dict.fromkeys(zip(t.cell_A, t.cell_B)))
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.2))
    width = 0.8 / len(methods)
    shades = ['#3d3c39', '#8a8984', '#2a78d6', '#7fb0ea']
    for k, m in enumerate(methods):
        sub = t[t.method == m].set_index(['cell_A', 'cell_B'])
        x = np.arange(len(pairs)) + (k - (len(methods) - 1) / 2) * width
        axes[0].bar(x, [sub.pearson_delta.get(p, np.nan) for p in pairs], width=width * 0.9, color=shades[k],
                    label=METHOD_LABEL.get(m, m))
        axes[1].bar(x, [sub.precision.get(p, np.nan) for p in pairs], width=width * 0.9, color=shades[k])
    exp = t.drop_duplicates(['cell_A', 'cell_B']).set_index(['cell_A', 'cell_B']).expected_precision
    axes[1].scatter(np.arange(len(pairs)), [exp.get(p) for p in pairs], marker='_', s=300, color='#c8443f', zorder=3,
                    label='random switches')
    short = {'GM12878': 'GM', 'K562': 'K562', 'H1': 'H1'}
    for ax in axes:
        ax.set_xticks(np.arange(len(pairs)))
        ax.set_xticklabels(['{} vs {}'.format(short[a], short[b]) for a, b in pairs])
    axes[0].set_ylabel('Pearson r, estimate vs Hi-C difference')
    axes[1].set_ylabel('Switches confirmed by Hi-C')
    axes[0].set_ylim(0, 0.65)
    axes[1].set_ylim(0, 0.6)
    axes[0].legend(frameon=False, fontsize=6, loc='upper left')
    axes[1].legend(frameon=False, fontsize=6, loc='upper right')
    save(fig, 'cell_differences.pdf')


def example_track(cell='H1', set_name='H1_encode_hg38', graph='H1_jung', chrom='chr2', lo=0, hi=120000000):
    path = os.path.join(cc.TABLES, 'chromatin_bins_{}.tsv.gz'.format(set_name))
    cap = os.path.join(cc.TABLES, 'capture_bins_{}_250kb.tsv.gz'.format(graph))
    if not os.path.exists(path):
        return
    bins = pd.read_csv(path, sep='\t')
    sub = bins[(bins.chrom == chrom) & (bins.start >= lo) & (bins.end <= hi)]
    tracks = [('Hi-C eigenvector', sub['ref_hic_{}_100kb'.format(cell)].values, sub.start.values),
              ('Marks PC1', sub.marks_pc1.values - np.nanmedian(sub.marks_pc1.values), sub.start.values)]
    if 'access_level' in sub:
        tracks.append(('DNase level', sub.access_level.values - np.nanmedian(sub.access_level.values),
                       sub.start.values))
    if os.path.exists(cap):
        c = pd.read_csv(cap, sep='\t')
        c = c[(c.chrom == chrom) & (c.start >= lo) & (c.end <= hi)]
        tracks.append(('Capture eigenvector (250 kb)', c.capture_eig_score.values, c.start.values))
    fig, axes = plt.subplots(len(tracks), 1, figsize=(6.5, 0.75 * len(tracks)), sharex=True)
    for ax, (label, values, starts) in zip(axes, tracks):
        sd = np.nanstd(values)
        v = np.clip(np.where(np.isfinite(values), values / (sd if sd > 0 else 1), 0), -3, 3)
        ax.fill_between(starts / 1e6, 0, v, where=v > 0, color='#c8443f', lw=0, step='post')
        ax.fill_between(starts / 1e6, 0, v, where=v <= 0, color='#2a78d6', lw=0, step='post')
        ax.axhline(0, color=GREY, lw=0.5)
        ax.set_ylim(-3, 3)
        ax.set_ylabel(label, rotation=0, ha='right', va='center')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
    axes[-1].set_xlabel('{} position (Mb), {}; values in SD units, clipped at 3'.format(chrom, cell))
    save(fig, 'example_track_{}_{}.pdf'.format(cell, chrom))


def capture_panels():
    path = os.path.join(cc.TABLES, 'capture_eval_summary.tsv')
    if not os.path.exists(path):
        return
    s = pd.read_csv(path, sep='\t')
    s = s[(s.chrom == 'genome') & s.reference.str.startswith('hic')]
    graphs = ['GM_chicago', 'GM_mifsud_counts', 'K562', 'H1_jung', 'hESC_H9']
    colour = {'GM_chicago': CELL_COLOUR['GM12878'], 'GM_mifsud_counts': '#7fb0ea', 'K562': CELL_COLOUR['K562'],
              'H1_jung': CELL_COLOUR['H1'], 'hESC_H9': '#7dd3b2'}
    methods = [('capture_eig_score', 'Capture eigenvector'), ('capture_density', 'Capture contact density'),
               ('chromatin_pc1', 'Marks PC1'), ('combined', 'Marks + capture')]
    for scope in ('covered_bins', 'all_bins'):
        fig, axes = plt.subplots(1, len(methods), figsize=(7.0, 2.0), sharey=True)
        for ax, (method, title) in zip(axes, methods):
            for g in graphs:
                sub = s[(s.graph == g) & (s.method == method) & (s.scope == scope)].sort_values('resolution')
                if len(sub):
                    ax.plot(sub.resolution / 1e6, sub.auc, marker='o', ms=4, lw=1.5, color=colour[g], label=g)
            ax.set_title(title, fontsize=8)
            ax.set_xscale('log')
            ax.set_xticks([0.25, 0.5, 1.0])
            ax.set_xticklabels(['0.25', '0.5', '1'])
            ax.axhline(0.5, color=GREY, lw=0.8, ls='--')
            ax.set_xlabel('Bin size (Mb)')
        axes[0].set_ylabel('AUC for Hi-C compartment A')
        axes[-1].legend(frameon=False, fontsize=6.5, loc='lower right')
        save(fig, 'capture_auc_{}.pdf'.format(scope))
    fig, ax = plt.subplots(figsize=(2.6, 2.0))
    for g in graphs:
        sub = s[(s.graph == g) & (s.method == 'capture_eig_score') & (s.scope == 'covered_bins')].sort_values(
            'resolution')
        if len(sub):
            ax.plot(sub.resolution / 1e6, sub.coverage, marker='o', ms=4, lw=1.5, color=colour[g], label=g)
    ax.set_xscale('log')
    ax.set_xticks([0.25, 0.5, 1.0])
    ax.set_xticklabels(['0.25', '0.5', '1'])
    ax.set_ylim(0, 1)
    ax.set_xlabel('Bin size (Mb)')
    ax.set_ylabel('Bins with a capture estimate')
    ax.legend(frameon=False, fontsize=6.5)
    save(fig, 'capture_coverage.pdf')


def domain_panels():
    path_b = os.path.join(cc.TABLES, 'domains_boundaries.tsv')
    if not os.path.exists(path_b):
        return
    b = pd.read_csv(path_b, sep='\t')
    labels = ['{} / {}'.format(g, t) for g, t in zip(b.graph, b.tads)]
    fig, ax = plt.subplots(figsize=(3.6, 0.3 * len(b) + 0.6))
    y = np.arange(len(b))
    ax.barh(y - 0.18, b.shift_enrichment, height=0.34, color='#2a78d6', label='vs circular shift')
    ax.barh(y + 0.18, b.matched_enrichment, height=0.34, color='#eb6834', label='vs node-matched gaps')
    ax.axvline(1, color=GREY, lw=0.8, ls='--')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel('Community boundaries near TAD boundaries\n(observed / expected)')
    ax.legend(frameon=False, fontsize=6.5, loc='lower right')
    save(fig, 'domains_boundary_enrichment.pdf')
    s = pd.read_csv(os.path.join(cc.TABLES, 'domains_spans.tsv'), sep='\t')
    for unit, name in (('community', 'domains_communities_within_tad.pdf'),
                       ('promoter_hub', 'domains_hubs_within_tad.pdf')):
        sub = s[s.unit == unit]
        fig, ax = plt.subplots(figsize=(3.6, 0.3 * len(sub) + 0.6))
        y = np.arange(len(sub))
        ax.barh(y - 0.18, sub.observed_fraction, height=0.34, color=DARK, label='observed')
        ax.barh(y + 0.18, sub.shift_fraction, height=0.34, color=GREY, label='shifted (expected)')
        ax.set_yticks(y)
        ax.set_yticklabels(['{} / {}'.format(g, t) for g, t in zip(sub.graph, sub.tads)])
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel('{} spans without a TAD boundary'.format('Community' if unit == 'community' else 'Hub'))
        ax.legend(frameon=False, fontsize=6.5, loc='lower right')
        save(fig, name)
    e = pd.read_csv(os.path.join(cc.TABLES, 'domains_edges.tsv'), sep='\t')
    fig, ax = plt.subplots(figsize=(3.6, 0.3 * len(e) + 0.6))
    y = np.arange(len(e))
    ax.barh(y - 0.18, e.observed_fraction, height=0.34, color=DARK, label='observed')
    ax.barh(y + 0.18, e.shift_fraction, height=0.34, color=GREY, label='shifted, same distances')
    ax.set_yticks(y)
    ax.set_yticklabels(['{} / {}'.format(g, t) for g, t in zip(e.graph, e.tads)])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel('Edges crossing a TAD boundary')
    ax.legend(frameon=False, fontsize=6.5, loc='lower right')
    save(fig, 'domains_edges_crossing.pdf')


def explanation_panel():
    path = os.path.join(cc.TABLES, 'domains_explanations.tsv')
    if not os.path.exists(path):
        return
    t = pd.read_csv(path, sep='\t')
    t = t[(t.roles == 'declared') & (t.comparison == 'hESC_vs_hNEC') & (t.tads == 'hicFindTADs_H1_50kb')]
    levels = [lv for lv in ('step1_marks', 'step1_state', 'step2', 'step3') if lv in set(t.level)]
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    x = np.arange(len(levels))
    for k, (cls, colour) in enumerate((('explaining', DARK), ('not_explaining', GREY))):
        sub = t[t['class'] == cls].set_index('level')
        obs = [sub.same_tad_fraction.get(lv, np.nan) for lv in levels]
        exp = [sub.expected_fraction.get(lv, np.nan) for lv in levels]
        ax.bar(x + (k - 0.5) * 0.38, obs, width=0.36, color=colour, label='{} (observed)'.format(cls.replace('_', ' ')))
        ax.scatter(x + (k - 0.5) * 0.38, exp, marker='_', s=160, linewidths=2, color='#c8443f', zorder=3,
                   label='distance-matched expectation' if k == 0 else None)
    ax.set_xticks(x)
    ax.set_xticklabels([lv.replace('_', ' ') for lv in levels])
    ax.set_ylim(0, 1)
    ax.set_ylabel('Strongest node in the promoter TAD')
    ax.set_title('hESC vs hNEC, H1 Micro-C TADs', fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    ax.set_ylim(0, 1.05)
    save(fig, 'explanations_within_tad.pdf')


def switch_panels():
    path = os.path.join(cc.OUT, 'switch', 'switch_summary.tsv')
    if not os.path.exists(path):
        return
    s = pd.read_csv(path, sep='\t')
    bins_path = os.path.join(cc.OUT, 'switch', 'hESC_vs_hNEC_quantile_bins.tsv.gz')
    if os.path.exists(bins_path):
        b = pd.read_csv(bins_path, sep='\t')
        fig, ax = plt.subplots(figsize=(2.6, 2.6))
        other = b[~b.switch.isin(['B->A', 'A->B'])]
        ax.scatter(other.score_A, other.score_B, s=1, color=GREY, alpha=0.3, lw=0, rasterized=True)
        for lab, colour in (('B->A', '#c8443f'), ('A->B', '#2a78d6')):
            sub = b[b.switch == lab]
            ax.scatter(sub.score_A, sub.score_B, s=2, color=colour, lw=0, label='{} ({})'.format(lab, len(sub)),
                       rasterized=True)
        ax.axhline(0, color=DARK, lw=0.5)
        ax.axvline(0, color=DARK, lw=0.5)
        ax.set_xlabel('Compartment score hESC (E008)')
        ax.set_ylabel('Compartment score hNEC (E009)')
        ax.legend(frameon=False, fontsize=6.5, markerscale=3)
        save(fig, 'switch_scatter_hESC_hNEC.pdf')
    sub = s[s.normalisation == 'quantile']
    fig, ax = plt.subplots(figsize=(2.8, 2.0))
    ax.bar(np.arange(len(sub)), sub.switch_fraction * 100, color=[DARK if 'control' not in c else GREY
                                                                  for c in sub.comparison])
    ax.set_xticks(np.arange(len(sub)))
    ax.set_xticklabels([c.replace('_control', '\n(control)').replace('_vs_', ' vs ') for c in sub.comparison])
    ax.set_ylabel('Bins switching (%)')
    save(fig, 'switch_rates.pdf')
    sub = s[(s.normalisation == 'quantile') & s.concordance.notna()]
    if len(sub):
        fig, ax = plt.subplots(figsize=(2.8, 2.0))
        x = np.arange(len(sub))
        ax.bar(x - 0.27, sub.concordance, width=0.26, color=DARK, label='switch bins')
        ax.bar(x, sub.matched_nonswitch_concordance, width=0.26, color=GREY, label='matched, no switch')
        ax.bar(x + 0.27, sub.matched_nonswitch_large_delta_concordance, width=0.26, color='#b9b8b2',
               label='matched, large change no switch')
        ax.axhline(0.5, color='#c8443f', lw=0.8, ls='--')
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace('_vs_', ' vs ') for c in sub.comparison])
        ax.set_ylim(0, 1)
        ax.set_ylabel('Expression change concordant')
        ax.legend(frameon=False, fontsize=6, loc='lower right')
        save(fig, 'switch_concordance.pdf')


def main():
    summary_path = os.path.join(cc.TABLES, 'chromatin_eval_summary.tsv')
    if os.path.exists(summary_path):
        summary = pd.read_csv(summary_path, sep='\t')
        chromatin_bars(summary, 'pearson', 'Pearson r with Hi-C eigenvector (genome)', 'chromatin_pearson.pdf')
        chromatin_bars(summary, 'sign_agreement', 'A/B sign agreement', 'chromatin_sign_agreement.pdf')
        chromatin_bars(summary, 'auc', 'AUC for Hi-C compartment A', 'chromatin_auc.pdf')
        chromatin_crosscell(summary)
        chromatin_sets(summary)
    cell_difference_panel()
    example_track()
    capture_panels()
    domain_panels()
    explanation_panel()
    switch_panels()


if __name__ == '__main__':
    main()
