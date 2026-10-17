"""Figures of the candidate trigger factors (egaTriggers, analysis/18_triggers.sh), one PDF per panel.

  enrichment_<comparison>.pdf   motifs enriched at the explaining nodes per gene group: the five motifs with the
                                smallest q per group, q the larger of the two backgrounds (matched random marked
                                nodes, non-explaining nodes of the same genes); filled and coloured by -log10 q when
                                q < 0.05, open otherwise; dot size the log2 odds ratio against random marked nodes.
  candidates_<comparison>_<gene>.pdf  the ten best ranked candidate factors of one gene: expression change of the
                                factor (bar), relation to the gene (colour: activator- / repressor-consistent),
                                number of the gene's scanned nodes with the motif and the steps with a hit (label).

Usage: python figure_triggers.py <triggers dir> <figure dir> [comparison:gene ...]
Without genes: Hoxa3 and Hoxa2 (Md_vs_PA2), the three hNEC bivalent -> active genes with the largest increase
(hESC_vs_hNEC) and the three Kasumi-1 genes with the largest absolute change explained beyond the promoter
(siMM_vs_siRE).
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COMPARISONS = ['hESC_vs_hNEC', 'Md_vs_PA2', 'siMM_vs_siRE', 't821_vs_FLT3ITD']
COLOURS = {'activator-consistent': '#2a78d6', 'repressor-consistent': '#eb6834'}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})


def combined(enrichment):
    """One row per group and motif: q = the larger q of the two backgrounds, effect against random marked nodes."""
    wide = enrichment.pivot_table(index=['group', 'motif_id', 'motif'], columns='background',
                                  values=['q', 'p', 'log2_odds', 'fg_fraction', 'bg_fraction'], aggfunc='first')
    out = pd.DataFrame(index=wide.index)
    q = wide['q']
    out['q'] = q.max(axis=1, skipna=False) if q.shape[1] == 2 else np.nan
    out['p'] = wide['p'].max(axis=1, skipna=False) if q.shape[1] == 2 else np.nan
    out['q_random'] = q.get('random_marked')
    out['q_same'] = q.get('same_genes')
    out['log2_odds'] = wide['log2_odds'].get('random_marked')
    out['fg_fraction'] = wide['fg_fraction'].get('random_marked')
    out['bg_fraction'] = wide['bg_fraction'].get('random_marked')
    return out.reset_index()


def enrichment_panel(enrichment, comparison, path, top=5, alpha=0.05):
    table = combined(enrichment).dropna(subset=['q'])
    hits = table.sort_values(['q', 'p']).groupby('group').head(top)
    groups = sorted(set(hits['group']))
    motifs = list(dict.fromkeys(hits.sort_values(['group', 'q'])['motif']))
    sub = table[table['group'].isin(groups) & table['motif'].isin(motifs)]
    sub = sub.sort_values('q').drop_duplicates(['group', 'motif'])
    x = sub['group'].map({g: i for i, g in enumerate(groups)})
    y = sub['motif'].map({m: i for i, m in enumerate(motifs)})
    significant = sub['q'] < alpha
    top_cell = sub.set_index(['group', 'motif']).index.isin(hits.set_index(['group', 'motif']).index)
    logq = -np.log10(sub['q'].clip(lower=1e-30))
    size = 12 + 40 * sub['log2_odds'].clip(lower=0, upper=4)
    fig, ax = plt.subplots(figsize=(0.55 * len(groups) + 2.6, 0.17 * len(motifs) + 1.6))
    weak = top_cell & ~significant
    ax.scatter(x[weak], y[weak], s=size[weak], facecolors='none', edgecolors='#86b6ef', linewidths=0.8)
    rest = ~top_cell & ~significant
    ax.scatter(x[rest], y[rest], s=4, color='#c3c2b7', linewidths=0)
    vmax = max(2.0, float(logq[significant].max()) if significant.any() else 2.0)
    points = ax.scatter(x[significant], y[significant], s=size[significant], c=logq[significant], cmap='Blues',
                        vmin=0, vmax=vmax, edgecolors='#fcfcfb', linewidths=0.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=60, ha='right')
    ax.set_yticks(range(len(motifs)))
    ax.set_yticklabels(motifs)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_title('{}: best {} motifs per group (filled: q < {} against both backgrounds; open: '
                 'not significant)'.format(comparison.replace('_vs_', ' vs '), top, alpha), fontsize=8)
    bar = fig.colorbar(points, ax=ax, shrink=0.5, pad=0.02)
    bar.set_label('$-\\log_{10} q$')
    for value in (1, 2, 4):
        ax.scatter([], [], s=12 + 40 * value, color='#86b6ef', label='log2 OR {}'.format(value))
    ax.legend(loc='upper left', bbox_to_anchor=(1.25, 1.0), frameon=False, labelspacing=1.2)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def candidate_panel(triggers, comparison, gene, path, top=10):
    rows = triggers[triggers['gene'].astype(str).str.split(',').apply(lambda gs: gene in gs)]
    if rows.empty:
        return False
    first = rows.iloc[0]
    # a gene with several promoter fragments: the first listed fragment
    rows = rows[rows['gene_node'] == first['gene_node']]
    rows = rows.sort_values('rank').drop_duplicates('tf').head(top)
    fig, ax = plt.subplots(figsize=(4.2, 0.28 * len(rows) + 1.0))
    y = np.arange(len(rows))
    ax.barh(y, rows['tf_log2FC'], height=0.6, color=[COLOURS[r] for r in rows['relation']])
    ax.axvline(0, color='#52514e', linewidth=0.6)
    labels = ['{} ({}{})'.format(tf, motif, '' if ev == 'named factor' else ', family')
              for tf, motif, ev in zip(rows['tf'], rows['motif'], rows['evidence'])]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    span = float(rows['tf_log2FC'].abs().max()) or 1.0
    for yi, (lfc, n, scanned, steps) in enumerate(zip(rows['tf_log2FC'], rows['n_nodes_with_motif'],
                                                      rows['n_nodes_scanned'], rows['steps_with_motif'])):
        ax.text(max(lfc, 0) + 0.04 * span, yi, '{}/{} {}'.format(n, scanned, steps),
                va='center', ha='left', fontsize=6.5, color='#52514e')
    ax.set_xlim(-1.1 * span, 2.2 * span)
    ax.set_xlabel('factor expression log2FC')
    ax.set_title('{} {} ({}, log2FC {:+.1f}; explained at {})\nlabel: scanned nodes with the motif, steps'.format(
        gene, first['gene_node'], comparison.replace('_vs_', ' vs '), first['gene_log2FC'], first['explained_at']),
        fontsize=7.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLOURS.values()]
    ax.legend(handles, list(COLOURS), loc='lower left', bbox_to_anchor=(0, 1.08), ncol=2, frameon=False,
              fontsize=6.5)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return True


def default_genes(directory):
    picks = [('Md_vs_PA2', 'Hoxa3'), ('Md_vs_PA2', 'Hoxa2')]
    path = os.path.join(directory, 'hESC_vs_hNEC_gene_triggers.tsv')
    if os.path.exists(path):
        t = pd.read_csv(path, sep='\t')
        biv = t[t['promoter_transition'].astype(str).str.startswith('bivalent -> active')]
        biv = biv.drop_duplicates('gene').sort_values('gene_log2FC', ascending=False)
        picks += [('hESC_vs_hNEC', g.split(',')[0]) for g in biv['gene'].head(3)]
    path = os.path.join(directory, 'siMM_vs_siRE_gene_triggers.tsv')
    if os.path.exists(path):
        t = pd.read_csv(path, sep='\t')
        t = t[t['explained_at'].astype(str).str.startswith('step')].drop_duplicates('gene')
        t = t.reindex(t['gene_log2FC'].abs().sort_values(ascending=False).index)
        picks += [('siMM_vs_siRE', g.split(',')[0]) for g in t['gene'].head(3)]
    return picks


def main(argv):
    directory, figures = argv[1], argv[2]
    os.makedirs(figures, exist_ok=True)
    for comparison in COMPARISONS:
        path = os.path.join(directory, comparison + '_motif_enrichment.tsv')
        if os.path.exists(path) and os.path.getsize(path) > 1:
            enrichment_panel(pd.read_csv(path, sep='\t'), comparison,
                             os.path.join(figures, 'enrichment_{}.pdf'.format(comparison)))
    picks = [tuple(a.split(':', 1)) for a in argv[3:]] or default_genes(directory)
    cache = {}
    for comparison, gene in picks:
        if comparison not in cache:
            path = os.path.join(directory, comparison + '_gene_triggers.tsv')
            cache[comparison] = pd.read_csv(path, sep='\t') if os.path.exists(path) else pd.DataFrame()
        if cache[comparison].empty or not candidate_panel(
                cache[comparison], comparison, gene,
                os.path.join(figures, 'candidates_{}_{}.pdf'.format(comparison, gene))):
            print('no candidates for {} in {}'.format(gene, comparison))


if __name__ == '__main__':
    main(sys.argv)
