"""Figures of the keratinocyte time course (egaTimecourse, analysis/22_timecourse.sh), one PDF per panel.

  trajectory_<gene>.pdf          one gene: log2 change from day 0 of expression and of each explanation level
                                 (promoter H3K27ac, promoter contacts, contact-weighted H3K27ac at steps 1-3); the
                                 levels explaining the gene in each pairwise run are listed under the plot
  ordering_<subset>_<layers>.pdf order categories of the event intervals: observed units (bar, 95% bootstrap interval
                                 over genes), the within-gene and the layer-marginal permutation null (mean and
                                 2.5-97.5 percentiles); p values of both nulls above the bars (within-gene first);
                                 ordering_donor2_<layers>.pdf the same with H3K27ac of donor 2 only
  knockdown_<factor>_candidate.pdf  H3K27ac change siRNA knockdown vs control at the explaining nodes of genes with the
                                 factor as candidate trigger and of the other explained genes (day 3 vs day 0 run)
  knockdown_<factor>_motif.pdf   the same for explaining nodes with and without a motif of the factor
  knockdown_effects.pdf          difference of medians (candidate - other) with 95% bootstrap intervals, per factor
                                 and gene direction

Usage: python figure_timecourse.py <timecourse dir> <figure dir> [gene ...]
Without genes: up to six keratinocyte differentiation genes explained in at least one pairwise run.
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
DAYS = ['D0', 'D3', 'D6']
DAY_VALUES = [0, 3, 6]
PAIRS = ['D0_vs_D3', 'D3_vs_D6', 'D0_vs_D6']
INK = '#2b2b2b'
MUTED = '#8a8a86'
LEVELS = [('expression', 'expression', INK), ('H3K27ac_promoter', 'promoter H3K27ac', '#2a78d6'),
          ('contacts', 'promoter contacts', '#eb6834'), ('H3K27ac_step1', 'H3K27ac step 1', '#1baf7a'),
          ('H3K27ac_step2', 'H3K27ac step 2', '#eda100'), ('H3K27ac_step3', 'H3K27ac step 3', '#e87ba4')]
CANDIDATES = ['KRT1', 'KRT10', 'IVL', 'LOR', 'FLG', 'SBSN', 'KRTDAP', 'DSC1', 'DSG1', 'GRHL3', 'ZNF750', 'KLF4',
              'EHF', 'HOPX', 'PRDM1', 'SPRR1A', 'LCE3D', 'CALML5', 'S100A8', 'TGM1', 'MAF', 'MAFB']
LABEL = {'forward': 'contact / mark\nbefore expression', 'simultaneous': 'same\ninterval',
         'reversed': 'expression\nfirst', 'other': 'other', 'before': 'first layer\nearlier',
         'same': 'same\ninterval', 'after': 'first layer\nlater'}
LAYER_TITLE = {'contact_mark_expression': 'contact, H3K27ac, expression', 'contact_expression': 'contact, expression',
               'mark_expression': 'H3K27ac, expression', 'contact_mark': 'contact, H3K27ac'}


def default_genes(traj, n=6):
    explained = traj[[c for c in traj.columns if c.startswith('explained_')]].fillna('').astype(str)
    has = explained.apply(lambda r: any(v for v in r), axis=1)
    pool = traj[has & (traj['pattern'] != '00')]
    genes = [g for g in CANDIDATES if g in set(pool['gene'])]
    return genes[:n]


def trajectory_panel(traj, gene, path):
    row = traj[traj['gene'] == gene].iloc[0]
    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    ends = []
    for key, label, colour in LEVELS:
        values = np.array([row['{}_{}'.format(key, d)] for d in DAYS], dtype=float)
        if np.isnan(values).any():
            continue
        values = values - values[0]
        ax.plot(DAY_VALUES, values, color=colour, lw=2.4 if key == 'expression' else 1.6, marker='o', ms=4,
                label=label, zorder=3 if key == 'expression' else 2)
        ends.append((values[-1], label))
    ax.axhline(0, color=MUTED, lw=0.6, zorder=1)
    ax.set_xticks(DAY_VALUES)
    ax.set_xticklabels(['day 0', 'day 3', 'day 6'])
    ax.set_ylabel('log2 change from day 0')
    ax.set_title('{} ({})'.format(gene, row['node']), fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc='center left', bbox_to_anchor=(1.0, 0.5))
    lines = []
    for pair in PAIRS:
        explained = str(row.get('explained_' + pair, '') or '')
        explained = explained if explained not in ('', 'nan') else 'none'
        lines.append('{}: {}'.format(pair.replace('_vs_', ' vs '), explained.replace(',', ', ')))
    fig.text(0.02, -0.02, '\n'.join(lines), fontsize=5.5, va='top', color=INK)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def ordering_panel(ordering, subset, layers, path, threshold=1.0, label=None):
    sub = ordering[(ordering['subset'] == subset) & (ordering['layers'] == layers) &
                   (ordering['chromatin_threshold'] == threshold)]
    if sub.empty or sub['observed'].sum() == 0:
        return False
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(3.0, 2.4))
    obs = sub['observed'].to_numpy(dtype=float)
    ax.bar(x, obs, width=0.6, color='#2a78d6', label='observed', zorder=2)
    ax.errorbar(x, obs, yerr=[obs - sub['boot_low'], sub['boot_high'] - obs], fmt='none', ecolor=INK, lw=0.8,
                capsize=2, zorder=3)
    ax.errorbar(x - 0.12, sub['null_mean'], yerr=[sub['null_mean'] - sub['null_low'], sub['null_high'] -
                                                    sub['null_mean']],
                fmt='D', color='#eb6834', ms=4, lw=1.0, capsize=2, label='within-gene null', zorder=4)
    if 'marginal_mean' in sub:
        ax.errorbar(x + 0.12, sub['marginal_mean'], yerr=[sub['marginal_mean'] - sub['marginal_low'],
                                                          sub['marginal_high'] - sub['marginal_mean']],
                    fmt='s', color='#1baf7a', ms=4, lw=1.0, capsize=2, label='layer-marginal null', zorder=4)
    for xi, (o, p) in enumerate(zip(obs, sub['p'])):
        top = max(sub['boot_high'].iloc[xi], sub['null_high'].iloc[xi],
                  sub.get('marginal_high', sub['null_high']).iloc[xi])
        text = 'p={:.2g}'.format(p)
        if 'marginal_p' in sub:
            text += '\np={:.2g}'.format(sub['marginal_p'].iloc[xi])
        ax.text(xi, top * 1.03 + 1, text, ha='center', fontsize=5, color=INK)
    highest = max(sub['boot_high'].max(), sub['null_high'].max(), sub.get('marginal_high', sub['null_high']).max())
    ax.set_ylim(0, highest * 1.3 + 2)
    ax.set_xticks(x)
    ax.set_xticklabels([LABEL[c] for c in sub['category']], fontsize=6)
    ax.set_ylabel('units (gene, explaining node)')
    ax.set_title('{}: {}\n{} units, {} genes'.format(label or subset, LAYER_TITLE[layers],
                                                     int(sub['units'].iloc[0]), int(sub['genes'].iloc[0])),
                 fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return True


def knockdown_panel(nodes, test, factor, evidence, path):
    row = test[(test['factor'] == factor) & (test['evidence'] == evidence) & (test['genes'] == 'all') &
               (test['measure'] == 'mean')]
    if row.empty or not row['nodes_candidate'].iloc[0] or pd.isna(row['median_difference'].iloc[0]):
        return False
    row = row.iloc[0]
    if evidence == 'any':
        cand = nodes['candidate_' + factor].astype(bool)
        ctrl = ~cand
        labels = ('{} candidate trigger'.format(factor), 'other explained')
    else:
        cand = nodes['motif_' + factor].astype(bool)
        ctrl = nodes['scanned'].astype(bool) & ~cand
        labels = ('{} motif'.format(factor), 'no {} motif'.format(factor))
    values = [nodes.loc[mask, 'log2_' + factor].dropna().to_numpy() for mask in (cand, ctrl)]
    fig, ax = plt.subplots(figsize=(2.4, 2.4))
    parts = ax.violinplot(values, positions=[0, 1], showextrema=False, widths=0.7)
    for body, colour in zip(parts['bodies'], ['#2a78d6', MUTED]):
        body.set_facecolor(colour)
        body.set_alpha(0.5)
    ax.boxplot(values, positions=[0, 1], widths=0.18, showfliers=False, medianprops={'color': INK})
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['{}\n({} nodes)'.format(labels[0], len(values[0])),
                        '{}\n({} nodes)'.format(labels[1], len(values[1]))], fontsize=6)
    ax.set_ylabel('H3K27ac log2 si{} / siControl'.format(factor))
    ax.set_title('median diff {:+.3f} [{:+.3f}, {:+.3f}]\nCliff delta {:+.3f}, p={:.2g}'.format(
        row['median_difference'], row['md_low'], row['md_high'], row['cliffs_delta'], row['p_two_sided']),
        fontsize=6.5)
    low, high = np.percentile(np.concatenate(values), [1, 99])
    ax.set_ylim(low - 0.2, high + 0.2)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return True


def knockdown_effects(test, path):
    sub = test[test['evidence'].isin(['any', 'motif at node']) & (test['measure'] == 'mean')].dropna(
        subset=['median_difference'])
    if sub.empty:
        return False
    fig, ax = plt.subplots(figsize=(3.2, 0.3 * len(sub) + 0.8))
    y = np.arange(len(sub))[::-1]
    ax.errorbar(sub['median_difference'], y, xerr=[sub['median_difference'] - sub['md_low'],
                                                   sub['md_high'] - sub['median_difference']],
                fmt='o', color='#2a78d6', ms=4, capsize=2, lw=1)
    ax.axvline(0, color=MUTED, lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(['si{} {}: {} genes ({} / {} nodes)'.format(
        r.factor, 'candidate' if r.evidence == 'any' else 'motif', r.genes, r.nodes_candidate, r.nodes_control)
        for r in sub.itertuples()],
                       fontsize=6)
    ax.set_xlabel('difference of medians, H3K27ac log2 KD / control\n(candidate-trigger genes - other genes)')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return True


def main():
    base, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    prefix = os.path.join(base, 'timecourse')
    traj = pd.read_csv(prefix + '_trajectories.tsv', sep='\t')
    genes = sys.argv[3:] or default_genes(traj)
    for gene in genes:
        if gene in set(traj['gene']):
            trajectory_panel(traj, gene, os.path.join(out, 'trajectory_{}.pdf'.format(gene)))
    ordering = pd.read_csv(prefix + '_ordering.tsv', sep='\t')
    for subset in ('all', 'promoter', 'step1', 'step2', 'step3'):
        ordering_panel(ordering, subset, 'contact_mark_expression',
                       os.path.join(out, 'ordering_{}_contact_mark_expression.pdf'.format(subset)))
    for layers in ('contact_expression', 'mark_expression', 'contact_mark'):
        ordering_panel(ordering, 'all', layers, os.path.join(out, 'ordering_all_{}.pdf'.format(layers)))
    if os.path.exists(os.path.join(base, 'timecourse_donor2_ordering.tsv')):
        donor2 = pd.read_csv(os.path.join(base, 'timecourse_donor2_ordering.tsv'), sep='\t')
        for layers in ('contact_mark_expression', 'mark_expression'):
            ordering_panel(donor2, 'all', layers, os.path.join(out, 'ordering_donor2_{}.pdf'.format(layers)),
                           label='H3K27ac donor 2 only')
    if os.path.exists(prefix + '_knockdown_test.tsv'):
        test = pd.read_csv(prefix + '_knockdown_test.tsv', sep='\t')
        nodes = pd.read_csv(prefix + '_knockdown_nodes.tsv', sep='\t')
        for factor in test['factor'].unique():
            knockdown_panel(nodes, test, factor, 'any', os.path.join(out, 'knockdown_{}_candidate.pdf'.format(factor)))
            if 'motif_' + factor in nodes:
                knockdown_panel(nodes, test, factor, 'motif at node',
                                os.path.join(out, 'knockdown_{}_motif.pdf'.format(factor)))
        knockdown_effects(test, os.path.join(out, 'knockdown_effects.pdf'))
    print('figures in', out, sorted(os.listdir(out)))


if __name__ == '__main__':
    main()
