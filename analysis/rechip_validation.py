"""Validation of the bivalency calls of EGA against the sequential ChIP (re-ChIP) of Minoux et al. 2017.

The re-ChIP (H3K27me3 then H3K4me2, GEO GSE89435) measures the co-occupancy of the two marks on the same nucleosomes;
its fold enrichment over the pooled input (analysis/17_minoux_fold.sh) is the reference. EGA calls a promoter
fragment (bait node of the Md / PA2 union graph) bivalent when an active mark (H3K27ac or H3K4me2) and the repressive
mark H3K27me3 are both present (mean fold enrichment over the fragment >= threshold, 1.5 / 2 / 3), active-only,
repressive-only or none (states.fragment_states with the declared roles), separately in Md and PA2.

Output (<out>/rechip_validation.tsv, one table with a `section` column):
  class           per condition and threshold: re-ChIP fold enrichment of the fragments of every call (n, median,
                  quartiles, fraction >= 2 and >= 3)
  call            per condition, threshold and re-ChIP cutoff (2, 3): precision and recall of the bivalent call;
                  rows with threshold 'score' give the AUC of the two-mark co-presence score min(max active FE,
                  repressive FE), whose cutoff is the presence threshold of the call
  transition      per threshold, for the promoters of the differentially expressed genes (<prefix>_gene_classes.tsv)
                  and for all bait fragments: state in Md -> state in PA2, the re-ChIP fold enrichment in both, the
                  median log2 change (PA2 + 1) / (Md + 1), the fraction losing re-ChIP signal, the fractions with
                  re-ChIP >= 2 in Md and in PA2 and a Wilcoxon signed-rank P (Md versus PA2); since the re-ChIP
                  signal of the two libraries differs as a whole, the change is also given relative to the
                  median change of all promoter fragments bivalent in both (change_minus_stable, Mann-Whitney
                  p_vs_stable) and as the change of the percentile rank among all promoter fragments of each
                  condition (median_rank_change, fraction_rank_down, Wilcoxon rank_wilcoxon_p)
  poised_in_both  differentially expressed genes explained as poised in both conditions by the cascade, at the
                  promoter (bivalent in both) and at step 1 (a contacting fragment poised in both): how many are
                  confirmed by re-ChIP >= 2 (>= 3) in both conditions; also all DE genes with a promoter bivalent in both
Figure <out>/rechip_validation.pdf: re-ChIP fold enrichment of the promoter fragments per EGA call (threshold 2) in
Md and PA2.

Usage: python rechip_validation.py <results dir> <out dir> [prefix (default hops_n/Md_vs_PA2)]
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from epigeneticgraphanalyzer.annotation import bigwig_interval_signal  # noqa: E402
from epigeneticgraphanalyzer.graph import read_graph  # noqa: E402
from epigeneticgraphanalyzer.states import fragment_states  # noqa: E402

CONDITIONS = ('Md', 'PA2')
MARKS = ('H3K27ac', 'H3K4me2', 'H3K27me3')
ROLES = {'H3K27ac': 'active', 'H3K4me2': 'active', 'H3K27me3': 'repressive'}
THRESHOLDS = (1.5, 2.0, 3.0)
CUTOFFS = (2.0, 3.0)
CALLS = {'bivalent': 'bivalent', 'active': 'active only', 'repressed': 'repressive only', 'unmarked': 'none'}
ORDER = ('bivalent', 'active only', 'repressive only', 'none')
COLOURS = {'Md': '#2a78d6', 'PA2': '#eb6834'}
INK, MUTED, GRID = '#1f1f1e', '#6b6a66', '#d9d8d2'


def signals(fold_dir, nodes):
    """{(condition, track): fold enrichment of every node} for the three marks and the re-ChIP."""
    parts = [n.replace(':', '-').split('-') for n in nodes]
    chroms = [p[0] for p in parts]
    starts = np.array([int(p[1]) for p in parts])
    ends = np.array([int(p[2]) for p in parts])
    return {(c, t): bigwig_interval_signal(os.path.join(fold_dir, '{}_{}.fold.bw'.format(c, t)), chroms, starts, ends)
            for c in CONDITIONS for t in MARKS + ('reChIP',)}


def auc(score, positive):
    """Area under the ROC curve (Mann-Whitney), ties one half."""
    positive = np.asarray(positive, dtype=bool)
    n1, n0 = positive.sum(), (~positive).sum()
    if not n1 or not n0:
        return np.nan
    ranks = stats.rankdata(score)
    return float((ranks[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def calls(values, condition, threshold):
    signal = {m: values[(condition, m)] for m in MARKS}
    states = fragment_states(signal, ROLES, threshold, np.ones(len(signal['H3K27me3']), dtype=bool))
    return np.array([CALLS[s] for s in states], dtype=object)


def describe(rechip):
    rechip = np.asarray(rechip, dtype=float)
    if not len(rechip):
        return dict(n=0)
    return dict(n=len(rechip), median_rechip=float(np.median(rechip)), q25_rechip=float(np.percentile(rechip, 25)),
                q75_rechip=float(np.percentile(rechip, 75)), fraction_ge_2=float(np.mean(rechip >= 2)),
                fraction_ge_3=float(np.mean(rechip >= 3)))


def rechip_change(values, index):
    md, pa2 = values[('Md', 'reChIP')][index], values[('PA2', 'reChIP')][index]
    return md, pa2, np.log2((pa2 + 1) / (md + 1))


def transition_rows(values, index, scope, threshold, stable):
    """State transitions Md -> PA2 of the fragments `index` with their re-ChIP change. The re-ChIP signal differs
    between the two libraries as a whole, so every change is also compared with that of the promoter fragments
    bivalent in both conditions (`stable`, all baits; Mann-Whitney P)."""
    a, b = calls(values, 'Md', threshold)[index], calls(values, 'PA2', threshold)[index]
    md, pa2, change = rechip_change(values, index)
    rows = []
    for label in pd.unique(pd.Series(a + ' -> ' + b)):
        mask = (a + ' -> ' + b) == label
        row = dict(section='transition', scope=scope, threshold=threshold, group=label, n=int(mask.sum()),
                   median_rechip_Md=float(np.median(md[mask])), median_rechip_PA2=float(np.median(pa2[mask])),
                   median_log2_change=float(np.median(change[mask])),
                   fraction_losing=float(np.mean(change[mask] < 0)),
                   fraction_ge_2_Md=float(np.mean(md[mask] >= 2)), fraction_ge_2_PA2=float(np.mean(pa2[mask] >= 2)))
        row['change_minus_stable'] = row['median_log2_change'] - float(np.median(stable))
        # rank of the fragment among all promoter fragments of its condition (percent), free of library scaling
        shift = values[('PA2', 'rank')][index][mask] - values[('Md', 'rank')][index][mask]
        row['median_rank_change'] = float(np.median(shift))
        row['fraction_rank_down'] = float(np.mean(shift < 0))
        if mask.sum() >= 5 and np.any(shift != 0):
            row['rank_wilcoxon_p'] = float(stats.wilcoxon(shift).pvalue)
        if mask.sum() >= 5 and np.any(md[mask] != pa2[mask]):
            row['wilcoxon_p'] = float(stats.wilcoxon(md[mask], pa2[mask]).pvalue)
        if mask.sum() >= 5 and label != 'bivalent -> bivalent':
            row['p_vs_stable'] = float(stats.mannwhitneyu(change[mask], stable).pvalue)
        rows.append(row)
    return rows


def figure(values, path):
    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    width = 0.36
    counts = {c: [] for c in ORDER}
    for k, condition in enumerate(CONDITIONS):
        call = calls(values, condition, 2.0)
        data = [np.maximum(values[(condition, 'reChIP')][call == c], 0.05) for c in ORDER]
        positions = np.arange(len(ORDER)) + (k - 0.5) * width
        parts = ax.boxplot(data, positions=positions, widths=width * 0.85, patch_artist=True, showfliers=False,
                           medianprops=dict(color=INK, lw=1), whiskerprops=dict(color=MUTED, lw=0.6),
                           capprops=dict(color=MUTED, lw=0.6), boxprops=dict(lw=0))
        for box in parts['boxes']:
            box.set_facecolor(COLOURS[condition])
        for c, d in zip(ORDER, data):
            counts[c].append('{:,}'.format(len(d)))
    ax.axhline(2, color=MUTED, lw=0.6, ls='--')
    ax.set_yscale('log')
    ax.set_ylim(0.05, None)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(['{}\n{}'.format(c.replace(' only', '\nonly') if ' ' in c else c + '\n', ' / '.join(counts[c]))
                        for c in ORDER], fontsize=6)
    ax.set_ylabel('re-ChIP fold enrichment', fontsize=7)
    ax.set_xlabel('EGA call of the promoter fragment (FE $\\geq$ 2; fragments Md / PA2)', fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid(axis='y', color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOURS[c]) for c in CONDITIONS]
    ax.legend(handles, CONDITIONS, fontsize=6, frameon=False, loc='upper right')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)


def main():
    results, out = sys.argv[1:3]
    prefix = os.path.join(results, sys.argv[3] if len(sys.argv) > 3 else 'hops_n/Md_vs_PA2')
    fold_dir = os.path.join(results, 'minoux_fold')
    os.makedirs(out, exist_ok=True)
    graph = read_graph(prefix + '_union.graphml')
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t', low_memory=False)
    baits = sorted(n for n, b in graph.nodes(data='bait') if b)
    extra = sorted(set(classes.get('step1_state_fragment', pd.Series(dtype=str)).dropna().astype(str)) - {''}
                   - set(baits))
    nodes = baits + extra
    position = {n: i for i, n in enumerate(nodes)}
    values = signals(fold_dir, nodes)
    bait_index = np.arange(len(baits))
    rows = []
    for condition in CONDITIONS:
        rechip = values[(condition, 'reChIP')][bait_index]
        active = np.maximum(values[(condition, 'H3K27ac')], values[(condition, 'H3K4me2')])[bait_index]
        score = np.minimum(active, values[(condition, 'H3K27me3')][bait_index])
        for cutoff in CUTOFFS:
            rows.append(dict(section='call', condition=condition, threshold='score', cutoff=cutoff,
                             n=len(baits), n_positive=int((rechip >= cutoff).sum()),
                             auc=auc(score, rechip >= cutoff)))
        for threshold in THRESHOLDS:
            call = calls(values, condition, threshold)[bait_index]
            for group in ORDER:
                rows.append(dict(section='class', condition=condition, threshold=threshold, group=group,
                                 **describe(rechip[call == group])))
            bivalent = call == 'bivalent'
            for cutoff in CUTOFFS:
                positive = rechip >= cutoff
                tp = int((bivalent & positive).sum())
                rows.append(dict(section='call', condition=condition, threshold=threshold, cutoff=cutoff,
                                 n=len(baits), n_called=int(bivalent.sum()), n_positive=int(positive.sum()),
                                 true_positives=tp, precision=tp / bivalent.sum() if bivalent.any() else np.nan,
                                 recall=tp / positive.sum() if positive.any() else np.nan,
                                 auc_binary=auc(bivalent.astype(float), positive)))
    de_index = np.array([position[n] for n in classes['node'] if n in position], dtype=np.int64)
    for condition in CONDITIONS:
        # percentile rank of the re-ChIP signal among the promoter fragments (bait ranks; other fragments 0)
        rank = np.zeros(len(nodes))
        rank[bait_index] = 100 * (stats.rankdata(values[(condition, 'reChIP')][bait_index]) - 0.5) / len(baits)
        values[(condition, 'rank')] = rank
    for threshold in THRESHOLDS:
        both = (calls(values, 'Md', threshold) == 'bivalent') & (calls(values, 'PA2', threshold) == 'bivalent')
        stable = rechip_change(values, bait_index[both[bait_index]])[2]
        rows += transition_rows(values, de_index, 'DE genes', threshold, stable)
        rows += transition_rows(values, bait_index, 'all baits', threshold, stable)

    def confirmed(fragments, cutoff):
        index = np.array([position[f] for f in fragments], dtype=np.int64)
        return int(((values[('Md', 'reChIP')][index] >= cutoff) & (values[('PA2', 'reChIP')][index] >= cutoff)).sum())

    poised = classes[classes.get('explained_type', pd.Series(dtype=str)) == 'poised unresolved']
    groups = [('promoter, explained as poised in both', poised.loc[poised['explained_at'] == 'gene', 'node']),
              ('step 1 fragment, explained as poised in both',
               poised.loc[poised['explained_at'] == 'step1', 'step1_state_fragment'].dropna().astype(str))]
    if 'promoter_unresolved' in classes:
        flag = classes['promoter_unresolved'].astype(str).str.lower() == 'true'
        groups.append(('promoter bivalent in both, all DE genes', classes.loc[flag, 'node']))
    for name, fragments in groups:
        fragments = [f for f in fragments if f in position]
        rows.append(dict(section='poised_in_both', threshold=2.0, group=name, n=len(fragments),
                         confirmed_ge_2=confirmed(fragments, 2.0) if fragments else 0,
                         confirmed_ge_3=confirmed(fragments, 3.0) if fragments else 0))
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(out, 'rechip_validation.tsv'), sep='\t', index=False)
    # re-ChIP of every bait in both conditions with the calls at threshold 2 (for gene look-ups)
    per_bait = pd.DataFrame({'node': baits, 'genes': [graph.nodes[n].get('genes', '') for n in baits]})
    for condition in CONDITIONS:
        for track in MARKS + ('reChIP',):
            per_bait['{}_{}'.format(track, condition)] = values[(condition, track)][bait_index]
        per_bait['call_' + condition] = calls(values, condition, 2.0)[bait_index]
    per_bait.to_csv(os.path.join(out, 'rechip_baits.tsv.gz'), sep='\t', index=False)
    figure({k: v[bait_index] for k, v in values.items()}, os.path.join(out, 'rechip_validation.pdf'))
    for section, part in table.groupby('section', sort=False):
        print(section)
        print(part.dropna(axis=1, how='all').drop(columns='section').to_string(index=False))


if __name__ == '__main__':
    main()
