"""egaExplain: an integrated report of a differential analysis, linking interaction changes, chromatin changes and
expression changes.

Input: the output prefix of egaDifferential (<prefix>_union.graphml, _gene_classes.tsv, _indirect_only_genes.tsv,
_gene_community_changes.tsv, _community_expression.tsv). Output:

  <out>_edge_chromatin.tsv  per chromatin mark: the log2 change of the mark at the two ends of gained, lost and shared
                            interactions, and whether gained (lost) interactions end on fragments whose mark changes
                            more (less) than at shared interactions of the same genomic distance (concordance, 0.5
                            without relation; Mann-Whitney P). Active marks are expected to rise at gained and to fall
                            at lost interactions, repressive marks (--repressive) the opposite.
  <out>_genes.tsv           per differentially expressed gene: expression change, interactions of its promoter
                            gained / lost / shared, community fate, the level (promoter, direct, indirect neighbours)
                            and the marks that explain the change, the nodes that explain it through indirect
                            contacts, and a one-line summary.
  <out>_communities.tsv     per community with a coordinated expression change (FDR < --fdr): rewiring, expression
                            change and the mean mark changes at its promoters and at their partners, with a summary.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from ._cli import base_parser, ensure_parent, setup_logging
from .graph import read_graph

log = logging.getLogger(__name__)

DISTANCE_BINS = np.array([0, 50e3, 100e3, 250e3, 500e3, 1e6, 2e6, np.inf])


def parse_arguments(args=None):
    parser = base_parser(__doc__.split('\n\n')[0])
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential')
    parser.add_argument('--names', nargs=2, default=None, help='Condition names A B (default: from the graph)')
    parser.add_argument('--repressive', nargs='*', default=['H3K27me3', 'H3K9me3'],
                        help='Repressive marks, expected to change opposite to activity')
    parser.add_argument('--fdr', type=float, default=0.1, help='Communities with a coordinated change at this FDR')
    parser.add_argument('--pseudo', type=float, default=1.0, help='Pseudo count of the log2 mark change')
    parser.add_argument('--outputPrefix', '-o', required=True)
    return parser.parse_args(args)


def condition_names(graph):
    """Condition names from the node attributes <mark>_<condition>."""
    keys = set()
    for _, data in graph.nodes(data=True):
        keys |= {k for k in data if '_' in k and k not in ('gene_ids',)}
        break
    suffixes = sorted({k.rsplit('_', 1)[1] for k in keys})
    return suffixes[:2]


def mark_changes(graph, names, pseudo):
    """Node x mark table of log2 (B + pseudo) / (A + pseudo)."""
    marks = sorted({k[:-len(names[0]) - 1] for _, d in graph.nodes(data=True) for k in d
                    if k.endswith('_' + names[0]) and k[:-len(names[0]) - 1] + '_' + names[1] in d})
    rows = {}
    for n, d in graph.nodes(data=True):
        rows[n] = [np.log2((float(d['{}_{}'.format(m, names[1])]) + pseudo) /
                           (float(d['{}_{}'.format(m, names[0])]) + pseudo))
                   if d.get('{}_{}'.format(m, names[0])) not in (None, '') else np.nan for m in marks]
    return pd.DataFrame.from_dict(rows, orient='index', columns=marks), marks


def edge_table(graph, changes):
    """One row per interaction: status, distance and the mean mark change of its two ends."""
    middle = {n: (int(d['start']) + int(d['end'])) // 2 for n, d in graph.nodes(data=True)}
    rows = []
    for u, v, d in graph.edges(data=True):
        rows.append((u, v, d.get('status', 'shared'), abs(middle[u] - middle[v])))
    edges = pd.DataFrame(rows, columns=['node_a', 'node_b', 'status', 'distance'])
    ends = (changes.loc[edges['node_a']].to_numpy() + changes.loc[edges['node_b']].to_numpy()) / 2
    for k, mark in enumerate(changes.columns):
        edges[mark] = ends[:, k]
    edges['bin'] = np.digitize(edges['distance'], DISTANCE_BINS)
    return edges


def matched_concordance(test, reference, bins_test, bins_reference):
    """P(test value > reference value of the same distance bin), ties one half."""
    wins, total = 0.0, 0
    for b in np.unique(bins_test):
        x = test[bins_test == b]
        y = np.sort(reference[bins_reference == b])
        if len(x) == 0 or len(y) == 0:
            continue
        below = np.searchsorted(y, x, side='left')
        equal = np.searchsorted(y, x, side='right') - below
        wins += (below + 0.5 * equal).sum()
        total += len(x) * len(y)
    return wins / total if total else np.nan


def edge_chromatin(edges, marks, repressive):
    rows = []
    shared = edges[edges['status'] == 'shared']
    for mark in marks:
        sign = -1.0 if mark in repressive else 1.0
        for status in ('gained', 'lost'):
            sub = edges[edges['status'] == status]
            x, y = sub[mark].to_numpy(), shared[mark].to_numpy()
            keep_x, keep_y = np.isfinite(x), np.isfinite(y)
            concordance = matched_concordance(x[keep_x], y[keep_y], sub['bin'].to_numpy()[keep_x],
                                              shared['bin'].to_numpy()[keep_y])
            # expected direction: activity up at gained, down at lost interactions (repressive marks opposite)
            expected_up = (status == 'gained') == (sign > 0)
            p = stats.mannwhitneyu(x[keep_x], y[keep_y], alternative='greater' if expected_up else 'less').pvalue \
                if keep_x.sum() and keep_y.sum() else np.nan
            rows.append(dict(mark=mark, repressive=mark in repressive, interactions=status, n=int(keep_x.sum()),
                             median_change=float(np.median(x[keep_x])) if keep_x.any() else np.nan,
                             shared_median_change=float(np.median(y[keep_y])) if keep_y.any() else np.nan,
                             concordance_vs_shared=concordance,
                             concordant=(concordance > 0.5) == expected_up if np.isfinite(concordance) else False,
                             mannwhitney_p_expected_direction=p))
    return pd.DataFrame(rows)


def describe_marks(values, marks, level, explaining=None):
    """Mark changes of at least 0.5 at one level, restricted to the explaining marks if given."""
    if explaining:
        marks = [m for m in marks if m in explaining]
    parts = []
    for mark in marks:
        value = values.get('{}_{}'.format(mark, level))
        if value is not None and np.isfinite(value) and abs(value) >= 0.5:
            parts.append('{} {}{:.1f}'.format(mark, '+' if value > 0 else '', value))
    return ', '.join(parts)


def gene_report(graph, prefix, names, marks, edges):
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    try:
        switches = pd.read_csv(prefix + '_gene_community_changes.tsv', sep='\t').set_index('node')
    except FileNotFoundError:
        switches = pd.DataFrame()
    try:
        indirect = pd.read_csv(prefix + '_indirect_only_genes.tsv', sep='\t').set_index('node')
    except FileNotFoundError:
        indirect = pd.DataFrame()
    status = {}
    for u, v, s in edges[['node_a', 'node_b', 'status']].itertuples(index=False):
        for n in (u, v):
            status.setdefault(n, {'gained': 0, 'lost': 0, 'shared': 0})[s] += 1
    rows = []
    for _, r in classes.iterrows():
        node = r['node']
        counts = status.get(node, {'gained': 0, 'lost': 0, 'shared': 0})
        gene = str(r['genes']).split(',')[0]
        level = r.get('explanation', '')
        fate = switches.loc[node, 'status'] if node in switches.index else ''
        transition = switches.loc[node, 'community_a_transition'] if node in switches.index else ''
        fate = '' if pd.isna(fate) else fate
        transition = '' if pd.isna(transition) else transition
        drivers = ''
        if node in indirect.index:
            items = []
            for item in str(indirect.loc[node, 'drivers']).split('; '):
                mark = item.split(':')[0]
                target = item.split('[')[-1].split(']')[0] if '[' in item else ''
                target = 'other end' if target == 'OE' else target
                if target and '{} at {}'.format(mark, target) not in items:
                    items.append('{} at {}'.format(mark, target))
            drivers = '; '.join(items[:6])
        if level in ('promoter', 'direct', 'indirect'):
            change = describe_marks(r, marks, level, str(r.get('explaining_marks', '')).split(','))
            where = {'promoter': 'at the promoter', 'direct': 'at direct partners',
                     'indirect': 'at partners two steps away'}[level]
            explanation = '{} {}'.format(change or r.get('explaining_marks', ''), where)
            if level == 'indirect' and drivers:
                explanation += ' ({})'.format(drivers)
        else:
            explanation = 'no concordant chromatin change'
        summary = '{} {}{:.1f}; interactions {} lost, {} gained, {} shared{}; {}'.format(
            gene, '+' if r['rna_log2FC'] > 0 else '', r['rna_log2FC'], counts['lost'], counts['gained'],
            counts['shared'], '; community {}'.format(transition) if transition else '', explanation)
        rows.append(dict(node=node, gene=gene, rna_log2FC=r['rna_log2FC'], interactions_lost=counts['lost'],
                         interactions_gained=counts['gained'], interactions_shared=counts['shared'],
                         community_status=fate, community_transition=transition, explanation_level=level,
                         explaining_marks=r.get('explaining_marks', ''), indirect_drivers=drivers,
                         summary=summary))
    return pd.DataFrame(rows).sort_values('rna_log2FC', key=np.abs, ascending=False)


def community_report(graph, prefix, changes, marks, fdr):
    try:
        coherent = pd.read_csv(prefix + '_community_expression.tsv', sep='\t')
        switches = pd.read_csv(prefix + '_gene_community_changes.tsv', sep='\t')
    except FileNotFoundError:
        return pd.DataFrame()
    coherent = coherent[coherent['coherence_q'] < fdr]
    rows = []
    for _, r in coherent.iterrows():
        side = 'community_b' if r['community_b'] >= 0 else 'community_a'
        key = r['community_b'] if r['community_b'] >= 0 else r['community_a']
        promoters = switches.loc[switches[side] == key, 'node'].tolist()
        promoters = [n for n in promoters if n in graph]
        partners = sorted({w for n in promoters for w in graph.neighbors(n)} - set(promoters))
        genes = switches.loc[switches[side] == key].dropna(subset=['rna_log2FC'])
        genes = genes.reindex(genes['rna_log2FC'].abs().sort_values(ascending=False).index)
        names = [g for g in ','.join(genes['genes'].astype(str)).split(',') if g and not g.startswith(('RP11-', 'AC0'))]
        row = dict(community=int(key), transition=r['transition'], n_genes=int(r['n_genes']),
                   mean_rna_log2FC=r['mean_rna_log2FC'], coherence_q=r['coherence_q'],
                   edges_shared=int(r['edges_shared']), edges_lost=int(r['edges_lost']),
                   edges_gained=int(r['edges_gained']), genes=','.join(dict.fromkeys(names))[:200])
        parts = []
        for mark in marks:
            at_promoters = float(changes.loc[promoters, mark].mean()) if promoters else np.nan
            at_partners = float(changes.loc[partners, mark].mean()) if partners else np.nan
            row['{}_promoters'.format(mark)] = at_promoters
            row['{}_partners'.format(mark)] = at_partners
            for value, where in ((at_promoters, 'promoters'), (at_partners, 'partners')):
                if np.isfinite(value) and abs(value) >= 0.3:
                    parts.append('{} {}{:.1f} at {}'.format(mark, '+' if value > 0 else '', value, where))
        row['summary'] = '{} ({}): {}, expression {}{:.1f}; interactions {} lost, {} gained, {} shared; {}'.format(
            row['genes'].split(',')[0] if row['genes'] else key, r['n_genes'], r['transition'],
            '+' if r['mean_rna_log2FC'] > 0 else '', r['mean_rna_log2FC'], row['edges_lost'], row['edges_gained'],
            row['edges_shared'], '; '.join(parts) if parts else 'no mean mark change of at least 0.3')
        rows.append(row)
    return pd.DataFrame(rows)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    graph = read_graph(args.prefix + '_union.graphml')
    names = args.names or condition_names(graph)
    changes, marks = mark_changes(graph, names, args.pseudo)
    log.info('conditions %s, marks %s', names, marks)
    edges = edge_table(graph, changes)
    association = edge_chromatin(edges, marks, set(args.repressive))
    association.to_csv(ensure_parent(args.outputPrefix + '_edge_chromatin.tsv'), sep='\t', index=False)
    genes = gene_report(graph, args.prefix, names, marks, edges)
    genes.to_csv(args.outputPrefix + '_genes.tsv', sep='\t', index=False)
    communities = community_report(graph, args.prefix, changes, marks, args.fdr)
    communities.to_csv(args.outputPrefix + '_communities.tsv', sep='\t', index=False)
    log.info('%d genes, %d communities; edge-chromatin association:\n%s', len(genes), len(communities),
             association[['mark', 'interactions', 'n', 'median_change', 'shared_median_change',
                          'concordance_vs_shared', 'mannwhitney_p_expected_direction']].to_string(index=False))


if __name__ == '__main__':
    main()
