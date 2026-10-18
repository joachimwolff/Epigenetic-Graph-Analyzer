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
                            gained / lost / shared, community fate, the level (promoter, direct neighbours, nodes at
                            2, 3, ... hops) and the marks that explain the change, the nodes that explain it through
                            indirect contacts, and a one-line summary. With egaDifferential --proximitySupport the
                            summary states how often graph neighbours at the explaining level are close in the
                            nucleus (support decreases with the number of hops). Every gene is described by the
                            level of the explanation cascade that explains it first (explanation_detail of
                            egaDifferential): at the promoter by marks with direction and value and / or by the
                            contact gain or loss of the promoter (interactions gained, lost, log2 change), at the
                            enhancers (step 1) by marks, the number of concordant fragments and the strongest, or at
                            hop k by marks and the number of nodes with the 3D proximity support of hop k and the
                            nodes carrying the change; a chromatin state transition of the promoter
                            (e.g. 'promoter bivalent in hESC (H3K4me3 4.1, H3K27me3 3.2), resolved to active in hNEC
                            by loss of H3K27me3') or of a contacting fragment explaining the gene is added.
                            With <prefix>_explanations_all_levels.tsv (egaDifferential), all_levels lists every
                            level in rank order with its result under the declared roles (all_levels_learned: the
                            roles learned at each level), also when an earlier level explains the gene: promoter
                            marks, contacts, promoter state, (with egaDifferential --structuralVariants /
                            --copyNumber: structural variant, copy number,) step 1 marks, step 1 states, step 2 ..
                            step n; the
                            state levels under both presence statistics (mean | max1kb), followed by the explaining
                            levels, e.g. 'Hoxa3 (up, log2 10.6): promoter marks: no concordant change (...);
                            contacts: ...; promoter state: bivalent in Md and PA2 (H3K4me2 5.1/5.2, H3K27me3 67/80;
                            mean) - poised, trigger not among the measured layers | ...; step 2: H3K27ac +3.4 at 12
                            nodes (5 concordant), strongest Hoxa2 promoter H3K27ac log2 +3.1, P 8.5e-05, q 0.004,
                            beyond mirror, 3D support n/a; explaining levels (mean): promoter_state;step2'.
                            With --tads (TAD / boundary calls, see egaTadCrossing): the gained / lost / shared
                            contacts of the promoter that cross a boundary (gained_crossing, lost_crossing,
                            shared_crossing) and a sentence (tad_crossing) naming the boundary, its strength and the
                            partner region, appended to the summary.
                            With the support columns of the all-level table (egaDifferential), appended after the
                            other columns (declared roles, fragment mean): first_level_support (stable,
                            threshold-sensitive, not testable: the contacts that carry the first explaining level),
                            first_level_threshold_sensitive (it rests only on threshold-sensitive contacts),
                            only_threshold_sensitive (every explaining level does), threshold_sensitive_levels,
                            support_note (per explaining level with contacts, e.g. 'step 2: rests on 2
                            threshold-sensitive contacts') and provenance (tags of the first explaining level: M
                            measured, C called, D derived, I interpretation, e.g. 'MC>D>I').
                            With egaDifferential --stepFeature decomposed: step_terms (per step level with a
                            concordant feature the combined feature, the mark term and the contact term of the marks
                            and which term drives it: marks, contacts, both, neither) and first_step_driver (the
                            driver of the first step level explaining the gene, declared roles, fragment mean).
  <out>_communities.tsv     per community with a coordinated expression change (FDR < --fdr): rewiring, expression
                            change and the mean mark changes at its promoters and at their partners, with a summary.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from ._cli import base_parser, ensure_parent, setup_logging
from .differential import STATE_LEVELS, explanation_level_label, level_hops
from .egaTadCrossing import add_tad_arguments, tad_kwargs
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
    add_tad_arguments(parser)
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


STEP_WORDS = {2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}


def level_location(level):
    """Where a level lies, as in 'H3K27ac +1.2 at partners three steps away'."""
    hops = level_hops(level)
    if hops == 0:
        return 'at the promoter'
    if hops == 1:
        return 'at direct partners'
    return 'at partners {} steps away'.format(STEP_WORDS.get(hops, hops))


def proximity_sentence(hops, support, interval):
    """'explained at 3 hops, where graph neighbours are close in the nucleus with probability 0.54 (0.52-0.56)
    against 0.5 for matched pairs'; empty without a support value."""
    if pd.isna(support) or pd.isna(hops) or int(hops) == 0:
        return ''
    text = 'explained at {} hop{}, where graph neighbours are close in the nucleus with probability {:.2f}'.format(
        int(hops), '' if int(hops) == 1 else 's', float(support))
    if isinstance(interval, str) and interval and interval != 'NA':
        text += ' ({})'.format(interval)
    return text + ' against 0.5 for matched pairs'


def all_level_text(gene, rna, levels, roles):
    """Every level of one gene in rank order with its result under the role set `roles` (rows of
    <prefix>_explanations_all_levels.tsv); state levels under both presence statistics."""
    parts = []
    levels = levels.sort_values('rank')
    presences = [p for p in ('mean', 'max1kb') if 'detail_{}_{}'.format(roles, p) in levels]
    for _, r in levels.iterrows():
        texts = [str(r['detail_{}_{}'.format(roles, p)]) for p in presences]
        if r['level'] not in STATE_LEVELS:
            texts = texts[:1]
        parts.append('{}: {}'.format(explanation_level_label(r['level']), ' | '.join(texts)))
    explaining = []
    for p in presences:
        hit = levels.loc[levels['explains_{}_{}'.format(roles, p)].astype(str).str.lower() == 'true', 'level']
        explaining.append('{} ({})'.format(';'.join(hit) or 'none', p))
    return '{} ({}, log2 {:.1f}): {}; explaining levels: {}'.format(
        gene, 'up' if rna > 0 else 'down', rna, '; '.join(parts), ', '.join(explaining))


def term_columns(levels, roles='declared', presence='mean'):
    """Mark term and contact term of the step levels of one gene (egaDifferential --stepFeature decomposed; rows of
    the all-level table): step_terms, per step with a concordant feature 'step 2: H3K27ac combined +1.3, mark term
    +0.2, contact term +1.1 (contacts)', and first_step_driver, the driver (marks, contacts, both, neither) of the
    first step that explains the gene. {} without the term columns."""
    driver = 'driver_' + roles
    if driver not in levels:
        return {}
    levels = levels.sort_values('rank')
    parts, first = [], ''
    for _, r in levels.iterrows():
        value = r.get(driver)
        if not isinstance(value, str) or not value:
            continue
        values = {}
        for column, name in (('combined_changes', 'combined'), ('mark_term_changes', 'mark term'),
                             ('contact_term_changes', 'contact term')):
            for item in str(r.get(column, '')).split(';'):
                if ':' in item:
                    mark, number = item.split(':')
                    values.setdefault(mark, []).append('{} {}'.format(name, number))
        concordant = set()
        for column in ('concordant_combined_', 'concordant_mark_term_', 'concordant_contact_term_'):
            text = r.get(column + roles, '')
            concordant |= set(text.split(',')) if isinstance(text, str) and text else set()
        shown = [m for m in values if m in concordant] or list(values)
        parts.append('{}: {} ({})'.format(explanation_level_label(r['level']),
                                          '; '.join('{} {}'.format(m, ', '.join(values[m])) for m in shown), value))
        explains = str(r.get('explains_{}_{}'.format(roles, presence), '')).lower() == 'true'
        if explains and not first:
            first = value
    return dict(step_terms=' | '.join(parts), first_step_driver=first)


def support_columns(levels, roles='declared', presence='mean'):
    """Support of the explaining levels of one gene (rows of the all-level table); {} without support columns."""
    key = '{}_{}'.format(roles, presence)
    if 'support_' + key not in levels:
        return {}
    levels = levels.sort_values('rank')
    hit = levels[levels['explains_' + key].astype(str).str.lower() == 'true']
    flagged = hit['threshold_sensitive_only_' + key].astype(str).str.lower() == 'true'
    notes = ['{}: {}'.format(explanation_level_label(r['level']), r['support_detail_' + key])
             for _, r in hit.iterrows() if isinstance(r['support_detail_' + key], str) and r['support_detail_' + key]]
    first = hit.iloc[0] if len(hit) else None
    support = first['support_' + key] if first is not None else ''
    return dict(first_level_support=support if isinstance(support, str) else '',
                first_level_threshold_sensitive=bool(flagged.iloc[0]) if len(hit) else False,
                only_threshold_sensitive=bool(len(hit) and flagged.all()),
                threshold_sensitive_levels=';'.join(hit.loc[flagged.to_numpy(), 'level']),
                support_note='; '.join(notes),
                provenance=str(first['provenance_chain']) if first is not None and 'provenance_chain' in hit and
                isinstance(first['provenance_chain'], str) else '')


def gene_report(graph, prefix, names, marks, edges):
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    try:
        all_levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', keep_default_na=False,
                                 na_values=[''])
        all_levels = dict(tuple(all_levels.groupby('node', sort=False)))
    except FileNotFoundError:
        all_levels = {}
    try:
        switches = pd.read_csv(prefix + '_gene_community_changes.tsv', sep='\t').set_index('node')
    except FileNotFoundError:
        switches = pd.DataFrame()
    try:
        indirect = pd.read_csv(prefix + '_indirect_only_genes.tsv', sep='\t').set_index('node')
    except (FileNotFoundError, pd.errors.EmptyDataError):
        # no file, or no gene first explained at two or more steps (empty table)
        indirect = pd.DataFrame()
    status = {}
    for u, v, s in edges[['node_a', 'node_b', 'status']].itertuples(index=False):
        for n in (u, v):
            status.setdefault(n, {'gained': 0, 'lost': 0, 'shared': 0})[s] += 1
    rows = []
    supported, extra = {}, {}
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
        explained = isinstance(level, str) and level not in ('', 'unexplained')
        proximity = ''
        if explained:
            change = describe_marks(r, marks, level, str(r.get('explaining_marks', '')).split(','))
            explanation = '{} {}'.format(change or r.get('explaining_marks', ''), level_location(level))
            if level_hops(level) >= 2 and drivers:
                explanation += ' ({})'.format(drivers)
            if 'proximity_support' in r:
                proximity = proximity_sentence(level_hops(level), r['proximity_support'],
                                               r.get('proximity_support_ci', ''))
                if proximity:
                    explanation += '; ' + proximity
        else:
            explanation = 'no concordant chromatin change'
        # explanation cascade of egaDifferential: the layers and marks at the gene and at step 1, chromatin states
        detail = r.get('explanation_detail', '')
        if isinstance(detail, str) and detail:
            explanation = detail
            at = str(r.get('explained_at', ''))
            if at.startswith('step') and at[4:].isdigit() and int(at[4:]) >= 2 and drivers and \
                    r.get('explained_type', '') == 'marks':
                explanation += ' (nodes: {})'.format(drivers)
        states = []
        for column in ('promoter_state_summary', 'step1_state_summary'):
            text = r.get(column, '')
            if column == 'step1_state_summary' and r.get('explained_type', '') not in ('state', 'poised unresolved'):
                continue
            if isinstance(text, str) and text and text not in explanation:
                states.append(text)
        if states:
            explanation += '; ' + '; '.join(states)
        summary ='{} {}{:.1f}; interactions {} lost, {} gained, {} shared{}; {}'.format(
            gene, '+' if r['rna_log2FC'] > 0 else '', r['rna_log2FC'], counts['lost'], counts['gained'],
            counts['shared'], '; community {}'.format(transition) if transition else '', explanation)
        row = dict(node=node, gene=gene, rna_log2FC=r['rna_log2FC'], interactions_lost=counts['lost'],
                   interactions_gained=counts['gained'], interactions_shared=counts['shared'],
                   community_status=fate, community_transition=transition, explanation_level=level,
                   explaining_marks=r.get('explaining_marks', ''), indirect_drivers=drivers,
                   summary=summary)
        for column in ('explained_at', 'explained_type', 'explained_by', 'promoter_transition', 'step1_transition'):
            if column in r:
                row[column] = r[column]
        if 'explaining_hops' in r:
            row['explaining_hops'] = level_hops(level) if explained else np.nan
        if 'proximity_support' in r:
            row['proximity_support'] = r['proximity_support']
            row['proximity_support_ci'] = r.get('proximity_support_ci', '')
        if node in all_levels:
            for roles, column in (('declared', 'all_levels'), ('learned', 'all_levels_learned')):
                if 'detail_{}_mean'.format(roles) in all_levels[node]:
                    row[column] = all_level_text(gene, r['rna_log2FC'], all_levels[node], roles)
            row.update(term_columns(all_levels[node]))
            supported.update({column: None for column in support_columns(all_levels[node])})
            extra[node] = support_columns(all_levels[node])
        rows.append(row)
    table = pd.DataFrame(rows)
    # support of the contacts that carry the explanation: appended after the other columns
    for column in supported:
        table[column] = [extra.get(n, {}).get(column, '' if column in ('first_level_support', 'support_note',
                                                                        'threshold_sensitive_levels', 'provenance')
                                              else False) for n in table['node']]
    return table.sort_values('rna_log2FC', key=np.abs, ascending=False)


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


def add_tad_crossing(genes, graph, args, names):
    """Columns gained / lost / shared_crossing and tad_crossing (sentence), appended to the summary."""
    from .tads import Analysis
    analysis = Analysis(args.prefix, args.tads, graph=graph, names=names, **tad_kwargs(args))
    rows = {}
    for node in genes['node']:
        if node in analysis.genes and node not in rows:
            rows[node] = analysis.crossing.gene_row(node, '')
    for s in ('gained', 'lost', 'shared'):
        genes['{}_crossing'.format(s)] = [rows[n]['{}_crossing'.format(s)] if n in rows else np.nan
                                          for n in genes['node']]
    sentences = []
    for node, gene in zip(genes['node'], genes['gene']):
        text = rows[node]['sentence'] if node in rows else ''
        sentences.append(text.replace('{}: '.format(node), '{}: '.format(gene), 1) if text else '')
    genes['tad_crossing'] = sentences
    genes['summary'] = [s if not t else '{}; boundaries: {}'.format(s, t.split(': ', 1)[1].rstrip('.'))
                        for s, t in zip(genes['summary'], sentences)]
    return genes


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
    if args.tads:
        genes = add_tad_crossing(genes, graph, args, names)
    genes.to_csv(args.outputPrefix + '_genes.tsv', sep='\t', index=False)
    communities = community_report(graph, args.prefix, changes, marks, args.fdr)
    communities.to_csv(args.outputPrefix + '_communities.tsv', sep='\t', index=False)
    log.info('%d genes, %d communities; edge-chromatin association:\n%s', len(genes), len(communities),
             association[['mark', 'interactions', 'n', 'median_change', 'shared_median_change',
                          'concordance_vs_shared', 'mannwhitney_p_expected_direction']].to_string(index=False))


if __name__ == '__main__':
    main()
