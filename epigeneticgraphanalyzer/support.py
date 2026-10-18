"""Threshold sensitivity of the union edges and provenance of the statements of a differential run.

Support of a union edge (egaDifferential, column and graph attribute 'support')
    stable               called in both conditions (shared), or gained / lost while the other condition holds the
                         pair below the floor: its score in every replicate is below --supportFloor, or the pair is
                         absent from replicate files that list scores down to the floor (absence then means a score
                         below the lowest score of the file);
    threshold-sensitive  gained / lost while the other condition lists the pair at a score of at least the floor in some
                         replicate (a sub-threshold score, e.g. 3 <= score < 5, or a call in only some replicates that
                         the replicate rule did not keep), or the pair is listed in only some replicates of its own
                         condition;
    not testable         gained / lost and absent from a replicate file of the other condition that holds only scores
                         above the floor (e.g. published CHiCAGO exports with score >= 5): absence there cannot be told
                         apart from a sub-threshold score.

    The lowest score of every input file decides whether absence is informative (<prefix>_support_inputs.tsv).

Support of an explanation (explanation_levels)
    A level that explains a gene is evaluated again with the threshold-sensitive gained / lost contacts at weight 0 in
    both conditions, and with those and the untested ones at weight 0 (graph distances kept): stable when it still
    explains without both, not testable when it explains without the threshold-sensitive contacts only,
    threshold-sensitive when it no longer explains without them (it rests on threshold-sensitive contacts). A step 1
    state takes the class of the best path to its contacting fragment (a path is as weak as its weakest contact).
    The contacts counted per class are the edges of the shortest paths from the promoter to the nodes that carry the
    level (the nodes of a step contributing in the direction of the expression change; the promoter's own gained or
    lost contacts in the direction of the contact change; the contacting fragment of a step 1 state).

Provenance of a statement
    measured        signal values: mark fold enrichment, read counts, expression change
    called          interaction calls, chromatin state calls, element classes, structural variant joins
    derived         graph steps and paths, contact-weighted mark terms, communities, hubs, support classes
    interpretation  the explanation of an expression change by a level (explains, result, detail)
"""

import numpy as np
import pandas as pd

STABLE, SENSITIVE, UNTESTED = 'stable', 'threshold-sensitive', 'not testable'
SUPPORT_CLASSES = (STABLE, SENSITIVE, UNTESTED)
# order of the weakest link along a path (larger is weaker) and of the best path to a node (smaller is better)
SUPPORT_CODE = {STABLE: 0, UNTESTED: 1, SENSITIVE: 2}
CODE_SUPPORT = {v: k for k, v in SUPPORT_CODE.items()}

PROVENANCE_LEVELS = ('measured', 'called', 'derived', 'interpretation')
PROVENANCE_TAGS = {'measured': 'M', 'called': 'C', 'derived': 'D', 'interpretation': 'I'}
PROVENANCE_TEXT = {
    'measured': 'signal value (mark fold enrichment, read count, expression change)',
    'called': 'call on a threshold (interaction call, chromatin state, element class, structural variant join)',
    'derived': 'computed from calls and signals (graph steps and paths, contact-weighted mark terms, communities, '
               'hubs, support of the contacts)',
    'interpretation': 'attribution of the expression change to a level (explains, result, explanation text)'}
# provenance of the values of every level of the all-level table and the chain of a statement at that level
LEVEL_PROVENANCE = {'promoter_marks': ('measured', 'M>I'), 'promoter_contacts': ('called', 'C>I'),
                    'promoter_state': ('called', 'M>C>I'), 'structural_variant': ('called', 'C>I'),
                    'copy_number': ('measured', 'M>I'), 'step1_marks': ('derived', 'MC>D>I'),
                    'step1_state': ('called', 'M>C>D>I')}
STEP_PROVENANCE = ('derived', 'MC>D>I')


def level_provenance(level):
    """(provenance of the values, chain of tags) of a level of explanation_level_names."""
    return LEVEL_PROVENANCE.get(level, STEP_PROVENANCE)


def default_floor(threshold):
    """Floor of a sub-threshold score that counts as threshold-sensitive: 0.6 x the threshold (3 at CHiCAGO 5)."""
    return 0.6 * float(threshold)


def _pair_key(table):
    return (table['chrom'] + ':' + table['bait_start'].astype(np.int64).astype(str) + '-' +
            table['bait_end'].astype(np.int64).astype(str) + '|' + table['chrom'] + ':' +
            table['other_start'].astype(np.int64).astype(str) + '-' + table['other_end'].astype(np.int64).astype(str))


def input_detection(tables, paths, names, threshold, floor):
    """One row per input file: condition, replicate, file, pairs, min_score, max_score, pairs_below_threshold,
    holds_subthreshold (the file lists pairs below the calling threshold), resolves_floor (its lowest score is at most
    the floor, so a pair absent from it scored below the floor) and absence ('below floor' / 'not testable').
    `tables` and `paths` map a condition name to its replicate tables and files."""
    rows = []
    for name in names:
        for k, (table, path) in enumerate(zip(tables[name], paths[name]), start=1):
            score = table['score'].astype(float).to_numpy() if len(table) else np.array([])
            score = score[np.isfinite(score)]
            low = float(score.min()) if len(score) else np.nan
            resolves = bool(len(score)) and low <= floor
            rows.append(dict(condition=name, replicate=k, file=str(path), pairs=int(len(table)), min_score=low,
                             max_score=float(score.max()) if len(score) else np.nan,
                             pairs_below_threshold=int((score < threshold).sum()),
                             holds_subthreshold=bool(len(score)) and low < threshold, resolves_floor=resolves,
                             absence='below floor' if resolves else 'not testable'))
    return pd.DataFrame(rows)


def _replicate_scores(tables, keys):
    """keys x replicates matrix of the score of every pair (NaN where a replicate does not list it)."""
    out = np.full((len(keys), len(tables)), np.nan)
    index = pd.Index(keys)
    for j, table in enumerate(tables):
        if not len(table):
            continue
        scores = pd.Series(table['score'].astype(float).to_numpy(), index=_pair_key(table).to_numpy())
        scores = scores.groupby(level=0).max()
        out[:, j] = scores.reindex(index).to_numpy()
    return out


def edge_support(edges, tables, detection, names, threshold, floor):
    """Support class of every union edge (see the module documentation). `tables` maps a condition name to its
    replicate tables (as read, before the replicate rule), `detection` is input_detection. Returns a frame with the
    columns support, support_reason, other_max_score (highest score of the pair in the replicates of the condition in
    which it is not called; NaN for shared edges and pairs listed nowhere there) and replicates_floor_<name> (replicates
    of each condition listing the pair at a score of at least the floor), index as `edges`."""
    name_a, name_b = names
    keys = (edges['node_a'] + '|' + edges['node_b']).to_numpy()
    scores = {n: _replicate_scores(tables[n], keys) for n in names}
    listed = {n: np.isfinite(scores[n]) for n in names}
    at_floor = {n: np.nan_to_num(scores[n], nan=-np.inf) >= floor for n in names}
    resolves = {n: detection.loc[detection['condition'] == n].sort_values('replicate')['resolves_floor']
                .to_numpy(dtype=bool) for n in names}
    lowest = {n: detection.loc[detection['condition'] == n, 'min_score'].min() for n in names}
    status = edges['status'].to_numpy()
    support = np.empty(len(edges), dtype=object)
    reason = np.empty(len(edges), dtype=object)
    other_max = np.full(len(edges), np.nan)
    for i in range(len(edges)):
        if status[i] == 'shared':
            support[i], reason[i] = STABLE, 'called in both'
            continue
        own, other = (name_b, name_a) if status[i] == 'gained' else (name_a, name_b)
        s_other = scores[other][i]
        if np.isfinite(s_other).any():
            other_max[i] = float(np.nanmax(s_other))
        n_other, n_own = len(s_other), listed[own].shape[1]
        floor_other = int(at_floor[other][i].sum())
        listed_own = int(listed[own][i].sum())
        if floor_other:
            if np.nanmax(s_other) >= threshold:
                support[i] = SENSITIVE
                reason[i] = 'called in {} of {} replicates of {} (max score {:.3g}), not kept by the replicate ' \
                            'rule'.format(int((np.nan_to_num(s_other, nan=-np.inf) >= threshold).sum()), n_other,
                                          other, other_max[i])
            else:
                support[i] = SENSITIVE
                reason[i] = 'sub-threshold score {:.3g} in {} ({:g} <= score < {:g})'.format(
                    other_max[i], other, floor, threshold)
        elif listed_own < n_own:
            support[i] = SENSITIVE
            reason[i] = 'listed in {} of {} replicates of {}'.format(listed_own, n_own, own)
        else:
            known = listed[other][i] | resolves[other]
            if known.all():
                support[i] = STABLE
                reason[i] = 'score below {:g} in {}'.format(floor, other) if listed[other][i].any() else \
                    'absent from {} (files list scores down to {:.3g})'.format(other, lowest[other])
            else:
                support[i] = UNTESTED
                reason[i] = 'absent from {} input holding scores >= {:.3g} only'.format(
                    other, float(np.nanmin(detection.loc[(detection['condition'] == other) &
                                                         ~detection['resolves_floor'], 'min_score'])))
    out = pd.DataFrame({'support': support, 'support_reason': reason, 'other_max_score': other_max},
                       index=edges.index)
    for n in names:
        out['replicates_floor_' + n] = at_floor[n].sum(axis=1).astype(int)
    return out


def support_summary(edges):
    """Edges per status and support class with the fraction of the status: status, support, edges, fraction."""
    rows = []
    for status in ('gained', 'lost', 'shared'):
        part = edges[edges['status'] == status]
        for cls in SUPPORT_CLASSES:
            n = int((part['support'] == cls).sum())
            rows.append(dict(status=status, support=cls, edges=n, fraction=n / len(part) if len(part) else np.nan))
    changed = edges[edges['status'].isin(['gained', 'lost'])]
    for cls in SUPPORT_CLASSES:
        n = int((changed['support'] == cls).sum())
        rows.append(dict(status='gained+lost', support=cls, edges=n,
                         fraction=n / len(changed) if len(changed) else np.nan))
    return pd.DataFrame(rows)


def edge_code(graph, u, v):
    """Support code of a graph edge (SUPPORT_CODE); an edge without the attribute counts as stable when shared and
    as not testable otherwise."""
    data = graph.edges[u, v]
    value = data.get('support')
    if value in SUPPORT_CODE:
        return SUPPORT_CODE[value]
    return 0 if data.get('status', 'shared') == 'shared' else 1


def graph_has_support(graph):
    return any('support' in d for _, _, d in graph.edges(data=True))


class PathSupport:
    """Shortest paths from one promoter up to `max_hops` steps: the predecessors of every node on its shortest paths
    and the best support class over these paths (weakest edge of a path, best path of a node)."""

    def __init__(self, graph, source, max_hops):
        self.graph = graph
        self.source = source
        self.dist = {source: 0}
        self.preds = {}
        self.best = {source: 0}
        if source not in graph:
            return
        frontier = [source]
        for k in range(1, max_hops + 1):
            nxt = []
            for u in frontier:
                for w in graph.adj[u]:
                    d = self.dist.get(w)
                    if d is None:
                        self.dist[w] = k
                        self.preds[w] = [u]
                        nxt.append(w)
                    elif d == k:
                        self.preds[w].append(u)
            for w in nxt:
                self.best[w] = min(max(self.best[p], edge_code(graph, p, w)) for p in self.preds[w])
            frontier = nxt

    def node_class(self, node):
        code = self.best.get(node)
        return CODE_SUPPORT[code] if code is not None else ''

    def carrying_edges(self, targets):
        """Edges of the shortest paths from the source to `targets`, as (u, v) with the support code."""
        seen, edges, stack = set(), {}, [t for t in targets if t in self.dist]
        while stack:
            m = stack.pop()
            if m in seen:
                continue
            seen.add(m)
            for p in self.preds.get(m, ()):
                key = (p, m) if p <= m else (m, p)
                edges[key] = edge_code(self.graph, p, m)
                stack.append(p)
        return edges

    def summary(self, targets):
        """Support of a level carried by the nodes `targets`: dict(support, edges per class, nodes per class)."""
        targets = [t for t in dict.fromkeys(targets) if t in self.dist and t != self.source]
        if not targets:
            return None
        edges = self.carrying_edges(targets)
        codes = [self.best[t] for t in targets]
        return _summary(list(edges.values()), codes)


def _summary(edge_codes, node_codes):
    counts = {cls: int(sum(1 for c in edge_codes if c == SUPPORT_CODE[cls])) for cls in SUPPORT_CLASSES}
    nodes = {cls: int(sum(1 for c in node_codes if c == SUPPORT_CODE[cls])) for cls in SUPPORT_CLASSES}
    if nodes[STABLE]:
        support = STABLE
    elif nodes[SENSITIVE] == len(node_codes):
        support = SENSITIVE
    else:
        support = UNTESTED
    return dict(support=support, edges=counts, nodes=nodes, n_nodes=len(node_codes))


def contact_summary(graph, node, contact_log2):
    """Support of the promoter's own contacts: the gained contacts when the summed strength rises, the lost ones when
    it falls (both when unchanged); each contact is its own path."""
    if node not in graph:
        return None
    wanted = ('gained',) if contact_log2 > 0 else ('lost',) if contact_log2 < 0 else ('gained', 'lost')
    codes = [edge_code(graph, node, w) for w in graph.adj[node]
             if graph.edges[node, w].get('status', 'shared') in wanted]
    if not codes:
        return None
    return _summary(codes, codes)


def support_text(summary, explains, cls=''):
    """Text of the support of a level: 'rests on 2 threshold-sensitive contacts' (explaining level of class
    threshold-sensitive), else the carrying contacts per class and what remains without the threshold-sensitive
    (and untested) ones; '' without carrying contacts."""
    if not summary:
        return ''
    edges = summary['edges']
    if explains and cls == SENSITIVE:
        n = edges[SENSITIVE]
        return 'rests on {} threshold-sensitive contact{}'.format(n, '' if n == 1 else 's') if n else \
            'rests on threshold-sensitive contacts'
    parts = ['{} {}'.format(edges[c], c) for c in SUPPORT_CLASSES if edges[c]]
    total = sum(edges.values())
    text = '{} contact{}'.format(', '.join(parts), '' if total == 1 else 's')
    if explains and cls == STABLE and (edges[SENSITIVE] or edges[UNTESTED]):
        text += '; explains without the threshold-sensitive and untested ones'
    elif explains and cls == UNTESTED:
        text += '; explains without the threshold-sensitive ones, not without the untested ones'
    return text


def counterfactual_values(graph, operators, mark_vectors, drop, pseudo_count=1.0):
    """Level values of the sources of `operators` when the edges of the support classes `drop` get weight 0 in both
    conditions; the graph distances (nodes of a step) stay those of `graph`. Returns {level: {mark: array over the
    sources}} (the formula of gene_features: log2 (H_B v_B + pseudo) / (H_A v_A + pseudo)) and 'contact' (the change
    of the summed strength of the own contacts, contact_direct)."""
    import networkx as nx
    from .differential import ActivityOperators
    reduced = nx.Graph()
    reduced.add_nodes_from(operators.nodes)
    for u, v, data in graph.edges(data=True):
        zero = data.get('support') in drop
        reduced.add_edge(u, v, r_a=0.0 if zero else data['r_a'], r_b=0.0 if zero else data['r_b'])
    reduced_ops = ActivityOperators(reduced, operators.sources, max_hops=operators.max_hops)
    out = {}
    for level in reduced_ops.levels:
        H_a, H_b, _ = reduced_ops.level(level)
        out[level] = {m: np.log2((H_b @ vb + pseudo_count) / (H_a @ va + pseudo_count))
                      for m, (va, vb) in mark_vectors.items()}
    ones = np.ones(len(reduced_ops.nodes))
    H_a, H_b, _ = reduced_ops.level('direct')
    out['contact'] = np.log2((H_b @ ones + pseudo_count) / (H_a @ ones + pseudo_count))
    return out


def all_levels_provenance(columns):
    """Provenance of every column of the all-level table (explanation_levels): column, provenance, tag."""
    rows = []
    for column in columns:
        base = column.split('_')[0]
        if column in ('node', 'genes', 'rank', 'level', 'hops', 'provenance', 'provenance_chain'):
            kind = 'identifier'
        elif column == 'rna_log2FC' or column.startswith(('mark_changes', 'mirror_changes', 'copy_number')):
            kind = 'measured'
        elif column.startswith(('contacts_', 'transition_', 'state_fragment_', 'sv_')):
            kind = 'called'
        elif column.startswith(('explains_', 'result_', 'detail_', 'concordant_marks_', 'signs_')):
            kind = 'interpretation'
        elif base in ('p', 'q') or column.startswith(('nodes_tested', 'contact_log2', 'proximity_support', 'n_concordant',
                                                      'strongest', 'beyond_mirror', 'support', 'edges_',
                                                      'threshold_sensitive_only_')):
            kind = 'derived'
        else:
            kind = 'derived'
        rows.append(dict(column=column, provenance=kind, tag=PROVENANCE_TAGS.get(kind, '')))
    return pd.DataFrame(rows)
