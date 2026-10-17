"""Mark roles and chromatin states of the differential mode (egaDifferential).

Roles
    Every mark has a role: 'active', 'repressive', 'enhancer' (active-type at fragments that are not baited
    promoters), 'architectural', 'neutral' or 'unknown'. Roles are declared (--markRoles, or --markSign with +1 =
    active, -1 = repressive) and can be learned from the data (learn_mark_roles): the association of the change of a
    mark with the expression change at the gene (step 0) and at its directly contacting fragments (step 1), with
    chromosome block-bootstrap intervals. The sign of a role (+1 active / enhancer, -1 repressive, 0 otherwise) is
    the expected direction of a mark change on expression in the explanation cascade.

States
    A mark is present on a fragment in a condition if its mean signal over the fragment reaches a threshold (fold
    enrichment, default 2). From the roles: a fragment carries an active-type mark (at a bait: an active-role mark;
    at another fragment: an active- or enhancer-role mark) and / or a repressive-role mark. Baits: 'bivalent' (both),
    'active', 'repressed', 'unmarked'; other fragments: 'poised enhancer' (both), 'active enhancer', 'repressed',
    'unmarked'. No mark name enters these rules; only the roles do.

Learned states
    learned_state_transitions describes every fragment by the combination of present marks (all marks, whatever
    their role) in each condition, keeps the combinations seen at least `min_fragments` times and rates each
    combination and each transition between combinations by the expression change of the genes whose promoter
    (step 0) or contacting fragments (step 1) make it.
"""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from .validation import block_resample_items, block_resamples

log = logging.getLogger(__name__)

MARK_ROLES = ('active', 'repressive', 'enhancer', 'architectural', 'neutral', 'unknown')
ROLE_SIGN = {'active': 1, 'enhancer': 1, 'repressive': -1, 'architectural': 0, 'neutral': 0, 'unknown': 0}
NOT_CALLABLE = 'states not callable'
POISED_STATES = ('bivalent', 'poised enhancer')
# (active-type mark present, repressive mark present) of every state
STATE_FLAGS = {'bivalent': (1, 1), 'active': (1, 0), 'repressed': (0, 1), 'unmarked': (0, 0),
               'poised enhancer': (1, 1), 'active enhancer': (1, 0)}


def bigwig_genome_mean(path):
    """Mean signal over the covered bases of a bigWig file (header sumData / nBasesCovered)."""
    import pyBigWig
    handle = pyBigWig.open(path)
    try:
        header = handle.header()
    finally:
        handle.close()
    return float(header['sumData']) / max(1, int(header['nBasesCovered']))


def parse_mark_roles(entries):
    """mark -> role from 'mark=role' entries and / or TSV files with the columns mark and role."""
    roles = {}
    for entry in entries or []:
        if '=' in entry:
            mark, role = entry.split('=', 1)
            roles[mark] = role
        else:
            table = pd.read_csv(entry, sep='\t', comment='#')
            if not {'mark', 'role'} <= set(table.columns):
                raise ValueError('{} needs the columns mark and role'.format(entry))
            roles.update(dict(zip(table['mark'].astype(str), table['role'].astype(str))))
    for mark, role in roles.items():
        if role not in MARK_ROLES:
            raise ValueError('Unknown role {} of {} (roles: {})'.format(role, mark, ', '.join(MARK_ROLES)))
    return roles


def resolve_mark_roles(marks, mark_sign=(), mark_roles=None):
    """Role of every mark. Marks without a declaration are 'active' (sign +1, as with --markSign); --markSign
    entries 'mark=+1|-1|0' give active / repressive / neutral; `mark_roles` (parse_mark_roles) take precedence."""
    roles = {m: 'active' for m in marks}
    for entry in mark_sign or ():
        mark, sign = entry.split('=')
        if mark in roles:
            sign = int(sign)
            roles[mark] = 'active' if sign > 0 else 'repressive' if sign < 0 else 'neutral'
    for mark, role in (mark_roles or {}).items():
        if mark in roles:
            roles[mark] = role
        else:
            log.warning('role of %s declared, but %s is not among the marks', mark, mark)
    return roles


def role_signs(roles):
    """mark -> expected direction of its change on expression (+1, -1, 0)."""
    return {m: ROLE_SIGN[r] for m, r in roles.items()}


def role_marks(roles):
    return {kind: [m for m, r in roles.items() if r == kind] for kind in ('active', 'enhancer', 'repressive')}


def states_callable(roles):
    """{'promoter': bivalent promoters callable (active and repressive marks), 'step1': poised enhancers callable
    (active or enhancer marks, and repressive marks)}."""
    marks = role_marks(roles)
    return {'promoter': bool(marks['active']) and bool(marks['repressive']),
            'step1': bool(marks['active'] or marks['enhancer']) and bool(marks['repressive'])}


def _any_present(signal, marks, threshold, n):
    if not marks:
        return np.zeros(n, dtype=bool)
    return np.any([np.asarray(signal[m], dtype=float) >= threshold for m in marks], axis=0)


def fragment_states(signal, roles, threshold, is_bait):
    """State of every fragment in one condition (`signal`: mark -> array, `is_bait`: bool array)."""
    is_bait = np.asarray(is_bait, dtype=bool)
    marks = role_marks(roles)
    n = len(is_bait)
    active = np.where(is_bait, _any_present(signal, marks['active'], threshold, n),
                      _any_present(signal, marks['active'] + marks['enhancer'], threshold, n))
    repressive = _any_present(signal, marks['repressive'], threshold, n)
    bait = np.select([active & repressive, active, repressive], ['bivalent', 'active', 'repressed'], 'unmarked')
    other = np.select([active & repressive, active, repressive], ['poised enhancer', 'active enhancer', 'repressed'],
                      'unmarked')
    return np.where(is_bait, bait, other).astype(object)


def transition_direction(state_a, state_b):
    """+1 towards activity (active-type mark gained or repressive mark lost), -1 towards repression, 0 for no
    change or opposite changes (e.g. unmarked -> bivalent)."""
    a, b = STATE_FLAGS[state_a], STATE_FLAGS[state_b]
    return int(np.sign((b[0] - a[0]) - (b[1] - a[1])))


def transition_label(state_a, state_b):
    """'bivalent -> active: repressive mark lost', named from the roles."""
    label = '{} -> {}'.format(state_a, state_b)
    if state_a == state_b:
        return label
    a, b = STATE_FLAGS[state_a], STATE_FLAGS[state_b]
    changes = ['{} mark {}'.format(kind, 'gained' if b[k] > a[k] else 'lost')
               for k, kind in enumerate(('active', 'repressive')) if a[k] != b[k]]
    return label + ': ' + ', '.join(changes)


def label_states(label):
    """(state A, state B) of a transition label."""
    a, b = label.split(':')[0].split(' -> ')
    return a, b


def mark_crossings(values_a, values_b, threshold, min_change=0.5, pseudo_count=1.0):
    """Marks whose presence differs between the conditions: (gained, lost, supported). `values_a` / `values_b` map
    mark -> signal of one fragment; supported: every crossing mark changes by at least `min_change`
    (log2 (B + pseudo) / (A + pseudo)) in the direction of the crossing, so that a signal moving from just below to
    just above the threshold does not make a transition."""
    gained, lost, supported = [], [], True
    for mark in values_a:
        a, b = float(values_a[mark]), float(values_b[mark])
        if (a >= threshold) == (b >= threshold):
            continue
        change = np.log2((b + pseudo_count) / (a + pseudo_count))
        if b >= threshold:
            gained.append(mark)
        else:
            lost.append(mark)
            change = -change
        if change < min_change:
            supported = False
    return gained, lost, supported


def _state_marks(roles, bait):
    """Marks whose presence makes the state of a bait (active, repressive) or of another fragment (active,
    enhancer, repressive), grouped by flag: ([active-type marks], [repressive marks])."""
    marks = role_marks(roles)
    return (marks['active'] if bait else marks['active'] + marks['enhancer']), marks['repressive']


def _transition(values_a, values_b, state_a, state_b, roles, bait, threshold, min_change, pseudo_count):
    """Crossing marks of the flags that change (gained, lost, supported)."""
    active, repressive = _state_marks(roles, bait)
    fa, fb = STATE_FLAGS[state_a], STATE_FLAGS[state_b]
    used = (active if fa[0] != fb[0] else []) + (repressive if fa[1] != fb[1] else [])
    return mark_crossings({m: values_a[m] for m in used}, {m: values_b[m] for m in used}, threshold, min_change,
                          pseudo_count)


def _signal_text(values):
    return ', '.join('{} {:.1f}'.format(m, v) for m, v in values.items())


def state_sentence(where, state_a, state_b, values_a, values_b, gained, lost, names=('A', 'B')):
    """One-line description of a state transition naming the marks, e.g. 'promoter bivalent in hESC (H3K4me3 4.1,
    H3K27me3 3.2), resolved to active in hNEC by loss of H3K27me3'."""
    name_a, name_b = names
    if state_a == state_b:
        text = '{} {} in {} ({}) and in {} ({})'.format(where, state_a, name_a, _signal_text(values_a), name_b,
                                                         _signal_text(values_b))
        return text + (': poised, trigger not measured' if state_a in POISED_STATES else '')
    by = ' and '.join(['loss of ' + m for m in lost] + ['gain of ' + m for m in gained])
    if state_a in POISED_STATES:
        text = '{} {} in {} ({}), resolved to {} in {}'.format(where, state_a, name_a, _signal_text(values_a),
                                                               state_b, name_b)
    else:
        text = '{} {} in {} ({}), {} in {} ({})'.format(where, state_a, name_a, _signal_text(values_a), state_b,
                                                        name_b, _signal_text(values_b))
    return text + (' by ' + by if by else '')


def gene_states(classes, graph, state_signals, roles, threshold=2.0, min_change=0.5, pseudo_count=1.0,
                names=('A', 'B')):
    """Chromatin state transitions of every differential gene at its promoter fragment (step 0) and at its directly
    contacting fragments (step 1).

    `state_signals` maps mark -> (node -> signal in A, node -> signal in B) on the scale of the presence threshold
    (fold enrichment); `roles` maps mark -> role. Promoter: fragment_states in both conditions, label
    transition_label (e.g. 'bivalent -> active: repressive mark lost'), direction (transition_direction) and support
    (mark_crossings of the marks whose role flag changes). promoter_explains: a supported transition from or to
    'bivalent' whose direction agrees with the expression change (bivalent -> active with up, bivalent ->
    repressed with down, active -> bivalent with down, repressed -> bivalent with up); promoter_unresolved:
    bivalent in both conditions.

    Step 1: the states of every fragment in contact with the promoter in either condition (union graph; baits among
    them take the promoter states, other fragments the enhancer states). step1_explains: at least one contacting
    fragment with a supported transition from or to a poised state (bivalent, poised enhancer) concordant with the
    expression change; step1_transition is then the most frequent such label (fragments counted in
    step1_transition_fragments, the one with the largest change in step1_state_fragment). step1_unresolved:
    otherwise a contacting fragment poised in both conditions. step1_transitions lists all labels of changed or
    poised fragments with their counts. Levels without the roles of their states get 'states not callable'.
    promoter_state_summary and step1_state_summary describe the transition in words, naming the marks.
    """
    callable_ = states_callable(roles)
    nodes = classes['node'].tolist()
    direction = np.sign(classes['rna_log2FC'].to_numpy(dtype=float))
    out = pd.DataFrame(index=pd.Index(nodes, name='node'))
    used_marks = [m for m, r in roles.items() if r in ('active', 'enhancer', 'repressive')]

    def values(node, marks, condition):
        return {m: float(state_signals[m][condition].get(node, 0.0)) for m in marks}

    def states_of(fragments, bait):
        signal_a = {m: np.array([state_signals[m][0].get(f, 0.0) for f in fragments]) for m in used_marks}
        signal_b = {m: np.array([state_signals[m][1].get(f, 0.0) for f in fragments]) for m in used_marks}
        return fragment_states(signal_a, roles, threshold, bait), fragment_states(signal_b, roles, threshold, bait)

    if callable_['promoter']:
        active, repressive = _state_marks(roles, True)
        promoter_marks = list(dict.fromkeys(active + repressive))
        state_a, state_b = states_of(nodes, np.ones(len(nodes), dtype=bool))
        for m in promoter_marks:
            out['promoter_{}_A'.format(m)] = [state_signals[m][0].get(n, 0.0) for n in nodes]
            out['promoter_{}_B'.format(m)] = [state_signals[m][1].get(n, 0.0) for n in nodes]
        labels, directions, supported, summaries = [], [], [], []
        for i, node in enumerate(nodes):
            va, vb = values(node, promoter_marks, 0), values(node, promoter_marks, 1)
            gained, lost, ok = _transition(va, vb, state_a[i], state_b[i], roles, True, threshold, min_change,
                                           pseudo_count)
            labels.append(transition_label(state_a[i], state_b[i]))
            directions.append(transition_direction(state_a[i], state_b[i]))
            supported.append(ok)
            interesting = state_a[i] != state_b[i] or state_a[i] in POISED_STATES
            summaries.append(state_sentence('promoter', state_a[i], state_b[i], va, vb, gained, lost, names)
                             if interesting else '')
        directions = np.array(directions, dtype=int)
        supported = np.array(supported, dtype=bool)
        poised = np.isin(state_a, POISED_STATES) | np.isin(state_b, POISED_STATES)
        out['promoter_state_A'] = state_a
        out['promoter_state_B'] = state_b
        out['promoter_transition'] = labels
        out['promoter_direction'] = directions
        out['promoter_supported'] = supported
        out['promoter_explains'] = poised & (state_a != state_b) & supported & (directions * direction > 0)
        out['promoter_unresolved'] = (state_a == 'bivalent') & (state_b == 'bivalent')
        out['promoter_state_summary'] = summaries
    else:
        out['promoter_transition'] = NOT_CALLABLE
        out['promoter_explains'] = False
        out['promoter_unresolved'] = False
        out['promoter_state_summary'] = ''
    if callable_['step1']:
        neighbours = {n: sorted(graph.neighbors(n)) if n in graph else [] for n in nodes}
        fragments = sorted({m for ns in neighbours.values() for m in ns})
        position = {m: i for i, m in enumerate(fragments)}
        is_bait = np.array([bool(graph.nodes[f].get('bait', False)) for f in fragments], dtype=bool)
        state_a, state_b = states_of(fragments, is_bait) if fragments else (np.zeros(0, object), np.zeros(0, object))
        columns = {key: [] for key in ('step1_poised_A', 'step1_poised_B', 'step1_transitions', 'step1_transition',
                                       'step1_transition_fragments', 'step1_explains', 'step1_unresolved',
                                       'step1_state_fragment', 'step1_state_summary')}
        for i, node in enumerate(nodes):
            entries = []
            for fragment in neighbours[node]:
                j = position[fragment]
                sa, sb = state_a[j], state_b[j]
                if sa == sb and sa not in POISED_STATES:
                    continue
                active, repressive = _state_marks(roles, is_bait[j])
                shown = list(dict.fromkeys(active + repressive))
                va, vb = values(fragment, shown, 0), values(fragment, shown, 1)
                gained, lost, ok = _transition(va, vb, sa, sb, roles, is_bait[j], threshold, min_change, pseudo_count)
                d = transition_direction(sa, sb)
                strength = max([abs(np.log2((vb[m] + pseudo_count) / (va[m] + pseudo_count))) for m in gained + lost]
                               + [0.0])
                entries.append(dict(fragment=fragment, label=transition_label(sa, sb), a=sa, b=sb,
                                    poised=sa in POISED_STATES or sb in POISED_STATES,
                                    concordant=sa != sb and ok and d * direction[i] > 0, strength=strength, va=va,
                                    vb=vb, gained=gained, lost=lost))
            counts = pd.Series([e['label'] for e in entries], dtype=object).value_counts()
            columns['step1_poised_A'].append(sum(1 for f in neighbours[node] if state_a[position[f]] in POISED_STATES))
            columns['step1_poised_B'].append(sum(1 for f in neighbours[node] if state_b[position[f]] in POISED_STATES))
            columns['step1_transitions'].append(';'.join('{}={}'.format(k, v) for k, v in counts.items()))
            concordant = [e for e in entries if e['poised'] and e['concordant']]
            unresolved = [e for e in entries if e['a'] == e['b'] and e['a'] in POISED_STATES]
            other = [e for e in entries if e['poised'] and e['a'] != e['b']]
            chosen, explains, open_poised = [], False, False
            for group, flag in ((concordant, 'explains'), (unresolved, 'unresolved'), (other, '')):
                if group:
                    top = pd.Series([e['label'] for e in group]).value_counts().index[0]
                    chosen = sorted([e for e in group if e['label'] == top], key=lambda e: -e['strength'])
                    explains, open_poised = flag == 'explains', flag == 'unresolved'
                    break
            columns['step1_explains'].append(explains)
            columns['step1_unresolved'].append(open_poised)
            columns['step1_transition'].append(chosen[0]['label'] if chosen else 'no poised fragment')
            columns['step1_transition_fragments'].append(len(chosen))
            columns['step1_state_fragment'].append(chosen[0]['fragment'] if chosen else '')
            if chosen:
                e = chosen[0]
                where = 'contacting fragment {}{}'.format(
                    e['fragment'], ' ({} of {})'.format(len(chosen), len(neighbours[node])) if len(chosen) > 1 else '')
                columns['step1_state_summary'].append(state_sentence(where, e['a'], e['b'], e['va'], e['vb'],
                                                                     e['gained'], e['lost'], names))
            else:
                columns['step1_state_summary'].append('')
        for key, value in columns.items():
            out[key] = value
    else:
        out['step1_transition'] = NOT_CALLABLE
        out['step1_explains'] = False
        out['step1_unresolved'] = False
        out['step1_state_summary'] = ''
    return out


def state_transition_table(classes, states, level='promoter'):
    """Number of differential genes per state transition at the promoter ('promoter') or at step 1 ('step1', the
    label chosen per gene by gene_states): genes, genes_up, genes_down, direction of the transition (+1 towards
    activity, -1 towards repression, 0 unchanged or undetermined), concordant (direction agrees with the expression
    change), supported_concordant (also supported by the mark changes; promoter only) and fraction_concordant."""
    column = level + '_transition'
    table = states.reindex(classes['node'].to_numpy())
    rna = np.sign(classes['rna_log2FC'].to_numpy(dtype=float))
    labels = table[column].astype(str).to_numpy()
    rows = []
    for label in pd.unique(labels):
        mask = labels == label
        d = 0 if label in (NOT_CALLABLE, 'no poised fragment') else transition_direction(*label_states(label))
        concordant = mask & (d * rna > 0)
        row = dict(level=level, transition=label, direction=d, genes=int(mask.sum()),
                   genes_up=int((mask & (rna > 0)).sum()), genes_down=int((mask & (rna < 0)).sum()),
                   concordant=int(concordant.sum()) if d else np.nan)
        if level == 'promoter' and 'promoter_supported' in table:
            row['supported_concordant'] = int((concordant & table['promoter_supported'].fillna(False)
                                               .to_numpy(dtype=bool)).sum()) if d else np.nan
        row['fraction_concordant'] = row['concordant'] / row['genes'] if d and row['genes'] else np.nan
        rows.append(row)
    result = pd.DataFrame(rows)
    if len(result):
        result = result.sort_values(['genes', 'transition'], ascending=[False, True]).reset_index(drop=True)
    return result


def _chromosomes(nodes):
    return np.array([str(n).split(':')[0] for n in nodes], dtype=object)


def _spearman(x, y):
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3 or np.all(x[keep] == x[keep][0]) or np.all(y[keep] == y[keep][0]):
        return np.nan
    return float(stats.spearmanr(x[keep], y[keep])[0])


def _concordance(x, y, x_threshold, y_threshold):
    """Fraction of items with |x| >= x_threshold and |y| >= y_threshold whose signs agree (0.5 without relation)."""
    keep = (np.abs(x) >= x_threshold) & (np.abs(y) >= y_threshold)
    return float(np.mean(np.sign(x[keep]) == np.sign(y[keep]))) if keep.any() else np.nan


def learn_mark_roles(features, marks, declared=None, chromatin_threshold=1.0, rna_threshold=1.0, n_boot=1000,
                     seed=0):
    """Role of every mark learned from the data.

    Per mark, over the genes of `features` (gene_features with contact_direct): Spearman rho of the mark change at
    the promoter (step 0) and summed over the direct partners (step 1, genes with partners) with the expression
    change, and of the promoter mark change with the change of the promoter's own interactions (contact_direct);
    concordance at step 0 and step 1: the fraction of genes with |expression change| >= `rna_threshold` and |mark
    change| >= `chromatin_threshold` whose signs agree (0.5 without relation; a repressive mark goes below 0.5).
    95% intervals from `n_boot` chromosome block-bootstrap resamples (validation.block_resample_items). learned_role:
    'active' if the interval of the expression rho at step 0 lies above 0, 'repressive' if below; where it includes 0
    the interval at step 1 decides; otherwise 'no consistent effect'. steps_disagree: the step-0 and step-1
    intervals exclude 0 on opposite sides. disagrees: declared and learned roles have opposite signs.
    """
    data = features.dropna(subset=['rna_log2FC']).reset_index(drop=True)
    rna = data['rna_log2FC'].to_numpy(dtype=float)
    contact = np.nan_to_num(data['contact_direct'].to_numpy(dtype=float)) if 'contact_direct' in data else None
    has_partners = data['n_direct'].to_numpy() > 0 if 'n_direct' in data else np.ones(len(data), dtype=bool)
    chroms = _chromosomes(data['node'])
    resamples = list(block_resample_items(chroms, n_boot, seed)) if len(data) else []
    rows = []
    for mark in marks:
        step0 = np.nan_to_num(data[mark + '_promoter'].to_numpy(dtype=float))
        step1 = np.nan_to_num(data[mark + '_direct'].to_numpy(dtype=float)) if mark + '_direct' in data else None

        def statistics(index):
            values = {'rho_step0': _spearman(step0[index], rna[index]),
                      'concordance_step0': _concordance(step0[index], rna[index], chromatin_threshold, rna_threshold)}
            if step1 is not None:
                sub = index[has_partners[index]]
                values['rho_step1'] = _spearman(step1[sub], rna[sub])
                values['concordance_step1'] = _concordance(step1[sub], rna[sub], chromatin_threshold, rna_threshold)
            if contact is not None:
                values['rho_contact'] = _spearman(step0[index], contact[index])
            return values

        observed = statistics(np.arange(len(data)))
        boot = [statistics(items) for items in resamples]
        row = {'mark': mark, 'declared_role': (declared or {}).get(mark, ''), 'n_genes': len(data)}
        for key, value in observed.items():
            values = np.array([b[key] for b in boot], dtype=float)
            row[key] = value
            row[key + '_ci_low'] = float(np.nanpercentile(values, 2.5)) if np.isfinite(values).any() else np.nan
            row[key + '_ci_high'] = float(np.nanpercentile(values, 97.5)) if np.isfinite(values).any() else np.nan
        # the gene's own fragment decides; step 1 where step 0 is undecided
        signs = {k: 1 if row[k + '_ci_low'] > 0 else -1 if row[k + '_ci_high'] < 0 else 0
                 for k in ('rho_step0', 'rho_step1') if k in observed}
        decided = signs.get('rho_step0', 0) or signs.get('rho_step1', 0)
        row['learned_role'] = {1: 'active', -1: 'repressive'}.get(decided, 'no consistent effect')
        row['steps_disagree'] = bool(signs.get('rho_step0', 0) * signs.get('rho_step1', 0) < 0)
        declared_sign = ROLE_SIGN.get(row['declared_role'], 0)
        learned_sign = {'active': 1, 'repressive': -1}.get(row['learned_role'], 0)
        row['learned_sign'] = learned_sign
        row['disagrees'] = bool(declared_sign and learned_sign and declared_sign != learned_sign)
        rows.append(row)
    return pd.DataFrame(rows)


def learned_roles(table, declared):
    """Roles for the cascade from learn_mark_roles: a learned active mark keeps a declared enhancer role (else
    'active'), a learned repressive mark is 'repressive', a mark without a consistent effect 'neutral'."""
    roles = {}
    for row in table.itertuples():
        if row.learned_sign > 0:
            roles[row.mark] = 'enhancer' if declared.get(row.mark) == 'enhancer' else 'active'
        elif row.learned_sign < 0:
            roles[row.mark] = 'repressive'
        else:
            roles[row.mark] = 'neutral'
    return roles


def mark_combinations(present, marks):
    """Combination label of every fragment from a dict mark -> bool array: present marks joined by '+', 'none'."""
    n = len(next(iter(present.values()))) if present else 0
    labels = np.full(n, '', dtype=object)
    for mark in marks:
        labels = np.where(present[mark], np.where(labels == '', mark, labels + '+' + mark), labels)
    return np.where(labels == '', 'none', labels).astype(object)


def learned_state_transitions(features, graph, state_signals, learned_signs, threshold=2.0, rna_threshold=1.0,
                              min_fragments=50, n_boot=1000, seed=0, min_genes=5):
    """Recurring mark combinations of fragments and their transitions, rated by the expression changes that follow.

    Every fragment of the union graph is described in each condition by the combination of marks present
    (state_signals >= `threshold`; all marks, whatever their role). Combinations seen at least `min_fragments` times
    over both conditions are kept, the rest are pooled as 'other'. A transition is the pair (combination in A,
    combination in B) of a fragment. Rows of kind 'transition' (level step0: the promoter of the gene makes it;
    step1: at least one fragment in contact with the promoter makes it) give the genes of `features` with that
    transition, their mean and median expression change, a 95% chromosome block-bootstrap interval of the mean
    (validation.block_resamples; groups of at least `min_genes` genes), the numbers of genes changing at least
    `rna_threshold` up and down, and concordance (fraction of those in the direction of the mean). Rows of kind 'combination' give the fragments
    per condition and, per level, the mean expression change of genes entering and leaving the combination, with
    intervals, and a label: 'bivalent-like' (holds marks of opposite `learned_signs` and a transition leaving it
    has an interval excluding 0, i.e. its resolution predicts the direction), 'active-like' (entering it goes with
    up or leaving it with down: interval above / below 0), 'repressive-like' (the opposite), 'neutral' otherwise
    (step 0 first, step 1 where step 0 is undecided).
    """
    marks = list(state_signals)
    nodes = list(graph.nodes)
    position = {n: i for i, n in enumerate(nodes)}
    present = [{m: np.array([state_signals[m][c].get(n, 0.0) for n in nodes]) >= threshold for m in marks}
               for c in (0, 1)]
    combo = [mark_combinations(p, marks) for p in present]
    counts = pd.Series(np.concatenate(combo)).value_counts()
    kept = set(counts.index[counts >= min_fragments])
    combo = [np.where(np.isin(c, list(kept)), c, 'other').astype(object) for c in combo]
    data = features.dropna(subset=['rna_log2FC'])
    data = data[data['node'].isin(position)].reset_index(drop=True)
    rna = data['rna_log2FC'].to_numpy(dtype=float)
    labels, resamples = block_resamples(_chromosomes(data['node']), n_boot, seed)
    chrom_index = np.searchsorted(labels, _chromosomes(data['node']))
    weights = np.array([np.bincount(r, minlength=len(labels)) for r in resamples], dtype=float) if len(labels) \
        else np.zeros((0, 0))

    def summary(index):
        index = np.unique(np.asarray(index, dtype=np.int64))
        values = rna[index]
        up = int((values >= rna_threshold).sum())
        down = int((values <= -rna_threshold).sum())
        out = dict(n_genes=len(index), mean_rna_log2FC=float(values.mean()) if len(index) else np.nan,
                   median_rna_log2FC=float(np.median(values)) if len(index) else np.nan, n_up=up, n_down=down,
                   ci_low=np.nan, ci_high=np.nan, concordance=np.nan)
        if len(index) >= min_genes:
            sums = np.bincount(chrom_index[index], weights=values, minlength=len(labels))
            sizes = np.bincount(chrom_index[index], minlength=len(labels)).astype(float)
            with np.errstate(invalid='ignore', divide='ignore'):
                boot = (weights @ sums) / (weights @ sizes)
            if np.isfinite(boot).any():
                out['ci_low'], out['ci_high'] = float(np.nanpercentile(boot, 2.5)), float(np.nanpercentile(boot, 97.5))
        if up + down:
            out['concordance'] = (up if out['mean_rna_log2FC'] > 0 else down) / (up + down)
        return out

    groups = {('step0', None): {}, ('step1', None): {}}
    for g, node in enumerate(data['node']):
        i = position[node]
        groups[('step0', None)].setdefault((combo[0][i], combo[1][i]), []).append(g)
        for pair in {(combo[0][position[f]], combo[1][position[f]]) for f in graph.neighbors(node)}:
            groups[('step1', None)].setdefault(pair, []).append(g)
    fragments = pd.Series(list(zip(combo[0], combo[1]))).value_counts()
    rows = []
    rated = {}
    for (level, _), members in groups.items():
        for (a, b), index in sorted(members.items()):
            s = summary(index)
            rated[(level, a, b)] = s
            rows.append(dict(kind='transition', level=level, combination_a=a, combination_b=b,
                             transition='{} -> {}'.format(a, b), changed=a != b,
                             n_fragments=int(fragments.get((a, b), 0)), **s))
    sign = {m: learned_signs.get(m, 0) for m in marks}

    def effect(entering, leaving):
        up = entering['ci_low'] > 0 or leaving['ci_high'] < 0
        down = entering['ci_high'] < 0 or leaving['ci_low'] > 0
        return 1 if up and not down else -1 if down and not up else 0

    for c in sorted(set(combo[0]) | set(combo[1])):
        row = dict(kind='combination', level='', combination_a=c, combination_b=c, transition=c,
                   n_fragments=int((combo[0] == c).sum() + (combo[1] == c).sum()),
                   fragments_A=int((combo[0] == c).sum()), fragments_B=int((combo[1] == c).sum()))
        members = c.split('+') if c not in ('none', 'other') else []
        opposite = any(sign[m] > 0 for m in members) and any(sign[m] < 0 for m in members)
        decided, resolution = 0, False
        for level in ('step0', 'step1'):
            entering = summary([g for (a, b), idx in groups[(level, None)].items() if b == c and a != c for g in idx])
            leaving = summary([g for (a, b), idx in groups[(level, None)].items() if a == c and b != c for g in idx])
            for name, s in (('entering', entering), ('leaving', leaving)):
                row['{}_{}_genes'.format(level, name)] = s['n_genes']
                row['{}_{}_mean'.format(level, name)] = s['mean_rna_log2FC']
                row['{}_{}_ci_low'.format(level, name)] = s['ci_low']
                row['{}_{}_ci_high'.format(level, name)] = s['ci_high']
            if not decided:
                decided = effect(entering, leaving)
            resolution |= any(rated[(level, a, b)]['ci_low'] > 0 or rated[(level, a, b)]['ci_high'] < 0
                              for (lv, a, b) in rated if lv == level and a == c and b != c)
        if c == 'other':
            row['label'] = 'pooled rare combinations'
        elif opposite and resolution:
            row['label'] = 'bivalent-like'
        else:
            row['label'] = {1: 'active-like', -1: 'repressive-like'}.get(decided, 'neutral')
        rows.append(row)
    return pd.DataFrame(rows)
