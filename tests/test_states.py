"""Mark roles (declared and learned), chromatin state calls and transitions on toy data."""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer.differential import cascade_table, classify_genes, explanation_cascade
from epigeneticgraphanalyzer.states import (NOT_CALLABLE, fragment_states, gene_states, learn_mark_roles,
                                            learned_level_signs, learned_roles, learned_state_transitions,
                                            mark_crossings,
                                            parse_mark_roles, resolve_mark_roles, role_signs, states_callable,
                                            transition_direction, transition_label)


def test_mark_sign_and_mark_roles_are_equivalent(tmp_path):
    marks = ['H3K27ac', 'H3K27me3', 'CTCF']
    by_sign = resolve_mark_roles(marks, ['H3K27me3=-1', 'H3K9me3=-1'])
    by_roles = resolve_mark_roles(marks, [], parse_mark_roles(['H3K27me3=repressive']))
    table = tmp_path / 'roles.tsv'
    pd.DataFrame({'mark': ['H3K27me3', 'H3K27ac'], 'role': ['repressive', 'active']}).to_csv(table, sep='\t',
                                                                                         index=False)
    by_file = resolve_mark_roles(marks, [], parse_mark_roles([str(table)]))
    assert role_signs(by_sign) == role_signs(by_roles) == role_signs(by_file) == {'H3K27ac': 1, 'H3K27me3': -1,
                                                                                   'CTCF': 1}
    # declared roles take precedence over --markSign; architectural marks have no expected direction
    roles = resolve_mark_roles(marks, ['H3K27me3=-1'], parse_mark_roles(['CTCF=architectural']))
    assert role_signs(roles) == {'H3K27ac': 1, 'H3K27me3': -1, 'CTCF': 0}
    with pytest.raises(ValueError):
        parse_mark_roles(['CTCF=insulator'])
    # the cascade is the same under both declarations
    from test_hops import features_table
    features = features_table()
    features['contact_direct'] = [1.2, -1.0, -1.5, 0.5, 0.0]
    signs_a = role_signs(resolve_mark_roles(['H3K27ac', 'H3K27me3'], ['H3K27me3=-1']))
    signs_b = role_signs(resolve_mark_roles(['H3K27ac', 'H3K27me3'], [], {'H3K27me3': 'repressive'}))
    a = explanation_cascade(classify_genes(features, ['H3K27ac', 'H3K27me3'], signs_a), ['H3K27ac', 'H3K27me3'],
                            signs_a)
    b = explanation_cascade(classify_genes(features, ['H3K27ac', 'H3K27me3'], signs_b), ['H3K27ac', 'H3K27me3'],
                            signs_b)
    pd.testing.assert_frame_equal(a, b)


def toy_features(n=400, seed=0):
    """Genes on eight chromosomes: the loss of mark R raises expression and the promoter's interactions; mark Z is
    noise."""
    rng = np.random.default_rng(seed)
    r = rng.normal(0, 1.5, n)
    return pd.DataFrame({
        'node': ['chr{}:{}-{}'.format(i % 8 + 1, 1000 * i, 1000 * i + 500) for i in range(n)],
        'rna_log2FC': -r + rng.normal(0, 0.7, n), 'contact_direct': -0.8 * r + rng.normal(0, 0.7, n),
        'n_direct': np.full(n, 2), 'R_promoter': r, 'R_direct': rng.normal(0, 1, n),
        'Z_promoter': rng.normal(0, 1.5, n), 'Z_direct': rng.normal(0, 1.5, n)})


def test_learned_roles():
    features = toy_features()
    table = learn_mark_roles(features, ['R', 'Z'], {'R': 'active', 'Z': 'unknown'}, n_boot=200, seed=1,
                             levels=['direct'])
    row = table.set_index(['mark', 'level'])
    # declared active, learned repressive at the promoter: the loss of R goes with expression up and interactions up
    assert row.loc[('R', 'promoter'), 'learned_role'] == 'repressive'
    assert row.loc[('R', 'promoter'), 'learned_sign'] == -1 and row.loc[('R', 'promoter'), 'rho_ci_high'] < 0
    assert row.loc[('R', 'promoter'), 'concordance'] < 0.5 and bool(row.loc[('R', 'promoter'), 'disagrees'])
    assert row.loc[('R', 'promoter_vs_contact'), 'rho_ci_high'] < 0
    assert row.loc[('R', 'promoter_vs_contact'), 'learned_role'] == 'falls with contact'
    # no level decides for another: R at step 1 is noise
    assert row.loc[('R', 'step1'), 'learned_role'] == 'no consistent effect'
    assert row.loc[('Z', 'promoter'), 'learned_role'] == 'no consistent effect'
    assert not bool(row.loc[('Z', 'promoter'), 'disagrees'])
    # the own interactions versus expression
    assert row.loc[('contact', 'promoter'), 'learned_role'] == 'active'
    assert learned_roles(table, {'R': 'active'}) == {'promoter': {'R': 'repressive', 'Z': 'neutral'},
                                                     'direct': {'R': 'neutral', 'Z': 'neutral'}}
    signs = learned_level_signs(table)
    assert signs.sign('R', 'promoter') == -1 and signs.sign('R', 'direct') == 0
    assert signs.sign('contact', 'promoter') == 1


def per_level_features(n=600, seed=2):
    """Mark M lost at the promoter and gained at step 1 of up-regulated genes (repressive at the promoter, activating
    at step 1); at step 2 it is noise."""
    rng = np.random.default_rng(seed)
    rna = rng.normal(0, 1.5, n)
    return pd.DataFrame({
        'node': ['chr{}:{}-{}'.format(i % 10 + 1, 1000 * i, 1000 * i + 500) for i in range(n)],
        'genes': ['g{}'.format(i) for i in range(n)], 'rna_log2FC': rna,
        'contact_direct': rng.normal(0, 0.3, n), 'n_direct': np.full(n, 3), 'n_indirect': np.full(n, 8),
        'M_promoter': -rna + rng.normal(0, 0.8, n), 'M_direct': rna + rng.normal(0, 0.8, n),
        'M_indirect': rng.normal(0, 1.5, n)})


def test_roles_learned_per_level_and_cascade():
    features = per_level_features()
    table = learn_mark_roles(features, ['M'], {'M': 'active'}, n_boot=200, seed=0)
    row = table.set_index(['mark', 'level'])
    assert row.loc[('M', 'promoter'), 'learned_role'] == 'repressive'
    assert row.loc[('M', 'step1'), 'learned_role'] == 'active'
    assert row.loc[('M', 'step2'), 'learned_role'] == 'no consistent effect'
    assert bool(row.loc[('M', 'promoter'), 'disagrees']) and not bool(row.loc[('M', 'step1'), 'disagrees'])
    assert list(row.loc[('M', 'promoter'), ['feature', 'hops']]) == ['promoter', 0]
    assert list(row.loc[('M', 'step2'), ['feature', 'hops']]) == ['indirect', 2]
    signs = learned_level_signs(table)
    # two up-regulated genes: M gained at the promoter only (explained by an active role, not by the learned
    # repressive one; gained at step 2, where M has no consistent effect), and M gained at step 1 only
    genes = pd.DataFrame({'node': ['chr1:1-2', 'chr1:3-4'], 'rna_log2FC': [2.0, 2.0], 'contact_direct': [0.0, 0.0],
                          'n_direct': [3, 3], 'n_indirect': [8, 8], 'M_promoter': [2.0, 0.0],
                          'M_direct': [0.0, 2.0], 'M_indirect': [2.0, 0.0]})
    declared = explanation_cascade(genes, ['M'], {'M': 1})
    learned = explanation_cascade(genes, ['M'], signs)
    assert list(declared['explained_at']) == ['gene', 'step1']
    assert list(learned['explained_at']) == ['unexplained', 'step1']
    # M lost at the promoter of an up-regulated gene explains it with the learned repressive role
    genes.loc[0, 'M_promoter'] = -2.0
    assert explanation_cascade(genes, ['M'], signs)['explained_at'].iloc[0] == 'gene'
    classes = classify_genes(genes, ['M'], signs)
    assert list(classes['explanation']) == ['promoter', 'direct']


def test_fragment_states_from_roles():
    roles = {'P': 'active', 'E': 'enhancer', 'R': 'repressive', 'C': 'architectural'}
    signal = {'P': np.array([3, 3, 0, 0, 0, 0]), 'E': np.array([0, 0, 3, 3, 0, 3]), 'R': np.array([3, 0, 3, 0, 3, 0]),
              'C': np.array([9, 9, 9, 9, 9, 9])}
    is_bait = np.array([True, True, True, False, False, True])
    states = fragment_states(signal, roles, 2.0, is_bait)
    # an enhancer-role mark is active-type only away from baits
    assert list(states) == ['bivalent', 'active', 'repressed', 'active enhancer', 'repressed', 'unmarked']
    other = fragment_states(signal, roles, 2.0, np.zeros(6, dtype=bool))
    assert other[2] == 'poised enhancer'
    assert states_callable(roles) == {'promoter': True, 'step1': True}
    assert states_callable({'P': 'active', 'C': 'architectural'}) == {'promoter': False, 'step1': False}
    assert states_callable({'E': 'enhancer', 'R': 'repressive'}) == {'promoter': False, 'step1': True}
    assert transition_direction('bivalent', 'active') == 1 and transition_direction('bivalent', 'repressed') == -1
    assert transition_direction('active', 'bivalent') == -1 and transition_direction('unmarked', 'bivalent') == 0
    assert transition_label('bivalent', 'active') == 'bivalent -> active: repressive mark lost'
    assert transition_label('poised enhancer', 'active enhancer') == \
        'poised enhancer -> active enhancer: repressive mark lost'
    # a signal moving just across the threshold makes no supported transition
    assert mark_crossings({'R': 2.1}, {'R': 1.9}, 2.0) == ([], ['R'], False)
    assert mark_crossings({'R': 3.2}, {'R': 0.5}, 2.0) == ([], ['R'], True)


def toy_bivalent():
    """Promoters p..s and u, other end e of gene t. P active mark, R repressive mark."""
    signal = {  # node: (P in A, R in A, P in B, R in B)
        'chr1:0-10': (4.1, 3.2, 4.5, 0.4),    # p: bivalent -> active, expression up
        'chr1:20-30': (4.0, 3.0, 3.9, 3.1),   # q: bivalent in both, expression down
        'chr1:40-50': (4.0, 3.0, 4.0, 0.5),   # r: bivalent -> active, expression down (discordant)
        'chr1:60-70': (4.0, 2.1, 4.0, 1.9),   # s: R just below the threshold in B (not supported)
        'chr1:80-90': (0.0, 0.0, 0.0, 0.0),   # t: unmarked promoter
        'chr1:100-110': (3.0, 3.0, 3.0, 0.2),  # e: poised enhancer -> active enhancer, contacts t
        'chr1:120-130': (4.0, 0.0, 4.0, 0.0),  # u: active in both
    }
    graph = nx.Graph()
    for node in signal:
        graph.add_node(node, bait=node != 'chr1:100-110')
    graph.add_edge('chr1:80-90', 'chr1:100-110', status='gained')
    graph.add_edge('chr1:0-10', 'chr1:120-130', status='shared')
    state_signals = {'P': ({n: v[0] for n, v in signal.items()}, {n: v[2] for n, v in signal.items()}),
                     'R': ({n: v[1] for n, v in signal.items()}, {n: v[3] for n, v in signal.items()})}
    nodes = ['chr1:0-10', 'chr1:20-30', 'chr1:40-50', 'chr1:60-70', 'chr1:80-90', 'chr1:120-130']
    classes = pd.DataFrame({'node': nodes, 'genes': ['P', 'Q', 'R', 'S', 'T', 'U'],
                            'rna_log2FC': [2.0, -2.0, -2.0, 2.0, 2.0, -2.0],
                            'n_direct': [1, 0, 0, 0, 1, 1], 'n_indirect': [0] * 6, 'contact_direct': [0.0] * 6,
                            'P_promoter': [0.0] * 6, 'P_direct': [0.0] * 6, 'P_indirect': [0.0] * 6,
                            'R_promoter': [0.0] * 6, 'R_direct': [0.0] * 6, 'R_indirect': [0.0] * 6})
    return graph, state_signals, classes


def test_bivalent_resolution_and_cascade():
    graph, state_signals, classes = toy_bivalent()
    roles = {'P': 'active', 'R': 'repressive'}
    states = gene_states(classes, graph, state_signals, roles, 2.0, names=('hESC', 'hNEC'))
    assert states.loc['chr1:0-10', 'promoter_transition'] == 'bivalent -> active: repressive mark lost'
    assert list(states['promoter_explains']) == [True, False, False, False, False, False]
    assert list(states['promoter_unresolved']) == [False, True, False, False, False, False]
    assert not states.loc['chr1:60-70', 'promoter_supported']
    assert states.loc['chr1:0-10', 'promoter_state_summary'] == \
        'promoter bivalent in hESC (P 4.1, R 3.2), resolved to active in hNEC by loss of R'
    assert states.loc['chr1:20-30', 'promoter_state_summary'].endswith('poised, trigger not measured')
    # step 1: the other end of t is a poised enhancer resolved to an active enhancer, concordant with up
    assert states.loc['chr1:80-90', 'step1_explains']
    assert states.loc['chr1:80-90', 'step1_transition'] == 'poised enhancer -> active enhancer: repressive mark lost'
    assert states.loc['chr1:80-90', 'step1_state_fragment'] == 'chr1:100-110'
    # a contacting promoter active in both conditions is no poised fragment
    assert states.loc['chr1:120-130', 'step1_transition'] == 'bivalent -> active: repressive mark lost'
    assert not states.loc['chr1:120-130', 'step1_explains']   # u goes down
    signs = role_signs(roles)
    cascade = explanation_cascade(classes, ['P', 'R'], signs, states=states)
    assert list(cascade['explained_at']) == ['gene', 'gene', 'unexplained', 'unexplained', 'step1', 'unexplained']
    assert list(cascade['explained_type'][:2]) == ['state', 'poised unresolved']
    assert cascade.loc['chr1:80-90', 'explained_type'] == 'state'
    table = cascade_table(cascade, max_hops=2, states_callable=states_callable(roles)).set_index('level')
    assert table.loc['gene_state', 'genes'] == 1 and table.loc['gene_poised_unresolved', 'genes'] == 1
    assert table.loc['step1_state', 'genes'] == 1 and table.loc['unexplained', 'genes'] == 3
    # without a repressive role no states: the cascade says so
    flat = {'P': 'active', 'R': 'neutral'}
    none = gene_states(classes, graph, state_signals, flat, 2.0)
    assert set(none['promoter_transition']) == {NOT_CALLABLE} and not none['promoter_explains'].any()
    table = cascade_table(explanation_cascade(classes, ['P', 'R'], role_signs(flat), states=none), max_hops=2,
                          states_callable=states_callable(flat)).set_index('level')
    assert table.loc['gene_state', 'note'] == NOT_CALLABLE and table.loc['gene_state', 'genes'] == 0


def test_learned_state_transitions():
    rng = np.random.default_rng(0)
    graph = nx.Graph()
    signal = {'P': ({}, {}), 'R': ({}, {})}
    rows = []
    for i in range(240):
        node = 'chr{}:{}-{}'.format(i % 6 + 1, 1000 * i, 1000 * i + 500)
        partner = 'chr{}:{}-{}'.format(i % 6 + 1, 1000 * i + 500, 1000 * i + 900)
        graph.add_edge(node, partner)
        resolved = i % 2 == 0
        signal['P'][0][node] = signal['P'][1][node] = 4.0
        signal['R'][0][node] = 4.0 if resolved else 0.0
        signal['R'][1][node] = 0.0
        for mark in 'PR':
            signal[mark][0][partner] = signal[mark][1][partner] = 0.0
        rows.append({'node': node, 'rna_log2FC': (2.0 if resolved else 0.0) + rng.normal(0, 0.3)})
    table = learned_state_transitions(pd.DataFrame(rows), graph, signal, {'P': 1, 'R': -1}, 2.0, min_fragments=50,
                                      n_boot=200)
    transitions = table[table['kind'] == 'transition'].set_index(['level', 'transition'])
    resolved = transitions.loc[('step0', 'P+R -> P')]
    assert resolved['n_genes'] == 120 and resolved['ci_low'] > 1.5 and resolved['concordance'] == 1.0
    combos = table[table['kind'] == 'combination'].set_index('combination_a')
    assert combos.loc['P+R', 'label'] == 'bivalent-like' and combos.loc['P', 'label'] == 'active-like'
    # without opposite learned roles the combination is only repressive-like (leaving it goes with up)
    table = learned_state_transitions(pd.DataFrame(rows), graph, signal, {'P': 1, 'R': 0}, 2.0, n_boot=200)
    combos = table[table['kind'] == 'combination'].set_index('combination_a')
    assert combos.loc['P+R', 'label'] == 'repressive-like'
