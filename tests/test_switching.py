import os

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer import switching


def toy_graph(path):
    """chr1 with three baited promoters PA (down), PB (up), PC (up) and four other ends:
    E loses PA and gains PB (planted switch); F keeps PA and PC (no switch); G gains PB and keeps PC (gain only);
    H keeps PA and PB (no switch, opposite expression)."""
    g = nx.Graph()

    def node(name, start, bait, genes='', rna=None, k27ac=(1.0, 1.0)):
        data = dict(chrom='chr1', start=start, end=start + 1000, bait=int(bait), genes=genes, gene_ids=genes,
                    H3K27ac_X=k27ac[0], H3K27ac_Y=k27ac[1])
        if rna is not None:
            data['rna_log2FC'] = rna
        g.add_node('chr1:{}-{}'.format(start, start + 1000), **data)
        return 'chr1:{}-{}'.format(start, start + 1000)

    pa = node('PA', 100000, True, 'GA', -2.0)
    pb = node('PB', 300000, True, 'GB', 2.0)
    pc = node('PC', 500000, True, 'GC', 1.0)
    e = node('E', 200000, False, k27ac=(1.0, 8.0))
    f = node('F', 400000, False)
    gg = node('G', 450000, False)
    h = node('H', 250000, False)

    def edge(u, v, status, ra, rb):
        g.add_edge(u, v, r_a=ra, r_b=rb, log2FC=float(np.log2((rb * 5 + 1) / (ra * 5 + 1))), qvalue=np.nan,
                   status=status)

    edge(e, pa, 'lost', 3.0, 0.0)
    edge(e, pb, 'gained', 0.0, 3.0)
    edge(f, pa, 'shared', 2.0, 2.1)
    edge(f, pc, 'shared', 2.0, 2.2)
    edge(gg, pb, 'gained', 0.0, 2.0)
    edge(gg, pc, 'shared', 2.0, 2.0)
    edge(h, pa, 'shared', 2.0, 2.0)
    edge(h, pb, 'shared', 2.0, 2.0)
    nx.write_graphml(g, path)
    return dict(PA=pa, PB=pb, PC=pc, E=e, F=f, G=gg, H=h)


def load(tmp_path):
    prefix = os.path.join(str(tmp_path), 'X_vs_Y')
    ids = toy_graph(prefix + '_union.graphml')
    nodes, edges = switching.read_union_graphml(prefix + '_union.graphml')
    return prefix, ids, nodes, edges


def test_planted_switch_found_and_non_switch_not(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    names = switching.condition_names(prefix, nodes)
    assert names == ('X', 'Y')
    events = switching.switch_events(nodes, edges, names)
    assert len(events) == 1
    event = events.iloc[0]
    assert event['element'] == ids['E']
    assert event['promoter_lost'] == ids['PA'] and event['promoter_gained'] == ids['PB']
    assert event['element_type'] == 'non-promoter element'
    assert bool(event['concordant']) and bool(event['concordant_strict'])
    assert bool(event['element_chromatin_changed'])
    assert event['distance_lost'] == 100000 and event['distance_gained'] == 100000
    # F (keeps both), G (gain without loss) and H (keeps both) are no switches
    assert not set(events['element']) & {ids['F'], ids['G'], ids['H']}


def test_weakened_contacts_count_as_switch(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    # H: contact with PA weakened (shared, log2 -1.5), contact with PB strengthened (shared, log2 +1.2)
    for k, row in edges.iterrows():
        pair = {row['u'], row['v']}
        if pair == {ids['H'], ids['PA']}:
            edges.loc[k, 'log2FC'] = -1.5
        if pair == {ids['H'], ids['PB']}:
            edges.loc[k, 'log2FC'] = 1.2
    events = switching.switch_events(nodes, edges, ('X', 'Y'))
    assert set(events['element']) == {ids['E'], ids['H']}
    # below the threshold no switch
    events = switching.switch_events(nodes, edges, ('X', 'Y'), contact_threshold=2.0)
    assert set(events['element']) == {ids['E']}


def test_same_gene_promoters_are_not_a_switch(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    nodes.loc[ids['PB'], 'genes'] = 'GA'
    assert switching.switch_events(nodes, edges, ('X', 'Y')).empty


def test_hub_switch(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    hubs, left, joined, p_nodes = switching.hub_switches(nodes, edges, ('X', 'Y'))
    assert list(hubs['element']) == [ids['E']]
    assert hubs.iloc[0]['genes_left'] == 'GA' and hubs.iloc[0]['genes_joined'] == 'GB'
    assert bool(hubs.iloc[0]['concordant'])


def test_competition(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    table, _ = switching.competition_pairs(nodes, edges)
    tested = table[table['tested']].set_index('promoter_2')
    # PB gains E shared with PA (PB up, PA down: competitive) and G shared with PC (both up: not competitive)
    assert len(tested) == 2
    row = tested.loc[ids['PA']]
    assert row['promoter_1'] == ids['PB'] and row['n_shared'] == 2
    assert row['delta'] == 2 and bool(row['competitive'])
    row = tested.loc[ids['PC']]
    assert row['promoter_1'] == ids['PB'] and row['delta'] == 1 and not bool(row['competitive'])
    # PA / PC share F with kept contacts only: not tested
    pair = table[(table['promoter_1'].isin([ids['PA'], ids['PC']])) & (table['promoter_2'].isin([ids['PA'], ids['PC']]))]
    assert len(pair) == 1 and not bool(pair.iloc[0]['tested'])


def test_null_draws_keep_values():
    rng = np.random.default_rng(1)
    values = np.arange(10, dtype=float)
    strata = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2])
    draws = switching.stratified_draws(values, strata, 50, rng)
    for s in np.unique(strata):
        idx = strata == s
        assert np.all(np.sort(draws[:, idx], axis=1) == np.sort(values[idx]))
    chroms = np.array(['a'] * 6 + ['b'] * 4)
    draws = switching.rotation_draws(values, chroms, np.arange(10), 50, rng)
    for c in ('a', 'b'):
        idx = chroms == c
        assert np.all(np.sort(draws[:, idx], axis=1) == np.sort(values[idx]))
    # a rotation keeps neighbours next to each other
    a = draws[:, :6]
    assert np.all(((a[:, 1:] - a[:, :-1]) % 6 == 1) | ((a[:, 1:] - a[:, :-1]) == -5))
    assert switching.permutation_p(5.0, np.arange(10)) == (1 + 5) / 11


def test_role_swap_test():
    # 30 elements with 3 events each: all concordant -> significant; balanced -> not
    clusters = np.repeat(np.arange(30), 3)
    concordant = np.ones(90, dtype=bool)
    discordant = np.zeros(90, dtype=bool)
    mean, p, null = switching.role_swap_test(concordant, discordant, clusters, n_perm=2000, seed=0)
    assert abs(mean - 0.5) < 0.02 and p < 0.001 and len(null) == 2000
    concordant = np.tile([True, False], 45)
    mean, p, _ = switching.role_swap_test(concordant, ~concordant, clusters, n_perm=2000, seed=0)
    assert p > 0.2


def test_null_test_detects_planted_concordance():
    rng = np.random.default_rng(0)
    n = 400
    index = ['p{}'.format(i) for i in range(n)]
    rna = rng.normal(size=n)
    promoters = pd.DataFrame({'chrom': 'chr1', 'start': np.arange(n) * 10000, 'end': np.arange(n) * 10000 + 1000,
                              'rna_log2FC': rna}, index=index)
    down = [index[i] for i in np.flatnonzero(rna < -0.5)]
    up = [index[i] for i in np.flatnonzero(rna > 0.5)]
    m = min(len(down), len(up))
    first, second = np.array(down[:m]), np.array(up[:m])
    pool = pd.Index(np.r_[first, second])
    levels = pd.Series(rng.normal(size=n), index=index)
    distances = pd.Series(50000.0, index=index)
    observed, perm, rot = switching.null_test(first, second, pool, levels, distances, promoters, n_perm=200, seed=0)
    assert observed['concordant'] == 1.0
    assert switching.permutation_p(observed['concordant'], perm['concordant']) < 0.01
    assert switching.permutation_p(observed['concordant'], rot['concordant']) < 0.01


def test_cli(tmp_path):
    prefix, ids, nodes, edges = load(tmp_path)
    out = os.path.join(str(tmp_path), 'out', 'X_vs_Y')
    switching.main(['-p', prefix, '-o', out, '--permutations', '20', '--geneSets', 'up=GB:up'])
    for suffix in ('_switch_events.tsv', '_switch_summary.tsv', '_competition.tsv', '_switch_hubs.tsv',
                   '_switch_biology.tsv'):
        assert os.path.exists(out + suffix)
    summary = pd.read_csv(out + '_switch_summary.tsv', sep='\t')
    assert summary.loc[summary['set'] == 'non-promoter element', 'n_events'].iloc[0] == 1
    figures = os.listdir(os.path.join(str(tmp_path), 'out', 'figures'))
    assert 'X_vs_Y_switch_top1_X.pdf' in figures and 'X_vs_Y_switch_top1_Y.pdf' in figures
    biology = pd.read_csv(out + '_switch_biology.tsv', sep='\t')
    assert biology.iloc[0]['events_as_gained'] == 1
