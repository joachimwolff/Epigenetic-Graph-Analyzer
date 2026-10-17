"""Motif scanning, enrichment and candidate trigger factors (triggers, egaTriggers) on toy sequences."""

import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest

pytest.importorskip('MOODS')

from epigeneticgraphanalyzer import egaTriggers  # noqa: E402
from epigeneticgraphanalyzer.triggers import (best_window, explaining_nodes, gene_triggers,  # noqa: E402
                                              matched_background, motif_enrichment, motif_genes, motif_thresholds,
                                              read_jaspar, scan_sequences, strata, tf_family)

PLANT = 'TGACGCCATG'
JASPAR = """>MA9999.1\tToy1::Pax3
A  [ 0 0 20 0 0 0 0 20 0 0 ]
C  [ 0 0 0 20 0 20 20 0 0 0 ]
G  [ 0 20 0 0 20 0 0 0 0 20 ]
T  [ 20 0 0 0 0 0 0 0 20 0 ]
>MA9998.1\tCTRL(var.2)
A  [ 20 20 0 20 0 20 0 0 ]
C  [ 0 0 0 0 20 0 0 20 ]
G  [ 0 0 20 0 0 0 20 0 ]
T  [ 0 0 0 0 0 0 0 0 ]
"""


def random_seq(rng, n):
    return ''.join(rng.choice(list('ACGT'), n))


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


@pytest.fixture
def motifs(tmp_path):
    path = tmp_path / 'motifs.jaspar'
    path.write_text(JASPAR)
    return read_jaspar(str(path))


def test_read_jaspar_and_names(motifs):
    assert [m['id'] for m in motifs] == ['MA9999.1', 'MA9998.1']
    assert motifs[0]['counts'].shape == (4, 10)
    assert motif_genes('Pax3::Foxo1') == ['PAX3', 'FOXO1']
    assert motif_genes('NR2F1(var.2)') == ['NR2F1']
    assert tf_family('HOXA2') == tf_family('Hoxb3') == 'HOX'
    assert tf_family('SOX10') == 'SOX' and tf_family('CEBPA') == 'CEBP' and tf_family('RUNX1') == 'RUNX'
    assert tf_family('ZNF263') is None


def test_planted_motif_found_on_both_strands(motifs):
    rng = np.random.default_rng(1)
    background = [0.25] * 4
    plain = [random_seq(rng, 200) for _ in range(30)]
    forward = [s[:90] + PLANT + s[100:] for s in plain]
    reverse = [s[:50] + revcomp(PLANT) + s[60:] for s in plain]
    counts, best = scan_sequences(plain + forward + reverse, motifs, background, pvalue=1e-4, threads=2, chunk=7)
    assert counts.shape == (90, 2)
    assert (counts[30:, 0] > 0).all()
    assert (counts[:30, 0] > 0).mean() < 0.1
    _, maxima = motif_thresholds(motifs, background)
    assert np.allclose(best[30:, 0], maxima[0], atol=1e-3)


def test_enrichment_detects_planted_motif(motifs):
    rng = np.random.default_rng(2)
    fg = [(lambda s: s[:200] + PLANT + s[210:])(random_seq(rng, 500)) for _ in range(40)]
    bg = [random_seq(rng, 500) for _ in range(200)]
    counts, _ = scan_sequences(fg + bg, motifs, [0.25] * 4, threads=1)
    table = motif_enrichment(counts > 0, np.arange(40), np.arange(40, 240))
    assert table.loc[0, 'q'] < 1e-10 and table.loc[0, 'log2_odds'] > 3
    assert table.loc[1, 'q'] > 0.05


def test_matched_background_follows_gc_strata():
    rng = np.random.default_rng(3)
    gc = np.concatenate([np.full(20, 0.62), rng.uniform(0.3, 0.7, 2000)])
    stratum = strata(gc, np.full(len(gc), 1000))
    bg = matched_background(np.arange(20), np.arange(20, len(gc)), stratum, ratio=3, rng=rng)
    assert len(bg) == 60
    assert set(stratum[bg]) == {stratum[0]}


def test_best_window_takes_the_peak():
    profile = np.zeros(80)
    profile[50:60] = 5.0
    offset = best_window(profile, 50, 500)
    assert offset == 50 * 50
    assert best_window(np.zeros(80), 50, 500) == 35 * 50


def toy_differential():
    """Gene G (up) at chr1:0-1000 with partner m1 (unchanged) and m2 two steps away gaining the active mark P;
    the promoter of the factor PAX3 is m1."""
    graph = nx.Graph()
    g, m1, m2, m3 = 'chr1:0-1000', 'chr1:2000-3000', 'chr1:4000-5000', 'chr1:6000-7000'
    for u, v in ((g, m1), (m1, m2), (m1, m3)):
        graph.add_edge(u, v, r_a=1.0, r_b=1.0, status='shared', log2FC=0.0, qvalue=np.nan)
    values = {g: (4, 4), m1: (2, 2), m2: (1, 15), m3: (2, 2.5)}
    for n in graph.nodes:
        chrom, rest = n.split(':')
        start, end = rest.split('-')
        graph.nodes[n].update(chrom=chrom, start=int(start), end=int(end), bait=n == g,
                              genes={g: 'G', m1: 'PAX3'}.get(n, ''), P_A=values[n][0], P_B=values[n][1])
    rows = []
    for rank, level in enumerate(['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks',
                                  'step1_state', 'step2'], start=1):
        explains = level == 'step2'
        rows.append(dict(node=g, genes='G', rna_log2FC=3.0, rank=rank, level=level, signs_declared='P:+1',
                         concordant_marks_declared='P' if explains else np.nan, explains_declared_mean=explains,
                         result_declared_mean='concordant' if explains else 'below threshold',
                         state_fragment_declared_mean=np.nan))
    levels = pd.DataFrame(rows)
    classes = pd.DataFrame(dict(node=[g], genes=['G'], rna_log2FC=[3.0], explained_at=['step2'],
                                explained_type=['marks'], promoter_transition=['active -> active']))
    return graph, levels, classes, (g, m1, m2, m3)


def test_explaining_nodes_toy():
    graph, levels, _, (g, m1, m2, m3) = toy_differential()
    nodes = explaining_nodes(graph, levels, ('A', 'B'), max_hops=2)
    flags = {(r.level, r.node): r.explaining for r in nodes.itertuples()}
    assert flags[('promoter', g)] is False
    assert flags[('step1', m1)] is False
    assert flags[('step2', m2)] is True and flags[('step2', m3)] is False


def test_gene_triggers_toy(motifs):
    graph, levels, classes, (g, m1, m2, m3) = toy_differential()
    nodes = explaining_nodes(graph, levels, ('A', 'B'), max_hops=2)
    index = {g: 0, m1: 1, m2: 2, m3: 3}
    counts = np.array([[0, 0], [0, 0], [2, 1], [0, 0]], dtype=np.uint16)
    best = np.where(counts > 0, 10.0, np.nan).astype(np.float32)
    expression = {'PAX3': (2.0, 1e-5, 100.0), 'TOY1': (-1.5, 0.01, 50.0), 'CTRL': (0.1, 0.9, 10.0),
                  'PAX7': (1.0, 0.001, 5.0)}
    table = gene_triggers(classes, nodes, index, counts, best, motifs, np.array([12.0, 10.0]), expression, graph)
    assert set(table['tf']) == {'PAX3', 'TOY1', 'PAX7'}
    pax3 = table.set_index('tf').loc['PAX3']
    assert pax3['relation'] == 'activator-consistent' and pax3['tf_hops'] == 1 and pax3['evidence'] == 'named factor'
    assert table.set_index('tf').loc['TOY1', 'relation'] == 'repressor-consistent'
    assert table.set_index('tf').loc['PAX7', 'evidence'] == 'family of PAX3'
    assert table['rank'].tolist() == list(range(1, len(table) + 1))
    assert table.iloc[-1]['tf'] == 'PAX7'
    # motif 0 hits 1 of 4 windows overall and 1 of the gene's 2 scanned windows (promoter, m2): P = 1 - 0.75^2
    assert np.isclose(pax3['motif_background_fraction'], 0.25)
    assert np.isclose(pax3['gene_motif_p'], 1 - 0.75 ** 2)


def test_egatriggers_end_to_end(tmp_path, motifs):
    graph, levels, classes, (g, m1, m2, m3) = toy_differential()
    rng = np.random.default_rng(4)
    genome = random_seq(rng, 8000)
    genome = genome[:4400] + PLANT + genome[4410:]
    d = str(tmp_path)
    (tmp_path / 'genome.fa').write_text('>chr1\n' + '\n'.join(genome[i:i + 60] for i in range(0, 8000, 60)) + '\n')
    (tmp_path / 'motifs.jaspar').write_text(JASPAR)
    prefix = os.path.join(d, 'A_vs_B')
    for n in graph.nodes:
        graph.nodes[n]['bait'] = int(graph.nodes[n]['bait'])
    nx.write_graphml(graph, prefix + '_union.graphml')
    levels.to_csv(prefix + '_explanations_all_levels.tsv', sep='\t', index=False)
    classes.to_csv(prefix + '_gene_classes.tsv', sep='\t', index=False)
    pd.DataFrame(dict(gene_name=['PAX3', 'TOY1', 'G'], log2FoldChange=[2.0, -1.0, 3.0], padj=[1e-4, 0.5, 1e-6],
                      baseMean=[10, 10, 10])).to_csv(os.path.join(d, 'expr.tsv'), sep='\t', index=False)
    out = os.path.join(d, 'out', 'A_vs_B')
    egaTriggers.main(['-p', prefix, '--genome', os.path.join(d, 'genome.fa'), '--motifs',
                      os.path.join(d, 'motifs.jaspar'), '--expression', os.path.join(d, 'expr.tsv'),
                      '--maxHops', '2', '--minGenes', '1', '--minWindows', '1', '--threads', '1', '-o', out])
    triggers = pd.read_csv(out + '_gene_triggers.tsv', sep='\t')
    assert triggers['tf'].tolist() == ['PAX3']
    assert triggers.loc[0, 'steps_with_motif'] == 'step2' and triggers.loc[0, 'n_nodes_with_motif'] == 1
    windows = pd.read_csv(out + '_windows.tsv.gz', sep='\t')
    assert (windows['win_end'] - windows['win_start']).max() == 1000
    assert os.path.exists(out + '_motif_enrichment.tsv')
