"""Mark presence from peak calls (--statePeaksA / --statePeaksB): peak files of several callers, overlap with the
fragments, states with peaks for all marks and with peaks for one mark and the threshold rule for the other."""

import gzip

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.elements import node_elements, tss_status
from epigeneticgraphanalyzer.states import (PeakPresence, fragment_states, gene_states, learned_state_transitions,
                                            mark_crossings, mark_present, peak_nodes, read_peaks)

# fragments of a toy union graph: name -> (chrom, start, end); four baited promoters and one other end
FRAGMENTS = {'chr1:1000-2000': ('chr1', 1000, 2000),   # p
             'chr1:3000-4000': ('chr1', 3000, 4000),   # q
             'chr1:5000-6000': ('chr1', 5000, 6000),   # r
             'chr1:7000-8000': ('chr1', 7000, 8000),   # t
             'chr1:9000-10000': ('chr1', 9000, 10000)}  # e, other end of t
P, Q, R, T, E = FRAGMENTS


def toy_graph():
    graph = nx.Graph()
    for node, (chrom, start, end) in FRAGMENTS.items():
        graph.add_node(node, chrom=chrom, start=start, end=end, bait=node != E)
    graph.add_edge(T, E, status='gained')
    return graph


def write_bed(path, rows, header=None, compress=False):
    text = (header + '\n' if header else '') + ''.join('\t'.join(str(v) for v in row) + '\n' for row in rows)
    if compress:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return str(path)


def presence(graph, path, min_overlap=1):
    nodes = list(graph.nodes)
    return peak_nodes(read_peaks([path]), nodes, [graph.nodes[n]['chrom'] for n in nodes],
                      [graph.nodes[n]['start'] for n in nodes], [graph.nodes[n]['end'] for n in nodes], min_overlap)


def test_peak_files_and_overlap(tmp_path):
    graph = toy_graph()
    # MACS2 narrowPeak (10 columns), with a track line
    narrow = write_bed(tmp_path / 'a.narrowPeak', [('chr1', 1900, 2100, 'p1', 50, '.', 5.0, 8.0, 6.0, 100),
                                                   ('chr1', 4000, 4500, 'p2', 50, '.', 5.0, 8.0, 6.0, 100)],
                       header='track type=narrowPeak')
    # the peak ending at 4000 does not reach q (half-open), the one at 1900-2100 overlaps p by 100 bp
    assert presence(graph, narrow) == {P}
    assert presence(graph, narrow, min_overlap=100) == {P}
    assert presence(graph, narrow, min_overlap=101) == set()
    # MACS broadPeak (9 columns), gzip; a broad domain spanning r and t, '1' prefix-less chromosome names
    broad = write_bed(tmp_path / 'a.broadPeak.gz', [('1', 4800, 8200, 'b1', 50, '.', 3.0, 4.0, 3.5)], compress=True)
    assert presence(graph, broad) == {R, T}
    # SEACR (6 columns) and a generic three-column BED with a comment line; several files are pooled
    seacr = write_bed(tmp_path / 'a.seacr.bed', [('chr1', 9500, 9600, 120.0, 30.0, 'chr1:9550-9551')])
    plain = write_bed(tmp_path / 'a.bed', [('chr1', 0, 1001)], header='# peaks')
    assert presence(graph, seacr) == {E}
    assert presence(graph, plain) == {P}
    nodes = list(graph.nodes)
    pooled = peak_nodes(read_peaks([seacr, plain]), nodes, [graph.nodes[n]['chrom'] for n in nodes],
                        [graph.nodes[n]['start'] for n in nodes], [graph.nodes[n]['end'] for n in nodes])
    assert pooled == {P, E}
    # a long peak far before a fragment does not hide a short one inside it
    mixed = write_bed(tmp_path / 'm.bed', [('chr1', 0, 900), ('chr1', 9100, 9200), ('chr1', 20000, 90000)])
    assert presence(graph, mixed) == {E}


def toy_states(peaks_for=('P', 'R')):
    """Promoters p, q, r, t and other end e. Mark P active, R repressive; signal (fold enrichment) and peaks
    disagree on purpose: the signal never reaches the threshold 2 where peaks call the mark present."""
    graph = toy_graph()
    # node: (P in A, R in A, P in B, R in B) signal
    signal = {P: (1.0, 1.5, 1.2, 0.2),    # peaks: bivalent -> active (R lost, change log2 2.5/1.2 > 0.5)
              Q: (1.0, 1.0, 1.0, 1.0),    # peaks: bivalent in both
              R: (0.0, 0.0, 0.0, 0.0),    # no peaks: unmarked
              T: (0.0, 0.0, 0.0, 0.0),    # unmarked promoter
              E: (1.5, 1.5, 1.5, 0.1)}    # peaks: poised enhancer -> active enhancer, contacts t
    peaks = {'P': ({P, Q, E}, {P, Q, E}), 'R': ({P, Q, E}, {Q})}
    state_signals = {}
    for k, mark in enumerate(('P', 'R')):
        values = ({n: v[k] for n, v in signal.items()}, {n: v[k + 2] for n, v in signal.items()})
        state_signals[mark] = PeakPresence(values, peaks[mark]) if mark in peaks_for else values
    classes = pd.DataFrame({'node': [P, Q, R, T], 'genes': ['GP', 'GQ', 'GR', 'GT'],
                            'rna_log2FC': [2.0, -2.0, 1.0, 2.0]})
    return graph, state_signals, classes


def test_presence_and_states_from_peaks():
    graph, state_signals, classes = toy_states()
    roles = {'P': 'active', 'R': 'repressive'}
    nodes = list(graph.nodes)
    # presence: peaks, not the signal (which stays below the threshold everywhere)
    assert list(mark_present(state_signals['R'], 0, nodes, 2.0)) == [True, True, False, False, True]
    assert list(mark_present(state_signals['R'], 1, nodes, 2.0)) == [False, True, False, False, False]
    states = gene_states(classes, graph, state_signals, roles, 2.0, names=('A', 'B'))
    assert list(states['promoter_state_A']) == ['bivalent', 'bivalent', 'unmarked', 'unmarked']
    assert list(states['promoter_state_B']) == ['active', 'bivalent', 'unmarked', 'unmarked']
    assert states.loc[P, 'promoter_transition'] == 'bivalent -> active: repressive mark lost'
    assert states.loc[P, 'promoter_supported'] and states.loc[P, 'promoter_explains']
    assert states.loc[Q, 'promoter_unresolved']
    # the state summary shows the signal of the marks
    assert states.loc[P, 'promoter_state_summary'] == \
        'promoter bivalent in A (P 1.0, R 1.5), resolved to active in B by loss of R'
    # step 1: the other end of t is a poised enhancer resolved to an active enhancer
    assert states.loc[T, 'step1_explains']
    assert states.loc[T, 'step1_transition'] == 'poised enhancer -> active enhancer: repressive mark lost'
    # the same states at every threshold: the peaks do not depend on it
    for threshold in (0.5, 1.5, 3.0):
        other = gene_states(classes, graph, state_signals, roles, threshold)
        pd.testing.assert_series_equal(other['promoter_state_A'], states['promoter_state_A'])
    # a peak lost without a signal change of at least --stateMinChange is no supported transition
    assert mark_crossings({'R': 1.0}, {'R': 1.0}, 2.0, present_a={'R': True}, present_b={'R': False}) == \
        ([], ['R'], False)
    assert mark_crossings({'R': 3.0}, {'R': 0.0}, 2.0, present_a={'R': True}, present_b={'R': False}) == \
        ([], ['R'], True)
    # element classes and learned combinations follow the peaks
    status = tss_status(graph, {n: (['G'], ['G']) for n in (P, Q, R, T)})
    elements = node_elements(graph, {'mean': state_signals}, 2.0, {'declared': roles}, status)
    a = elements[elements['condition'] == 'A'].set_index('node')
    b = elements[elements['condition'] == 'B'].set_index('node')
    assert a.loc[P, 'present_marks_mean'] == 'P+R' and b.loc[P, 'present_marks_mean'] == 'P'
    assert a.loc[P, 'class_declared_mean'] == 'bivalent promoter'
    assert b.loc[P, 'class_declared_mean'] == 'active promoter'
    assert a.loc[R, 'present_marks_mean'] == 'none'
    features = classes.assign(rna_log2FC=classes['rna_log2FC'])
    table = learned_state_transitions(features, graph, state_signals, {'P': 1, 'R': -1}, 2.0, min_fragments=1,
                                      n_boot=10)
    transitions = table[(table['kind'] == 'transition') & (table['level'] == 'step0')].set_index('transition')
    assert transitions.loc['P+R -> P', 'n_genes'] == 1 and transitions.loc['P+R -> P+R', 'n_genes'] == 1


def test_mixed_peaks_and_fold_enrichment():
    graph, state_signals, classes = toy_states(peaks_for=('P',))
    roles = {'P': 'active', 'R': 'repressive'}
    # R keeps the threshold rule: its signal (at most 1.5) is below 2, so no fragment carries R; P from peaks
    states = gene_states(classes, graph, state_signals, roles, 2.0)
    assert list(states['promoter_state_A']) == ['active', 'active', 'unmarked', 'unmarked']
    assert list(states['promoter_state_B']) == ['active', 'active', 'unmarked', 'unmarked']
    # at threshold 1 the signal of R calls it present at p, q (A; q also in B) and e: bivalent where P has peaks
    states = gene_states(classes, graph, state_signals, roles, 1.0)
    assert list(states['promoter_state_A']) == ['bivalent', 'bivalent', 'unmarked', 'unmarked']
    assert list(states['promoter_state_B']) == ['active', 'bivalent', 'unmarked', 'unmarked']
    nodes = list(graph.nodes)
    assert list(mark_present(state_signals['R'], 0, nodes, 1.0)) == [True, True, False, False, True]
    assert list(mark_present(state_signals['P'], 0, nodes, 1.0)) == [True, True, False, False, True]


def test_threshold_path_unchanged():
    """Plain tuples (no peaks) and peaks equal to signal >= threshold give the same states: the default path keeps
    the threshold rule."""
    graph, state_signals, classes = toy_states(peaks_for=())
    roles = {'P': 'active', 'R': 'repressive'}
    threshold = 1.0
    nodes = list(graph.nodes)
    as_peaks = {m: PeakPresence(v, tuple({n for n in nodes if v[c][n] >= threshold} for c in (0, 1)))
                for m, v in state_signals.items()}
    pd.testing.assert_frame_equal(gene_states(classes, graph, state_signals, roles, threshold),
                                  gene_states(classes, graph, as_peaks, roles, threshold))
    signal = {m: np.array([v[0][n] for n in nodes]) for m, v in state_signals.items()}
    is_bait = np.array([graph.nodes[n]['bait'] for n in nodes])
    present = {m: mark_present(as_peaks[m], 0, nodes, threshold) for m in as_peaks}
    assert list(fragment_states(signal, roles, threshold, is_bait)) == \
        list(fragment_states(signal, roles, threshold, is_bait, present))


def test_differential_cli_with_peaks(tmp_path):
    """egaDifferential with --statePeaksA / --statePeaksB for the repressive mark: the states follow the peaks, the
    mark changes of the explanation levels stay those of the signal tracks."""
    import os

    import pytest

    from epigeneticgraphanalyzer import egaDifferential
    from synthetic import fragments, simulate_capture, write_bigwig, write_tss_bed

    d = str(tmp_path)
    ibed_a, planted, baits = simulate_capture(d, seed=7, n_noise=0, name='a')
    ibed_b, _, _ = simulate_capture(d, seed=8, planted=planted, n_noise=0, name='b')
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    active = write_bigwig(os.path.join(d, 'ac.bw'), {b: 4.0 for b in baits})
    repressive = write_bigwig(os.path.join(d, 'me3.bw'), {})   # 1.0 everywhere: below the threshold 2
    # H3K27me3 peaks on every bait in A (the first one overlapping its fragment by 50 bp only), none in B
    rows = []
    for k, (chrom, i) in enumerate(baits):
        start, end = fragments(chrom)[i][1], fragments(chrom)[i][2]
        rows.append((chrom, start + 100, start + 600, 'peak{}'.format(k), 100, '.', 6.0, 10.0, 8.0, 250)
                    if k else (chrom, end - 50, end + 1000, 'peak0', 100, '.', 6.0, 10.0, 8.0, 250))
    peaks_a = write_bed(tmp_path / 'me3_A.narrowPeak', rows)
    peaks_b = write_bed(tmp_path / 'me3_B.bed', [])
    rng = np.random.default_rng(0)
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)), 'B': rng.uniform(5, 50, len(baits))})
    expression_path = os.path.join(d, 'expression.tsv')
    expression.to_csv(expression_path, sep='\t', index=False)
    common = ['--interactionsA', ibed_a, '--interactionsB', ibed_b, '-f', 'ibed', '--names', 'A', 'B', '--tss', tss,
              '--tracksA', 'H3K27ac=' + active, 'H3K27me3=' + repressive, '--tracksB', 'H3K27ac=' + active,
              'H3K27me3=' + repressive, '--markSign', 'H3K27me3=-1', '--expression', expression_path,
              '--expressionColumns', 'A', 'B', '--geneKey', 'genes', '--permutations', '2', '--bootstrap', '5',
              '--roleBootstrap', '5', '--genePermutations', '0', '--communityPermutations', '10', '--threads', '1']
    plain = os.path.join(d, 'out', 'plain')
    egaDifferential.main(common + ['-o', plain])
    peaks = os.path.join(d, 'out', 'peaks')
    egaDifferential.main(common + ['--statePeaksA', 'H3K27me3=' + peaks_a, '--statePeaksB',
                                   'H3K27me3=' + peaks_b, '-o', peaks])
    bait_nodes = {'{}:{}-{}'.format(*fragments(c)[i]) for c, i in baits}
    for prefix, bivalent_in_a in ((plain, False), (peaks, True)):
        elements = pd.read_csv(prefix + '_elements.tsv', sep='\t')
        a = elements[(elements['condition'] == 'A') & elements['node'].isin(bait_nodes)]
        b = elements[(elements['condition'] == 'B') & elements['node'].isin(bait_nodes)]
        assert len(a) == len(baits)
        assert set(b['class_declared_mean']) == {'active promoter'}
        assert set(a['class_declared_mean']) == {'bivalent promoter' if bivalent_in_a else 'active promoter'}
        assert set(a['class_declared_max1kb']) == set(a['class_declared_mean'])
    # with --statePeakMinOverlap 51 the peak overlapping the first bait by 50 bp no longer counts
    strict = os.path.join(d, 'out', 'strict')
    egaDifferential.main(common + ['--statePeaksA', 'H3K27me3=' + peaks_a, '--statePeaksB',
                                   'H3K27me3=' + peaks_b, '--statePeakMinOverlap', '51', '-o', strict])
    elements = pd.read_csv(strict + '_elements.tsv', sep='\t').set_index(['condition', 'node'])
    first = '{}:{}-{}'.format(*fragments(baits[0][0])[baits[0][1]])
    second = '{}:{}-{}'.format(*fragments(baits[1][0])[baits[1][1]])
    assert elements.loc[('A', first), 'class_declared_mean'] == 'active promoter'
    assert elements.loc[('A', second), 'class_declared_mean'] == 'bivalent promoter'
    marks = pd.read_csv(peaks + '_state_marks.tsv', sep='\t').set_index(['mark', 'condition', 'presence'])
    assert marks.loc[('H3K27me3', 'A', 'mean'), 'scale'] == 'peaks'
    assert marks.loc[('H3K27me3', 'B', 'mean'), 'fraction_ge_2'] == 0.0
    assert marks.loc[('H3K27me3', 'A', 'mean'), 'fraction_ge_3'] > 0
    assert str(marks.loc[('H3K27ac', 'A', 'mean'), 'scale']) == '1'
    # bivalent -> active transitions at the promoters, unsupported: the H3K27me3 signal does not change
    classes = pd.read_csv(peaks + '_gene_classes.tsv', sep='\t')
    assert set(classes['promoter_transition']) == {'bivalent -> active: repressive mark lost'}
    assert not classes['promoter_supported'].any()
    # the mark changes of the levels come from the signal tracks, whatever the peaks
    pd.testing.assert_frame_equal(pd.read_csv(plain + '_gene_features.tsv', sep='\t'),
                                  pd.read_csv(peaks + '_gene_features.tsv', sep='\t'))
    # peaks of a mark in one condition only, or of a mark without tracks, are refused
    with pytest.raises(ValueError):
        egaDifferential.main(common + ['--statePeaksA', 'H3K27me3=' + peaks_a, '-o', peaks])
    with pytest.raises(ValueError):
        egaDifferential.main(common + ['--statePeaksA', 'H3K4me3=' + peaks_a, '--statePeaksB',
                                       'H3K4me3=' + peaks_b, '-o', peaks])
