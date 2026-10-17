"""Element classes of fragments on toy data: TSS-based promoter calls of other ends, class labels from the roles of the
present marks (no mark name used), roles learned per element type, transitions, external annotation."""

import networkx as nx
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.differential import ActivityOperators, explanation_levels, gene_features, \
    classify_genes
from epigeneticgraphanalyzer.elements import (assign_elements, element_agreement, element_class_counts,
                                              element_classes, element_transitions, node_elements, node_types,
                                              roles_per_element_type, stratum_masks, tss_status)
from epigeneticgraphanalyzer.graph import annotate_genes


def add_node(graph, node, bait, genes=''):
    chrom, rest = node.split(':')
    start, end = rest.split('-')
    graph.add_node(node, chrom=chrom, start=int(start), end=int(end), bait=bait, genes=genes)


def test_promoter_calls_of_other_ends():
    graph = nx.Graph()
    add_node(graph, 'chr1:0-1000', True, 'A')          # baited promoter
    add_node(graph, 'chr1:5000-6000', False, 'B')      # other end with a TSS: non-baited promoter
    add_node(graph, 'chr1:9000-10000', True, '')       # bait without a TSS
    add_node(graph, 'chr1:20000-21000', False, '')     # enhancer candidate
    status = tss_status(graph)
    assert status['chr1:0-1000'][:2] == ('baited promoter', 'promoter')
    assert status['chr1:5000-6000'] == ('non-baited promoter', 'promoter', 'B')
    assert status['chr1:9000-10000'][:2] == ('baited non-promoter', 'non-promoter element')
    assert status['chr1:20000-21000'][:2] == ('non-promoter element', 'non-promoter element')
    # an other annotation (--annotation): a TSS at 20500 makes the last fragment a non-baited promoter
    nodes = pd.DataFrame([(n, d['chrom'], d['start'], d['end']) for n, d in graph.nodes(data=True)],
                         columns=['node', 'chrom', 'start', 'end'])
    tss = pd.DataFrame({'chrom': ['chr1', 'chr1'], 'tss': [500, 20500], 'strand': ['+', '-'],
                        'gene_id': ['A', 'C'], 'gene_name': ['A', 'C'], 'gene_type': ['', '']})
    other = tss_status(graph, annotate_genes(nodes, tss))
    assert other['chr1:20000-21000'] == ('non-baited promoter', 'promoter', 'C')
    assert other['chr1:5000-6000'][0] == 'non-promoter element'


def present(*rows):
    """mark -> bool array from rows of present mark sets."""
    marks = sorted({m for r in rows for m in r})
    return {m: np.array([m in r for r in rows]) for m in marks}


def test_class_labels_from_roles():
    # rows: fragments; the same rules for any mark names
    rows = [{'X', 'Z'}, {'X'}, {'Y', 'Z'}, {'Y'}, {'Z'}, {'W'}, {'C'}, {'U'}, set()]
    for names in ({'X': 'X', 'Y': 'Y', 'Z': 'Z'}, {'X': 'H3K4me3', 'Y': 'H3K4me1', 'Z': 'H3K27me3'}):
        rename = dict(names, W='W', C='C', U='U')
        p = present(*[{rename[m] for m in r} for r in rows])
        for m in rename.values():
            p.setdefault(m, np.zeros(len(rows), dtype=bool))
        roles = {rename['X']: 'active', rename['Y']: 'enhancer', rename['Z']: 'repressive', 'W': 'repressive',
                 'C': 'architectural', 'U': 'unknown'}
        promoter = element_classes(p, roles, np.ones(len(rows), dtype=bool), primary=[rename['Z']])
        element = element_classes(p, roles, np.zeros(len(rows), dtype=bool), primary=[rename['Z']])
        assert list(promoter) == ['bivalent promoter', 'active promoter', 'poised promoter', 'primed promoter',
                                  'repressed promoter', 'heterochromatin-like promoter', 'architectural promoter',
                                  'promoter with marks of no role', 'unmarked promoter']
        assert list(element) == ['poised enhancer', 'active enhancer', 'poised enhancer', 'primed enhancer',
                                 'repressed element', 'heterochromatin-like element', 'architectural element',
                                 'element with marks of no role', 'unmarked element']
    # without a primary mark every repressive mark is primary; without an active-role mark no 'primed' call
    p = present({'W'}, {'Y'})
    roles = {'W': 'repressive', 'Y': 'enhancer'}
    assert list(element_classes(p, roles, np.zeros(2, dtype=bool))) == ['repressed element',
                                                                        'enhancer-marked element']


def toy_run(n_genes=240, seed=0):
    """Genes on 12 chromosomes, each with a non-baited promoter partner and an element partner. The mark M rises with
    expression at the partner promoters and falls with it at the elements."""
    rng = np.random.default_rng(seed)
    graph = nx.Graph()
    expression, m_a, m_b, presence_a, presence_b = {}, {}, {}, {}, {}
    for i in range(n_genes):
        chrom = 'chr{}'.format(i % 12 + 1)
        base = (i // 12) * 100000
        g, p, e = ['{}:{}-{}'.format(chrom, base + o, base + o + 1000) for o in (0, 20000, 40000)]
        add_node(graph, g, True, 'G{}'.format(i))
        add_node(graph, p, False, 'P{}'.format(i))
        add_node(graph, e, False, '')
        for u, v in ((g, p), (g, e)):
            graph.add_edge(u, v, r_a=1.0, r_b=1.0, status='shared', log2FC=0.0, qvalue=np.nan)
        rna = rng.normal(0, 2)
        expression[g] = rna
        m_a.update({g: 4.0, p: 4.0, e: 4.0})
        m_b.update({g: 4.0, p: 4.0 * 2 ** (rna + rng.normal(0, 0.3)), e: 4.0 * 2 ** (-rna + rng.normal(0, 0.3))})
        presence_a.update({g: 3.0, p: 3.0, e: 3.0 if i % 2 else 0.0})
        presence_b.update({g: 3.0, p: 3.0, e: 0.0 if i % 2 else 3.0})
    return graph, {'M': (m_a, m_b)}, expression, {'mean': {'M': (presence_a, presence_b)}}


def test_roles_learned_per_element_type_and_transitions():
    graph, marks, expression, presence = toy_run()
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    roles = {'M': 'active'}
    elements = node_elements(graph, presence, 2.0, {'declared': roles}, tss_status(graph), names=('A', 'B'))
    assert set(elements['class_declared_mean']) == {'active promoter', 'active enhancer', 'unmarked element'}
    counts = element_class_counts(elements)
    non_baited = counts[(counts['tss_status'] == 'non-baited promoter') & (counts['condition'] == 'A')]
    assert int(non_baited['fragments'].sum()) == 240
    masks = stratum_masks(elements, operators)
    assert {'promoter', 'non-promoter element', 'non-baited promoter', 'class: active enhancer'} <= set(masks)
    table = roles_per_element_type(features, operators, vectors, ['M'], roles, masks, n_boot=50, min_genes=50)
    step1 = table[(table['level'] == 'step1')].set_index('element_type')
    assert step1.loc['promoter', 'learned_role'] == 'active'
    assert step1.loc['non-baited promoter', 'learned_role'] == 'active'
    assert step1.loc['non-promoter element', 'learned_role'] == 'repressive'
    # the class strata hold the genes with a node of that class (in either condition) at that level
    assert step1.loc['class: active enhancer', 'n_genes'] == 240
    # the promoter level of the element stratum has no gene (own fragments are promoters)
    assert not ((table['element_type'] == 'non-promoter element') & (table['level'] == 'promoter')).any()
    # class transitions: element switching off (active -> unmarked) and on, at step 1
    transitions = element_transitions(elements, operators, features, n_boot=50)
    step1 = transitions[(transitions['step'] == 1)].set_index('transition')
    assert step1.loc['active enhancer -> unmarked element', 'n_genes'] == 120
    assert step1.loc['unmarked element -> active enhancer', 'n_fragments'] == 120
    assert step1.loc['active promoter -> active promoter', 'n_genes'] == 240


def test_all_levels_split_promoters_and_elements():
    graph, marks, expression, presence = toy_run(n_genes=24)
    operators = ActivityOperators(graph, list(expression), max_hops=2)
    vectors = {m: (operators.vector(a), operators.vector(b)) for m, (a, b) in marks.items()}
    features = gene_features(graph, marks, expression, operators=operators, mark_vectors=vectors)
    classes = classify_genes(features[features['rna_log2FC'].abs() > 1], ['M'], {'M': 1})
    elements = node_elements(graph, presence, 2.0, {'declared': {'M': 'active'}}, tss_status(graph))
    levels = explanation_levels(classes, graph, operators, vectors, ['M'], {'declared': {'M': 1}},
                                node_types=node_types(elements))
    step1 = levels[levels['level'] == 'step1_marks']
    explained = step1[step1['explains_declared_mean']]
    assert len(explained)
    # M follows expression at the partner promoters only: every concordant node there is a promoter
    assert (explained['n_concordant_promoters_declared'] >= 1).all()
    assert (explained['n_concordant_elements_declared'] == 0).all()
    assert (explained['strongest_type_declared'] == 'promoter').all()
    assert explained['strongest_promoter_declared'].str.contains('P').all()


def test_external_annotation_agreement():
    graph = nx.Graph()
    add_node(graph, 'chr1:0-1000', True, 'A')
    add_node(graph, 'chr1:5000-6000', False, '')
    add_node(graph, 'chr1:9000-10000', False, '')
    annotation = pd.DataFrame({'chrom': ['chr1', 'chr1', 'chr1'], 'start': [100, 5100, 5800],
                               'end': [400, 5300, 9100], 'name': ['PLS', 'dELS', 'CTCF-only']})
    assigned = assign_elements(graph, annotation)
    # largest overlap wins: 5100-5300 (200 bp) against 5800-6000 (200 bp, tie: first) on the second fragment
    assert assigned == {'chr1:0-1000': 'PLS', 'chr1:5000-6000': 'dELS', 'chr1:9000-10000': 'CTCF-only'}
    presence = {'mean': {'X': ({'chr1:0-1000': 3.0, 'chr1:5000-6000': 3.0, 'chr1:9000-10000': 0.0},) * 2}}
    elements = node_elements(graph, presence, 2.0, {'declared': {'X': 'active'}}, tss_status(graph))
    table = element_agreement(elements, {'A': assigned})
    hits = table.set_index(['ega_class', 'external'])['fragments'].to_dict()
    assert hits == {('active promoter', 'PLS'): 1, ('active enhancer', 'dELS'): 1,
                    ('unmarked element', 'CTCF-only'): 1}
