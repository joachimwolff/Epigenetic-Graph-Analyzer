"""egaReport on a small synthetic egaDifferential run: sections, gene reports, all-level table, no network access."""

import os
import re

import numpy as np
import pandas as pd
import pytest

from epigeneticgraphanalyzer import egaDifferential, egaExplain, egaReport
from synthetic import simulate_capture, write_bigwig, write_tss_bed


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    d = str(tmp_path_factory.mktemp('report'))
    da, db = os.path.join(d, 'a'), os.path.join(d, 'b')
    os.makedirs(da)
    os.makedirs(db)
    planted = sorted(simulate_capture(da, seed=7)[1])
    gained = {(c, i, j + 3) for c, i, j in planted[:20]}
    ibed_a = [simulate_capture(da, seed=s, planted=planted, n_noise=0, name='a{}'.format(s))[0] for s in (7, 9)]
    ibed_b, _, baits = simulate_capture(db, seed=8, planted=set(planted[20:]) | gained, n_noise=0)
    tss = write_tss_bed(os.path.join(d, 'tss.bed'), baits)
    rng = np.random.default_rng(0)
    bait_set = set(baits)
    up = {(c, i) for c, i in baits[::3]}
    ac_a = write_bigwig(os.path.join(d, 'ac_a.bw'), {})
    ac_b_values = {(c, j): 4.0 for c, _, j in gained}
    ac_b_values.update({b: 6.0 for b in up})
    ac_b = write_bigwig(os.path.join(d, 'ac_b.bw'), ac_b_values)
    me_a = write_bigwig(os.path.join(d, 'me_a.bw'), {b: 4.0 for b in bait_set})
    me_b = write_bigwig(os.path.join(d, 'me_b.bw'), {b: (1.0 if b in up else 4.0) for b in bait_set})
    expression = pd.DataFrame({'gene': ['GENE{}'.format(k) for k in range(len(baits))],
                               'A': rng.uniform(5, 50, len(baits)),
                               'B': [rng.uniform(150, 300) if b in up else rng.uniform(5, 50) for b in baits]})
    expression_path = os.path.join(d, 'expression.tsv')
    expression.to_csv(expression_path, sep='\t', index=False)
    prefix = os.path.join(d, 'out', 'cond1_vs_cond2')
    egaDifferential.main(['--interactionsA'] + ibed_a + ['--interactionsB', ibed_b, '-f', 'ibed', '--names', 'cond1',
                          'cond2', '--tss', tss, '--tracksA', 'H3K27ac=' + ac_a, 'H3K27me3=' + me_a,
                          '--tracksB', 'H3K27ac=' + ac_b, 'H3K27me3=' + me_b,
                          '--markRoles', 'H3K27ac=active', 'H3K27me3=repressive',
                          '--expression', expression_path, '--expressionColumns', 'A', 'B', '--geneKey', 'genes',
                          '--permutations', '5', '--bootstrap', '20', '--roleBootstrap', '20', '--roles', 'both',
                          '--genePermutations', '200', '--communityPermutations', '200', '--threads', '1',
                          '-o', prefix])
    explain = os.path.join(d, 'out', 'explain', 'cond1_vs_cond2')
    egaExplain.main(['-p', prefix, '-o', explain])
    classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
    return dict(dir=d, prefix=prefix, classes=classes, baits=baits, gained=gained)


def external_resources(text):
    """URLs the browser would load: src / href attributes, CSS url() and @import (not plain text)."""
    found = re.findall(r'(?:src|href)\s*=\s*["\']\s*(?:https?:)?//', text, flags=re.I)
    found += re.findall(r'url\(\s*["\']?\s*(?:https?:)?//', text, flags=re.I)
    found += re.findall(r'@import', text)
    return found


FORBIDDEN_IN_BODY = ('/mnt/', 'step1_state', 'promoter_marks', 'promoter_state', 'nan', 'explains, rank',
                     'egaDifferential (declared roles)', '_gene_classes.tsv', 'chr1:', '\u2014')


def check_body(text, run):
    """The report body (outside the appendix) reads as text: no paths, internal names, reprs or raw coordinates."""
    body = egaReport.body_outside_appendix(text)
    # the text of the body without the inline figures
    prose = re.sub(r'<svg.*?</svg>', '', body, flags=re.S)
    for word in FORBIDDEN_IN_BODY:
        assert word not in prose, word
    assert run['dir'] not in prose
    # no Python reprs (the word in a sentence, 'None of ...', is fine)
    assert not re.search(r'\bNone\b(?! of)', prose)
    assert not re.search(r"\[\'|\{\'|\bTrue\b|\bFalse\b", prose)
    # comma-joined lists keep a space after the comma
    assert not re.search(r'GENE\d+,GENE\d+', re.sub(r'<[^>]+>', '', prose))
    return body


def test_main_report(run):
    out = os.path.join(run['dir'], 'reports')
    outputs = egaReport.main(['-p', run['prefix'], '--top', '2', '--pdf', 'none', '-o', out])
    html_path = outputs[0][0]
    assert os.path.basename(html_path) == 'cond1_vs_cond2_report.html'
    text = open(html_path).read()
    for anchor in ('findings', 'cascade', 'genes', 'communities', 'appendix'):
        assert 'id="{}"'.format(anchor) in text, anchor
    assert 'cond1' in text and 'cond2' in text
    # figures are inline SVG with captions below; no external resources, no em-dashes
    assert text.count('<svg') >= 3
    assert re.search(r'</div><figcaption><span class="label">Figure 1\.</span>', text)
    assert external_resources(text) == []
    assert 'http' not in re.sub(r'<p[^>]*>.*?</p>', '', text, flags=re.S)
    body = check_body(text, run)
    # key findings and gene cards with the at-a-glance paragraph and the ladder
    assert '<ul class="findings">' in body and 'class="card"' in body and 'class="glance"' in body
    assert 'Explanation ladder' in body
    # the numbers of the summary come from the run outputs (appendix) and the key findings (body)
    edges = pd.read_csv(run['prefix'] + '_edges.tsv.gz', sep='\t')
    assert '{:,}'.format(int((edges['status'] == 'gained').sum())) in text
    assert '{:,} genes differ in expression'.format(len(run['classes'])) in body
    # the appendix is a collapsible section with the paths
    assert '<details><summary>' in text and run['prefix'] in text


def test_gene_reports(run):
    out = os.path.join(run['dir'], 'genes')
    classes = run['classes']
    explained = classes[classes['explained_at'] != 'unexplained']
    genes = [str(g).split(',')[0] for g in explained['genes'].head(2)]
    outputs = egaReport.main(['-p', run['prefix'], '--gene'] + genes + ['--pdf', 'none', '-o', out])
    pages = [o[0] for o in outputs]
    assert os.path.exists(os.path.join(out, 'cond1_vs_cond2_genes_index.html'))
    for gene in genes:
        path = os.path.join(out, 'cond1_vs_cond2_gene_{}.html'.format(gene))
        assert path in pages
        text = open(path).read()
        for anchor in ('glance', 'interactions', 'evidence', 'states', 'hub', 'appendix'):
            assert 'id="{}"'.format(anchor) in text, anchor
        body = check_body(text, run)
        # readable level labels in the evidence table, the expression change as a factor
        assert 'Promoter: histone marks' in body and 'Promoter: own contacts' in body
        assert re.search(r'Expression of <span class="gene">{}</span> is about [0-9.,]+-fold (higher|lower)'
                         .format(gene), body)
        assert 'Distance control' in body
        assert external_resources(text) == []
    index = open(os.path.join(out, 'cond1_vs_cond2_genes_index.html')).read()
    assert all('cond1_vs_cond2_gene_{}.html'.format(g) in index for g in genes)


def test_formatting():
    assert egaReport.fmt_pq(8.4999e-05) == '8.5 \u00d7 10<sup>\u22125</sup>'
    assert egaReport.fmt_pq(8.4999e-05, markup=False) == '8.5 \u00d7 10^-5'
    assert egaReport.fmt_pq(0.0116849) == '0.012'
    assert egaReport.fmt_pq(0.004975) == '0.005'
    assert egaReport.fold(10.604116) == '1,560'
    assert egaReport.fold(-1.73499) == '3.3'
    assert egaReport.fmt_log2(-0.04) == '0.0'
    assert egaReport.fmt_log2(-1.37) == '\u22121.4'
    assert egaReport.join_words(['a', 'b', 'c']) == 'a, b and c'
    assert egaReport.level_label('step1_state') == 'Direct partners: chromatin state'
    assert egaReport.level_label('step3') == '3 steps away'
    assert egaReport.readable_note('repressive mark gained') == 'gain of the repressive mark'
    assert egaReport.state_text('repressed -> bivalent: active mark gained', ('Md', 'PA2')) == \
        'repressed in Md, bivalent in PA2'
    assert egaReport.fmt_distance(-120400) == '120 kb'
    assert egaReport.fmt_distance(1037000) == '1.0 Mb'


def test_places(run):
    data = egaReport.Run(run['prefix'])
    graph = data.graph
    baits = [n for n, d in graph.nodes(data=True) if d.get('bait') and d.get('genes')]
    others = [n for n, d in graph.nodes(data=True) if not d.get('bait') and d.get('chrom') ==
              graph.nodes[baits[0]]['chrom']]
    gene = egaReport.gene_names(graph.nodes[baits[0]]['genes'])[0]
    places = egaReport.Places(graph)
    assert places.label(baits[0], markup=False) == '{} promoter'.format(gene)
    label = places.label(others[0], baits[0], gene, markup=False)
    assert re.fullmatch(r'element [0-9.]+ (bp|kb|Mb) from {}( \(near \S+\))?'.format(gene), label), label
    places.strands = {gene: '+'}
    assert re.search(r'(upstream|downstream) of {}'.format(gene), places.label(others[0], baits[0], gene,
                                                                               markup=False))


def test_places_name_element_classes(run):
    data = egaReport.Run(run['prefix'])
    graph = data.graph
    lookup = data.element_lookup()
    assert lookup and set(lookup) == set(graph.nodes)
    places = egaReport.Places(graph, elements=lookup, names=data.names)
    baits = [n for n, d in graph.nodes(data=True) if d.get('bait') and d.get('genes')]
    stable = next(n for n in baits if lookup[n][0] == lookup[n][1])
    gene = egaReport.gene_names(graph.nodes[stable]['genes'])[0]
    adjective = egaReport.PROMOTER_ADJECTIVE[lookup[stable][0]]
    assert places.label(stable, markup=False, article=True) == 'the {} promoter of {}'.format(adjective, gene)
    changed = next(n for n in baits if lookup[n][0] != lookup[n][1])
    assert re.search(r'promoter \(\w+ in cond1, \w+ in cond2\)$', places.label(changed, markup=False))
    # an other end gaining H3K27ac: an element named by its class in each condition
    gained = next(n for n, i in lookup.items() if i[0] == 'unmarked element' and i[1] == 'active enhancer')
    label = places.label(gained, stable, gene, markup=False, article=True)
    assert label.startswith('an element ') and 'unmarked element in cond1, active enhancer in cond2' in label, label
    other = next(n for n, i in lookup.items() if i[0] == i[1] == 'unmarked element' and
                 graph.nodes[n]['chrom'] == graph.nodes[stable]['chrom'])
    assert places.label(other, stable, gene, markup=False, article=True).startswith('an unmarked element ')
    text = egaReport.Report(data, None, False, places).class_counts_text(['active enhancer', 'active enhancer',
                                                                         'unmarked element'])
    assert text == 'two active enhancers and one unmarked element'


def test_levels_fallback_matches_cascade(run, tmp_path):
    # without the all-level table the report derives the levels from the other outputs
    import glob
    import shutil
    prefix = str(tmp_path / os.path.basename(run['prefix']))
    for path in glob.glob(run['prefix'] + '*'):
        if not path.endswith('_explanations_all_levels.tsv'):
            shutil.copy(path, str(tmp_path / os.path.basename(path)))
    data = egaReport.Run(prefix)
    assert data.levels_source.endswith('(derived)')
    first = data.first_level()
    # every gene explained by egaDifferential has an explaining level at its first cascade level
    for node, level in first.items():
        rows = data.levels[(data.levels['node'] == node) & data.levels['explains']]
        if level == 'unexplained':
            continue
        assert level in set(rows['cascade_level']), (node, level)
        assert rows['rank'].min() == 1


def test_all_level_table_is_read(run, tmp_path):
    """With <prefix>_explanations_all_levels.tsv the report uses it (gene x level rows with rank, P, q, mirror)."""
    import shutil
    prefix = str(tmp_path / 'copy')
    for path in os.listdir(os.path.dirname(run['prefix'])):
        if path.startswith('cond1_vs_cond2'):
            shutil.copy(os.path.join(os.path.dirname(run['prefix']), path),
                        str(tmp_path / path.replace('cond1_vs_cond2', 'copy')))
    classes = run['classes']
    node, genes = classes['node'].iloc[0], classes['genes'].iloc[0]
    common = dict(node=node, genes=genes, rna_log2FC=2.0)
    table = pd.DataFrame([
        dict(common, rank=1, level='promoter_marks', hops=0, nodes_tested=1, mark_changes='H3K27ac:+2.10',
             proximity_support=1.0, proximity_support_ci='NA', concordant_marks_declared='H3K27ac',
             explains_declared_mean=True, result_declared_mean='concordant', detail_declared_mean='H3K27ac +2.1',
             explains_learned_mean=True, result_learned_mean='concordant'),
        dict(common, rank=4, level='step1_marks', hops=1, nodes_tested=3, mark_changes='H3K27ac:+0.20',
             mirror_changes='H3K27ac:+0.10', proximity_support=0.66, proximity_support_ci='0.64-0.68',
             p_declared=0.5, q_declared=0.9, explains_declared_mean=False, result_declared_mean='below threshold',
             explains_learned_mean=False, result_learned_mean='below threshold'),
        dict(common, rank=6, level='step2', hops=2, nodes_tested=12, mark_changes='H3K27ac:+1.40',
             mirror_changes='H3K27ac:+0.20', proximity_support=0.58, proximity_support_ci='0.56-0.60',
             concordant_marks_declared='H3K27ac', n_concordant_nodes_declared=4, p_declared=0.004, q_declared=0.03,
             beyond_mirror_declared=True, explains_declared_mean=True, result_declared_mean='concordant',
             explains_learned_mean=False, result_learned_mean='no consistent role')])
    table.to_csv(prefix + '_explanations_all_levels.tsv', sep='\t', index=False)
    data = egaReport.Run(prefix, explain_prefix='')
    assert data.levels_source == 'copy_explanations_all_levels.tsv'
    rows = data.levels[data.levels['node'] == node].set_index('level')
    assert rows.loc['step2', 'rank'] == 2 and rows.loc['step2', 'q'] == 0.03
    assert rows.loc['step2', 'learned'].startswith('does not explain')
    assert bool(rows.loc['step2', 'beyond_mirror'])
    assert not rows.loc['step1', 'explains']
    gene = str(genes).split(',')[0]
    out = str(tmp_path / 'out')
    outputs = egaReport.main(['-p', prefix, '--explainPrefix', '', '--gene', gene, '--pdf', 'none', '-o', out])
    text = open(outputs[0][0]).read()
    assert 'copy_explanations_all_levels.tsv' in text
    # P and q of the supported step in the at-a-glance paragraph and the evidence table
    body = egaReport.body_outside_appendix(text)
    assert 'q = 0.03' in body and '2 steps away' in body
    assert 'copy_explanations_all_levels.tsv' not in body


def test_matplotlib_pdf(run):
    out = os.path.join(run['dir'], 'pdf')
    classes = run['classes']
    gene = str(classes.loc[classes['explained_at'] != 'unexplained', 'genes'].iloc[0]).split(',')[0]
    outputs = egaReport.main(['-p', run['prefix'], '--gene', gene, '--pdf', 'matplotlib', '-o', out])
    html_path, pdf_path, engine = outputs[0]
    assert engine == 'matplotlib'
    with open(pdf_path, 'rb') as handle:
        assert handle.read(4) == b'%PDF'


def test_community_report(run):
    table = pd.read_csv(run['prefix'] + '_gene_community_changes.tsv', sep='\t')
    community = int(table.loc[table['community_b'] >= 0, 'community_b'].iloc[0])
    out = os.path.join(run['dir'], 'community')
    outputs = egaReport.main(['-p', run['prefix'], '--community', str(community), '--pdf', 'none', '-o', out])
    text = open(outputs[0][0]).read()
    assert 'Community {}'.format(community) in text
    assert external_resources(text) == []
