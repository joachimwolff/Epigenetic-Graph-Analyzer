"""egaExplorer on the synthetic egaDifferential run of test_report: one offline HTML file with the search index and
the embedded data."""

import os
import re

from epigeneticgraphanalyzer import explorer
from test_report import external_resources, run  # noqa: F401  (module fixture: the synthetic run)


def build(run, tmp_path, *extra):  # noqa: F811
    out = str(tmp_path / 'explorer')
    path = explorer.main(['-p', run['prefix'], '-o', out] + list(extra))
    with open(path) as handle:
        return path, handle.read()


def test_explorer_page(run, tmp_path):  # noqa: F811
    path, text = build(run, tmp_path)
    assert os.path.basename(path) == 'cond1_vs_cond2_explorer.html'
    # search box first, the URL hash and the inline script and style
    assert '<input id="q" type="search"' in text
    assert text.index('id="q"') < text.index('id="tabs"')
    assert '<style>' in text and 'prefers-color-scheme: dark' in text
    assert 'location.hash' in text and 'DecompressionStream' in text
    # no external resources, no em dash, no placeholder left
    assert external_resources(text) == []
    assert not re.search(r'<(?:script|link|img|iframe)[^>]+(?:src|href)=', text)
    assert 'https://' not in text and 'http://' not in text.replace('http://www.w3.org/2000/svg', '')
    assert '—' not in text
    assert not re.search(r'\{\{\w+\}\}', text)


def test_search_index_and_data(run, tmp_path):  # noqa: F811
    _, text = build(run, tmp_path)
    index = explorer.search_index(text)
    names = [g[0] for g in index['genes']]
    differential = {g[0] for g in index['genes'] if g[1]}
    classes = run['classes']
    expected = {g for genes in classes['genes'] for g in str(genes).split(',')}
    assert expected <= set(names)
    assert expected == differential
    assert index['meta']['names'] == ['cond1', 'cond2']
    data = explorer.decode(text)
    assert set(data) >= {'meta', 'graph', 'genes', 'table', 'communities'}
    graph = data['graph']
    n = len(graph['nodes']['s'])
    assert n == data['meta']['nNodes'] and len(graph['edges']['u']) == data['meta']['nEdges']
    assert all(0 <= u < n for u in graph['edges']['u'] + graph['edges']['v'])
    assert set(graph['nodes']['m']) == set(data['meta']['marks'])
    assert len(data['table']['rows']) == len(classes)
    # every differential gene has the text of egaReport and the ladder; gene names are tagged
    gene = sorted(differential)[0]
    entry = data['genes'][gene]
    assert entry['d'] == 1 and entry['glance'].startswith('Expression of <x-g>{}</x-g>'.format(gene))
    assert entry['ladder']['levels'][0] == 'promoter_marks'
    assert entry['evidence'] and entry['evidence'][0]['cols'][0] == 'Level'
    assert 'hub' in entry and '<span class="gene">' not in entry['hub']['h']
    assert all(0 <= i < n for i in entry['n'])
    # nodes that carry an explanation point into the graph
    for item in (x for e in data['genes'].values() for x in e.get('x', [])):
        assert 0 <= item[0] < n and item[1] >= 1
    # shared method notes are kept once
    assert len(data['meta']['notes']) == len(set(data['meta']['notes']))


def test_differential_only(run, tmp_path):  # noqa: F811
    _, text = build(run, tmp_path, '--differentialOnly')
    index = explorer.search_index(text)
    assert all(g[1] for g in index['genes'])
    data = explorer.decode(text)
    assert 'switching' not in data and 'triggers' not in data
