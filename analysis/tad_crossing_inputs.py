"""Inputs of egaBoundaries for the two conditions of a differential run, from its union edges and union graph.

usage: tad_crossing_inputs.py <egaDifferential prefix> <bigWig with the chromosome sizes> <hicFindTADs boundaries.bed
       ...> <output dir>

Writes <out>/<condition>_nodes.tsv and _edges.tsv (the edges called in the condition, weight = its CHiCAGO score;
format of egaCreateGraph), <out>/chrom.sizes (from the bigWig header) and, per hicFindTADs boundary file,
<out>/<name>_chr.bed with UCSC chromosome names (reference of egaBoundaries --reference).
"""

import os
import sys

import pandas as pd
import pyBigWig

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import tads as td  # noqa: E402
from epigeneticgraphanalyzer.genome import normalize_chrom, parse_node_id  # noqa: E402
from epigeneticgraphanalyzer.graph import read_graph  # noqa: E402


def main(prefix, bigwig, *rest):
    *boundary_files, out = rest
    os.makedirs(out, exist_ok=True)
    edges = td.read_union_edges(prefix)
    names = td.condition_names(edges)
    graph = read_graph(prefix + '_union.graphml')
    in_a, in_b = td.presence(edges, names)
    for name, flag in zip(names, (in_a, in_b)):
        part = edges[flag]
        e = pd.DataFrame({'node_a': part['node_a'].values, 'node_b': part['node_b'].values,
                          'weight': part['score_{}'.format(name)].values,
                          'score': part['score_{}'.format(name)].values, 'distance': part['distance'].values})
        nodes = sorted(set(e['node_a']) | set(e['node_b']))
        rows = []
        for n in nodes:
            c, s, t = parse_node_id(n)
            d = graph.nodes[n] if n in graph else {}
            rows.append(dict(node=n, chrom=c, start=s, end=t, bait=bool(d.get('bait', False)),
                             genes=d.get('genes', '')))
        pd.DataFrame(rows).to_csv(os.path.join(out, '{}_nodes.tsv'.format(name)), sep='\t', index=False)
        e.to_csv(os.path.join(out, '{}_edges.tsv'.format(name)), sep='\t', index=False)
        print('{}: {} nodes, {} edges'.format(name, len(rows), len(e)))
    handle = pyBigWig.open(bigwig)
    with open(os.path.join(out, 'chrom.sizes'), 'w') as sizes:
        for chrom, size in handle.chroms().items():
            sizes.write('{}\t{}\n'.format(chrom, size))
    for path in boundary_files:
        bed = pd.read_csv(path, sep='\t', header=None)
        bed[0] = bed[0].astype(str).map(normalize_chrom)
        name = os.path.basename(path).replace('.bed', '')
        bed.to_csv(os.path.join(out, '{}_chr.bed'.format(name)), sep='\t', header=False, index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
