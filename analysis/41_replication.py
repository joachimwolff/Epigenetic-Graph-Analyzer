"""Replication of the explanations beyond the promoter and its direct contacts between two comparisons of the same
transition (here hESC / hNEC of Freire-Pritchett et al. and H1 / NPC of Jung et al.), per gene.

Genes: differentially expressed in both comparisons in the same direction, and in both without an explanation at the
promoter (marks, interactions) or step 1 (marks). For every such gene the script reports whether step 2 or step 3
explains it in each comparison (explains_declared_mean of the step levels, with P and q), and, for the genes explained
in both, the explaining elements of each comparison: the nodes at the explaining step whose explaining marks change
concordantly (sign of the mark x sign of the expression change x log2 change >= --minChange). An element of the
reference is replicated when an explaining element of the replication lies within --window bp of it (distance
between the fragments); control: the explaining elements of the replication mirrored around the promoter of the gene.

Outputs: <out>_genes.tsv (one row per gene), <out>_elements.tsv (one row per explaining element of the reference with
its distance to the nearest replication element and to the nearest mirrored one) and a summary on stdout.

Usage: python 41_replication.py <reference prefix> <reference deseq2> <replication prefix> <replication deseq2>
           <out prefix> [--window 10000] [--minChange 1]
"""

import argparse

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

EARLY = ['promoter_marks', 'promoter_contacts', 'step1_marks']


def log2_change(attrs, mark, a, b):
    return float(np.log2((float(attrs.get('{}_{}'.format(mark, b), 0) or 0) + 1) /
                         (float(attrs.get('{}_{}'.format(mark, a), 0) or 0) + 1)))


def load(prefix, deseq2):
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False)
    levels['explains'] = levels['explains_declared_mean'].astype(str) == 'True'
    graph = nx.read_graphml(prefix + '_union.graphml')
    names = [v for v in graph.graph.get('names', '').split(',') if v] or None
    expression = pd.read_csv(deseq2, sep='\t').dropna(subset=['padj']).drop_duplicates('gene_name')
    return levels, graph, names, expression.set_index('gene_name')


def conditions_of(graph):
    """Condition names from the mark attributes of the nodes (<mark>_<condition>), in attribute order."""
    attrs = next(iter(graph.nodes(data=True)))[1]
    found = []
    for key in attrs:
        if '_' in key and key.split('_', 1)[0].startswith('H3'):
            c = key.split('_', 1)[1]
            if c not in found:
                found.append(c)
    return found[:2]


def per_gene(levels, expression):
    """gene -> (node, sign, explained step or None, P, q) for genes without an early explanation."""
    early = set(levels[levels.explains & levels.level.isin(EARLY)].node)
    out = {}
    for node, sub in levels[levels.level.isin(['step2', 'step3'])].groupby('node', sort=False):
        genes = [g for g in str(sub.genes.iloc[0]).split(',') if g in expression.index]
        if node in early or not genes:
            continue
        fc = float(sub.rna_log2FC.iloc[0])
        gene = min(genes, key=lambda g: abs(expression.loc[g, 'log2FoldChange'] - fc))
        hit = sub[sub.explains].sort_values('level')
        step = hit.level.iloc[0] if len(hit) else None
        p = float(hit.p_declared.iloc[0]) if len(hit) else np.nan
        q = float(hit.q_declared.iloc[0]) if len(hit) else np.nan
        if gene not in out or (step and not out[gene][2]):
            out[gene] = (node, np.sign(fc), step, p, q)
    return out


def explaining_elements(graph, levels, node, step, sign, conditions, min_change):
    row = levels[(levels.node == node) & (levels.level == step)].iloc[0]
    signs = dict(x.split(':') for x in str(row.signs_declared).split(';') if ':' in x)
    marks = [m for m in str(row.concordant_marks_declared).split(',') if m in signs]
    k = int(step.replace('step', ''))
    distance = nx.single_source_shortest_path_length(graph, node, cutoff=k)
    a, b = conditions
    return [n for n, d in distance.items() if d == k and
            any(float(signs[m]) * sign * log2_change(graph.nodes[n], m, a, b) >= min_change for m in marks)]


def span(graph, n):
    return int(graph.nodes[n]['start']), int(graph.nodes[n]['end'])


def gap(x, y):
    return max(0, max(x[0], y[0]) - min(x[1], y[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('reference', 'reference_deseq2', 'replication', 'replication_deseq2', 'output'):
        parser.add_argument(name)
    parser.add_argument('--window', type=int, default=10000)
    parser.add_argument('--minChange', type=float, default=1.0)
    parser.add_argument('--padj', type=float, default=0.05)
    args = parser.parse_args()
    ref_levels, ref_graph, _, ref_expr = load(args.reference, args.reference_deseq2)
    rep_levels, rep_graph, _, rep_expr = load(args.replication, args.replication_deseq2)
    ref_cond, rep_cond = conditions_of(ref_graph), conditions_of(rep_graph)
    ref, rep = per_gene(ref_levels, ref_expr), per_gene(rep_levels, rep_expr)
    rows, elements = [], []
    for gene in sorted(set(ref) & set(rep)):
        if not (ref_expr.loc[gene, 'padj'] < args.padj and rep_expr.loc[gene, 'padj'] < args.padj):
            continue
        if ref[gene][1] != rep[gene][1]:
            continue
        row = dict(gene=gene, sign=int(ref[gene][1]), reference_node=ref[gene][0], replication_node=rep[gene][0],
                   reference_step=ref[gene][2], reference_P=ref[gene][3], reference_q=ref[gene][4],
                   replication_step=rep[gene][2], replication_P=rep[gene][3], replication_q=rep[gene][4])
        if ref[gene][2] and rep[gene][2]:
            ref_el = explaining_elements(ref_graph, ref_levels, ref[gene][0], ref[gene][2], ref[gene][1], ref_cond,
                                         args.minChange)
            rep_el = explaining_elements(rep_graph, rep_levels, rep[gene][0], rep[gene][2], rep[gene][1], rep_cond,
                                         args.minChange)
            centre = sum(span(rep_graph, rep[gene][0])) / 2
            rep_spans = [span(rep_graph, n) for n in rep_el]
            mirrored = [(int(2 * centre - e), int(2 * centre - s)) for s, e in rep_spans]
            for n in ref_el:
                s = span(ref_graph, n)
                d = min((gap(s, r) for r in rep_spans), default=np.nan)
                dm = min((gap(s, r) for r in mirrored), default=np.nan)
                elements.append(dict(gene=gene, element=n, distance=d, mirrored_distance=dm))
            row['reference_elements'] = len(ref_el)
            row['replication_elements'] = len(rep_el)
            within = [e for e in elements if e['gene'] == gene]
            row['replicated_elements'] = sum(e['distance'] <= args.window for e in within)
            row['mirrored_elements'] = sum(e['mirrored_distance'] <= args.window for e in within)
        rows.append(row)
    genes = pd.DataFrame(rows)
    genes.to_csv(args.output + '_genes.tsv', sep='\t', index=False)
    pd.DataFrame(elements).to_csv(args.output + '_elements.tsv', sep='\t', index=False)
    ref_x = genes.reference_step.notna()
    rep_x = genes.replication_step.notna()
    table = [[int((ref_x & rep_x).sum()), int((ref_x & ~rep_x).sum())],
             [int((~ref_x & rep_x).sum()), int((~ref_x & ~rep_x).sum())]]
    odds, p = fisher_exact(table)
    print('{} genes differentially expressed in both in the same direction without an early explanation in both'
          .format(len(genes)))
    print('explained beyond step 1: in both {}, reference only {}, replication only {}, neither {}; odds ratio {:.2f}, '
          'P {:.2g}'.format(table[0][0], table[0][1], table[1][0], table[1][1], odds, p))
    if elements:
        e = pd.DataFrame(elements)
        print('explaining elements of the reference within {} bp of a replication element: {} of {} ({:.0%}); '
              'mirrored control {} ({:.0%})'.format(args.window, int((e.distance <= args.window).sum()), len(e),
                                                     (e.distance <= args.window).mean(),
                                                     int((e.mirrored_distance <= args.window).sum()),
                                                     (e.mirrored_distance <= args.window).mean()))


if __name__ == '__main__':
    main()
