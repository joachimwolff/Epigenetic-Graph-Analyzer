"""Key results of the example run (example/run_example.sh), the files compared with example/expected/:
  graphs.tsv           nodes and edges of the Md, PA2 and union graphs
  cascade.tsv          explanation cascade of the differentially expressed genes (<prefix>_explanation_cascade.tsv)
  hoxa3_summary.tsv    every level of the two Hoxa3 promoter fragments (declared mark roles, fragment mean)
  reports.tsv          the reports written by egaReport

Usage: python summarize_example.py <example output dir> <summary dir>
"""

import os
import sys

import networkx as nx
import pandas as pd

HOXA3_COLUMNS = ['node', 'genes', 'rna_log2FC', 'level', 'hops', 'nodes_tested', 'explains_declared_mean',
                 'result_declared_mean', 'strongest_declared', 'p_declared', 'q_declared', 'beyond_mirror_declared',
                 'detail_declared_mean']
REPORTS = ['Md_vs_PA2_report.html', 'Md_vs_PA2_gene_Hoxa3.html', 'Md_vs_PA2_gene_Hoxa2.html',
           'Md_vs_PA2_genes_index.html']


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    run, out = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    rows = []
    for name, path in (('Md', 'graphs/Md.graphml'), ('PA2', 'graphs/PA2.graphml'),
                       ('union', 'Md_vs_PA2_union.graphml')):
        g = nx.read_graphml(os.path.join(run, path))
        rows.append({'graph': name, 'nodes': g.number_of_nodes(), 'edges': g.number_of_edges()})
    pd.DataFrame(rows).to_csv(os.path.join(out, 'graphs.tsv'), sep='\t', index=False)

    cascade = pd.read_csv(os.path.join(run, 'Md_vs_PA2_explanation_cascade.tsv'), sep='\t')
    cascade.to_csv(os.path.join(out, 'cascade.tsv'), sep='\t', index=False, float_format='%.6g')

    levels = pd.read_csv(os.path.join(run, 'Md_vs_PA2_explanations_all_levels.tsv'), sep='\t')
    hoxa3 = levels[levels['genes'].astype(str).str.split(',').apply(lambda genes: 'Hoxa3' in genes)]
    hoxa3[HOXA3_COLUMNS].to_csv(os.path.join(out, 'hoxa3_summary.tsv'), sep='\t', index=False, float_format='%.6g')

    present = [{'report': name, 'written': os.path.getsize(os.path.join(run, 'reports', name)) > 0
                if os.path.exists(os.path.join(run, 'reports', name)) else False} for name in REPORTS]
    pd.DataFrame(present).to_csv(os.path.join(out, 'reports.tsv'), sep='\t', index=False)
    explained = hoxa3[hoxa3['explains_declared_mean'].astype(str) == 'True']
    print('Hoxa3: levels with an observation:',', '.join(sorted(set(explained['level']))) or 'none')


if __name__ == '__main__':
    main()
