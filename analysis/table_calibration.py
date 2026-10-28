"""LaTeX tables of the calibration with unchanged genes (42_calibration.sh, 42b_calibration_hesc_md.sh), of the links of
the explaining elements (44_explaining_elements_eqtl.py) and of the threshold sensitivity (46_threshold_overlap.py).

Writes <out>/table_calibration_rows.tex (share of promoter fragments with a concordant observation per level, real
run against the mean and range of the unchanged-gene runs) and <out>/table_links_rows.tex (explaining elements against
window, distance and graph controls, per gene and per element).

Usage: python table_calibration.py <results dir> <out dir>
"""

import glob
import os
import sys

import pandas as pd

LEVELS = [('promoter_marks', 'Promoter marks'), ('promoter_contacts', 'Promoter interactions'),
          ('promoter_state', 'Promoter state'), ('step1_marks', 'Step 1 marks'), ('step1_state', 'Step 1 states'),
          ('step2', 'Step 2 marks'), ('step3', 'Step 3 marks')]
COMPARISONS = [('hESC / hNEC', 'revision_bioinf/decomposed/hESC_vs_hNEC', 'calibration_hesc_md/hESC_unchanged_*',
                'hESC_vs_hNEC'),
               ('Md / PA2', 'revision_bioinf/decomposed/Md_vs_PA2', 'calibration_hesc_md/Md_unchanged_*', 'Md_vs_PA2'),
               ('naive B / GM12878', 'gm12878_naiveB/naiveB_vs_GM12878', 'calibration/unchanged_*',
                'naiveB_vs_GM12878')]
LINKS = [('ENCODE-rE2G', 're2g'), ('GTEx LCL eQTL', 'gtex_lcl'), ('GEUVADIS fine-mapped eQTL', 'geuvadis_susie')]


def levels(prefix):
    t = pd.read_csv(prefix + '_explanation_levels.tsv', sep='\t')
    t = t[(t.roles == 'declared') & (t.presence == 'mean')].set_index('level')
    return t.fraction_any, int(t.differential_genes.iloc[0])


def pct(x):
    return '{:.1f}'.format(100 * x).replace('-', '$-$')


def main():
    results, out = sys.argv[1], sys.argv[2]
    columns, header = [], []
    for name, real, nulls, run in COMPARISONS:
        observed, n = levels(os.path.join(results, real))
        runs = [levels(p[:-len('_explanation_levels.tsv')])[0]
                for p in sorted(glob.glob(os.path.join(results, nulls, run + '_explanation_levels.tsv')))]
        null = pd.concat(runs, axis=1)
        columns.append((observed, null))
        header.append((name, n, len(runs)))
    lines = []
    for key, label in LEVELS:
        cells = []
        for observed, null in columns:
            cells.append('{} & {} ({}--{})'.format(pct(observed[key]), pct(null.loc[key].mean()),
                                                   pct(null.loc[key].min()), pct(null.loc[key].max())))
        lines.append('{} & {} \\\\'.format(label, ' & '.join(cells)))
    lines.append('\\midrule')
    lines.append('Fragments & {} \\\\'.format(' & '.join('\\multicolumn{{2}}{{c}}{{{:,} / {} runs}}'.format(n, k)
                                                         for _, n, k in header)))
    with open(os.path.join(out, 'table_calibration_rows.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')

    lines = []
    for label, key in LINKS:
        s = pd.read_csv(os.path.join(results, 'explaining_eqtl', 'gm12878_{}_p1_summary.tsv'.format(key)),
                        sep='\t').set_index('set')

        def ci(name):
            return '{} ({} to {})'.format(pct(s.loc[name, 'share']), pct(s.loc[name, 'low']), pct(s.loc[name, 'high']))
        lines.append('{} & {} & {} & {} & {} & {} & {:g} / {:g} & {} / {} & {} / {} \\\\'.format(
            label, int(s.genes.iloc[0]), pct(s.loc['eqtl', 'share']), pct(s.loc['window', 'share']),
            pct(s.loc['distance_share', 'share']), pct(s.loc['graph', 'share']),
            s.loc['median_elements_per_gene', 'share'], s.loc['median_window_per_gene', 'share'],
            pct(s.loc['pooled_elements_linked', 'share']), pct(s.loc['pooled_window_linked', 'share']),
            pct(s.loc['genes_with_linked_element', 'share']), pct(s.loc['genes_with_linked_window', 'share'])))
        lines.append(' & \\multicolumn{{9}}{{l}}{{difference to window {}, to distance {}, to graph {} '
                     'percentage points (95\\% interval)}} \\\\'.format(
                         ci('eqtl_minus_window'), ci('eqtl_minus_distance_share'), ci('eqtl_minus_graph')))
    with open(os.path.join(out, 'table_links_rows.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    print(open(os.path.join(out, 'table_calibration_rows.tex')).read())
    print(open(os.path.join(out, 'table_links_rows.tex')).read())


if __name__ == '__main__':
    main()
