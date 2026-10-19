"""Weightings of the step feature (egaDifferential --stepWeighting) compared on the mark term.

Usage: python step_weighting_summary.py <out_dir> <product_dir> <weighting> [<weighting> ...]

<product_dir> holds the runs of 31_step_decomposition.sh with the product weighting, <out_dir>/<weighting> those of
33_step_weighting.sh (each with step_decomposition_genes.tsv of step_decomposition_summary.py). Per comparison,
weighting and role set (declared, learned), on the mark term:
  first     differential genes by first explaining level, grouped as promoter (marks, contacts, state), step1
            (marks, state), step2, step3, none, and by level;
  beyond    per step level: genes explained (any) / first explained there, with explaining marks beyond the mirror,
            with per-gene permutation P < 0.05, and both;
  models    nested-model delta R2 of direct, step2 (indirect), step3 (hop3), the delta R2 of the mirrored feature,
            their difference and its bootstrap interval (permutation.json mirror_comparison);
  hoxa3     the step-2 row of the Hoxa3 fragment chr6:52207514-52213060 (Md / PA2);
  jaccard   overlap of the sets of genes first explained at step 2 between the weightings.
Writes <out_dir>/step_weighting_{first,beyond,models,hoxa3,jaccard}.tsv.
"""

import itertools
import json
import os
import sys

import pandas as pd

COMPARISONS = ('hESC_vs_hNEC', 'Md_vs_PA2')
GROUPS = {'promoter_marks': 'promoter', 'promoter_contacts': 'promoter', 'promoter_state': 'promoter',
          'step1_marks': 'step1', 'step1_state': 'step1', 'step2': 'step2', 'step3': 'step3', 'none': 'none'}
STEP_LEVELS = {'step1_marks': 'direct', 'step2': 'indirect', 'step3': 'hop3'}
HOXA3 = 'chr6:52207514-52213060'
ROLES = ('declared', 'learned')


def main():
    out, product = sys.argv[1], sys.argv[2]
    directories = {'product': product, **{w: os.path.join(out, w) for w in sys.argv[3:]}}
    first_rows, beyond_rows, model_rows, hoxa_rows, sets = [], [], [], [], {}
    for weighting, directory in directories.items():
        genes = pd.read_csv(os.path.join(directory, 'step_decomposition_genes.tsv'), sep='\t')
        genes = genes[genes['feature'] == 'mark term']
        for comparison in COMPARISONS:
            part_c = genes[genes['comparison'] == comparison]
            for roles in ROLES:
                part = part_c[part_c['roles'] == roles]
                grouped = part['first'].map(GROUPS).value_counts()
                detailed = part['first'].value_counts()
                first_rows.append(dict(comparison=comparison, weighting=weighting, roles=roles,
                                       differential_genes=len(part),
                                       **{g: int(grouped.get(g, 0)) for g in dict.fromkeys(GROUPS.values())},
                                       **{'first_' + lv: int(detailed.get(lv, 0)) for lv in GROUPS}))
                for level in STEP_LEVELS:
                    any_ = part['any'].fillna('').str.split(';').apply(lambda v, lv=level: lv in v)
                    first = part['first'] == level
                    beyond = part['beyond_' + level].astype(float) == 1
                    p = part['p_mark_term_' + level].astype(float) < 0.05
                    beyond_rows.append(dict(comparison=comparison, weighting=weighting, roles=roles, level=level,
                                            any=int(any_.sum()), any_beyond_mirror=int((any_ & beyond).sum()),
                                            any_p_lt_0_05=int((any_ & p).sum()),
                                            any_beyond_mirror_p_lt_0_05=int((any_ & beyond & p).sum()),
                                            first=int(first.sum()), first_beyond_mirror=int((first & beyond).sum()),
                                            first_beyond_mirror_p_lt_0_05=int((first & beyond & p).sum())))
                sets[(comparison, roles, weighting)] = set(part.loc[part['first'] == 'step2', 'node'])
            prefix = os.path.join(directory, comparison)
            nested = pd.read_csv(prefix + '_nested_models.tsv', sep='\t')
            with open(prefix + '_permutation.json') as handle:
                mirror = json.load(handle)['mirror_comparison']
            for row, (level, key) in enumerate(STEP_LEVELS.items(), start=1):
                model_rows.append(dict(comparison=comparison, weighting=weighting, level=level,
                                       delta_r2=float(nested['delta_r2'].iloc[row]),
                                       delta_r2_mirror=mirror.get('gain_{}_mirror'.format(key)),
                                       graph_minus_mirror=mirror.get('{}_graph_minus_mirror'.format(key)),
                                       ci_low=mirror.get('{}_ci_low'.format(key)),
                                       ci_high=mirror.get('{}_ci_high'.format(key)),
                                       fraction_bootstrap_above_0=mirror.get(
                                           '{}_fraction_bootstrap_above_0'.format(key))))
            levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', keep_default_na=False,
                                 na_values=[''], low_memory=False)
            hit = levels[(levels['node'] == HOXA3) & (levels['level'] == 'step2')]
            for _, r in hit.iterrows():
                row = dict(comparison=comparison, weighting=weighting, node=HOXA3, genes=r['genes'],
                           nodes_tested=r['nodes_tested'], mark_changes=r['mark_changes'],
                           mirror_changes=r['mirror_changes'])
                for roles in ROLES:
                    for column in ('concordant_marks', 'p', 'q', 'n_concordant_nodes', 'beyond_mirror'):
                        row['{}_{}'.format(column, roles)] = r.get('{}_{}'.format(column, roles))
                hoxa_rows.append(row)
    jaccard_rows = []
    for comparison in COMPARISONS:
        for roles in ROLES:
            for a, b in itertools.combinations(directories, 2):
                s_a, s_b = sets[(comparison, roles, a)], sets[(comparison, roles, b)]
                union = s_a | s_b
                jaccard_rows.append(dict(comparison=comparison, roles=roles, weighting_a=a, weighting_b=b,
                                         n_a=len(s_a), n_b=len(s_b), shared=len(s_a & s_b), union=len(union),
                                         jaccard=len(s_a & s_b) / len(union) if union else float('nan')))
    tables = {'first': first_rows, 'beyond': beyond_rows, 'models': model_rows, 'hoxa3': hoxa_rows,
              'jaccard': jaccard_rows}
    with pd.option_context('display.width', 250, 'display.max_columns', 60, 'display.max_rows', 500):
        for name, rows in tables.items():
            table = pd.DataFrame(rows)
            table.to_csv(os.path.join(out, 'step_weighting_{}.tsv'.format(name)), sep='\t', index=False)
            print('==', name)
            print(table.to_string(index=False))


if __name__ == '__main__':
    main()
