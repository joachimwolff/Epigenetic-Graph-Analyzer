"""Combined step feature, mark term and contact term per explanation level (egaDifferential --stepFeature decomposed).

Usage: python step_decomposition_summary.py <decomposed_dir> <combined_hops_n_dir>

<decomposed_dir> holds the runs of 31_step_decomposition.sh (<comparison>_explanations_all_levels.tsv,
_gene_classes.tsv, _mark_roles.tsv), <combined_hops_n_dir> the runs of 16_hops_n.sh with the combined feature
(the same inputs). Per comparison, role set (declared; learned: the roles learned in each run) and feature
(combined, mark term, contact term): the differential genes first explained at, and explained at (any), step 1
marks, step 2 and step 3; among them those whose explaining marks do not all change concordantly at the mirrored
loci (beyond mirror, the mirror of the same feature); and for the mark term the per-gene permutation P < 0.05.
The non-step levels (promoter marks, contacts, promoter state, step 1 states) are those of the run. Drivers: per step
level the genes explained by the combined feature (any / first) split by the term that changes concordantly
(marks, contacts, both, neither). Learned roles per level (element type all and non-promoter element) learned on the
combined feature and on the mark term. Writes <decomposed_dir>/step_decomposition_{levels,drivers,roles,genes}.tsv.
"""

import os
import sys

import numpy as np
import pandas as pd

COMPARISONS = ('hESC_vs_hNEC', 'Md_vs_PA2')
STEP_LEVELS = {'step1_marks': 'direct', 'step2': 'indirect', 'step3': 'hop3'}
FEATURES = {'combined': ('concordant_combined_', '_combined_mirror'),
            'mark term': ('concordant_mark_term_', '_mirror'),
            'contact term': ('concordant_contact_term_', '_contacts_mirror')}


def is_true(series):
    return series.astype(str).str.lower() == 'true'


def signs(text):
    out = {}
    for item in str(text).split(';'):
        if ':' in item:
            mark, value = item.split(':')
            out[mark] = int(value)
    return out


def marks_of(text):
    return [m for m in str(text).split(',') if m and m != 'nan']


def per_gene(levels, classes, roles):
    """One row per gene and feature: first explaining level, any-explained step levels, beyond mirror."""
    classes = classes.set_index('node')
    rows = []
    for node, gene in levels.groupby('node', sort=False):
        gene = gene.sort_values('rank')
        direction = np.sign(float(gene['rna_log2FC'].iloc[0]))
        for feature, (prefix, mirror_suffix) in FEATURES.items():
            first, any_levels, beyond = None, [], {}
            for _, r in gene.iterrows():
                level = r['level']
                if level in STEP_LEVELS:
                    concordant = marks_of(r[prefix + roles])
                    explains = bool(concordant)
                    if explains:
                        sign = signs(r['signs_' + roles])
                        fl = STEP_LEVELS[level]
                        mirrored = 0
                        for m in concordant:
                            value = classes.loc[node].get('{}_{}{}'.format(m, fl, mirror_suffix), np.nan)
                            if pd.notna(value) and np.sign(value) * sign[m] == direction and abs(value) >= 1:
                                mirrored += 1
                        beyond[level] = mirrored < len(concordant)
                else:
                    explains = str(r['explains_{}_mean'.format(roles)]).lower() == 'true'
                if explains:
                    if first is None:
                        first = level
                    any_levels.append(level)
            p = {lv: float(gene.loc[gene['level'] == lv, 'p_' + roles].iloc[0]) for lv in STEP_LEVELS
                 if (gene['level'] == lv).any()}
            rows.append(dict(node=node, genes=gene['genes'].iloc[0], rna_log2FC=float(gene['rna_log2FC'].iloc[0]),
                             roles=roles, feature=feature, first=first or 'none', any=';'.join(any_levels),
                             **{'beyond_' + lv: beyond.get(lv, np.nan) for lv in STEP_LEVELS},
                             **{'p_mark_term_' + lv: p.get(lv, np.nan) for lv in STEP_LEVELS}))
    return pd.DataFrame(rows)


def level_counts(genes):
    rows = []
    for (roles, feature), part in genes.groupby(['roles', 'feature'], sort=False):
        for level in STEP_LEVELS:
            first = part['first'] == level
            any_ = part['any'].str.split(';').apply(lambda v: level in v)
            beyond = part['beyond_' + level].astype(float) == 1
            row = dict(roles=roles, feature=feature, level=level, differential_genes=len(part),
                       first=int(first.sum()), any=int(any_.sum()),
                       first_beyond_mirror=int((first & beyond).sum()), any_beyond_mirror=int((any_ & beyond).sum()))
            if feature == 'mark term':
                p = part['p_mark_term_' + level].astype(float) < 0.05
                row.update(first_p_lt_0_05=int((first & p).sum()), first_beyond_mirror_p_lt_0_05=int(
                    (first & beyond & p).sum()), any_p_lt_0_05=int((any_ & p).sum()))
            rows.append(row)
    return pd.DataFrame(rows)


def drivers(levels, genes, roles):
    combined = genes[(genes['roles'] == roles) & (genes['feature'] == 'combined')].set_index('node')
    rows = []
    for level in STEP_LEVELS:
        part = levels[levels['level'] == level].set_index('node')
        values = part['driver_' + roles].fillna('').astype(str)
        concordant = part['concordant_combined_' + roles].fillna('').astype(str) != ''
        for subset, nodes in (('any combined', part.index[concordant]),
                              ('first combined', combined.index[combined['first'] == level]),
                              ('any mark term', part.index[values.isin(['marks', 'both'])]),
                              ('any contact term', part.index[values.isin(['contacts', 'both'])])):
            counts = values.reindex(nodes).value_counts()
            rows.append(dict(roles=roles, level=level, genes=subset, n=len(nodes),
                             **{d: int(counts.get(d, 0)) for d in ('marks', 'contacts', 'both', 'neither')}))
    return pd.DataFrame(rows)


def roles_table(decomposed, combined, comparison):
    out = []
    for feature, directory in (('combined', combined), ('mark term', decomposed)):
        path = os.path.join(directory, comparison + '_mark_roles.tsv')
        if not os.path.exists(path):
            continue
        table = pd.read_csv(path, sep='\t')
        table = table[table['element_type'].isin(['all', 'non-promoter element']) &
                      (table['versus'] == 'expression')]
        out.append(table.assign(feature=feature, comparison=comparison)[
            ['comparison', 'feature', 'element_type', 'mark', 'level', 'n_genes', 'rho', 'rho_ci_low', 'rho_ci_high',
             'declared_role', 'learned_role']])
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    decomposed, combined = sys.argv[1], sys.argv[2]
    all_counts, all_drivers, all_roles, all_genes = [], [], [], []
    for comparison in COMPARISONS:
        prefix = os.path.join(decomposed, comparison)
        levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', keep_default_na=False,
                             na_values=[''])
        for column in [c for c in levels if c.startswith(('concordant_', 'driver_', 'signs_'))]:
            levels[column] = levels[column].fillna('')
        classes = pd.read_csv(prefix + '_gene_classes.tsv', sep='\t')
        genes = pd.concat([per_gene(levels, classes, roles) for roles in ('declared', 'learned')], ignore_index=True)
        counts = level_counts(genes).assign(comparison=comparison)
        # cross-check: the combined feature of this run against the 16_hops_n run (declared roles)
        base_path = os.path.join(combined, comparison + '_explanations_all_levels.tsv')
        if os.path.exists(base_path):
            base = pd.read_csv(base_path, sep='\t', keep_default_na=False, na_values=[''])
            for level in STEP_LEVELS:
                n = int(is_true(base.loc[base['level'] == level, 'explains_declared_mean']).sum())
                counts.loc[(counts['roles'] == 'declared') & (counts['feature'] == 'combined') &
                           (counts['level'] == level), 'any_in_hops_n_run'] = n
        all_counts.append(counts)
        all_drivers.append(pd.concat([drivers(levels, genes, roles) for roles in ('declared', 'learned')])
                           .assign(comparison=comparison))
        all_roles.append(roles_table(decomposed, combined, comparison))
        all_genes.append(genes.assign(comparison=comparison))
        # Hoxa3 (and any gene named on the command line of interest) at every step level
        hit = levels[levels['genes'].astype(str).str.split(',').apply(lambda g: 'Hoxa3' in g)]
        if len(hit):
            columns = ['genes', 'level', 'combined_changes', 'mark_term_changes', 'contact_term_changes',
                       'concordant_combined_declared', 'concordant_mark_term_declared',
                       'concordant_contact_term_declared', 'driver_declared', 'p_declared', 'q_declared',
                       'beyond_mirror_declared', 'p_learned', 'q_learned', 'driver_learned']
            print('== {} Hoxa3'.format(comparison))
            print(hit[hit['level'].isin(STEP_LEVELS)][columns].to_string(index=False))
            if os.path.exists(base_path):
                b = base[base['genes'].astype(str).str.split(',').apply(lambda g: 'Hoxa3' in g)]
                print('combined (16_hops_n):')
                print(b[b['level'].isin(STEP_LEVELS)][['genes', 'level', 'mark_changes', 'concordant_marks_declared',
                                                      'p_declared', 'q_declared', 'beyond_mirror_declared',
                                                      'p_learned']].to_string(index=False))
    counts = pd.concat(all_counts, ignore_index=True)
    counts = counts[['comparison'] + [c for c in counts if c != 'comparison']]
    table_drivers = pd.concat(all_drivers, ignore_index=True)
    table_drivers = table_drivers[['comparison'] + [c for c in table_drivers if c != 'comparison']]
    roles = pd.concat(all_roles, ignore_index=True)
    counts.to_csv(os.path.join(decomposed, 'step_decomposition_levels.tsv'), sep='\t', index=False)
    table_drivers.to_csv(os.path.join(decomposed, 'step_decomposition_drivers.tsv'), sep='\t', index=False)
    roles.to_csv(os.path.join(decomposed, 'step_decomposition_roles.tsv'), sep='\t', index=False)
    pd.concat(all_genes, ignore_index=True).to_csv(os.path.join(decomposed, 'step_decomposition_genes.tsv'),
                                                   sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 500):
        print(counts.to_string(index=False))
        print(table_drivers.to_string(index=False))
        print(roles.to_string(index=False))


if __name__ == '__main__':
    main()
