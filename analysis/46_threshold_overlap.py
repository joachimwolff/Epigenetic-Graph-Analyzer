"""Overlap of the explanations beyond step 1 between two runs of the same comparison (for example CHiCAGO score 5 and
8, 45_threshold_sensitivity.sh), per promoter fragment.

For each run, the promoter fragments without an explanation at the promoter (marks, contacts, state) or at step 1
(marks, state) and with a candidate at step 2 or 3 are listed, all and with a per-gene P < 0.05 beyond the mirrored
loci. Reported: the fragments of the first run with a candidate, and how many of them are unexplained before step 2 in
the second run, receive a candidate there (at any of the two steps, and at the same step), or are explained earlier
in the second run.

Usage: python 46_threshold_overlap.py <run prefix 1> <run prefix 2> <out.tsv>
"""

import sys

import pandas as pd

EARLY = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state']


def candidates(prefix):
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False,
                         usecols=['node', 'genes', 'level', 'explains_declared_mean', 'p_declared',
                                  'beyond_mirror_declared'])
    levels['explains'] = levels.explains_declared_mean.astype(str) == 'True'
    levels['supported'] = levels.explains & (levels.p_declared < 0.05) & \
        (levels.beyond_mirror_declared.astype(str) == 'True')
    early = set(levels[levels.explains & levels.level.isin(EARLY)].node)
    later = levels[~levels.node.isin(early) & levels.level.isin(['step2', 'step3'])]
    first = later[later.explains].sort_values('level').drop_duplicates('node').set_index('node').level
    supported = later[later.supported].sort_values('level').drop_duplicates('node').set_index('node').level
    return set(levels.node), early, first, supported


def main():
    one, two, output = sys.argv[1:4]
    nodes_1, early_1, first_1, supported_1 = candidates(one)
    nodes_2, early_2, first_2, supported_2 = candidates(two)
    rows = []
    for name, set_1, set_2 in (('any', first_1, first_2), ('supported', supported_1, supported_2)):
        tested = [n for n in set_1.index if n in nodes_2]
        rows.append(dict(
            candidates=name, run1=len(set_1), run2=len(set_2), run1_tested_in_run2=len(tested),
            explained_earlier_in_run2=sum(n in early_2 for n in tested),
            candidate_in_run2=sum(n in set_2.index for n in tested),
            same_step_in_run2=sum(n in set_2.index and set_2[n] == set_1[n] for n in tested),
            any_candidate_in_run2=sum(n in first_2.index for n in tested)))
    table = pd.DataFrame(rows)
    table.to_csv(output, sep='\t', index=False)
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
