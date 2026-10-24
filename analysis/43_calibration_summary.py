"""Explanations beyond step 1 in the real run and in the calibration nulls of 42_calibration.sh.

For every run, the promoter fragments of differentially expressed genes (or of unchanged genes posing as such) without
an explanation at the promoter (marks, contacts, state) or at step 1 (marks, state) are counted, and among them the
share with a candidate explanation at step 2 or, failing that, at step 3 (explains_declared_mean): any, with a per-gene
P < 0.05 beyond the mirrored loci, and with q < 0.05. Rows: real, then each null run (<out>/<null>_<seed> or
<out>/bin<width>/rewired_<seed>), with the rewiring statistics where present.

Usage: python 43_calibration_summary.py <real run prefix> <calibration dir>
"""

import glob
import os
import sys

import pandas as pd

EARLY = ['promoter_marks', 'promoter_contacts', 'promoter_state', 'step1_marks', 'step1_state']


def summary(prefix):
    levels = pd.read_csv(prefix + '_explanations_all_levels.tsv', sep='\t', low_memory=False,
                         usecols=['node', 'level', 'explains_declared_mean', 'p_declared', 'q_declared',
                                  'beyond_mirror_declared'])
    levels['explains'] = levels.explains_declared_mean.astype(str) == 'True'
    early = set(levels[levels.explains & levels.level.isin(EARLY)].node)
    later = levels[~levels.node.isin(early) & levels.level.isin(['step2', 'step3'])].copy()
    later['supported'] = later.explains & (later.p_declared < 0.05) & \
        (later.beyond_mirror_declared.astype(str) == 'True')
    later['q05'] = later.explains & (later.q_declared < 0.05)
    row = dict(fragments=levels.node.nunique(), unexplained_early=later.node.nunique())
    for name, sub in later.groupby('level'):
        row[name] = int(sub.explains.sum())
        row[name + '_supported'] = int(sub.supported.sum())
        row[name + '_q05'] = int(sub.q05.sum())
    per_node = later.groupby('node')[['explains', 'supported', 'q05']].any()
    for column in per_node:
        row['beyond_' + column] = int(per_node[column].sum())
        row['share_' + column] = per_node[column].mean()
    return row


def main():
    real, out = sys.argv[1], sys.argv[2]
    rows = [dict(run='real', **summary(real))]
    for prefix in sorted(glob.glob(os.path.join(out, '*', 'naiveB_vs_GM12878_explanations_all_levels.tsv')) +
                         glob.glob(os.path.join(out, 'bin*', '*', 'naiveB_vs_GM12878_explanations_all_levels.tsv'))):
        prefix = prefix[:-len('_explanations_all_levels.tsv')]
        run = os.path.relpath(os.path.dirname(prefix), out)
        row = dict(run=run, **summary(prefix))
        log = os.path.join(os.path.dirname(prefix), 'rewiring.log')
        if os.path.exists(log):
            last = open(log).read().strip().splitlines()[-1]
            stats = dict(x.split('=') for x in last.split('\t') if '=' in x)
            row['edges_changed'] = float(stats.get('fraction_edges_changed', 'nan'))
        rows.append(row)
    table = pd.DataFrame(rows)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    print(table.to_string(index=False, float_format=lambda x: '{:.3f}'.format(x)))
    table.to_csv(os.path.join(out, 'calibration_summary.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
