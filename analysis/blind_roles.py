"""Blind test of the role discovery (analysis/18_blind_roles.sh): the roles learned for tracks given without a role,
per level and element type of the nodes, against the known biology of the marks.

Known biology (expected sign of the change of the mark on expression): H3K9me3 repressive (-1, constitutive and
facultative heterochromatin); H3K36me3 transcription-associated (+1, deposited over transcribed gene bodies; a baited
promoter fragment covers mostly the promoter, so the signal there may be weak); promoter DNA methylation repressive
(-1) and gene-body DNA methylation positive with expression (+1); H2A.Z marks active and poised promoters and
enhancers and has no single expected sign (0). Marks not listed have no expectation.

Usage: python blind_roles.py <egaDifferential prefix> <output TSV> <mark> [<mark> ...]
Writes one row per mark, level and element type: n_genes, rho with its 95% chromosome block-bootstrap interval,
learned_role, expected_role and verdict ('as expected', 'opposite', 'no consistent effect', 'no expectation').
"""

import sys

import pandas as pd

EXPECTED = {'H3K9me3': -1, 'H3K36me3': 1, 'H2A.Z': 0, 'DNAme_promoter': -1, 'DNAme_genebody': 1}
NAMES = {1: 'active', -1: 'repressive', 0: 'no single expectation'}


def verdict(expected, learned):
    if expected is None or expected == 0:
        return 'no expectation'
    if learned == 0:
        return 'no consistent effect'
    return 'as expected' if learned == expected else 'opposite'


def main():
    prefix, output = sys.argv[1:3]
    marks = sys.argv[3:]
    roles = pd.read_csv(prefix + '_mark_roles.tsv', sep='\t')
    if 'element_type' not in roles:
        roles['element_type'] = 'all'
    roles = roles[roles['mark'].isin(marks) & (roles['versus'] == 'expression')].copy()
    roles['expected_role'] = roles['mark'].map(lambda m: NAMES.get(EXPECTED.get(m), 'none'))
    roles['verdict'] = [verdict(EXPECTED.get(m), int(s)) for m, s in zip(roles['mark'], roles['learned_sign'])]
    columns = ['mark', 'element_type', 'level', 'n_genes', 'rho', 'rho_ci_low', 'rho_ci_high', 'concordance',
               'learned_role', 'expected_role', 'verdict']
    roles = roles[columns]
    roles.to_csv(output, sep='\t', index=False)
    main_types = roles[roles['element_type'].isin(['all', 'promoter', 'non-promoter element'])]
    with pd.option_context('display.width', 200, 'display.max_rows', 500):
        print(main_types.to_string(index=False))
        print()
        print(roles.groupby(['mark', 'verdict']).size().unstack(fill_value=0))


if __name__ == '__main__':
    main()
