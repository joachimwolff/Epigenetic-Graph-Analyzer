"""Supplementary Table S32: structural variants and copy number as explanation levels of the two AML comparisons
(egaDifferential --structuralVariants / --svContacts / --copyNumber, analysis/22_structural_variants.sh).

Part A: per comparison the strongest joins detected in the capture Hi-C read pairs (capture_sv_calls.py; at most
ten per comparison by read pairs, after the t(8;21)), with the breakends (100 kb bins), read
pairs and observed / expected per condition and the change between the conditions. Part B: the levels per
comparison: genes with the promoter within the window of a breakend (all, differential, enhancer-hijacking
candidates) and the differential genes explained by the structural variant and copy number levels (any level and
first level, declared roles, fragment mean). Part C: genes near the breakends of the t(8;21) and the differential
genes near any breakend, with the partner side, its active elements and the promoter's own contacts across the
junction.

Usage: python table_sv.py <results dir> <output dir> [joins per comparison]
"""

import os
import sys

import pandas as pd

COMPARISONS = (('siMM_vs_siRE', 'Kasumi-1 siMM vs siRE', ('siMM', 'siRE')),
               ('t821_vs_FLT3ITD', 't(8;21) vs FLT3-ITD', ('t821', 'FLT3ITD')))
FOCUS = ('RUNX1', 'RUNX1T1')


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('#', r'\#')


def num(x, digits=1):
    return '--' if pd.isna(x) else '{:,.{}f}'.format(x, digits).replace('-', '$-$')


def integer(x):
    return '--' if pd.isna(x) else '{:,}'.format(int(x))


def position(chrom, start):
    return '{}:{:.1f}'.format(chrom, (start + 50000) / 1e6)


def longtable(columns, header, rows, caption, label):
    lines = [r'{\footnotesize', r'\begin{longtable}{@{}' + columns + '@{}}', r'\caption{' + caption + r'}\label{'
             + label + r'}\\', r'\toprule', header + r' \\', r'\midrule', r'\endfirsthead', r'\toprule',
             header + r' \\', r'\midrule', r'\endhead']
    return lines + rows + [r'\bottomrule', r'\end{longtable}}']


def read(directory, comparison, suffix):
    path = os.path.join(directory, comparison + suffix)
    return pd.read_csv(path, sep='\t') if os.path.exists(path) else None


def part_a(directory, top):
    rows = []
    for comparison, label, (a, b) in COMPARISONS:
        joins = read(directory, comparison, '_sv_breakpoints.tsv')
        if joins is None:
            continue
        joins = joins.assign(score=joins[['pairs_' + a, 'pairs_' + b]].max(axis=1))
        focus = joins[joins['name'] == 't(8;21)']
        show = pd.concat([focus, joins.sort_values('score', ascending=False)]).drop_duplicates('sv_id').head(top)
        if rows:
            rows.append(r'\midrule')
        for i, r in enumerate(show.itertuples()):
            r = r._asdict()
            rows.append(' & '.join([label if i == 0 else '', tex(r['name']), position(r['chrom1'], r['start1']),
                                    position(r['chrom2'], r['start2']), integer(r['pairs_' + a]),
                                    integer(r['pairs_' + b]), num(r['obs_exp_' + a], 0), num(r['obs_exp_' + b], 0),
                                    tex(r['status'])]) + r' \\')
    header = 'Comparison & Join & Breakend 1 (Mb) & Breakend 2 (Mb) & Pairs A & Pairs B & O/E A & O/E B & Change'
    caption = ('Joins detected in the capture Hi-C read pairs of the two AML comparisons: inter-chromosomal and '
               'intra-chromosomal ($\\geq$5\\,Mb) 1\\,Mb bin pairs with at least 50 read pairs and an observed over '
               'expected of at least 20 in a condition, bin pairs within 2\\,Mb merged, breakends at 100\\,kb (the bin '
               'with the largest share of its inter-chromosomal coverage in pairs with the other side; resolution '
               'about 0.5\\,Mb). Pairs and O/E (observed / expected) summed over the 1\\,Mb bin pairs of the join; A / '
               'B: siMM / siRE (Kasumi-1, two libraries each) and t(8;21) / FLT3-ITD (first 30\\,M read pairs of one '
               'library per patient). Change: gained / lost (called in one condition, support changing at least '
               'twofold), present in both, similar support in both (called in one condition by the thresholds '
               'only). The t(8;21) first, then the joins with most read pairs.')
    return longtable('llllrrrrl', header, rows, caption, 'tab:s32')


def part_b(directory):
    rows = []
    for comparison, label, _ in COMPARISONS:
        genes = read(directory, comparison, '_sv_genes.tsv')
        counts = read(directory, comparison, '_explanation_levels.tsv')
        if genes is None or counts is None:
            continue
        counts = counts[(counts['roles'] == 'declared') & (counts['presence'] == 'mean')].set_index('level')

        def value(level, column):
            return integer(counts.loc[level, column]) if level in counts.index else '--'
        differential = genes[genes['differential'].astype(bool)]
        rows.append(' & '.join([label, integer(genes['node'].nunique()), integer(differential['node'].nunique()),
                                integer(genes.loc[genes['hijacking_candidate'].astype(bool), 'node'].nunique()),
                                integer(counts['differential_genes'].iloc[0]),
                                value('structural_variant', 'genes_any'), value('structural_variant', 'genes_first'),
                                value('copy_number', 'genes_any'), value('copy_number', 'genes_first')]) + r' \\')
    header = ('Comparison & Near breakend & Differential & Hijacking cand. & Differential genes & SV any & SV first & '
              'CN any & CN first')
    caption = ('Structural variant (SV) and copy number (CN) levels. Near breakend: promoter fragments with an '
               'expression value within 1\\,Mb of a breakend; differential among them; enhancer-hijacking '
               'candidates (join present in a condition, an active element on the partner side and contacts of the '
               'promoter across the junction). Differential genes of the comparison and those explained by the SV '
               'level (join gained with the expression rising or lost with it falling, an active partner element in '
               'the condition with the join and the junction contacts changing in the same direction) and by the CN '
               'level (relative copy number at the promoter changing by at least log2 0.5 with the expression), at '
               'any level and as the first level of the cascade (declared roles).')
    return longtable('lrrrrrrrr', header, rows, caption, 'tab:s32b')


def part_c(directory, top=25):
    rows = []
    for comparison, label, (a, b) in COMPARISONS:
        genes = read(directory, comparison, '_sv_genes.tsv')
        classes = read(directory, comparison, '_gene_classes.tsv')
        if genes is None:
            continue
        results = classes.set_index('node')['sv_result'] if classes is not None and 'sv_result' in classes else {}
        names = genes['genes'].fillna('').astype(str)
        focus = genes[(genes['sv_name'] == 't(8;21)') & names.map(lambda g: bool(set(g.split(',')) & set(FOCUS)))]
        differential = genes[genes['differential'].astype(bool)].assign(order=lambda d: -d['rna_log2FC'].abs())
        show = pd.concat([focus, differential.sort_values('order').drop(columns='order')]).drop_duplicates(
            ['node', 'sv_id']).head(top)
        if rows:
            rows.append(r'\midrule')
        for i, r in enumerate(show.itertuples()):
            r = r._asdict()
            gene = ','.join(str(r['genes']).split(',')[:2])
            result = results.get(r['node'], '') if r['differential'] else 'not differential'
            junction = '{} / {}'.format(integer(r.get('junction_pairs_' + a)), integer(r.get('junction_pairs_' + b)))
            rows.append(' & '.join([label if i == 0 else '', tex(gene), num(r['rna_log2FC'], 1), tex(r['sv_name']),
                                    num(r['distance'] / 1e6, 2), tex(r['sv_status']),
                                    tex(','.join(str(r['partner_genes']).split(',')[:2])
                                        if isinstance(r['partner_genes'], str) else ''),
                                    '{} / {} of {}'.format(int(r['partner_active_' + a]), int(r['partner_active_' + b]),
                                                           int(r['partner_nodes'])),
                                    junction, 'yes' if r['hijacking_candidate'] else 'no',
                                    tex(result if isinstance(result, str) else '')]) + r' \\')
    header = ('Comparison & Gene & log2FC & Join & Distance (Mb) & Change & Partner genes & Active partner nodes A / '
              'B & Junction pairs A / B & Hijacking & SV level')
    caption = ('Genes near the breakends: RUNX1 and RUNX1T1 at the t(8;21), then the differential genes near any '
               'breakend by absolute expression change (at most 25 per comparison). Distance from the promoter '
               'fragment to the breakend; partner genes: genes of the graph nodes within 1\\,Mb of the other '
               'breakend; active partner nodes: nodes there with an active element class (DNase present) in A / B; '
               'junction pairs: read pairs of the promoter fragment with the partner side in A / B; SV level: result '
               'of the structural variant level for the differential genes.')
    return longtable('llrlrllllll', header, rows, caption, 'tab:s32c')


def main(argv):
    results, output = argv[1], argv[2]
    top = int(argv[3]) if len(argv) > 3 else 10
    directory = os.path.join(results, 'sv')
    os.makedirs(output, exist_ok=True)
    lines = part_a(directory, top) + [''] + part_b(directory) + [''] + part_c(directory)
    with open(os.path.join(output, 'table_s32.tex'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv)
