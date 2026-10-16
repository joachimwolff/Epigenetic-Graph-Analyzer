"""Supplementary Tables S1 to S13 of the manuscript as LaTeX, generated from the result directory.

S1 data sets and graphs; S2 multi-way support by graph distance and pair class; S3 per-contact confirmation;
S4 promoter hubs; S5 communities with a coordinated expression change; S6 nested models and mirrored loci;
S7 classification of differentially expressed genes; S8 candidates for regulation through indirect contacts;
S9 virtual capture benchmark; S10 CRISPRi benchmark; S11 inter-chromosomal rearrangements and copy number of
Kasumi-1; S12 chromatin assortativity (ChAseR) by graph distance; S13 promoter hubs against a polymer model.

Usage: python supplement_tables.py <results dir> <output dir>
"""

import json
import os
import re
import sys

import pandas as pd

R = None


def path(*parts):
    return os.path.join(R, *parts)


def tex(text):
    return str(text).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('#', r'\#')


def num(x, digits=3):
    if pd.isna(x):
        return '--'
    return '{:.{}f}'.format(x, digits).replace('-', '$-$')


def pval(p):
    if pd.isna(p):
        return '--'
    if p == 0:
        return '$<10^{-300}$'
    if p >= 0.01:
        return '{:.2f}'.format(p)
    mantissa, exponent = '{:.1e}'.format(p).split('e')
    return '${}\\times10^{{{}}}$'.format(mantissa, int(exponent))


def integer(x):
    return '{:,}'.format(int(x))


def table(name, columns, header, rows, caption, label, size=r'\small', long=False, colsep=None):
    """LaTeX table; normal tables shrink to the text width if wider (adjustbox), long tables use `colsep`."""
    spacing = r'\setlength{\tabcolsep}{%s}' % colsep if colsep else ''
    lines = []
    if long:
        lines += ['{' + size + spacing, r'\begin{longtable}{@{}' + columns + '@{}}', r'\caption{' + caption + r'}\label{'
                  + label + r'}\\', r'\toprule', header + r' \\', r'\midrule', r'\endfirsthead', r'\toprule',
                  header + r' \\', r'\midrule', r'\endhead']
        lines += rows
        lines += [r'\bottomrule', r'\end{longtable}}']
    else:
        lines += [r'\begin{table}[!h]', r'\centering' + size + spacing,
                  r'\begin{adjustbox}{max width=\textwidth}', r'\begin{tabular}{@{}' + columns + '@{}}', r'\toprule',
                  header + r' \\', r'\midrule']
        lines += rows
        lines += [r'\bottomrule', r'\end{tabular}', r'\end{adjustbox}', r'\caption{' + caption + '}',
                  r'\label{' + label + '}', r'\end{table}']
    with open(os.path.join(OUT, name), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def graph_stats(prefix):
    nodes = pd.read_csv(prefix + '_nodes.tsv', sep='\t', usecols=['bait'])
    edges = sum(1 for _ in open(prefix + '_edges.tsv')) - 1
    communities = sum(1 for _ in open(prefix + '_communities.tsv')) - 1
    return len(nodes), int(nodes['bait'].astype(str).eq('True').sum()), edges, communities


def condition_graphs(prefix):
    """Nodes, baits and edges of the graphs of both conditions, from the egaDifferential log."""
    stats = []
    for line in open(prefix + '.log'):
        match = re.search(r'graph: (\d+) nodes \((\d+) baits\), (\d+) edges', line)
        if match:
            stats.append(tuple(int(x) for x in match.groups()))
    return stats[-2:]


def s1():
    graphs = [
        ('GM12878', 'Cairns et al. 2016, GSE81503', 'hg19', 'CHiCAGO, score $\\geq$ 5', 'graphs/GM_chicago'),
        ('GM12878', 'Mifsud et al. 2015', 'hg19', 'strong core, $\\geq$ 20 reads', 'graphs/GM_mifsud_c20_r1'),
        ('H9 hESC', 'Freire-Pritchett et al. 2017, GSE86821', 'hg19', 'CHiCAGO, score $\\geq$ 12', 'graphs/hESC'),
        ('H9 hNEC', 'Freire-Pritchett et al. 2017, GSE86821', 'hg19', 'CHiCAGO, score $\\geq$ 12', 'graphs/hNPC'),
        ('H1-hESC', 'Jung et al. 2019, raw reads', 'hg19', 'HiCUP, CHiCAGO, score $\\geq$ 5', 'graphs/H1'),
        ('K562 (DMSO)', 'GSE250135', 'GRCh38', 'CHiCAGO, score $\\geq$ 5', 'k562/K562_dmso'),
        ('GM12878 Hi-C', 'virtual capture, \\citet{rao2014}', 'hg19', 'CHiCAGO (HiCExplorer), 5\\,kb',
         'virtual_capture/graphs/VC_GM12878'),
        ('H1-hESC Micro-C', 'virtual capture, \\citet{krietenstein2020}', 'GRCh38', 'CHiCAGO (HiCExplorer), 5\\,kb',
         'virtual_capture/graphs/VC_H1_microc'),
    ]
    rows = []
    for cell, source, assembly, calls, prefix in graphs:
        n, baits, e, c = graph_stats(path(prefix))
        rows.append(' & '.join([cell, source, assembly, calls, integer(n), integer(baits), integer(e), integer(c)])
                    + r' \\')
    rows.append(r'\midrule')
    pairs = [
        (('H9 hESC', 'H9 hNEC'), 'Freire-Pritchett et al. 2017', 'hg19', 'CHiCAGO, score $\\geq$ 12',
         'differential/hESC_vs_hNEC'),
        (('Md (mouse)', 'PA2 (mouse)'), 'Kessler et al. 2023, GSE211902', 'mm10', 'CHiCAGO, 2 replicates$^a$',
         'kessler/Md_vs_PA2'),
        (('Kasumi-1 siMM', 'Kasumi-1 siRE'), 'Ptasinska et al. 2019, raw reads', 'GRCh38',
         'CHiCAGO, 2 replicates$^a$', 'aml_kasumi1/siMM_vs_siRE'),
        (('AML t(8;21)', 'AML FLT3-ITD'), 'Assi et al. 2019, GSE108316', 'GRCh38', 'HOMER, FDR $\\leq$ 0.01',
         'aml_assi/t821_vs_FLT3ITD'),
    ]
    for cells, source, assembly, calls, prefix in pairs:
        for cell, (n, baits, e) in zip(cells, condition_graphs(path(prefix))):
            rows.append(' & '.join([cell, source, assembly, calls, integer(n), integer(baits), integer(e), '--'])
                        + r' \\')
    table('table_s1.tex', 'p{3.2cm}p{4.6cm}lp{3.6cm}rrrr',
          'Cell type & Source & Assembly & Interaction calls & Nodes & Baits & Edges & Communities', rows,
          'Promoter interaction graphs. Edges are cis interactions within 2\\,Mb; communities by Leiden '
          '(resolution 1). Bait status follows the capture design (hg19 HindIII design of \\citet{mifsud2015}; '
          'HindIII fragments holding a TSS for Kasumi-1; promoter bins for virtual capture). Lower part: graphs of '
          'the conditions of the differential comparisons, built by egaDifferential. $^a$Present in both '
          'replicates and score $\\geq$ 5 in at least one.', 'tab:s1')


def s2():
    hops = [('GM12878 (CHiCAGO)', 'GM12878 SPRITE', 'hops/GM_chicago_sprite_GM12878_hg19_HindIII_fragments.tsv'),
            ('H9 hESC', 'H1-hESC GAM', 'hops/hESC_gam_H1_hg38_30kb.tsv'),
            ('H1-hESC', 'H1-hESC GAM', 'hops/H1_gam_H1_hg38_30kb.tsv'),
            ('GM12878, virtual', 'GM12878 SPRITE', 'virtual_capture/graphs/VC_GM12878_sprite_hops.tsv'),
            ('H1-hESC, virtual', 'H1-hESC GAM', 'virtual_capture/graphs/VC_H1_microc_gam_hops.tsv')]
    rows = [r'\multicolumn{9}{@{}l}{\textit{By graph distance}} \\']
    for graph, assay, f in hops:
        t = pd.read_csv(path(f), sep='\t').set_index('steps')
        rows.append(' & '.join([graph, assay] + [num(t.loc[k, 'concordance']) for k in range(1, 6)]
                               + ['{:.2f}--{:.2f}'.format(t.loc[1, 'median_distance'] / 1e6,
                                                          t.loc[5, 'median_distance'] / 1e6), '']) + r' \\')
    rows += [r'\midrule', r'\multicolumn{9}{@{}l}{\textit{By pair class}} \\',
             r' & & Direct & Indirect & Other-other & Triplets & Communities$^a$ & Communities$^b$ & \\', r'\midrule']
    runs = [('GM12878 (CHiCAGO)', 'GM12878 SPRITE', 'GM_chicago_sprite'),
            ('GM12878 (Mifsud)', 'GM12878 SPRITE', 'GM_mifsud_sprite'),
            ('H9 hESC', 'H1-hESC GAM', 'hESC_gam'), ('H1-hESC', 'H1-hESC GAM', 'H1_gam'),
            ('H9 hESC, cross-cell', 'GM12878 SPRITE', 'hESC_sprite_crosscell'),
            ('GM12878, cross-cell', 'H1-hESC GAM', 'GM_chicago_gam_crosscell'),
            ('H9 hNEC, cross-state', 'H1-hESC GAM', 'hNEC_gam_crossstate')]
    for graph, assay, name in runs:
        t = pd.read_csv(path('validation', name + '_summary.tsv'), sep='\t')
        t = t[t['stratum'] == 'all']

        def get(s, cls, metric):
            hit = t[(t['set'] == s) & (t['pair_class'] == cls) & (t['metric'] == metric)]
            return hit['fraction_above_control'].iloc[0] if len(hit) else float('nan')
        rows.append(' & '.join([graph, assay, num(get('direct', 'all', 'npmi')), num(get('indirect', 'all', 'npmi')),
                                num(get('indirect', 'oe-oe', 'npmi')), num(get('triplets', 'all', 'c_uvw')),
                                num(get('communities', 'null_shift', 'groups_with_3_members')),
                                num(get('communities', 'null_locus', 'groups_with_3_members')), '']) + r' \\')
    table('table_s2.tex', 'llrrrrrrl',
          'Graph & Multi-way data & 1 step & 2 steps & 3 steps & 4 steps & 5 steps & Median distance (Mb) &', rows,
          'Multi-way support of the interaction graphs: probability that a graph item co-occurs more in the '
          'multi-way data than a matched control (0.5 expected without relation). Top, fragment pairs at graph '
          'distance 1 to 5 (50,000 per distance) against pairs at the same offset from the same pools of captured '
          'fragments and coverage classes; the median genomic distance of the pairs is given for 1 and 5 steps. '
          'Bottom, NPMI concordance for direct interactions, indirect contacts, indirect contacts between two other '
          'ends, mediated triplets, and the number of multi-way groups holding at least three loci of a community '
          'against the community shifted along the chromosome ($^a$any position; $^b$positions in captured '
          'chromatin). Cross-cell runs test a graph against multi-way data of another cell type.', 'tab:s2')


DISTANCES = ('0-100kb', '100-250kb', '250-500kb', '500-1000kb', '1000-2000kb')


def s3():
    runs = [('GM12878 (CHiCAGO)', 'SPRITE', 'GM_chicago_sprite'), ('H9 hESC', 'GAM', 'hESC_gam'),
            ('H1-hESC', 'GAM', 'H1_gam')]
    rows = []
    for graph, assay, name in runs:
        t = pd.read_csv(path('validation', name + '_confirm_summary.tsv'), sep='\t')
        for s in ('direct', 'indirect'):
            sub = t[t['set'] == s]
            overall = sub[(sub['pair_class'] == 'all') & (sub['stratum'] == 'all')].iloc[0]
            strata = [('all', overall['n'], overall['confirmed_p'], overall['control_confirmed_p'])]
            # distance strata are reported per pair class: pooled, weighted by the number of pairs
            classes = sub[(sub['pair_class'] != 'all') & sub['stratum'].isin(DISTANCES)]
            for stratum in DISTANCES:
                part = classes[classes['stratum'] == stratum]
                if len(part):
                    n = part['n'].sum()
                    strata.append((stratum.replace('-', '--').replace('kb', '\\,kb'), n,
                                   (part['n'] * part['confirmed_p']).sum() / n,
                                   (part['n'] * part['control_confirmed_p']).sum() / n))
            for i, (_, r) in enumerate(sub[sub['stratum'].str.startswith('weight_quartile')].iterrows()):
                strata.append(('score quartile {}'.format(r['stratum'][-1]), r['n'], r['confirmed_p'],
                               r['control_confirmed_p']))
            for stratum, n, confirmed, control in strata:
                rows.append(' & '.join([graph, assay, s, stratum, integer(n), '{:.1f}'.format(100 * confirmed),
                                        '{:.1f}'.format(100 * control)]) + r' \\')
        rows.append(r'\midrule')
    table('table_s3.tex', 'lllrrrr',
          'Graph & Multi-way data & Contacts & Stratum & Pairs & Exceeding (\\%) & Controls exceeding (\\%)', rows[:-1],
          'Per-contact support: fraction of direct interactions and indirect contacts whose distance-normalised '
          'co-occurrence exceeds the 95th percentile of their matched controls (5\\% expected), and the same '
          'fraction for held-out controls; all contacts, by genomic distance (pair classes pooled) and, for direct '
          'interactions, by quartile of the interaction score.', 'tab:s3', size=r'\scriptsize', long=True, colsep='4pt')


def s4():
    rows = []
    for graph, assay, f in (('GM12878 (CHiCAGO)', 'SPRITE', 'hubs/GM_chicago_sprite_hub_summary.tsv'),
                            ('H9 hESC', 'GAM', 'hubs/hESC_gam_hub_summary.tsv')):
        t = pd.read_csv(path(f), sep='\t')
        for _, r in t.iterrows():
            low, _, high = str(r['hub_size']).lstrip('>=').partition('-')
            size = '$\\geq$ ' + low if str(r['hub_size']).startswith('>=') else (
                low if high in ('', low) else low + '--' + high)
            rows.append(' & '.join([graph, assay, size, integer(r['hubs']),
                                    '{:.1f}'.format(100 * r['hubs_with_promoter_plus2_in_one_group']),
                                    '{:.1f}'.format(100 * r['shifted']), num(r['median_groups_promoter_plus2'], 1),
                                    num(r['shifted_median_groups'], 1), num(r['mean_largest_in_one_group'], 2),
                                    num(r['shifted_mean_largest'], 2)]) + r' \\')
        rows.append(r'\midrule')
    table('table_s4.tex', 'lllrrrrrrr',
          'Graph & Data & Hub size & Hubs & With group (\\%) & Shifted (\\%) & Groups & Shifted & Largest & Shifted',
          rows[:-1],
          'Promoter hubs (a promoter and its direct partners) in multi-way groups. With group: hubs with at least '
          'one SPRITE cluster or GAM tube that contains the promoter and at least two of its partners; groups: '
          'median number of such groups per hub; largest: mean largest number of hub elements in one group. '
          'Shifted: the same hubs shifted along the chromosome (ten shifts per hub). The 30\\,kb windows of GAM '
          'contain most hubs completely, so that the GAM data do not discriminate hubs from shifted hubs.',
          'tab:s4')


def community_genes(prefix):
    """Gene names per community of A and of B, ordered by absolute expression change."""
    changes = pd.read_csv(prefix + '_gene_community_changes.tsv', sep='\t',
                          usecols=['genes', 'community_a', 'community_b', 'rna_log2FC'])
    changes = changes.dropna(subset=['rna_log2FC'])
    changes = changes.reindex(changes['rna_log2FC'].abs().sort_values(ascending=False).index)
    result = {}
    for side in ('community_a', 'community_b'):
        result[side] = {}
        for community, group in changes.dropna(subset=[side]).groupby(side):
            names = []
            for genes in group['genes'].astype(str):
                for gene in genes.split(','):
                    if gene not in names and not gene.startswith(('RP11-', 'RP1-', 'RP4-', 'RP5-', 'CTD-', 'AC0',
                                                                  'AL')):
                        names.append(gene)
            result[side][int(community)] = names
    return result


def s5():
    rows = []
    for label, prefix, italic in (('hESC / hNEC', 'differential/hESC_vs_hNEC', False),
                                  ('Md / PA2', 'kessler/Md_vs_PA2', True),
                                  ('Kasumi-1 siMM / siRE', 'aml_kasumi1/siMM_vs_siRE', False),
                                  ('AML t(8;21) / FLT3-ITD', 'aml_assi/t821_vs_FLT3ITD', False)):
        t = pd.read_csv(path(prefix + '_community_expression.tsv'), sep='\t')
        t = t[t['coherence_q'] < 0.1].sort_values('coherence_q')
        genes = community_genes(path(prefix))
        for _, r in t.iterrows():
            if r['community_b'] >= 0:
                names = genes['community_b'].get(int(r['community_b']), [])[:8]
            else:
                names = genes['community_a'].get(int(r['community_a']), [])[:8]
            rows.append(' & '.join([label, str(int(r['community_b'])) if r['community_b'] >= 0 else '--',
                                    r['transition'], str(int(r['n_genes'])),
                                    '{}/{}/{}'.format(int(r['edges_shared']), int(r['edges_lost']),
                                                      int(r['edges_gained'])),
                                    num(r['rewired_fraction'], 2), num(r['mean_rna_log2FC'], 2),
                                    num(r['coherence_q'], 3), (('\\textit{' + tex(', '.join(names)) + '}') if italic else tex(', '.join(names))) or '--']) + r' \\')
        rows.append(r'\midrule')
    table('table_s5.tex', 'p{1.9cm}rp{1.4cm}rrrrr>{\\raggedright\\arraybackslash}p{4.4cm}',
          'Comparison & Community & Transition & Genes & Shared/lost/gained & Rewired & Mean $\\log_2$FC & $q$ & '
          'Genes with the largest change', rows[:-1],
          'Communities of the second condition with a coordinated expression change of their genes (FDR $<$ 0.1, '
          '10,000 random gene sets of the same size). Transition: fate of the matched community of the first '
          'condition (--, no matched community in the second condition; genes of the first); shared/lost/gained: interactions of the community; rewired: fraction of lost and gained '
          'interactions. The AML patient comparison '
          'uses an exploratory expression contrast of three patients per subtype.', 'tab:s5', size=r'\scriptsize',
          long=True, colsep='3pt')


COMPARISONS = (('hESC / hNEC', 'differential/hESC_vs_hNEC'), ('Md / PA2', 'kessler/Md_vs_PA2'),
               ('Kasumi-1 siMM / siRE', 'aml_kasumi1/siMM_vs_siRE'),
               ('AML t(8;21) / FLT3-ITD', 'aml_assi/t821_vs_FLT3ITD'))


def s6():
    rows = []
    for label, prefix in COMPARISONS:
        nested = pd.read_csv(path(prefix + '_nested_models.tsv'), sep='\t').set_index('model')
        perm = json.load(open(path(prefix + '_permutation.json')))
        m = perm['mirror_comparison']
        rows.append(' & '.join([
            label, integer(nested.loc['promoter', 'n_genes']), num(nested.loc['promoter', 'r2'], 4),
            num(nested.loc['promoter+direct', 'delta_r2'], 4), pval(nested.loc['promoter+direct', 'f_pvalue']),
            num(nested.loc['promoter+direct+indirect', 'delta_r2'], 4),
            pval(nested.loc['promoter+direct+indirect', 'f_pvalue']), pval(perm['empirical_p_delta_r2']),
            num(m['gain_direct_mirror'], 4), num(m['gain_indirect_mirror'], 4),
            '{} [{}, {}]'.format(num(m['indirect_graph_minus_mirror'], 4), num(m['indirect_ci_low'], 4),
                                 num(m['indirect_ci_high'], 4)),
            '{:.0f}'.format(100 * m['indirect_fraction_bootstrap_above_0'])]) + r' \\')
    table('table_s6.tex', 'lrrrrrrrrrlr',
          'Comparison & Genes & $R^2$ promoter & $\\Delta R^2$ direct & $P$ & $\\Delta R^2$ indirect & $P$ & '
          '$P$ perm. & Mirror direct & Mirror indirect & Indirect graph $-$ mirror [95\\% CI] & Above 0 (\\%)', rows,
          'Variance of RNA $\\log_2$ fold changes explained by chromatin changes at the promoter, and the '
          'additional variance explained by contact-weighted changes at direct and at indirect neighbours (nested '
          'linear models, $F$ test). $P$ perm.: indirect level against 200 permutations of the node signals. '
          'Mirror: the same gains for loci at the mirrored genomic offsets; the difference between graph and '
          'mirrored loci for the indirect level with its 95\\% bootstrap interval and the fraction of 1,000 '
          'bootstrap resamples above zero.', 'tab:s6', size=r'\scriptsize\setlength{\tabcolsep}{4pt}')


def s7():
    rows = []
    for label, prefix in COMPARISONS:
        perm = json.load(open(path(prefix + '_permutation.json')))
        counts = perm['class_counts']
        test = pd.read_csv(path(prefix + '_indirect_gene_test.tsv'), sep='\t')
        rows.append(' & '.join([label, integer(perm['n_differential_genes']), integer(sum(counts.values()))]
                               + [integer(counts.get(k, 0)) for k in ('promoter', 'direct', 'indirect', 'unexplained')]
                               + ['{:.1f} ({:.1f} $\\pm$ {:.1f})'.format(100 * perm['observed_fraction_indirect'],
                                                                         100 * perm['null_fraction_indirect_mean'],
                                                                         100 * perm['null_fraction_indirect_sd']),
                                  '{} / {}'.format(int((test['indirect_p'] < 0.05).sum()), len(test)),
                                  str(int((test['indirect_q'] < 0.1).sum()))]) + r' \\')
    table('table_s7.tex', 'lrrrrrrlrr',
          'Comparison & DE genes & DE promoter fragments & Promoter & Direct & Indirect only & Unexplained & '
          'Indirect only, \\% (null) & '
          '$P<0.05$ & $q<0.1$', rows,
          'Hierarchical classification of the baited promoter fragments of differentially expressed genes (several genes can share a fragment) by the first level (promoter, direct '
          'neighbours, indirect neighbours) with a concordant chromatin change, the fraction explained at the '
          'indirect level only against 200 permutations of the node signals (mean $\\pm$ s.d.), and the per-gene '
          'permutation test of the concordant change at indirect neighbours (1,000 permutations, '
          'Benjamini-Hochberg). Candidates are listed in Table~S8.', 'tab:s7')


def s8():
    rows = []
    for label, prefix in COMPARISONS:
        t = pd.read_csv(path(prefix + '_indirect_only_genes.tsv'), sep='\t')
        t = t[t['beyond_mirror'].astype(bool) & (t['indirect_p'] < 0.05)]
        t = t.reindex(t['rna_log2FC'].abs().sort_values(ascending=False).index)
        for _, r in t.iterrows():
            genes = [g for g in str(r['genes']).split(',') if g != 'nan'][:3]
            drivers = []
            for item in str(r['drivers']).split('; '):
                gene = item.split('[')[-1].split(']')[0] if '[' in item else ''
                if gene and gene not in drivers:
                    drivers.append(gene)
            rows.append(' & '.join([label, '\\textit{' + tex(', '.join(genes)) + '}', num(r['rna_log2FC'], 2),
                                    str(int(r['n_direct'])), str(int(r['n_indirect'])), tex(str(r['explaining_marks']).replace(',', ', ')),
                                    pval(r['indirect_p']), num(r['indirect_q'], 2),
                                    tex(', '.join(dict.fromkeys(g.strip() for d in drivers[:5] for g in d.split(','))))]) + r' \\')
        rows.append(r'\midrule')
    table('table_s8.tex', 'p{1.5cm}>{\\raggedright\\arraybackslash}p{2.1cm}rrr>{\\raggedright\\arraybackslash}p{1.7cm}rr'
          '>{\\raggedright\\arraybackslash}p{3.5cm}',
          'Comparison & Gene & RNA $\\log_2$FC & Direct & Indirect & Marks & $P$ & $q$ & Explaining two-step nodes',
          rows[:-1],
          'Candidates for regulation through indirect contacts: differentially expressed genes explained only at '
          'indirect neighbours, without a concordant change of the explaining marks at the mirrored loci and with '
          'a per-gene permutation $P<0.05$. With the measured marks, a change two steps away is the best '
          'available explanation of the expression change of these genes; perturbation experiments are required '
          'to establish regulation. Direct and indirect: numbers of direct and two-step neighbours; explaining '
          'two-step nodes: annotated genes of the nodes with the largest concordant contribution (OE, other end).',
          'tab:s8', size=r'\scriptsize', long=True, colsep='3pt')


def s9():
    rows = []
    for label, f in (('GM12878 Hi-C vs CHiCAGO', 'GM12878_benchmark.tsv'),
                     ('GM12878 CHiCAGO vs Mifsud (capture)', 'GM12878_cairns_vs_mifsud20_benchmark.tsv'),
                     ('K562 Hi-C vs CHiCAGO', 'K562_benchmark.tsv'),
                     ('H1-hESC Micro-C vs CHiCAGO, all calls', 'H1_microc_benchmark.tsv'),
                     ('H1-hESC Micro-C vs CHiCAGO, matched number', 'H1_microc_matched_benchmark.tsv')):
        t = pd.read_csv(path('virtual_capture', f), sep='\t')
        for stratum in t['stratum'].unique():
            exact = t[(t['stratum'] == stratum) & (t['tolerance_bins'] == 0)].iloc[0]
            tol = t[(t['stratum'] == stratum) & (t['tolerance_bins'] == 1)].iloc[0]
            rows.append(' & '.join([label if stratum == 'all' else '', tex(stratum), integer(exact['shared_bait_bins']),
                                    integer(exact['real_pairs']), integer(exact['virtual_pairs']),
                                    '{:.1f} ({:.1f})'.format(100 * exact['precision'], 100 * tol['precision']),
                                    '{:.1f} ({:.1f})'.format(100 * exact['recall'], 100 * tol['recall'])]) + r' \\')
        rows.append(r'\midrule')
    table('table_s9.tex', 'llrrrrr',
          'Comparison & Distance & Shared baits & Reference calls & Virtual calls & Precision (\\%) & Recall (\\%)',
          rows[:-1],
          'Virtual promoter capture: agreement of CHiCAGO calls from Hi-C or Micro-C contact matrices (5\\,kb) with '
          'promoter capture calls of the same cell type on shared baited bins, exact and, in brackets, with a '
          'tolerance of one bin. Two capture experiments of GM12878 serve as reference for the agreement between '
          'capture data sets.', 'tab:s9', size=r'\scriptsize', long=True, colsep='4pt')


def s10():
    names = {'direct': 'direct interaction', 'two_hop': 'two-step path', 'hops': 'graph distance',
             'community': 'same community', 'ppr': 'personalised PageRank', 'contact_weight': 'contact weight'}
    rows = []
    for label, f in (('Training', 'bench_training_added_value.tsv'), ('Held out', 'bench_heldout_added_value.tsv')):
        t = pd.read_csv(path('k562', f), sep='\t')
        t = t[t['set'] == 'scorable']
        for _, r in t.iterrows():
            rows.append(' & '.join([label, tex(names.get(r['feature'], r['feature'])), tex(r['base']),
                                    num(r['coefficient'], 3), num(r['lr_statistic'], 2), pval(r['lr_p']),
                                    integer(r['n']), integer(r['positives'])]) + r' \\')
            label = ''
        rows.append(r'\midrule')
    table('table_s10.tex', 'lllrrrrr',
          'Set & Feature & Base model & Coefficient & LR statistic & $P$ & Pairs & Positives', rows[:-1],
          'CRISPRi benchmark of K562 enhancer-gene pairs \\citep{fulco2019,gasperini2019,gschwind2026} with a '
          'promoter in the K562 graph: added value of each graph feature over genomic distance and enhancer degree '
          'in a logistic model (likelihood ratio test), on the training pairs and on held-out pairs.', 'tab:s10',
          size=r'\scriptsize')


def s11():
    reps = {}
    for rep in ('siMM_rep1', 'siRE_rep1'):
        t = pd.read_csv(path('aml_kasumi1', 'sv', rep + '_trans_hotspots.tsv'), sep='\t')
        reps[rep] = t.set_index(['chrom1', 'start1', 'chrom2', 'start2'])
    a = reps['siMM_rep1']
    a = a[a['pairs'] >= 100].sort_values('obs_exp', ascending=False)
    rows = []
    seen = set()
    for key, r in a.iterrows():
        region = (key[0], key[1] // 3000000, key[2], key[3] // 3000000)
        if region in seen:
            continue
        seen.add(region)
        other = reps['siRE_rep1'].loc[key] if key in reps['siRE_rep1'].index else None
        rows.append(' & '.join(['{}:{}'.format(key[0], key[1] // 1000000), '{}:{}'.format(key[2], key[3] // 1000000),
                                integer(r['pairs']), num(r['expected'], 1), num(r['obs_exp'], 0),
                                integer(other['pairs']) if other is not None else '--',
                                num(other['obs_exp'], 0) if other is not None else '--']) + r' \\')
        if len(rows) == 10:
            break
    rows.append(r'\midrule')
    cn = {rep: pd.read_csv(path('aml_kasumi1', 'sv', rep + '_copy_number_chromosomes.tsv'), sep='\t')
          .set_index('chrom')['relative_coverage'] for rep in ('siMM_rep1', 'siRE_rep1')}
    chroms = list(cn['siMM_rep1'].index)
    for i in range(0, len(chroms), 4):
        cells = []
        for c in chroms[i:i + 4]:
            cells.append('{} {} / {}'.format(c, num(cn['siMM_rep1'][c], 2), num(cn['siRE_rep1'][c], 2)))
        rows.append(r'\multicolumn{7}{@{}l}{' + ';\\quad '.join(cells) + r'} \\')
    table('table_s11.tex', 'llrrrrr',
          'Bin 1 (chrom:Mb) & Bin 2 (chrom:Mb) & Pairs & Expected & Obs/exp & Pairs (2nd) & Obs/exp (2nd)', rows,
          'Karyotype of Kasumi-1 in capture Hi-C. Top, the ten strongest inter-chromosomal 1\\,Mb bin pairs with at '
          'least 100 read pairs in the control library (siMM), one per 3\\,Mb region, with read pairs and observed '
          'over expected in the second library (siRE). Bottom, coverage of the other ends of captured '
          'inter-chromosomal pairs per chromosome relative to the autosomal median (siMM / siRE).', 'tab:s11',
          size=r'\scriptsize')


def s12():
    t = pd.read_csv(path('competitors', 'chaser', 'chaser_assortativity_GM12878.tsv'), sep='\t')
    marks = ['H3K27ac', 'H3K4me1', 'H3K27me3', 'H3K4me3', 'CTCF']
    labels = {'step 1': '1 step (direct)', 'step 1, distance-matched random': '1 step, distance-matched random$^a$',
              'step 2': '2 steps', 'step 3': '3 steps', 'step 4': '4 steps', 'step 5': '5 steps'}
    rows = []
    for network, label in labels.items():
        sub = t[t['network'] == network].set_index('mark')
        cells = [num(sub.loc[m, 'assortativity']) + (' ({})'.format(num(sub.loc[m, 'sd']))
                                                     if pd.notna(sub.loc[m, 'sd']) else '') for m in marks]
        rows.append(' & '.join([label] + cells) + r' \\')
    table('table_s12.tex', 'l' + 'r' * len(marks), 'Fragment pairs & ' + ' & '.join(marks), rows,
          'Chromatin assortativity (ChAseR; \\citealt{pancaldi2016,madridmencia2020}) of fragment pairs of the '
          'GM12878 graph by graph distance. For each distance from one to five steps, up to 50,000 cis pairs within '
          '2\\,Mb were sampled; features are the mean $\\log_2$ signal of each GM12878 mark over a fragment. '
          '$^a$Mean (standard deviation) of ten randomisations of the direct interactions that preserve the '
          'distance distribution.', 'tab:s12')


def s13():
    t = pd.read_csv(path('competitors', 'polymer', 'GM12878_polymer_vs_graph.tsv'), sep='\t')
    coords = {'TAD10': 'chr12:11.69--12.21\\,Mb', 'TAD19': 'chr18:60.68--61.12\\,Mb'}
    rows = []
    for _, r in t.iterrows():
        if r['viewpoint'].startswith('pooled'):
            rows.append(r'\midrule')
            rows.append(' & '.join(['Pooled$^a$', '', '', integer(r['pairs']), integer(r['hub_pairs']),
                                    num(r['spearman_polymer'], 2), num(r['hub_coefficient']), pval(r['hub_p'])])
                        + r' \\')
            continue
        domain, _, viewpoint = r['viewpoint'].split('.')[:3]
        rows.append(' & '.join(['{} ({})'.format(domain, coords[domain]), viewpoint[2:].replace('-', '--'),
                                integer(r['partner_bins']), integer(r['pairs']), integer(r['hub_pairs']),
                                num(r['spearman_polymer'], 2), num(r['hub_coefficient']), pval(r['hub_p'])]) + r' \\')
    table('table_s13.tex', 'llrrrrrr',
          'Domain & Viewpoint bins & Partner bins & Bin pairs & Hub pairs & $\\rho$ polymer & Hub coefficient & '
          '$P$ hub', rows,
          'Promoter hubs of the GM12878 graph against three-body contact probabilities predicted from Hi-C by the '
          'polymer model of \\citet{liu2021}, in the GM12878 SPRITE example of that study (5\\,kb bins, hg19). For '
          'each promoter viewpoint, bin pairs $(j,k)$ outside the viewpoint are compared between the predicted '
          'probability and the number of SPRITE clusters containing the viewpoint, $j$ and $k$. Hub pairs are pairs '
          'of bins that both contain direct graph partners of the viewpoint promoter. $\\rho$, Spearman correlation '
          'of prediction and observation; hub coefficient and $P$, from the linear model '
          '$\\log(\\mathrm{observed}+1)\\sim\\log(\\mathrm{predicted})+\\mathrm{hub}$, fitted for viewpoints with at '
          'least three hub pairs. $^a$All viewpoints, with the viewpoint as covariate.', 'tab:s13')


def s14():
    rows = []
    for label, name in (('hESC / hNEC', 'hESC_vs_hNEC'), ('Md / PA2', 'Md_vs_PA2'),
                        ('Kasumi-1 siMM / siRE', 'siMM_vs_siRE'), ('AML t(8;21) / FLT3-ITD', 't821_vs_FLT3ITD')):
        t = pd.read_csv(path('explain', name + '_edge_chromatin.tsv'), sep='\t')
        for _, r in t.iterrows():
            rows.append(' & '.join([label, tex(r['mark']) + (' (repressive)' if r['repressive'] else ''),
                                    r['interactions'], integer(r['n']), num(r['median_change'], 2),
                                    num(r['shared_median_change'], 2), num(r['concordance_vs_shared'], 3),
                                    'yes' if r['concordant'] else 'no',
                                    pval(r['mannwhitney_p_expected_direction'])]) + r' \\')
            label = ''
        rows.append(r'\midrule')
    table('table_s14.tex', 'lllrrrrlr',
          'Comparison & Mark & Interactions & n & Median change & Shared & Concordance & Expected & $P$',
          rows[:-1],
          'Interaction changes and chromatin changes at the interacting fragments (egaExplain). Median change: '
          '$\\log_2$ change of the mark (second over first condition), mean of the two ends of each interaction; '
          'concordance: probability that an interaction gained (lost) in the second condition has a larger mark '
          'change at its ends than a shared interaction of the same genomic distance (0.5 without relation); '
          'expected: the concordance lies on the side expected for the mark (active marks up at gained and down at '
          'lost interactions, repressive marks the opposite); $P$, one-sided Mann-Whitney test in the expected '
          'direction without distance matching.', 'tab:s14', size=r'\scriptsize')


def main():
    global R, OUT
    R, OUT = sys.argv[1:3]
    for build in (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14):
        build()
        print(build.__name__, flush=True)


if __name__ == '__main__':
    main()
