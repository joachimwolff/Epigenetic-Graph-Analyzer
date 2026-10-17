"""eQTL input tables of egaEqtl in hg19 and hg38: the tested pairs (leads) and all significant pairs per data set.

  eqtlgen         eQTLGen phase I whole-blood cis-eQTLs (FDR < 0.05, hg19); lead = smallest P per gene
  gtex_lcl        GTEx v8 Cells_EBV-transformed_lymphocytes signif_variant_gene_pairs (GRCh38); lead per gene
  gtex_wb         GTEx v8 Whole_Blood signif_variant_gene_pairs (GRCh38); lead per gene
  geuvadis_susie  eQTL Catalogue GEUVADIS LCL SuSiE credible sets (QTD000110, GRCh38); leads = variants with
                  PIP >= 0.5, significant = all credible-set variants

The other assembly is reached with the UCSC chains (pyliftover, rows that do not lift or change chromosome are
dropped). Writes <out>/<name>.<assembly>.{leads,significant}.tsv.gz and <out>/prepare_summary.tsv.

Usage: python eqtl_prepare.py <raw eqtl dir> <annotation dir> <output dir>
"""

import os
import sys

import pandas as pd

from epigeneticgraphanalyzer.eqtl import (fine_mapped, lead_per_gene, lift_positions, read_eqtlgen, read_gtex_pairs,
                                          read_susie)


def main():
    raw, annotation, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    chains = {('hg19', 'hg38'): os.path.join(annotation, 'hg19ToHg38.over.chain.gz'),
              ('hg38', 'hg19'): os.path.join(annotation, 'hg38ToHg19.over.chain.gz')}
    gtex = os.path.join(raw, 'gtex_v8', 'GTEx_Analysis_v8_eQTL.tar')
    sources = {
        'eqtlgen': ('hg19', lambda: read_eqtlgen(os.path.join(
            raw, 'eqtlgen', '2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'))),
        'gtex_lcl': ('hg38', lambda: read_gtex_pairs(gtex, 'Cells_EBV-transformed_lymphocytes')),
        'gtex_wb': ('hg38', lambda: read_gtex_pairs(gtex, 'Whole_Blood')),
        'geuvadis_susie': ('hg38', lambda: read_susie(os.path.join(raw, 'eqtlcatalogue_geuvadis_lcl',
                                                                   'QTD000110.credible_sets.tsv.gz'))),
    }
    summary = []
    for name, (assembly, reader) in sources.items():
        if len(sys.argv) > 4 and name not in sys.argv[4].split(','):
            continue
        df = reader()
        leads = fine_mapped(df) if name == 'geuvadis_susie' else lead_per_gene(df)
        significant = df[['chrom', 'pos', 'gene_id']].drop_duplicates()
        tables = {assembly: (leads, significant)}
        other = 'hg19' if assembly == 'hg38' else 'hg38'
        tables[other] = tuple(lift_positions(t, chains[(assembly, other)]) for t in (leads, significant))
        for build, (lead_table, sig_table) in tables.items():
            columns = [c for c in ('chrom', 'pos', 'gene_id', 'pvalue', 'pip') if c in lead_table]
            lead_table[columns].to_csv(os.path.join(out, '{}.{}.leads.tsv.gz'.format(name, build)), sep='\t',
                                       index=False)
            sig_table[['chrom', 'pos', 'gene_id']].to_csv(
                os.path.join(out, '{}.{}.significant.tsv.gz'.format(name, build)), sep='\t', index=False)
            summary.append({'dataset': name, 'assembly': build, 'native': build == assembly,
                            'input_rows': len(df), 'leads': len(lead_table),
                            'lead_genes': lead_table['gene_id'].nunique(), 'significant_pairs': len(sig_table),
                            'significant_genes': sig_table['gene_id'].nunique()})
        print(name, 'done', flush=True)
    path = os.path.join(out, 'prepare_summary.tsv')
    previous = pd.read_csv(path, sep='\t') if os.path.exists(path) else pd.DataFrame()
    table = pd.DataFrame(summary)
    if len(previous):
        table = pd.concat([previous[~previous['dataset'].isin(table['dataset'])], table], ignore_index=True)
    table.to_csv(path, sep='\t', index=False)


if __name__ == '__main__':
    main()
