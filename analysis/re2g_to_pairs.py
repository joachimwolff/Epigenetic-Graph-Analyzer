"""ENCODE-rE2G predictions as element-gene pairs in the layout of the eQTL inputs (chrom, pos, gene_id), so that
44_explaining_elements_eqtl.py tests whether an explaining element holds an element that ENCODE-rE2G links to the same
gene. pos is the midpoint of the predicted element, lifted with a UCSC chain (rows that do not lift or change
chromosome are dropped); self-promoter rows are dropped.

Usage: python re2g_to_pairs.py <rE2G thresholded bed.gz> <chain> <out.tsv.gz> [--minScore 0]
"""

import argparse

import pandas as pd
from pyliftover import LiftOver


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('predictions')
    parser.add_argument('chain')
    parser.add_argument('output')
    parser.add_argument('--minScore', type=float, default=0.0)
    args = parser.parse_args()
    table = pd.read_csv(args.predictions, sep='\t', usecols=['#chr', 'start', 'end', 'TargetGeneEnsemblID',
                                                             'isSelfPromoter', 'Score'])
    table = table[(table.Score >= args.minScore) & ~table.isSelfPromoter.astype(str).str.lower().eq('true')]
    lift = LiftOver(args.chain)
    rows = []
    for chrom, start, end, gene in zip(table['#chr'], table.start, table.end, table.TargetGeneEnsemblID):
        hit = lift.convert_coordinate(chrom, (start + end) // 2)
        if hit and hit[0][0] == chrom:
            rows.append((chrom, int(hit[0][1]), str(gene).split('.')[0]))
    out = pd.DataFrame(rows, columns=['chrom', 'pos', 'gene_id'])
    out.to_csv(args.output, sep='\t', index=False)
    print('{} predictions, {} lifted, {} genes'.format(len(table), len(out), out.gene_id.nunique()))


if __name__ == '__main__':
    main()
