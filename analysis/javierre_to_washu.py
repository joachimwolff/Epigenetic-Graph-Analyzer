"""One cell type of the Javierre et al. 2016 promoter capture Hi-C peak matrix as WashU text for EGA.

The peak matrix (PCHiC_peak_matrix_cutoff5, hg19, HindIII) holds every interaction with a CHiCAGO score >= 5 in at
least one of 17 blood cell types and the scores of all cell types. All cis rows with a score of the chosen cell type
are written ('chrom,start,end<TAB>chrom,start,end<TAB>score', bait first), so that scores below the calling threshold
remain available; egaDifferential --scoreThreshold selects the interactions.

Usage: python javierre_to_washu.py <PCHiC_peak_matrix_cutoff5.txt.gz> <cell type column, e.g. nB> <output.txt.gz>
"""

import sys

import pandas as pd


def main():
    matrix, column, output = sys.argv[1:4]
    table = pd.read_csv(matrix, sep='\t', usecols=['baitChr', 'baitStart', 'baitEnd', 'oeChr', 'oeStart', 'oeEnd',
                                                   column], dtype={'baitChr': str, 'oeChr': str})
    table = table[(table['baitChr'] == table['oeChr']) & (table[column] > 0)]
    out = pd.DataFrame({
        'bait': 'chr' + table['baitChr'] + ',' + table['baitStart'].astype(str) + ',' + table['baitEnd'].astype(str),
        'other': 'chr' + table['oeChr'] + ',' + table['oeStart'].astype(str) + ',' + table['oeEnd'].astype(str),
        'score': table[column].round(2)})
    out.to_csv(output, sep='\t', header=False, index=False, compression='gzip')
    print('{}: {} cis interactions, {} with score >= 5'.format(column, len(out), int((out['score'] >= 5).sum())))


if __name__ == '__main__':
    main()
