"""CHiCAGO design (rmap, baitmap) for a promoter capture of restriction fragments holding transcription start sites.

The rmap is the HiCUP digest (fragment ids in digest order); baited are the fragments containing a TSS of the given
GENCODE gene types, as in promoter capture designs that target fragments overlapping annotated promoters (e.g.
Ptasinska et al. 2019: protein-coding, non-coding, antisense, snRNA, miRNA and snoRNA transcripts).

Usage: python promoter_baitmap.py <HiCUP digest> <genes.gtf.gz> <output prefix> [gene types, comma separated]
"""

import sys

import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.io import read_tss_from_gtf

DEFAULT_TYPES = 'protein_coding,lincRNA,antisense,snRNA,miRNA,snoRNA'


def main():
    digest, gtf, prefix = sys.argv[1:4]
    types = (sys.argv[4] if len(sys.argv) > 4 else DEFAULT_TYPES).split(',')
    frags = pd.read_csv(digest, sep='\t', skiprows=2, header=None, usecols=[0, 1, 2],
                        names=['chrom', 'start', 'end'])
    frags = frags[frags['chrom'].str.match(r'^chr([0-9]+|X|Y)$')].reset_index(drop=True)
    frags['id'] = np.arange(1, len(frags) + 1)
    frags.to_csv(prefix + '.rmap', sep='\t', header=False, index=False)
    tss = read_tss_from_gtf(gtf, types)
    names = {}
    for chrom, sites in tss.groupby('chrom'):
        f = frags[frags['chrom'] == chrom]
        if len(f) == 0:
            continue
        index = np.searchsorted(f['start'].to_numpy(), sites['tss'].to_numpy() + 1, side='right') - 1
        ok = (index >= 0) & (sites['tss'].to_numpy() + 1 <= f['end'].to_numpy()[np.clip(index, 0, len(f) - 1)])
        for i, gene in zip(f['id'].to_numpy()[index[ok]], sites['gene_name'].to_numpy()[ok]):
            names.setdefault(i, set()).add(gene)
    baits = frags[frags['id'].isin(names)].copy()
    baits['name'] = [','.join(sorted(names[i])) for i in baits['id']]
    baits.to_csv(prefix + '.baitmap', sep='\t', header=False, index=False)
    baits.iloc[:, :4].to_csv(prefix + '.baitmap_4col.txt', sep='\t', header=False, index=False)
    print('{} fragments, {} baited ({} TSS of {})'.format(len(frags), len(baits), len(tss), ','.join(types)))


if __name__ == '__main__':
    main()
