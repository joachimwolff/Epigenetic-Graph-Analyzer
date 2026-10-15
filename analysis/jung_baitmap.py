"""CHiCAGO bait map of the promoter capture design of Jung et al. 2019 on the hg19 HindIII digest.

Jung et al. captured +-300 bp around the two HindIII sites closest to each promoter and report interactions
for fragment units that span one or more HindIII sites (midpoint to midpoint of the neighbouring sites). The
baited HindIII fragments are the two digest fragments flanking every HindIII site inside a baited unit; baited
units are the first fragment of the promoter-other tables and both fragments of the promoter-promoter tables.

Usage: python jung_baitmap.py <digest.rmap> <GSE86189_X.po.all.txt.bz2> <GSE86189_X.pp.all.txt.bz2> <out.baitmap>
"""

import sys

import numpy as np
import pandas as pd


def units(path, columns):
    ids = pd.concat([pd.read_csv(path, sep='\t', usecols=[c])[c] for c in columns]).drop_duplicates()
    parts = ids.str.split('.', expand=True)
    return pd.DataFrame({'chrom': parts[0], 'start': parts[1].astype(np.int64), 'end': parts[2].astype(np.int64)})


def main():
    rmap_path, po, pp, output = sys.argv[1:5]
    rmap = pd.read_csv(rmap_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'id'])
    baited = pd.concat([units(po, ['frag1']), units(pp, ['frag1', 'frag2'])]).drop_duplicates()
    selected = set()
    for chrom, group in baited.groupby('chrom'):
        frags = rmap[rmap['chrom'] == chrom].sort_values('start')
        # HindIII sites lie between consecutive digest fragments: fragment k ends at site k
        site = frags['end'].to_numpy()[:-1]
        left_id = frags['id'].to_numpy()[:-1]
        right_id = frags['id'].to_numpy()[1:]
        lo = np.searchsorted(site, group['start'].to_numpy(), side='left')
        hi = np.searchsorted(site, group['end'].to_numpy(), side='right')
        for a, b in zip(lo, hi):
            selected.update(left_id[a:b])
            selected.update(right_id[a:b])
    baitmap = rmap[rmap['id'].isin(selected)].copy()
    baitmap['name'] = 'jung2019_bait'
    baitmap.to_csv(output, sep='\t', header=False, index=False)
    print('{} baited units -> {} baited HindIII fragments'.format(len(baited), len(baitmap)))


if __name__ == '__main__':
    main()
