"""Virtual promoter capture: CHiCAGO input from a Hi-C or Micro-C contact matrix (.mcool).

The bins of one resolution act as restriction fragments (rmap); bins holding a transcription start site of a
selected gene type (GENCODE GTF) are the baits. Every pixel with a baited bin becomes one chinput row with the raw
read count (bait-to-bait pixels once per bait, as bam2chicago writes them); trans pixels keep distSign NA, so
CHiCAGO estimates its technical noise as for real capture data. The design (.rmap, .baitmap) and the chinput are
written to the output directory; the .npb/.nbpb/.poe files are made by makeDesignFiles_py3.py of chicagoTools.

Usage: python virtual_capture.py <file.mcool::resolutions/5000> <genes.gtf.gz> <output dir> <name>
       [gene types, comma separated, default protein_coding]
"""

import os
import sys

import h5py
import numpy as np
import pandas as pd

from epigeneticgraphanalyzer.io import read_tss_from_gtf

CHUNK = 50000000


def main():
    uri, gtf, out, name = sys.argv[1:5]
    gene_types = sys.argv[5].split(',') if len(sys.argv) > 5 else ['protein_coding']
    path, group_path = uri.split('::', 1) if '::' in uri else (uri, '/')
    os.makedirs(out, exist_ok=True)
    with h5py.File(path, 'r') as handle:
        group = handle[group_path]
        chroms = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in group['chroms/name'][:]])
        bin_chrom = chroms[group['bins/chrom'][:]]
        bin_start = group['bins/start'][:].astype(np.int64)
        bin_end = group['bins/end'][:].astype(np.int64)
        keep_chrom = np.array([c in {'chr{}'.format(i) for i in list(range(1, 23)) + ['X']} for c in chroms])
        usable = keep_chrom[group['bins/chrom'][:]]
        # rmap: 1-based closed coordinates, ids are the bin index + 1
        ids = np.arange(1, len(bin_start) + 1)
        rmap = pd.DataFrame({'chrom': bin_chrom, 'start': bin_start + 1, 'end': bin_end, 'id': ids})[usable]
        rmap.to_csv(os.path.join(out, name + '.rmap'), sep='\t', header=False, index=False)

        tss = read_tss_from_gtf(gtf, gene_types)
        is_bait = np.zeros(len(bin_start), dtype=bool)
        names = {}
        offsets = {c: np.nonzero(bin_chrom == c)[0] for c in chroms[keep_chrom]}
        for chrom, sites in tss.groupby('chrom'):
            if chrom not in offsets:
                continue
            index = offsets[chrom]
            hit = np.searchsorted(bin_start[index], sites['tss'].to_numpy(), side='right') - 1
            for b, gene in zip(index[np.clip(hit, 0, len(index) - 1)], sites['gene_name']):
                is_bait[b] = True
                names.setdefault(b, set()).add(gene)
        bait_index = np.nonzero(is_bait & usable)[0]
        baitmap = pd.DataFrame({'chrom': bin_chrom[bait_index], 'start': bin_start[bait_index] + 1,
                                'end': bin_end[bait_index], 'id': ids[bait_index],
                                'name': [','.join(sorted(names[b])) for b in bait_index]})
        baitmap.to_csv(os.path.join(out, name + '.baitmap'), sep='\t', header=False, index=False)
        baitmap.iloc[:, :4].to_csv(os.path.join(out, name + '.baitmap_4col.txt'), sep='\t', header=False,
                                   index=False)
        print('{} bins, {} baited bins ({} genes)'.format(len(rmap), len(bait_index),
                                                         tss['gene_name'].nunique()), flush=True)

        middle = (bin_start + bin_end) // 2
        length = bin_end - bin_start
        bait_mask = is_bait & usable
        pixels = group['pixels']
        n = pixels['bin1_id'].shape[0]
        chinput = os.path.join(out, name + '.chinput')
        rows = 0
        with open(chinput, 'w') as handle_out:
            handle_out.write('# samplename={}\tbamname={}\tbaitmapfile={}.baitmap\tdigestfile={}.rmap\n'.format(
                name, os.path.basename(path), name, name))
            handle_out.write('baitID\totherEndID\tN\totherEndLen\tdistSign\n')
            for start in range(0, n, CHUNK):
                b1 = pixels['bin1_id'][start:start + CHUNK]
                b2 = pixels['bin2_id'][start:start + CHUNK]
                count = pixels['count'][start:start + CHUNK]
                ok = (b1 != b2) & usable[b1] & usable[b2] & (count > 0)
                for bait, other in ((b1, b2), (b2, b1)):
                    sel = ok & bait_mask[bait]
                    if not sel.any():
                        continue
                    bb, oo, cc = bait[sel], other[sel], count[sel]
                    cis = bin_chrom[bb] == bin_chrom[oo]
                    dist = np.where(cis, middle[oo] - middle[bb], 0)
                    frame = pd.DataFrame({'baitID': ids[bb], 'otherEndID': ids[oo], 'N': cc,
                                          'otherEndLen': length[oo], 'distSign': dist})
                    frame['distSign'] = frame['distSign'].astype('Int64').mask(~cis)
                    frame.to_csv(handle_out, sep='\t', header=False, index=False, na_rep='NA')
                    rows += len(frame)
                print('{:.0%} of pixels, {} chinput rows'.format(min(1.0, (start + CHUNK) / n), rows), flush=True)
    print('done', chinput, rows)


if __name__ == '__main__':
    main()
