"""Capture Hi-C contacts of Rubin et al. 2017 (GSE84660, keratinocyte differentiation, hg19) as one ibed per donor
and day for egaCreateGraph / egaDifferential.

Usage: python rubin_chic_to_ibed.py <GSE84660_CHiC_contacts_b2b.bed.gz> <GSE84660_CHiC_contacts_b2g.bed.gz> <out_dir>

The deposited tables list the contacts called by CHiCAGO (FDR < 0.01) at any of day 0, 3 and 6 with the raw read
count of every donor and day; the b2g table (bait to non-bait) is in addition restricted to other ends overlapping an
H3K27ac peak at any day. No per-day CHiCAGO score is deposited, so the score of a contact in one sample is derived
from its read count:

    score = 5 * (count / size factor) / expected(distance)

size factor: median-of-ratios of the sample over all contacts (b2b and b2g); expected(distance): median depth-
normalised count over the six samples of the contacts in the same distance bin (50 quantile bins of log distance,
cis contacts). A score of 5 (the CHiCAGO calling threshold used elsewhere, --scoreThreshold 5 --referenceScore 5)
is the distance-matched median of the deposited significant contacts. A contact is listed in a sample when it has at
least one read there; with --replicates all a day keeps the contacts with reads in both donors and a score of at
least 5 in at least one donor.

Outputs: <out_dir>/D{0,3,6}_B{1,2}.ibed.gz, baits.bed, conversion_summary.tsv.
"""

import os
import sys

import numpy as np
import pandas as pd

SAMPLES = ['B1_D0', 'B1_D3', 'B1_D6', 'B2_D0', 'B2_D3', 'B2_D6']
COLUMNS = ['bait_chr', 'bait_start', 'bait_end', 'other_chr', 'other_start', 'other_end', 'contact_id', 'bait_id',
           'other_id'] + SAMPLES


def read_table(path, kind):
    table = pd.read_csv(path, sep='\t', comment='#', header=None, names=COLUMNS)
    table['kind'] = kind
    return table


def size_factors(counts):
    """DESeq2 median-of-ratios size factors over rows without zeros."""
    positive = (counts > 0).all(axis=1)
    logs = np.log(counts[positive])
    reference = logs.mean(axis=1)
    return np.exp((logs.sub(reference, axis=0)).median(axis=0))


def main():
    b2b_path, b2g_path, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    table = pd.concat([read_table(b2b_path, 'b2b'), read_table(b2g_path, 'b2g')], ignore_index=True)
    total = len(table)
    table = table[table['bait_chr'] == table['other_chr']].reset_index(drop=True)
    counts = table[SAMPLES].astype(float)
    factors = size_factors(counts)
    normalised = counts / factors
    distance = ((table['other_start'] + table['other_end']) // 2 - (table['bait_start'] + table['bait_end']) // 2).abs()
    bins = pd.qcut(np.log10(distance.clip(lower=1)), 50, labels=False, duplicates='drop')
    expected = normalised.mean(axis=1).groupby(bins).median()
    scores = 5.0 * normalised.div(expected.reindex(bins).to_numpy(), axis=0)
    baits = pd.concat([table[['bait_chr', 'bait_start', 'bait_end']].set_axis(['chrom', 'start', 'end'], axis=1),
                       table.loc[table['kind'] == 'b2b', ['other_chr', 'other_start', 'other_end']]
                       .set_axis(['chrom', 'start', 'end'], axis=1)]).drop_duplicates().sort_values(['chrom', 'start'])
    baits.to_csv(os.path.join(out_dir, 'baits.bed'), sep='\t', header=False, index=False)
    bait_name = table['bait_chr'] + ':' + table['bait_start'].astype(str)
    other_name = np.where(table['kind'] == 'b2b', table['other_chr'] + ':' + table['other_start'].astype(str), '.')
    rows = []
    for sample in SAMPLES:
        donor, day = sample.split('_')
        listed = counts[sample] > 0
        ibed = pd.DataFrame({
            'bait_chr': table['bait_chr'], 'bait_start': table['bait_start'], 'bait_end': table['bait_end'],
            'bait_name': bait_name, 'otherEnd_chr': table['other_chr'], 'otherEnd_start': table['other_start'],
            'otherEnd_end': table['other_end'], 'otherEnd_name': other_name,
            'N_reads': counts[sample].astype(int), 'score': scores[sample].round(3)})[listed]
        path = os.path.join(out_dir, '{}_{}.ibed.gz'.format(day, donor))
        ibed.to_csv(path, sep='\t', index=False)
        rows.append({'sample': sample, 'size_factor': factors[sample], 'contacts_total': total,
                     'contacts_cis': len(table), 'listed': int(listed.sum()),
                     'score_ge_5': int((listed & (scores[sample] >= 5)).sum())})
    for day in ('D0', 'D3', 'D6'):
        both = (counts['B1_' + day] > 0) & (counts['B2_' + day] > 0)
        called = (scores['B1_' + day] >= 5) | (scores['B2_' + day] >= 5)
        rows.append({'sample': day + '_rule_all', 'listed': int(both.sum()), 'score_ge_5': int((both & called).sum()),
                     'contacts_total': total, 'contacts_cis': len(table),
                     'b2b': int((both & called & (table['kind'] == 'b2b')).sum()),
                     'b2g': int((both & called & (table['kind'] == 'b2g')).sum())})
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(out_dir, 'conversion_summary.tsv'), sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
