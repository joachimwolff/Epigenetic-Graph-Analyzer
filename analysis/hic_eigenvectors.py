"""Hi-C compartment eigenvectors with the phasing rule of cooltools (reference for the compartment evaluation).

hicPCA (18_hic_reference.sh) returns the first eigenvector, which on several chromosomes follows the chromosome
arms or another leading pattern instead of the compartments (e.g. H1 chr9, chr13, chr21: r -0.1 to -0.9 with the
published 4DN eigenvector). Here, per chromosome: balanced cis matrix of the mcool (cooler 'weight'), observed /
expected over valid bins, Pearson correlation, the three leading eigenvectors; the first one is kept unless another
correlates at least 1.5 times more strongly with GC content (absolute r), and the vector is oriented positively with
GC (cooltools convention), via epigeneticgraphanalyzer.compartments.contact_eigenvector. (Taking the plain maximum
picked a wrong eigenvector on 6 GM12878 chromosomes, checked against the Rao et al. 2014 subcompartments.)

Needs cooler (run with the python of the HiCExplorer environment).
Usage: python hic_eigenvectors.py [resolutions ...]   (default 100000 250000)
Output: <results>/compartments_domains/hic_reference/<cell>_eig_<res>.bedgraph and hic_eigenvector_choice.tsv
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compartments_common as cc  # noqa: E402
from epigeneticgraphanalyzer.compartments import contact_eigenvector  # noqa: E402

MCOOL_DIR = os.environ.get('EGA_MCOOL', os.path.expanduser('~/data/ega_virtual_capture/mcool'))


def main():
    import cooler
    resolutions = [int(r) for r in sys.argv[1:]] or [100000, 250000]
    choices = []
    for res in resolutions:
        bins = cc.bins_for('hg38', res)
        gc = cc.gc_content(bins, 'hg38')
        for cell, accession in cc.MCOOL.items():
            clr = cooler.Cooler('{}/{}.mcool::/resolutions/{}'.format(MCOOL_DIR, accession, res))
            out = []
            for chrom in cc.CHROMS:
                m = clr.matrix(balance=True).fetch(chrom)
                weights = clr.bins().fetch(chrom)['weight'].values
                valid = np.isfinite(weights)
                idx = np.flatnonzero((bins.chrom == chrom).values)
                n = min(len(idx), m.shape[0])
                vec, index = contact_eigenvector(np.nan_to_num(m[:n, :n]), valid=valid[:n], orient_by=gc[idx[:n]],
                                                 min_contacts=1, n_candidates=3)
                b = bins.iloc[idx[:n]]
                out.append(pd.DataFrame({'chrom': chrom, 'start': b.start.values, 'end': b.end.values, 'value': vec}))
                choices.append({'cell': cell, 'resolution': res, 'chrom': chrom, 'eigenvector': index + 1})
                print(cell, res, chrom, 'eigenvector', index + 1, flush=True)
            table = pd.concat(out).dropna()
            table.to_csv(os.path.join(cc.REF, '{}_eig_{}kb.bedgraph'.format(cell, res // 1000)), sep='\t',
                         header=False, index=False, float_format='%.6g')
    pd.DataFrame(choices).to_csv(os.path.join(cc.REF, 'hic_eigenvector_choice.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
