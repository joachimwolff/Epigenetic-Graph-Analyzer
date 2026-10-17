"""Shared inputs of the compartment and domain analyses (compartments_*.py, domains_graph.py, compartment_switch.py).

Paths follow the layout of the analysis machine; override the roots with EGA_RAW, EGA_RESULTS and EGA_HOPS.
"""

import gzip
import os
import sys

import numpy as np
import pandas as pd

# the package of this checkout before any installed copy
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer.compartments import interval_density, make_bins
from epigeneticgraphanalyzer.genome import read_chrom_sizes

RAW = os.environ.get('EGA_RAW', '/mnt/processed/wolffjoa/ega/raw')
RESULTS = os.environ.get('EGA_RESULTS', '/mnt/processed/wolffjoa/ega/results')
OUT = os.path.join(RESULTS, 'compartments_domains')
# hops_n tables (cascade); a fixed snapshot while the differential outputs are being regenerated
HOPS = os.environ.get('EGA_HOPS', os.path.join(RESULTS, 'hops_n', 'before_elements'))
REF = os.path.join(OUT, 'hic_reference')
TABLES = os.path.join(OUT, 'tables')
FIGURES = os.path.join(OUT, 'figures')
HG38_FASTA = os.environ.get('EGA_HG38_FASTA', os.path.expanduser('~/data/hic-multiome/genome/hg38.fa'))
HG19_FASTA = os.path.join(RAW, 'annotation', 'hg19.fa.gz')
CHROMS = ['chr{}'.format(c) for c in list(range(1, 23)) + ['X']]

# ENCODE tracks (GRCh38) of the hic-multiome project: same assembly as the 4DN mcool files
ENCODE_HG38 = {
    'GM12878': os.path.expanduser('~/data/hic-multiome/hg38_gm12878'),
    'K562': '/mnt/processed/wolffjoa/hic-multiome/hg38_k562',
    'H1': '/mnt/processed/wolffjoa/hic-multiome/hg38_h1',
}
ROLES = {
    'H3K27ac': 'active', 'H3K9ac': 'active', 'H3K4me3': 'active', 'H3K4me2': 'active', 'H3K4me1': 'enhancer',
    'H3K36me3': 'transcription', 'H3K79me2': 'transcription', 'H4K20me1': 'transcription',
    'H3K27me3': 'repressive', 'H3K9me3': 'heterochromatin', 'CTCF': 'other', 'RAD21': 'other', 'SMC3': 'other',
    'H3K27AC': 'active', 'H3K4ME1': 'enhancer', 'H3K4ME3': 'active', 'H3K27ME3': 'repressive',
}
MCOOL = {'GM12878': '4DNFIIV4M7TF', 'K562': '4DNFI18UHVRO', 'H1': '4DNFI9GMP2J8'}


def chrom_sizes(assembly):
    sizes = read_chrom_sizes(os.path.join(RAW, 'annotation', '{}.chrom.sizes'.format(assembly)))
    return {c: sizes[c] for c in CHROMS if c in sizes}


def bins_for(assembly, bin_size):
    return make_bins(chrom_sizes(assembly), bin_size)


def read_bedgraph_on_bins(path, bins):
    """Values of a bedGraph (or 4DN bigWig) at the bin midpoints; NaN where missing; exact 0 kept as 0."""
    out = np.full(len(bins), np.nan)
    mids = ((bins.start.values + bins.end.values) // 2).astype(np.int64)
    if path.endswith(('.bw', '.bigWig', '.bigwig')):
        import pyBigWig
        handle = pyBigWig.open(path)
        sizes = handle.chroms()
        for i, (c, m) in enumerate(zip(bins.chrom.values, mids)):
            if c in sizes and m < sizes[c]:
                v = handle.values(c, int(m), int(m) + 1)[0]
                out[i] = v
        return out
    table = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'value'])
    table = table[~table.chrom.astype(str).str.startswith(('track', 'browser', '#'))]
    table = table.astype({'start': np.int64, 'end': np.int64, 'value': float})
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        sub = table[table.chrom == chrom].sort_values('start')
        if sub.empty:
            continue
        m = mids[idx]
        k = np.searchsorted(sub.start.values, m, side='right') - 1
        ok = (k >= 0) & (m < sub.end.values[np.clip(k, 0, None)])
        vals = np.full(len(m), np.nan)
        vals[ok] = sub.value.values[k[ok]]
        out[idx] = vals
    return out


def reference_on_positions(path, chroms, positions, bin_size=None):
    """Reference eigenvector value at arbitrary positions (e.g. lifted hg19 bin midpoints)."""
    fake = pd.DataFrame({'chrom': chroms, 'start': np.asarray(positions), 'end': np.asarray(positions) + 1})
    return read_bedgraph_on_bins(path, fake)


def lift_bin_midpoints(bins, chain):
    """hg19 -> hg38 (or back) midpoint liftover of bins; returns (chrom, position) arrays (None / -1 when unmapped
    or on another chromosome)."""
    from epigeneticgraphanalyzer.genome import LiftOver
    lift = LiftOver(chain)
    chroms = np.empty(len(bins), dtype=object)
    pos = np.full(len(bins), -1, dtype=np.int64)
    mids = (bins.start.values + bins.end.values) // 2
    for i, (c, m) in enumerate(zip(bins.chrom.values, mids)):
        hit = lift.position(c, int(m))
        if hit is not None and hit[0] == c:
            chroms[i], pos[i] = hit
    return chroms, pos


def gene_density(bins, assembly):
    genes = pd.read_csv(os.path.join(REF, 'genes_{}_protein_coding.bed'.format(assembly)), sep='\t', header=None)
    tss = np.where(genes[5] == '+', genes[1], genes[2] - 1)
    return interval_density(bins, genes[0].values, tss)


def gc_content(bins, assembly):
    """GC fraction per bin from the genome FASTA (cached as a table next to the Hi-C reference)."""
    bin_size = int((bins.end - bins.start).max())
    cache = os.path.join(REF, 'gc_{}_{}.tsv.gz'.format(assembly, bin_size))
    if os.path.exists(cache):
        table = pd.read_csv(cache, sep='\t')
        return bins.merge(table, on=['chrom', 'start'], how='left')['gc'].values
    fasta = HG38_FASTA if assembly == 'hg38' else HG19_FASTA
    opener = gzip.open if fasta.endswith('.gz') else open
    wanted = set(bins.chrom)
    rows = []

    def flush(chrom, chunks):
        if chrom not in wanted:
            return
        seq = ''.join(chunks).upper()
        for start in range(0, len(seq), bin_size):
            s = seq[start:start + bin_size]
            gc = s.count('G') + s.count('C')
            at = s.count('A') + s.count('T')
            rows.append((chrom, start, gc / (gc + at) if gc + at > 0.5 * len(s) else np.nan))

    chrom, chunks = None, []
    with opener(fasta, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if chrom is not None:
                    flush(chrom, chunks)
                chrom, chunks = line[1:].split()[0], []
            elif chrom in wanted:
                chunks.append(line.strip())
        if chrom is not None:
            flush(chrom, chunks)
    table = pd.DataFrame(rows, columns=['chrom', 'start', 'gc'])
    table.to_csv(cache, sep='\t', index=False)
    return bins.merge(table, on=['chrom', 'start'], how='left')['gc'].values
