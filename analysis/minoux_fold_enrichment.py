"""Fold enrichment over input for the ChIP-seq tracks of Minoux et al. 2017 (GEO GSE89435).

The GEO wig files hold alignments per 100 bp window per million reads (QuasR qExportWig, scaling 1e6). The inputs
(three libraries of the mandibular arch, E10.5) are averaged. Fold enrichment in windows of `--window` bp is the ChIP
signal divided by a local input level, the larger of the input mean in `--background` bp around the window and the
genome-wide input mean (as the local lambda of MACS), so that windows with little input do not inflate the ratio.
The inputs serve both arches: they reflect mappability, GC content and copy number, shared by the two tissues.

Usage: python minoux_fold_enrichment.py <chrom sizes or reference.bw> <output.bw> <chip.wig.gz>
                                        --input <input.wig.gz> [<input.wig.gz> ...] [--window 1000] [--background 10000]
"""

import argparse
import gzip

import numpy as np
import pyBigWig

STEP = 100


def read_wig(path, sizes):
    """{chrom: per-100-bp values} of a fixedStep wig with step = span = 100."""
    tracks = {c: np.zeros((n + STEP - 1) // STEP, dtype=np.float32) for c, n in sizes.items()}
    chrom, position, values = None, 0, None
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('track'):
                continue
            if line.startswith('fixedStep'):
                fields = dict(f.split('=') for f in line.split()[1:])
                if int(fields.get('step', STEP)) != STEP:
                    raise ValueError('{}: step {} not supported'.format(path, fields['step']))
                chrom = fields['chrom']
                position = (int(fields['start']) - 1) // STEP
                values = tracks.get(chrom)
                continue
            if values is not None and position < len(values):
                values[position] = float(line)
            position += 1
    return tracks


def chrom_sizes(reference):
    if reference.endswith(('.bw', '.bigwig', '.bigWig')):
        handle = pyBigWig.open(reference)
        sizes = dict(handle.chroms())
        handle.close()
    else:
        sizes = {c: int(n) for c, n in (line.split()[:2] for line in open(reference) if line.strip())}
    return {c: n for c, n in sizes.items() if '_' not in c}


def window_mean(values, width):
    """Mean of consecutive blocks of `width` steps (last block partial)."""
    n = (len(values) + width - 1) // width
    padded = np.zeros(n * width, dtype=np.float64)
    padded[:len(values)] = values
    counts = np.full(n, width, dtype=np.float64)
    counts[-1] = len(values) - (n - 1) * width
    return padded.reshape(n, width).sum(axis=1) / counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('sizes')
    parser.add_argument('output')
    parser.add_argument('chip')
    parser.add_argument('--input', nargs='+', required=True)
    parser.add_argument('--window', type=int, default=1000)
    parser.add_argument('--background', type=int, default=10000)
    args = parser.parse_args()
    sizes = chrom_sizes(args.sizes)
    chip = read_wig(args.chip, sizes)
    inputs = [read_wig(p, sizes) for p in args.input]
    w = args.window // STEP
    b = args.background // args.window
    control = {c: np.mean([track[c] for track in inputs], axis=0) for c in sizes}
    genome_mean = (sum(float(control[c].sum()) for c in sizes) / sum(len(control[c]) for c in sizes)) * w
    out = pyBigWig.open(args.output, 'w')
    out.addHeader(sorted(sizes.items()))
    for chrom in sorted(sizes):
        signal = window_mean(chip[chrom], w) * w
        local = window_mean(control[chrom], w) * w
        # local input level: mean over the `b` windows centred on each window
        kernel = np.ones(b) / b
        local = np.convolve(local, kernel, mode='same')
        fold = signal / np.maximum(local, genome_mean)
        starts = np.arange(len(fold), dtype=np.int64) * args.window
        ends = np.minimum(starts + args.window, sizes[chrom])
        out.addEntries([chrom] * len(fold), starts.tolist(), ends=ends.tolist(), values=fold.astype(float).tolist())
    out.close()
    print('{}: genome input level {:.4g} per {} bp window'.format(args.output, genome_mean, args.window))


if __name__ == '__main__':
    main()
