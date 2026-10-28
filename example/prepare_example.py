"""Cut the example data set (example/data/) out of the full Md / PA2 inputs: the Hoxa locus of mouse chromosome 6
(mm10), chr6:49,500,000-54,500,000, which holds the Hoxa cluster, HIRE1 and HIRE2 and the domain boundary crossed by
the gained Hoxa2 contacts.

Inputs (public data, see REPRODUCE.md and scripts/download_data.sh):
  <data_dir>/kessler/GSM650520[89]_*.txt.gz, GSM650521[01]_*.txt.gz  CHiCAGO calls of Kessler et al. 2023 (GSE211902),
                                                                    WashU text, two replicates per tissue
  <data_dir>/kessler/E10.5_{MD,PA2}_{H3K27ac,H3K4me2,H3K27me3}_R1.bigwig  ChIP-seq of Minoux et al. 2017 (GSE89435)
  <data_dir>/kessler/mm10.ncbiRefSeq.gtf                            UCSC RefSeq annotation of mm10
  <results_dir>/minoux_fold/{Md,PA2}_{mark}.fold.bw                  fold enrichment over the pooled input (17_minoux_fold.sh)
  <results_dir>/kessler/deseq2_PA2_vs_Md.tsv                         PyDESeq2 of the RNA-seq counts of GSE89434
  <results_dir>/tad_crossing/hicfindtads/PA2_25kb_*                  hicFindTADs on the PA2 Hi-C (28_tad_crossing.sh)
  <results_dir>/tad_crossing/kessler2023_HIRE_mm10.bed               HIRE1 / HIRE2 of Kessler et al. 2023

Interactions are kept when both ends lie in the window; bigWig files keep the window of chromosome 6 at their native
resolution; the expression table keeps the genes with a transcript in the window (statistics of the genome-wide
PyDESeq2 run, not recomputed).

Usage: python prepare_example.py <data_dir> <results_dir> <example/data>
"""

import glob
import gzip
import os
import sys

import pandas as pd
import pyBigWig

CHROM = 'chr6'
START = 49_500_000
END = 54_500_000
CHROM_LENGTH = 149_736_546          # mm10 chr6
MARKS = ('H3K27ac', 'H3K4me2', 'H3K27me3')


def inside(chrom, start, end):
    return chrom == CHROM and int(start) >= START and int(end) <= END


def cut_washu(source, target):
    kept = 0
    with gzip.open(source, 'rt') as handle, gzip.open(target, 'wt', compresslevel=9) as out:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            a = fields[0].split(',')
            b = fields[1].split(',')
            if inside(*a) and inside(*b):
                out.write(line)
                kept += 1
    return kept


def cut_bigwig(source, target):
    bw = pyBigWig.open(source)
    intervals = bw.intervals(CHROM, START, END) or []
    bw.close()
    out = pyBigWig.open(target, 'w')
    out.addHeader([(CHROM, CHROM_LENGTH)])
    if intervals:
        out.addEntries([CHROM] * len(intervals), [i[0] for i in intervals], ends=[i[1] for i in intervals],
                       values=[float(i[2]) for i in intervals])
    out.close()
    return len(intervals)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    data, results, out = sys.argv[1:4]
    k = os.path.join(data, 'kessler')
    os.makedirs(out, exist_ok=True)
    for gsm in ('GSM6505208', 'GSM6505209', 'GSM6505210', 'GSM6505211'):
        source = glob.glob(os.path.join(k, gsm + '_*.txt.gz'))[0]
        print(os.path.basename(source), cut_washu(source, os.path.join(out, os.path.basename(source))))
    for tissue, label in (('MD', 'Md'), ('PA2', 'PA2')):
        for mark in MARKS:
            name = 'E10.5_{}_{}_R1.bigwig'.format(tissue, mark)
            print(name, cut_bigwig(os.path.join(k, name), os.path.join(out, name)))
            name = '{}_{}.fold.bw'.format(label, mark)
            print(name, cut_bigwig(os.path.join(results, 'minoux_fold', name), os.path.join(out, name)))
    # annotation: every line of a gene with a feature in the window
    genes = set()
    with open(os.path.join(k, 'mm10.ncbiRefSeq.gtf')) as handle:
        lines = [line for line in handle if line.startswith(CHROM + '\t')]
    for line in lines:
        f = line.split('\t')
        if int(f[3]) < END and int(f[4]) > START:
            genes.add(f[8].split('gene_id "', 1)[1].split('"', 1)[0])
    with gzip.open(os.path.join(out, 'mm10.ncbiRefSeq.chr6_hoxa.gtf.gz'), 'wt', compresslevel=9) as handle:
        for line in lines:
            if line.split('gene_id "', 1)[1].split('"', 1)[0] in genes:
                handle.write(line)
    deseq = pd.read_csv(os.path.join(results, 'kessler', 'deseq2_PA2_vs_Md.tsv'), sep='\t')
    deseq = deseq[deseq['gene_name'].isin(genes)]
    deseq.to_csv(os.path.join(out, 'deseq2_PA2_vs_Md.chr6_hoxa.tsv'), sep='\t', index=False)
    print('genes', len(genes), 'in the expression table', len(deseq))
    # domains and boundaries of PA2 (hicFindTADs), HIRE regions
    tads = os.path.join(results, 'tad_crossing', 'hicfindtads', 'PA2_25kb')
    for suffix in ('_domains.bed', '_boundaries.bed', '_score.bedgraph'):
        table = pd.read_csv(tads + suffix, sep='\t', header=None, comment='#')
        # hicFindTADs names the chromosomes of these coolers without 'chr'
        table = table[table[0].astype(str).isin([CHROM, CHROM[3:]]) & (table[1] < END) & (table[2] > START)]
        table.to_csv(os.path.join(out, 'PA2_25kb' + suffix), sep='\t', header=False, index=False)
    hire = pd.read_csv(os.path.join(results, 'tad_crossing', 'kessler2023_HIRE_mm10.bed'), sep='\t', header=None)
    hire.to_csv(os.path.join(out, 'kessler2023_HIRE_mm10.bed'), sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
