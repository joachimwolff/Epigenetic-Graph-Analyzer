#!/usr/bin/env bash
# RNA-seq coverage of the mouse E10.5 mandibular arch (Md) and second pharyngeal arch (PA2) for the genome-browser
# panels: the reads of the three replicates per tissue (Minoux et al. 2017, GSE89434) are mapped together with STAR,
# and the unique alignments give one coverage track per tissue in reads per million (bigWig, unstranded).
#   samples.tsv  columns condition, replicate, gsm, srr, layout, fastq (files separated by ';', R1 before R2)
#   star_index   STAR mm10 genome directory
# STAR keeps its temporary files (and the FIFOs of --readFilesCommand) under $TMPDIR on a local disk, since network
# shares may not support FIFOs.
set -euo pipefail
samples=${1:?usage: minoux_rnaseq_coverage.sh <samples.tsv> <star_index> <out_dir> [threads]}
index=${2:?}
out=${3:?}
threads=${4:-16}
star=${STAR:-STAR}
bg2bw=${BEDGRAPHTOBIGWIG:-bedGraphToBigWig}
mkdir -p "$out"
for condition in Md PA2; do
    # R1 (and R2) files of all replicates of the tissue, comma-separated for --readFilesIn
    r1=$(awk -F'\t' -v c="$condition" -v d="$(dirname "$samples")" 'NR > 1 && $1 == c {
            n = split($6, f, ";"); printf "%s%s/%s", (s++ ? "," : ""), d, f[1] }' "$samples")
    r2=$(awk -F'\t' -v c="$condition" -v d="$(dirname "$samples")" 'NR > 1 && $1 == c {
            n = split($6, f, ";"); if (n > 1) printf "%s%s/%s", (s++ ? "," : ""), d, f[2] }' "$samples")
    prefix="$out/$condition."
    if [[ ! -s "${prefix}Aligned.sortedByCoord.out.bam" ]]; then
        "$star" --runThreadN "$threads" --genomeDir "$index" --readFilesIn "$r1" ${r2:+"$r2"} \
            --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 \
            --limitBAMsortRAM 20000000000 --outFileNamePrefix "$prefix" \
            --outTmpDir "${TMPDIR:-/tmp}/star_${condition}_$$"
    fi
    "$star" --runMode inputAlignmentsFromBAM --inputBAMfile "${prefix}Aligned.sortedByCoord.out.bam" \
        --outWigType bedGraph --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix "$prefix" \
        --outTmpDir "${TMPDIR:-/tmp}/star_${condition}_wig_$$"
    LC_ALL=C sort -k1,1 -k2,2n "${prefix}Signal.Unique.str1.out.bg" > "${prefix}unique.bg"
    "$bg2bw" "${prefix}unique.bg" "$index/chrNameLength.txt" "$out/E10.5_${condition}_RNA-seq.rpm.bw"
    rm -f "${prefix}unique.bg" "${prefix}"Signal.*.bg
done
