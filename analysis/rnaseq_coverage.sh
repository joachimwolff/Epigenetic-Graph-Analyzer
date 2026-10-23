#!/usr/bin/env bash
# RNA-seq coverage per condition for the genome-browser panels: the reads of all replicates of a condition are mapped
# together with STAR, and the unique alignments give one coverage track per condition in reads per million (bigWig,
# unstranded), <out_dir>/<name prefix><condition>_RNA-seq.rpm.bw.
# Used for the mouse E10.5 Md and PA2 (Minoux et al. 2017, GSE89434; prefix E10.5_, mm10) and for hESC and hNEC
# (Freire-Pritchett et al. 2017, GSE86821; hg19 with GENCODE v19 junctions).
# Coverage per million reads grows with the read length, so the reads of every condition are clipped at their 3' end
# to the shortest read length of all runs (taken from the first read of each file).
#   samples.tsv  columns condition, replicate, gsm, srr, layout, fastq (files separated by ';', R1 before R2)
#   star_index   STAR mm10 genome directory
# STAR keeps its temporary files (and the FIFOs of --readFilesCommand) under $TMPDIR on a local disk, since network
# shares may not support FIFOs.
set -euo pipefail
samples=${1:?usage: rnaseq_coverage.sh <samples.tsv> <star_index> <out_dir> [threads] [name prefix]}
index=${2:?}
out=${3:?}
threads=${4:-16}
name=${5:-}
star=${STAR:-STAR}
bg2bw=${BEDGRAPHTOBIGWIG:-bedGraphToBigWig}
mkdir -p "$out"
dir=$(dirname "$samples")
# length of the first read; zcat ends with SIGPIPE once awk has it
read_length() { zcat "$1" | awk 'NR == 2 {print length($0); exit}' || true; }
shortest=$(awk -F'\t' 'NR > 1 {split($6, f, ";"); print f[1]}' "$samples" | while read -r f; do
    read_length "$dir/$f"; done | sort -n | head -1)
for condition in $(awk -F'\t' 'NR > 1 {print $1}' "$samples" | sort -u); do
    # R1 (and R2) files of all replicates of the tissue, comma-separated for --readFilesIn
    r1=$(awk -F'\t' -v c="$condition" -v d="$(dirname "$samples")" 'NR > 1 && $1 == c {
            n = split($6, f, ";"); printf "%s%s/%s", (s++ ? "," : ""), d, f[1] }' "$samples")
    r2=$(awk -F'\t' -v c="$condition" -v d="$(dirname "$samples")" 'NR > 1 && $1 == c {
            n = split($6, f, ";"); if (n > 1) printf "%s%s/%s", (s++ ? "," : ""), d, f[2] }' "$samples")
    prefix="$out/$condition."
    first=$(awk -F'\t' -v c="$condition" 'NR > 1 && $1 == c {split($6, f, ";"); print f[1]; exit}' "$samples")
    clip=$(( $(read_length "$dir/$first") - shortest ))
    if [[ ! -s "${prefix}Aligned.sortedByCoord.out.bam" ]]; then
        "$star" --runThreadN "$threads" --genomeDir "$index" --readFilesIn "$r1" ${r2:+"$r2"} \
            --readFilesCommand zcat --clip3pNbases "$clip" --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --limitBAMsortRAM 20000000000 --outFileNamePrefix "$prefix" \
            --outTmpDir "${TMPDIR:-/tmp}/star_${condition}_$$"
    fi
    "$star" --runMode inputAlignmentsFromBAM --inputBAMfile "${prefix}Aligned.sortedByCoord.out.bam" \
        --outWigType bedGraph --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix "$prefix" \
        --outTmpDir "${TMPDIR:-/tmp}/star_${condition}_wig_$$"
    LC_ALL=C sort -k1,1 -k2,2n "${prefix}Signal.Unique.str1.out.bg" > "${prefix}unique.bg"
    "$bg2bw" "${prefix}unique.bg" "$index/chrNameLength.txt" "$out/${name}${condition}_RNA-seq.rpm.bw"
    rm -f "${prefix}unique.bg" "${prefix}"Signal.*.bg
done
