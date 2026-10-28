#!/usr/bin/env bash
# H1-hESC promoter capture Hi-C of Jung et al. 2019 (GSE86189, SRP083285) from raw reads with standard tools:
# HiCUP (HindIII, hg19, bowtie2) per sequencing run, then all runs of the seven libraries pooled into one sample (as
# in Jung et al.; a single library has a median of 90 reads per bait, below CHiCAGO's minNPerBait of 250),
# bam2chicago and CHiCAGO with default settings; interactions with score >= 5 exported as ibed and WashU text.
#
# Usage: 01_h1_chicago.sh <work dir with sra/*.fastq.gz and design/> <bowtie2 index prefix> <HiCUP digest>
#                         <chicagoTools dir> [threads]
# design/ holds hg19_HindIII.rmap, jung2019_H1.baitmap (analysis/jung_baitmap.py), its 4-column version and the
# .npb/.nbpb/.poe files of makeDesignFiles_py3.py (--minFragLen 150 --maxFragLen 40000).
set -euo pipefail
work=${1:?usage: 01_h1_chicago.sh <work dir> <bowtie2 index> <digest> <chicagoTools dir> [threads]}
index=${2:?}
digest=${3:?}
tools=${4:?}
threads=${5:-32}
here=$(cd "$(dirname "$0")" && pwd)
bowtie2=$(command -v bowtie2)
# library (GEO sample) -> sequencing runs
declare -A libraries=(
    [ADS127]="SRR4094708" [ADS160]="SRR4094709 SRR4094710" [ADS171]="SRR4094711 SRR4094712"
    [ADS180]="SRR4094713 SRR4094714" [ADS187]="SRR4094715 SRR4094716" [ADS389]="SRR4094717"
    [ADS393]="SRR4094718"
)
mkdir -p "$work/hicup" "$work/bam" "$work/chinput" "$work/chicago"
for library in "${!libraries[@]}"; do
    bams=()
    for run in ${libraries[$library]}; do
        out="$work/hicup/$run"
        if ! ls "$out"/*.hicup.bam > /dev/null 2>&1; then
            mkdir -p "$out"
            hicup --bowtie2 "$bowtie2" --index "$index" --digest "$digest" --format Sanger --threads "$threads" \
                --zip --outdir "$out" "$work/sra/${run}_1.fastq.gz" "$work/sra/${run}_2.fastq.gz" \
                > "$out/hicup.log" 2>&1
        fi
        bams+=("$(ls "$out"/*.hicup.bam)")
    done
    merged="$work/bam/H1_$library.bam"
    [[ -s "$merged" ]] || samtools merge -n -@ "$threads" -f "$merged" "${bams[@]}"
    if [[ ! -s "$work/chinput/H1_$library/H1_$library.chinput" ]]; then
        (cd "$work/chinput" && bash "$tools/bam2chicago.sh" "$merged" "$work/design/jung2019_H1.baitmap" \
            "$work/design/hg19_HindIII.rmap" "H1_$library" > "H1_$library.log" 2>&1)
    fi
done
pooled="$work/bam/H1_pooled.bam"
[[ -s "$pooled" ]] || samtools merge -n -@ "$threads" -f "$pooled" "$work"/bam/H1_ADS*.bam
if [[ ! -s "$work/chinput/H1_pooled/H1_pooled.chinput" ]]; then
    (cd "$work/chinput" && bash "$tools/bam2chicago.sh" "$pooled" "$work/design/jung2019_H1.baitmap" \
        "$work/design/hg19_HindIII.rmap" "H1_pooled" > "H1_pooled.log" 2>&1)
fi
Rscript "$here/run_chicago.R" "$work/design" "$work/chicago/H1" "$work/chinput/H1_pooled/H1_pooled.chinput" \
    > "$work/chicago/H1.log" 2>&1
