#!/usr/bin/env bash
# H1-derived neural progenitor (NPC) promoter capture Hi-C of Jung et al. 2019 (GSE86189, library ADS414, runs
# SRR4094802 and SRR4094803) from raw reads, with the commands of 01_h1_chicago.sh: HiCUP (HindIII, hg19, bowtie2) per
# run, the runs merged, bam2chicago and CHiCAGO with default settings on the H1 design (same capture design);
# interactions with score >= 5 exported as ibed and WashU text.
#
# Usage: 38_npc_chicago.sh <work dir with sra/*.fastq.gz> <H1 design dir> <bowtie2 index prefix> <HiCUP digest>
#                          <chicagoTools dir> [threads]
set -euo pipefail
work=${1:?usage: 38_npc_chicago.sh <work dir> <design dir> <bowtie2 index> <digest> <chicagoTools dir> [threads]}
design=${2:?}
index=${3:?}
digest=${4:?}
tools=${5:?}
threads=${6:-16}
here=$(cd "$(dirname "$0")" && pwd)
bowtie2=$(command -v bowtie2)
runs=(SRR4094802 SRR4094803)
mkdir -p "$work/hicup" "$work/bam" "$work/chinput" "$work/chicago"
[[ -e "$work/design" ]] || ln -s "$design" "$work/design"
bams=()
for run in "${runs[@]}"; do
    out="$work/hicup/$run"
    if ! ls "$out"/*.hicup.bam > /dev/null 2>&1; then
        mkdir -p "$out"
        hicup --bowtie2 "$bowtie2" --index "$index" --digest "$digest" --format Sanger --threads "$threads" \
            --zip --outdir "$out" "$work/sra/${run}_1.fastq.gz" "$work/sra/${run}_2.fastq.gz" > "$out/hicup.log" 2>&1
    fi
    bams+=("$(ls "$out"/*.hicup.bam)")
done
merged="$work/bam/NPC_ADS414.bam"
[[ -s "$merged" ]] || samtools merge -n -@ "$threads" -f "$merged" "${bams[@]}"
if [[ ! -s "$work/chinput/NPC_ADS414/NPC_ADS414.chinput" ]]; then
    (cd "$work/chinput" && bash "$tools/bam2chicago.sh" "$merged" "$work/design/jung2019_H1.baitmap" \
        "$work/design/hg19_HindIII.rmap" "NPC_ADS414" > "NPC_ADS414.log" 2>&1)
fi
Rscript "$here/run_chicago.R" "$work/design" "$work/chicago/NPC" "$work/chinput/NPC_ADS414/NPC_ADS414.chinput" \
    > "$work/chicago/NPC.log" 2>&1
