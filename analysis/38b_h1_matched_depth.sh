#!/usr/bin/env bash
# H1 promoter capture Hi-C at the depth of the NPC library (38_npc_chicago.sh): the pooled H1 HiCUP alignments of
# 01_h1_chicago.sh are subsampled by read name (samtools view -s, both mates of a pair kept together) to the fraction
# NPC captured reads / H1 captured reads (sums of the chinput read counts), then bam2chicago and CHiCAGO with default
# settings on the same design, so that H1 and NPC are called from the same number of captured reads.
#
# Usage: 38b_h1_matched_depth.sh <H1 work dir> <NPC work dir> <chicagoTools dir> [threads] [seed]
set -euo pipefail
h1=${1:?usage: 38b_h1_matched_depth.sh <H1 work dir> <NPC work dir> <chicagoTools dir> [threads] [seed]}
npc=${2:?}
tools=${3:?}
threads=${4:-16}
seed=${5:-1}
here=$(cd "$(dirname "$0")" && pwd)
captured() { awk 'NR > 2 {s += $3} END {print s}' "$1"; }
n_h1=$(captured "$h1/chinput/H1_pooled/H1_pooled.chinput")
n_npc=$(captured "$npc/chinput/NPC_ADS414/NPC_ADS414.chinput")
fraction=$(awk -v a="$n_npc" -v b="$n_h1" 'BEGIN {printf "%.4f", a / b}')
echo "captured reads H1 $n_h1, NPC $n_npc, fraction $fraction"
out="$h1/matched"
mkdir -p "$out/chinput" "$out/chicago"
bam="$out/H1_matched.bam"
# samtools -s takes seed.fraction (integer part the seed)
[[ -s "$bam" ]] || samtools view -@ "$threads" -b -s "$seed${fraction#0}" -o "$bam" "$h1/bam/H1_pooled.bam"
if [[ ! -s "$out/chinput/H1_matched/H1_matched.chinput" ]]; then
    (cd "$out/chinput" && bash "$tools/bam2chicago.sh" "$bam" "$h1/design/jung2019_H1.baitmap" \
        "$h1/design/hg19_HindIII.rmap" "H1_matched" > "H1_matched.log" 2>&1)
fi
echo "captured reads H1 matched $(captured "$out/chinput/H1_matched/H1_matched.chinput")"
Rscript "$here/run_chicago.R" "$h1/design" "$out/chicago/H1_matched" "$out/chinput/H1_matched/H1_matched.chinput" \
    > "$out/chicago/H1_matched.log" 2>&1
