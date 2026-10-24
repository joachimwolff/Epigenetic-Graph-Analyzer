#!/usr/bin/env bash
# Calibration with unchanged genes (pseudo_de.py, as in 42_calibration.sh) for hESC / hNEC and Md / PA2, with the
# inputs and settings of 31_step_decomposition.sh: genes whose expression does not change receive the real fold
# changes; 43_calibration_summary.py compares the explanation levels with the real runs.
#
# Usage: 42b_calibration_hesc_md.sh <data_dir> <results_dir> <out_dir> [seeds] [gene permutations] [threads]
set -euo pipefail
data=${1:?usage: 42b_calibration_hesc_md.sh <data_dir> <results_dir> <out_dir> [seeds] [gene permutations] [threads]}
results=${2:?}
out=${3:?}
seeds=${4:-3}
gene_permutations=${5:-200000}
threads=${6:-8}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
common=(--verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 --chromatinThreshold 1
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads"
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000
        --stepFeature decomposed)
study="$data/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC"
roadmap() { for mark in H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
k="$data/kessler"
kessler() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
fold() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$results/minoux_fold/$1_${mark}.fold.bw"; done; }

hesc() {  # <seed>
    local dir="$out/hESC_unchanged_$1"
    mkdir -p "$dir"
    python "$here/pseudo_de.py" "$results/differential/deseq2_hNEC_vs_hESC.tsv" "$dir/expression.tsv" --seed "$1" \
        > "$dir/expression.log"
    egaDifferential "${common[@]}" \
        --interactionsA "$study/GSE86821_hESC12_merged_washU.txt.gz" \
        --interactionsB "$study/GSE86821_hNPC12_merged_washU.txt.gz" \
        --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
        --scoreThreshold 12 --referenceScore 12 --names hESC hNEC \
        --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E008) --tracksB $(roadmap E009) \
        --markSign H3K27me3=-1 --markRoles H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
        --expression "$dir/expression.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
        --alpha 0.05 --geneKey genes \
        --proximitySupport "$results/reviewer_controls/hESC_gam_hops_coverage.tsv" \
        -o "$dir/hESC_vs_hNEC" > "$dir.log" 2>&1
}

md() {  # <seed>
    local dir="$out/Md_unchanged_$1"
    mkdir -p "$dir"
    python "$here/pseudo_de.py" "$results/revision_bioinf/decomposed/deseq2_PA2_vs_Md.tsv" "$dir/expression.tsv" \
        --seed "$1" > "$dir/expression.log"
    egaDifferential "${common[@]}" \
        --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
        --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
        --format washu --replicates intersect --names Md PA2 \
        --tss "$k/mm10.ncbiRefSeq.gtf" --geneTypes \
        --tracksA $(kessler MD) --tracksB $(kessler PA2) --markSign H3K27me3=-1 \
        --markRoles H3K27ac=active H3K4me2=active H3K27me3=repressive \
        --stateTracksA $(fold Md) --stateTracksB $(fold PA2) \
        --expression "$dir/expression.tsv" --log2FCColumn log2FoldChange \
        --padjColumn padj --alpha 0.05 --geneKey genes \
        -o "$dir/Md_vs_PA2" > "$dir.log" 2>&1
}

for seed in $(seq 0 $((seeds - 1))); do
    hesc "$seed" &
    md "$seed" &
done
wait
