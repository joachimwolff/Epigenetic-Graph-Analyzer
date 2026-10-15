#!/usr/bin/env bash
# Differential mode in human: H9 hESC against hESC-derived neuroectoderm (hNEC), CHiCAGO calls of
# Freire-Pritchett et al. 2017 (GSE86821, merged replicates, score >= 12), RNA-seq of the same study (PyDESeq2,
# analysis/hesc_deseq2.py) and Roadmap Epigenomics marks of H9 (E008) and H9-derived neuronal progenitors (E009);
# E009 has no H3K27ac, so the marks shared by both epigenomes are used.
set -euo pipefail
data=${1:?usage: 04_differential.sh <data_dir> <results_dir> [permutations]}
results=${2:?}
permutations=${3:-200}
gtf="$data/annotation/gencode.v19.annotation.gtf.gz"
study="$data/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC"
baits="$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap"
out="$results/differential"
mkdir -p "$out"
tracks() {
    local eid=$1
    for mark in H3K4me1 H3K4me3 H3K27me3; do
        printf '%s=%s ' "$mark" "$data/roadmap/$eid-$mark.fc.signal.bigwig"
    done
}
[[ -s "$out/deseq2_hNEC_vs_hESC.tsv" ]] || python "$(dirname "$0")/hesc_deseq2.py" \
    "$study/GSE86821_RNA-seq_rawcounts_hESChNPC.txt.gz" "$out/deseq2_hNEC_vs_hESC.tsv"
egaDifferential --verbose \
    --interactionsA "$study/GSE86821_hESC12_merged_washU.txt.gz" \
    --interactionsB "$study/GSE86821_hNPC12_merged_washU.txt.gz" \
    --format washu --baits "$baits" --scoreThreshold 12 --referenceScore 12 --names hESC hNEC \
    --tss "$gtf" --tracksA $(tracks E008) --tracksB $(tracks E009) --markSign H3K27me3=-1 \
    --expression "$out/deseq2_hNEC_vs_hESC.tsv" --log2FCColumn log2FoldChange --padjColumn padj --alpha 0.05 \
    --geneKey genes --rnaThreshold 1 --chromatinThreshold 1 --permutations "$permutations" --bootstrap 1000 \
    --seed 0 -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1
