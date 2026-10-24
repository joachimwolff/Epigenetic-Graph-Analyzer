#!/usr/bin/env bash
# Sensitivity of the explanations beyond step 1 between naive B cells and GM12878 to the calling threshold: the run of
# 35_gm12878_naiveB.sh repeated with CHiCAGO score >= <threshold> for both conditions (absence still tested against
# the sub-threshold scores of naive B); 46_threshold_overlap.py compares the explained genes with the run at score 5.
#
# Usage: 45_threshold_sensitivity.sh <data_dir> <real run dir> <out_dir> [threshold] [gene permutations] [threads]
set -euo pipefail
data=${1:?usage: 45_threshold_sensitivity.sh <data_dir> <real run dir> <out_dir> [threshold] [gene permutations] [threads]}
real=${2:?}
out=${3:?}
threshold=${4:-8}
gene_permutations=${5:-200000}
threads=${6:-10}
mkdir -p "$out"
nb="$data/naiveB"
roadmap() { for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
egaDifferential --verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 --chromatinThreshold 1 \
    --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --stepFeature decomposed \
    --interactionsA "$nb/javierre2016/naiveB_nB_washU.txt.gz" \
    --interactionsB "$data/chicago_calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz" \
    --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold "$threshold" --referenceScore 5 --names naiveB GM12878 \
    --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E032) --tracksB $(roadmap E116) \
    --markSign H3K27me3=-1 --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
    --expression "$real/deseq2_GM12878_vs_naiveB.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes -o "$out/naiveB_vs_GM12878" > "$out/naiveB_vs_GM12878.log" 2>&1
