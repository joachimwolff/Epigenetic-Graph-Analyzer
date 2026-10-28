#!/usr/bin/env bash
# Mark term and contact term of the step features (egaDifferential --stepFeature decomposed) for hESC / hNEC and
# Md / PA2, with the inputs and settings of 16_hops_n.sh. The step levels of the cascade, the learned roles, the
# mirror control and the per-gene permutation test use the mark term (marks of B versus A under the fixed strengths
# r_bar = (r_A + r_B) / 2); the combined feature and the contact term are written beside it. The Md / PA2 expression
# table is recomputed by kessler_deseq2.py (repeated GEO libraries used once). Outputs in <out>/ only;
# step_decomposition_summary.py compares the three features per level and the learned roles with the combined
# feature of <results>/hops_n. An optional sixth argument sets egaDifferential --stepWeighting (product, the default,
# unweighted or normalised; 33_step_weighting.sh).
set -euo pipefail
data=${1:?usage: 31_step_decomposition.sh <data_dir> <results_dir> <out_dir> [gene permutations] [threads] [step weighting]}
results=${2:?}
out=${3:?}
gene_permutations=${4:-200000}
threads=${5:-6}
weighting=${6:-product}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
common=(--verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 --chromatinThreshold 1
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads"
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000
        --stepFeature decomposed)
if [ "$weighting" != product ]; then
    common+=(--stepWeighting "$weighting")
fi

# hESC / hNEC
study="$data/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC"
roadmap() { for mark in H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
egaDifferential "${common[@]}" \
    --interactionsA "$study/GSE86821_hESC12_merged_washU.txt.gz" \
    --interactionsB "$study/GSE86821_hNPC12_merged_washU.txt.gz" \
    --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 12 --referenceScore 12 --names hESC hNEC \
    --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E008) --tracksB $(roadmap E009) \
    --markSign H3K27me3=-1 --markRoles H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
    --expression "$results/differential/deseq2_hNEC_vs_hESC.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/reviewer_controls/hESC_gam_hops_coverage.tsv" \
    -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1

# Md / PA2
k="$data/kessler"
python "$here/kessler_deseq2.py" "$k/GSE89434_Minoux_RNAseq_gene_counts.csv.gz" "$out/deseq2_PA2_vs_Md.tsv"
kessler() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
fold() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$results/minoux_fold/$1_${mark}.fold.bw"; done; }
egaDifferential "${common[@]}" \
    --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
    --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
    --format washu --replicates intersect --names Md PA2 \
    --tss "$k/mm10.ncbiRefSeq.gtf" --geneTypes \
    --tracksA $(kessler MD) --tracksB $(kessler PA2) --markSign H3K27me3=-1 \
    --markRoles H3K27ac=active H3K4me2=active H3K27me3=repressive \
    --stateTracksA $(fold Md) --stateTracksB $(fold PA2) \
    --expression "$out/deseq2_PA2_vs_Md.tsv" --log2FCColumn log2FoldChange \
    --padjColumn padj --alpha 0.05 --geneKey genes \
    -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1

for comparison in hESC_vs_hNEC Md_vs_PA2; do
    egaExplain -p "$out/$comparison" -o "$out/explain/$comparison" > "$out/explain_$comparison.log" 2>&1
done
python "$here/step_decomposition_summary.py" "$out" "$results/hops_n" > "$out/summary.log" 2>&1
