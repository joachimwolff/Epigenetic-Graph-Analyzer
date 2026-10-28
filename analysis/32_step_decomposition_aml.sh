#!/usr/bin/env bash
# egaDifferential --stepFeature decomposed for Kasumi-1 siMM / siRE and AML t(8;21) / FLT3-ITD, with the inputs
# and settings of 16_hops_n.sh (as 31_step_decomposition.sh for hESC / hNEC and Md / PA2). Outputs in <out>/ only.
set -euo pipefail
results=${1:?usage: 32_step_decomposition_aml.sh <results_dir> <kasumi1 dir> <assi2019 dir> <out_dir> [gene permutations] [threads]}
kasumi=${2:?}
assi=${3:?}
out=${4:?}
gene_permutations=${5:-200000}
threads=${6:-6}
data=$(dirname "$results")/raw
data=${EGA_DATA:-$data}
mkdir -p "$out"
common=(--verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 --chromatinThreshold 1
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads"
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000
        --stepFeature decomposed)

# Kasumi-1 siMM / siRE (08_aml_kasumi1.sh)
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
egaDifferential "${common[@]}" \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 5 --names siMM siRE \
    --markRoles DNase_d2=active DNase_d10=active --stateScale genomeMean \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    --expression "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1

# AML t(8;21) / FLT3-ITD (07_aml_assi.sh)
bw="$assi/bigwig"
join() { local IFS=,; echo "$*"; }
t821_dhs=$(join "$bw"/GSM289362[789]_t821*_DHS.bw "$bw"/GSM2893630_t821KIT-2_DHS.bw)
itd_dhs=$(join "$bw"/GSM289360[67]_ITD-*_DHS.bw)
egaDifferential "${common[@]}" \
    --interactionsA "$assi/GSM2895201_t821-1_CHiC_interactions.txt.gz" \
    --interactionsB "$assi/GSM2895200_FLT3-ITD_CHiC_interactions.txt.gz" \
    --format homer --baits "$data/chicago_calls/design_hg38_HindIII_lifted/baits_hg38_lifted.bed" --maxFdr 0.01 \
    --referenceScore 2 --names t821 FLT3ITD --markRoles DHS=active --stateScale genomeMean \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" --tracksA "DHS=$t821_dhs" --tracksB "DHS=$itd_dhs" \
    --expression "$results/aml_assi/expression_FLT3ITD_vs_t821.tsv" --log2FCColumn log2FoldChange \
    --padjColumn pvalue --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1

for comparison in siMM_vs_siRE t821_vs_FLT3ITD; do
    egaExplain -p "$out/$comparison" -o "$out/explain/$comparison" > "$out/explain_$comparison.log" 2>&1
done
