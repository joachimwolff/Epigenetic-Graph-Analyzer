#!/usr/bin/env bash
# Differential mode with explanations up to three graph steps (egaDifferential --maxHops 3): promoter, 1 hop
# (direct), 2 hops (indirect), 3 hops. The four comparisons of 04, 05, 07 and 08 with the same inputs and settings,
# written to <results>/hops_n/ (the outputs of 04-08 are not touched). The expression tables and bigWig files
# prepared by 04-08 are reused, so those scripts must have run before.
# --proximitySupport adds the 3D support of the explaining level of every gene (analysis/hop_distance.py step curves):
#   hESC / hNEC:          GAM of H1 hESC on the hESC graph (11_reviewer_controls.sh)
#   Md / PA2:             none (mouse, no multi-way data)
#   Kasumi-1 and AML:     HiPore-C of K562 (09_porec.sh), the closest myeloid cell line with multi-way data
set -euo pipefail
data=${1:?usage: 16_hops_n.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> [permutations] [threads]}
results=${2:?}
kasumi=${3:?}
assi=${4:?}
permutations=${5:-200}
threads=${6:-24}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/hops_n"
mkdir -p "$out"
common=(--verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1
        --permutations "$permutations" --bootstrap 1000 --seed 0 --threads "$threads")

# hESC / hNEC (04_differential.sh)
study="$data/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC"
roadmap() { for mark in H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
egaDifferential "${common[@]}" \
    --interactionsA "$study/GSE86821_hESC12_merged_washU.txt.gz" \
    --interactionsB "$study/GSE86821_hNPC12_merged_washU.txt.gz" \
    --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 12 --referenceScore 12 --names hESC hNEC \
    --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E008) --tracksB $(roadmap E009) \
    --markSign H3K27me3=-1 \
    --expression "$results/differential/deseq2_hNEC_vs_hESC.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/reviewer_controls/hESC_gam_hops_coverage.tsv" \
    -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1

# Md / PA2 (05_kessler_md_pa2.sh)
k="$data/kessler"
kessler() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
egaDifferential "${common[@]}" \
    --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
    --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
    --format washu --replicates intersect --names Md PA2 \
    --tss "$k/mm10.ncbiRefSeq.gtf" --geneTypes \
    --tracksA $(kessler MD) --tracksB $(kessler PA2) --markSign H3K27me3=-1 \
    --expression "$results/kessler/deseq2_PA2_vs_Md.tsv" --log2FCColumn log2FoldChange \
    --padjColumn padj --alpha 0.05 --geneKey genes \
    -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1

# Kasumi-1 siMM / siRE (08_aml_kasumi1.sh)
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
egaDifferential "${common[@]}" \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 5 --names siMM siRE \
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
    --referenceScore 2 --names t821 FLT3ITD \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" --tracksA "DHS=$t821_dhs" --tracksB "DHS=$itd_dhs" \
    --expression "$results/aml_assi/expression_FLT3ITD_vs_t821.tsv" --log2FCColumn log2FoldChange \
    --padjColumn pvalue --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1

# explanation cascade per comparison (one PDF each)
python "$here/figure_explanation_cascade.py" "$results" "$out/figures"
