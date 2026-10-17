#!/usr/bin/env bash
# Blind test of the role discovery: the hESC / hNEC comparison of 16_hops_n.sh with additional Roadmap tracks of E008
# (H9 hESC) and E009 (H9-derived neuronal progenitors) given with the role 'unknown' (no declared role, no sign):
# H3K9me3, H3K36me3 and H2A.Z (each when present for both epigenomes; a track of one epigenome only, e.g. the WGBS
# methylation of E008, cannot enter a comparison). The roles of the unknown tracks are learned per level and per
# element type of the nodes (<prefix>_mark_roles.tsv, column element_type) and compared with the known biology in
# analysis/blind_roles.py. Outputs: <results>/hops_n/blind/. Same settings as 16_hops_n.sh.
set -euo pipefail
data=${1:?usage: 18_blind_roles.sh <data_dir> <results_dir> [permutations] [threads]}
results=${2:?}
permutations=${3:-200}
threads=${4:-24}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/hops_n/blind"
mkdir -p "$out"
roadmap="$data/roadmap"
known=(H3K4me1 H3K4me3 H3K27me3)
unknown=()
for mark in H3K9me3 H3K36me3 H2A.Z; do
    if [[ -e "$roadmap/E008-$mark.fc.signal.bigwig" && -e "$roadmap/E009-$mark.fc.signal.bigwig" ]]; then
        unknown+=("$mark")
    fi
done
tracks() { for mark in "${known[@]}" "${unknown[@]}"; do printf '%s=%s ' "$mark" "$roadmap/$1-$mark.fc.signal.bigwig"; done; }
roles=(H3K4me3=active H3K4me1=enhancer H3K27me3=repressive)
for mark in "${unknown[@]}"; do roles+=("$mark=unknown"); done
study="$data/chicago_calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC"
egaDifferential --verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1 \
    --permutations "$permutations" --bootstrap 1000 --seed 0 --threads "$threads" \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --interactionsA "$study/GSE86821_hESC12_merged_washU.txt.gz" \
    --interactionsB "$study/GSE86821_hNPC12_merged_washU.txt.gz" \
    --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 12 --referenceScore 12 --names hESC hNEC \
    --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(tracks E008) --tracksB $(tracks E009) \
    --markSign H3K27me3=-1 --markRoles "${roles[@]}" --primaryRepressive H3K27me3 \
    --expression "$results/differential/deseq2_hNEC_vs_hESC.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/reviewer_controls/hESC_gam_hops_coverage.tsv" \
    -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1
python "$here/blind_roles.py" "$out/hESC_vs_hNEC" "$out/blind_roles.tsv" "${unknown[@]}" > "$out/blind_roles.log" 2>&1
# DNA methylation of hESC only (no WGBS of E009): a state annotation of the element classes and promoters of hESC
if [[ -e "$roadmap/E008_WGBS_FractionalMethylation.bigwig" ]]; then
    python "$here/blind_methylation.py" "$out/hESC_vs_hNEC" "$roadmap/E008_WGBS_FractionalMethylation.bigwig" hESC \
        "$out/methylation_E008" > "$out/methylation_E008.log" 2>&1
fi
