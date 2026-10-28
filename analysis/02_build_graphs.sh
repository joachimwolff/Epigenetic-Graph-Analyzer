#!/usr/bin/env bash
# Interaction graphs and communities from published interaction calls.
#  * GM12878, CHiCAGO calls of Cairns et al. 2016 (GSE81503), score >= 5
#  * GM12878, Mifsud et al. 2015 calls restricted to the strong core (>= 20 reads, <= 2 Mb)
#  * H9 hESC and hESC-derived neural progenitors (hNPC), CHiCAGO calls of Freire-Pritchett et al. 2017
#    (GSE86821, merged replicates, score >= 12 as published); these WashU tracks list pairs in both
#    orientations, so the bait status of all hg19 sets is taken from the capture design
#  * H1-hESC, CHiCAGO calls from the raw reads of Jung et al. 2019 (analysis/01_h1_chicago.sh), score >= 5
#  * K562 (DMSO), CHiCAGO calls of GSE250135 (GRCh38), score >= 5
set -euo pipefail
data=${1:?usage: 02_build_graphs.sh <data_dir> <results_dir>}
results=${2:?}
gtf="$data/annotation/gencode.v19.annotation.gtf.gz"
gtf38="$data/annotation/gencode.v29.annotation.gtf.gz"
calls="$data/chicago_calls"
# bait map of the hg19 HindIII promoter capture design shared by all hg19 data sets (Javierre et al. 2016, OSF)
baits="$calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap"
mkdir -p "$results/graphs"
graph() {
    local name=$1
    shift
    egaCreateGraph "$@" --maxDistance 2000000 --weight score --communityMethod leiden --resolution 1.0 --seed 42 \
        -o "$results/graphs/$name" > "$results/graphs/$name.log" 2>&1
}
graph GM_chicago -i "$calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz" \
    -f washu --scoreThreshold 5 --baits "$baits" --tss "$gtf" &
graph GM_mifsud_c20_r1 -i "$data/mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt" -f mifsud \
    --minCount 20 --baits "$baits" --tss "$gtf" &
graph hESC -i "$calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC/GSE86821_hESC12_merged_washU.txt.gz" \
    -f washu --scoreThreshold 12 --baits "$baits" --tss "$gtf" &
graph hNPC -i "$calls/GSE86821_freire_pritchett_2017_H9_hESC_hNPC/GSE86821_hNPC12_merged_washU.txt.gz" \
    -f washu --scoreThreshold 12 --baits "$baits" --tss "$gtf" &
graph H1 -i "$calls/H1_jung2019_reprocessed/H1.ibed.gz" -f ibed --scoreThreshold 5 --tss "$gtf" &
graph K562 -i "$data/k562/contacts/GSM7974410_k562_dmso_washU.txt.gz" -f washu --scoreThreshold 5 --tss "$gtf38" &
wait
for name in GM_chicago GM_mifsud_c20_r1 hESC hNPC H1; do
    egaIndirectContacts -g "$results/graphs/$name.graphml" --maxDistance 2000000 --triplets \
        -o "$results/graphs/$name"
done
