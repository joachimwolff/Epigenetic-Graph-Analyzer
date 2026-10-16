#!/usr/bin/env bash
# Robustness of the per-gene candidates for regulation through indirect contacts (reviewer requests), from the
# outputs of 04_differential.sh, 05_kessler_md_pa2.sh and 08_aml_kasumi1.sh:
#   gene_test     the per-gene permutation test of the indirect level recomputed from the union graphs with more
#                 permutations (analysis/indirect_gene_test.py; same seed, so the first 1000 permutations are those
#                 of the original run, which is reproduced exactly with 1000)
#   subthreshold  direct pairs below the calling threshold between each candidate promoter and its explaining
#                 two-hop nodes, and the candidates reclassified with these pairs as direct edges
#                 (analysis/subthreshold_contacts.py), on the lowest score each data set provides:
#                 Kasumi-1   CHiCAGO exports at score >= 1; pairs at score >= 2 under the replicate rule of the run
#                 Md/PA2     published WashU files hold score >= 5 only; pairs called in one of the two replicates
#                            (the run keeps pairs called in both)
#                 hESC/hNEC  published files hold score >= 12 only (the calling threshold); H1-hESC CHiCAGO calls
#                            at score >= 5 (Jung et al. 2019, analysis/01_h1_chicago.sh, same HindIII fragments) as a
#                            proxy for both conditions, no hNEC counterpart exists
set -euo pipefail
results=${1:?usage: 13_candidate_robustness.sh <results_dir> <data_dir> <kasumi1 dir> [gene_test|subthreshold|all]}
data=${2:?}
kasumi=${3:?}
step=${4:-all}
permutations=${PERMUTATIONS:-20000 200000}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/robustness_candidates"
mkdir -p "$out"
if [[ $step == gene_test || $step == all ]]; then
    for n in $permutations; do
        python "$here/indirect_gene_test.py" --prefix "$results/differential/hESC_vs_hNEC" --names hESC hNEC \
            --markSign H3K27me3=-1 --permutations "$n" -o "$out/hESC_vs_hNEC_perm$n" > "$out/hESC_vs_hNEC_perm$n.log" 2>&1
        python "$here/indirect_gene_test.py" --prefix "$results/kessler/Md_vs_PA2" --names Md PA2 \
            --markSign H3K27me3=-1 --permutations "$n" -o "$out/Md_vs_PA2_perm$n" > "$out/Md_vs_PA2_perm$n.log" 2>&1
        python "$here/indirect_gene_test.py" --prefix "$results/aml_kasumi1/siMM_vs_siRE" --names siMM siRE \
            --permutations "$n" -o "$out/siMM_vs_siRE_perm$n" > "$out/siMM_vs_siRE_perm$n.log" 2>&1
    done
fi
[[ $step == subthreshold || $step == all ]] || exit 0
# candidates from the largest gene test; the reclassified candidates are tested again with 20,000 permutations
c=${permutations##* }
n=20000
k="$data/kessler"
kessler_tracks() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
python "$here/subthreshold_contacts.py" --prefix "$results/kessler/Md_vs_PA2" \
    --candidates "$out/Md_vs_PA2_perm${c}_indirect_only_genes.tsv" --names Md PA2 --markSign H3K27me3=-1 \
    --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
    --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
    --format washu --replicates union --referenceScore 5 --bins 5 12 \
    --tracksA $(kessler_tracks MD) --tracksB $(kessler_tracks PA2) --permutations "$n" \
    -o "$out/Md_vs_PA2_subthreshold" > "$out/Md_vs_PA2_subthreshold.log" 2>&1
roadmap_tracks() { for mark in H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
h1="$data/chicago_calls/H1_jung2019_reprocessed/H1.ibed.gz"
python "$here/subthreshold_contacts.py" --prefix "$results/differential/hESC_vs_hNEC" \
    --candidates "$out/hESC_vs_hNEC_perm${c}_indirect_only_genes.tsv" --names hESC hNEC --markSign H3K27me3=-1 \
    --interactionsA "$h1" --interactionsB "$h1" --format ibed \
    --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 5 --referenceScore 12 --bins 5 12 \
    --tracksA $(roadmap_tracks E008) --tracksB $(roadmap_tracks E009) --permutations "$n" \
    -o "$out/hESC_vs_hNEC_subthreshold_H1" > "$out/hESC_vs_hNEC_subthreshold_H1.log" 2>&1
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
python "$here/subthreshold_contacts.py" --prefix "$results/aml_kasumi1/siMM_vs_siRE" \
    --candidates "$out/siMM_vs_siRE_perm${c}_indirect_only_genes.tsv" --names siMM siRE \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 2 --referenceScore 5 --bins 1 2 5 \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    --permutations "$n" -o "$out/siMM_vs_siRE_subthreshold" > "$out/siMM_vs_siRE_subthreshold.log" 2>&1
