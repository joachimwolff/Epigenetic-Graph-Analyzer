#!/usr/bin/env bash
# Stricter controls for the multi-way support (reviewer requests):
#   * step curve with coverage-matched controls (as described in the manuscript) and 95% chromosome block bootstrap
#   * the same with controls also matched for A/B compartments of both ends and for domain membership
#     (GM12878: Rao et al. 2014 subcompartments and Arrowhead domains, hg19; H1-hESC: 4DN Micro-C compartment
#     eigenvector and domains between insulation boundaries, GRCh38; analysis/prepare_domains.py)
#   * cell-type specificity: GM12878-only, hESC-only and shared interactions in GM12878 SPRITE and H1-hESC GAM
# Runs one after another (each run holds a multi-way matrix in memory).
set -euo pipefail
data=${1:?usage: 11_reviewer_controls.sh <data_dir> <results_dir>}
results=${2:?}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/reviewer_controls"
mkdir -p "$out"
sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
gam="$results/multiway/gam_H1_hg38_30kb.npz"
chain="$data/annotation/hg19ToHg38.over.chain.gz"
gm_comp="$data/domains/GM12878_rao2014_hg19/GSE63525_GM12878_subcompartments.bed.gz"
gm_dom="$data/domains/GM12878_rao2014_hg19/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
h1_dir="$data/domains/H1_4dn_hg38"
[[ -s "$h1_dir/H1_domains_hg38.bed" ]] || python "$here/prepare_domains.py" "$h1_dir/4DNFI475YIT8.bw" \
    "$h1_dir/4DNFIED5HLDC.bed.gz" "$h1_dir/H1_compartments_hg38.bedGraph" "$h1_dir/H1_domains_hg38.bed"
h1_comp="$h1_dir/H1_compartments_hg38.bedGraph"
h1_dom="$h1_dir/H1_domains_hg38.bed"

hops() {   # <graph> <multiway> <name> [chain] ; coverage-matched and domain-matched step curves
    local graph=$1 multiway=$2 name=$3 chain_option=() comp dom
    [[ $# -gt 3 ]] && chain_option=(--chain "$4")
    if [[ "$multiway" == "$sprite" ]]; then comp=$gm_comp; dom=$gm_dom; else comp=$h1_comp; dom=$h1_dom; fi
    [[ -s "$out/${name}_hops_coverage.tsv" ]] || python "$here/hop_distance.py" "$results/graphs/$graph.graphml" \
        "$multiway" "$out/${name}_hops_coverage.tsv" "${chain_option[@]}" > "$out/${name}_hops_coverage.log" 2>&1
    [[ -s "$out/${name}_hops_domains.tsv" ]] || python "$here/hop_distance.py" "$results/graphs/$graph.graphml" \
        "$multiway" "$out/${name}_hops_domains.tsv" "${chain_option[@]}" --compartments "$comp" --domains "$dom" \
        > "$out/${name}_hops_domains.log" 2>&1
}
hops GM_chicago "$sprite" GM_chicago_sprite
hops hESC "$gam" hESC_gam "$chain"
hops H1 "$gam" H1_gam "$chain"
[[ -s "$out/GM_vs_hESC_domains.tsv" ]] || python "$here/celltype_specific_edges.py" \
    "$results/graphs/GM_chicago.graphml" "$results/graphs/hESC.graphml" GM12878 hESC "$out/GM_vs_hESC_domains.tsv" \
    --multiway SPRITE_GM12878 "$sprite" --multiway GAM_H1 "$gam" "$chain" \
    --domainsFor SPRITE_GM12878 "$gm_comp" "$gm_dom" --domainsFor GAM_H1 "$h1_comp" "$h1_dom" \
    > "$out/GM_vs_hESC_domains.log" 2>&1
