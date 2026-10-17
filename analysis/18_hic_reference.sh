#!/usr/bin/env bash
# Hi-C reference for the compartment and domain analyses (compartments_domains.py): A/B compartment eigenvectors
# (hicPCA, 100 kb and 250 kb, first eigenvector, sign oriented by protein-coding gene density) and TADs (hicFindTADs,
# 50 kb) of the 4DN mcool files of GM12878 (Hi-C), K562 (Hi-C) and H1-hESC (Micro-C), all GRCh38.
# Output: <results>/compartments_domains/hic_reference/<cell>_pc1_<res>.bedgraph, <cell>_50kb_domains.bed,
#         <cell>_50kb_boundaries.bed, <cell>_50kb_tad_score.bm
# Requires HiCExplorer (hicPCA, hicFindTADs; set HICEXPLORER_BIN to the directory holding them).
set -euo pipefail
mcool=${1:?usage: 18_hic_reference.sh <mcool_dir> <data_dir> <results_dir> [threads]}
data=${2:?}
results=${3:?}
threads=${4:-8}
bin=${HICEXPLORER_BIN:-}
[ -n "$bin" ] && bin="$bin/"
out="$results/compartments_domains/hic_reference"
mkdir -p "$out"
genes="$out/genes_hg38_protein_coding.bed"
if [ ! -s "$genes" ]; then
    zcat "$data/annotation/gencode.v29.annotation.gtf.gz" \
        | awk -F'\t' '$3=="gene" && $9 ~ /gene_type "protein_coding"/ {match($9,/gene_name "[^"]+"/);
              print $1"\t"$4-1"\t"$5"\t"substr($9,RSTART+11,RLENGTH-12)"\t0\t"$7}' \
        | grep -E '^chr([0-9]+|X)\s' | sort -k1,1 -k2,2n > "$genes"
fi
chroms=$(for c in $(seq 1 22) X; do printf 'chr%s ' "$c"; done)
while read -r cell accession; do
    for res in 100000 250000; do
        label=$((res / 1000))kb
        [ -s "$out/${cell}_pc1_${label}.bedgraph" ] && continue
        "${bin}hicPCA" -m "$mcool/$accession.mcool::/resolutions/$res" --whichEigenvectors 1 \
            -o "$out/${cell}_pc1_${label}.bedgraph" --format bedgraph --extraTrack "$genes" \
            --chromosomes $chroms
    done
    if [ ! -s "$out/${cell}_50kb_domains.bed" ]; then
        "${bin}hicFindTADs" -m "$mcool/$accession.mcool::/resolutions/50000" --outPrefix "$out/${cell}_50kb" \
            --minDepth 150000 --maxDepth 500000 --step 50000 --thresholdComparisons 0.01 --delta 0.01 \
            --correctForMultipleTesting fdr --chromosomes $chroms -p "$threads"
    fi
done <<'EOF'
GM12878 4DNFIIV4M7TF
K562 4DNFI18UHVRO
H1 4DNFI9GMP2J8
EOF
