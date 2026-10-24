#!/usr/bin/env bash
# Replication of neural induction with independent data: H1 hESC versus H1-derived neural progenitors (NPC), with the
# settings of 31_step_decomposition.sh, for the comparison with hESC / hNEC (Freire-Pritchett et al. 2017, H9).
#   promoter capture Hi-C: Jung et al. 2019 (GSE86189), H1 subsampled to the captured reads of the NPC library
#     (38b_h1_matched_depth.sh) and NPC (38_npc_chicago.sh), both CHiCAGO on the same design, score >= 5
#   histone marks: Roadmap consolidated fold change over input, E003 (H1) and E007 (H1-derived NPC), H3K27ac,
#     H3K4me1, H3K4me3, H3K27me3
#   expression: RNA-seq of Xie et al. 2013 (GSE16256), two Illumina replicates per cell type, STAR gene counts
#     (rnaseq_coverage.sh) and PyDESeq2 (star_counts_deseq2.py)
# A second run uses the strongest H1 interactions in the number of NPC calls (--matchCalls), as a sensitivity check.
#
# Usage: 40_h1_npc.sh <data_dir> <H1 matched ibed> <NPC ibed> <RNA counts dir> <out_dir> [gene permutations] [threads]
set -euo pipefail
data=${1:?usage: 40_h1_npc.sh <data_dir> <H1 ibed> <NPC ibed> <counts dir> <out_dir> [gene permutations] [threads]}
h1=${2:?}
npc=${3:?}
counts=${4:?}
out=${5:?}
gene_permutations=${6:-200000}
threads=${7:-6}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
python "$here/star_counts_deseq2.py" "$data/h1_npc_rna/fastq/samples.tsv" "$counts" \
    "$data/annotation/gencode.v19.annotation.gtf.gz" H1 NPC "$out/deseq2_NPC_vs_H1.tsv"
# sensitivity: the strongest H1 calls in the number of NPC calls
n_npc=$(($(wc -l < "$npc") - 1))
{ head -1 "$h1"; tail -n +2 "$h1" | sort -t$'\t' -k10,10gr | head -n "$n_npc"; } > "$out/H1_matched_top.ibed"
roadmap() { for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
for run in matched top; do
    calls=$h1
    [[ $run == top ]] && calls="$out/H1_matched_top.ibed"
    egaDifferential --verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 \
        --chromatinThreshold 1 --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
        --stepFeature decomposed \
        --interactionsA "$calls" --interactionsB "$npc" --format ibed --scoreThreshold 5 --referenceScore 5 \
        --names H1 NPC --tss "$data/annotation/gencode.v19.annotation.gtf.gz" \
        --tracksA $(roadmap E003) --tracksB $(roadmap E007) \
        --markSign H3K27me3=-1 --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
        --expression "$out/deseq2_NPC_vs_H1.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
        --alpha 0.05 --geneKey genes -o "$out/H1_vs_NPC_$run" > "$out/H1_vs_NPC_$run.log" 2>&1
done
