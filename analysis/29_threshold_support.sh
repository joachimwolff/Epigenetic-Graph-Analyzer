#!/usr/bin/env bash
# Threshold sensitivity of the gained and lost contacts and provenance of the statements (egaDifferential support
# classes, egaExplain / egaReport / egaExplorer flags and tags, egaReport plain overlay view), on two data sets with
# different inputs:
#   Kasumi-1 siMM / siRE: own CHiCAGO exports down to score 1 (08_aml_kasumi1.sh), so a gained / lost pair can be
#     looked up below the threshold in the other condition: stable (score < floor) or threshold-sensitive
#     (floor <= score < 5, or called in one replicate only)
#   Md / PA2: published CHiCAGO calls (score >= 5) of two replicates per tissue (05_kessler_md_pa2.sh): presence in
#     one replicate of the other tissue is informative (threshold-sensitive), absence from both is not testable
# The settings of 16_hops_n.sh (maxHops 3, both role sets, both presence statistics); support floor 3 and calling
# score 5 (--supportFloor, --supportThreshold). Outputs: <results>/flags_demo/ (Md_vs_PA2_*, siMM_vs_siRE_*),
# explain/, reports/ (Hoxa2, CCNA1 with the overlay panels), explorer/, support_summary.tsv and
# support_genes.tsv (analysis/support_summary.py).
# Stage: all (default), differential (egaDifferential only), displays (egaExplain, egaReport, egaExplorer, summary).
set -euo pipefail
data=${1:?usage: 29_threshold_support.sh <data_dir> <results_dir> <kasumi1 dir> [stage] [threads]}
results=${2:?}
kasumi=${3:?}
stage=${4:-all}
threads=${5:-16}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=${PYTHON:-python}
out="$results/flags_demo"
mkdir -p "$out"/{explain,reports,explorer}
common=(--verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads"
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000
        --supportThreshold 5 --supportFloor 3)
k="$data/kessler"
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
if [ "$stage" = all ] || [ "$stage" = differential ]; then
    # Md / PA2: expression of kessler_deseq2.py (repeated GEO libraries dropped)
    "$py" "$here/kessler_deseq2.py" "$k/GSE89434_Minoux_RNAseq_gene_counts.csv.gz" "$out/deseq2_PA2_vs_Md.tsv" \
        > "$out/deseq2.log" 2>&1
    kessler() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
    fold() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$results/minoux_fold/$1_${mark}.fold.bw"; done; }
    "$py" -m epigeneticgraphanalyzer.egaDifferential "${common[@]}" \
        --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
        --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
        --format washu --replicates intersect --names Md PA2 \
        --tss "$k/mm10.ncbiRefSeq.gtf" --geneTypes \
        --tracksA $(kessler MD) --tracksB $(kessler PA2) --markSign H3K27me3=-1 \
        --markRoles H3K27ac=active H3K4me2=active H3K27me3=repressive \
        --stateTracksA $(fold Md) --stateTracksB $(fold PA2) \
        --expression "$out/deseq2_PA2_vs_Md.tsv" --log2FCColumn log2FoldChange \
        --padjColumn padj --alpha 0.05 --geneKey genes \
        -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1 &
    "$py" -m epigeneticgraphanalyzer.egaDifferential "${common[@]}" \
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
        -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1 &
    wait
fi
if [ "$stage" = all ] || [ "$stage" = displays ]; then
    for comparison in Md_vs_PA2 siMM_vs_siRE; do
        "$py" -m epigeneticgraphanalyzer.egaExplain -p "$out/$comparison" -o "$out/explain/$comparison" \
            > "$out/explain_$comparison.log" 2>&1
    done
    # TAD boundaries for the boundary crossings of the overlay view: hicFindTADs of PA2 (28_tad_crossing.sh) and of
    # K562 (24_compartments_domains.sh, the closest myeloid Hi-C), boundaries at or above the median strength
    md_tads=(--tads "$results/tad_crossing/hicfindtads/PA2_25kb" --tadMinStrength 0.5)
    k_tads=(--tads "$results/compartments_domains/hic_reference/K562_50kb" --tadMinStrength 0.5)
    "$py" -m epigeneticgraphanalyzer.egaReport -p "$out/Md_vs_PA2" --gene Hoxa2 --overlay "${md_tads[@]}" \
        --annotation "$k/mm10.ncbiRefSeq.gtf" -o "$out/reports" > "$out/reports/Md_vs_PA2.log" 2>&1
    "$py" -m epigeneticgraphanalyzer.egaReport -p "$out/siMM_vs_siRE" --gene CCNA1 --overlay "${k_tads[@]}" \
        --annotation "$data/annotation/gencode.v29.annotation.gtf.gz" -o "$out/reports" \
        > "$out/reports/siMM_vs_siRE.log" 2>&1
    "$py" -m epigeneticgraphanalyzer.explorer -p "$out/Md_vs_PA2" "${md_tads[@]}" --annotation "$k/mm10.ncbiRefSeq.gtf" \
        -o "$out/explorer/" > "$out/explorer/Md_vs_PA2.log" 2>&1
    "$py" -m epigeneticgraphanalyzer.explorer -p "$out/siMM_vs_siRE" "${k_tads[@]}" \
        --annotation "$data/annotation/gencode.v29.annotation.gtf.gz" -o "$out/explorer/" \
        > "$out/explorer/siMM_vs_siRE.log" 2>&1
    "$py" "$here/support_summary.py" "$out" Md_vs_PA2 siMM_vs_siRE --gene siMM_vs_siRE=CCNA1 Md_vs_PA2=Hoxa2 \
        > "$out/support_summary.log" 2>&1
fi
