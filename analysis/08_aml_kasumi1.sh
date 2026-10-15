#!/usr/bin/env bash
# Differential mode on Kasumi-1 (t(8;21) AML) after siRNA knockdown of RUNX1-ETO (siRE) against mismatch control
# (siMM), hg38 (Ptasinska et al. 2014, 2019):
#   capture Hi-C: HiCUP per replicate (kasumi1/map_runs.sh), bam2chicago and the CHiCAGO pipeline of HiCExplorer
#   (kasumi1/call_chicago.sh) on a design of the HindIII fragments holding a TSS (analysis/promoter_baitmap.py),
#   exported at score >= 1 and combined by the replicate rule (present in both replicates, score >= 5 in one)
#   chromatin: DNase-seq of siMM and siRE at day 2 and day 10 after transfection (deposited on hg19, lifted to hg38)
#   expression: RNA-seq day 4, three replicates, salmon and PyDESeq2 (analysis/kasumi_rna.sh)
set -euo pipefail
kasumi=${1:?usage: 08_aml_kasumi1.sh <kasumi1 dir> <data_dir> <results_dir> [permutations]}
data=${2:?}
results=${3:?}
permutations=${4:-200}
gtf="$data/annotation/gencode.v29.annotation.gtf.gz"
sizes="$data/annotation/hg38.chrom.sizes"
out="$results/aml_kasumi1"
bw="$kasumi/dnase/bigwig"
mkdir -p "$out" "$bw"
# the deposited DNase-seq tracks are on hg19 (GAPDH promoter signal at the hg19 position): lifted to hg38
chain="$data/annotation/hg19ToHg38.over.chain.gz"
for f in "$kasumi"/dnase/*.wig.gz; do
    output="$bw/$(basename "$f" .wig.gz).bw"
    [[ -s "$output" ]] && continue
    python "$(dirname "$0")/liftover_wig.py" "$f" "$chain" "$sizes" "$output.bedGraph"
done
for f in "$kasumi"/dnase/*.wig.gz; do
    output="$bw/$(basename "$f" .wig.gz).bw"
    [[ -s "$output" ]] && continue
    bedGraphToBigWig "$output.bedGraph" "$sizes" "$output"
    rm "$output.bedGraph"
done
calls="$kasumi/chicago"
egaDifferential --verbose \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 5 --names siMM siRE --tss "$gtf" \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    --expression "$out/deseq2_siRE_vs_siMM.tsv" --log2FCColumn log2FoldChange --padjColumn padj --alpha 0.05 \
    --geneKey genes --rnaThreshold 1 --chromatinThreshold 1 --permutations "$permutations" --bootstrap 1000 \
    --seed 0 -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1
