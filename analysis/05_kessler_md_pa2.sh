#!/usr/bin/env bash
# Differential mode on mouse E10.5 cranial neural crest: mandibular arch (Md) versus second pharyngeal arch (PA2).
#   capture Hi-C: Kessler et al. 2023 (GSE211902), CHiCAGO calls of two replicates per tissue
#   RNA-seq counts and histone marks: Minoux et al. 2017 (GSE89434, GSE89435)
set -euo pipefail
data=${1:?usage: 05_kessler_md_pa2.sh <data_dir> <results_dir>}
results=${2:?}
here=$(cd "$(dirname "$0")" && pwd)
k="$data/kessler"
mkdir -p "$k" "$results/kessler"
geo=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6505nnn
for f in GSM6505208/suppl/GSM6505208_Md_CapHiC_Merged_All-Tr_2ndBioRep_Chicago_output_washU_text.txt.gz \
         GSM6505209/suppl/GSM6505209_Md_CapHiC_3rdBioRep_Chicago_output_washU_text.txt.gz \
         GSM6505210/suppl/GSM6505210_PA2_CapHiC_2ndBioRep_Chicago_output_washU_text.txt.gz \
         GSM6505211/suppl/GSM6505211_PA2_E10.5wt_CapHiC_3rdBioRep_Tr1_Chicago_output_washU_text.txt.gz; do
    [[ -s "$k/$(basename "$f")" ]] || curl -sSf -o "$k/$(basename "$f")" "$geo/$f"
done
counts="$k/GSE89434_Minoux_RNAseq_gene_counts.csv.gz"
[[ -s "$counts" ]] || curl -sSf -o "$counts" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89434/suppl/GSE89434_Minoux_RNAseq_gene_counts.csv.gz
# histone mark bigWigs of E10.5 Md and PA2 (Minoux et al. 2017) and the mm10 RefSeq annotation are expected in $k:
#   E10.5_{MD,PA2}_{H3K27ac,H3K4me2,H3K27me3}_R1.bigwig, mm10.ncbiRefSeq.gtf
python "$here/kessler_deseq2.py" "$counts" "$results/kessler/deseq2_PA2_vs_Md.tsv"
tracks() {
    for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done
}
egaDifferential --verbose --genePermutations 200000 \
    --interactionsA "$k"/GSM6505208_*.txt.gz "$k"/GSM6505209_*.txt.gz \
    --interactionsB "$k"/GSM6505210_*.txt.gz "$k"/GSM6505211_*.txt.gz \
    --format washu --replicates intersect --names Md PA2 \
    --tss "$k/mm10.ncbiRefSeq.gtf" --geneTypes \
    --tracksA $(tracks MD) --tracksB $(tracks PA2) --markSign H3K27me3=-1 \
    --expression "$results/kessler/deseq2_PA2_vs_Md.tsv" --log2FCColumn log2FoldChange \
    --padjColumn padj --alpha 0.05 --geneKey genes --rnaThreshold 1 --chromatinThreshold 1 \
    --permutations 200 --bootstrap 1000 --seed 0 -o "$results/kessler/Md_vs_PA2" \
    > "$results/kessler/Md_vs_PA2.log" 2>&1
