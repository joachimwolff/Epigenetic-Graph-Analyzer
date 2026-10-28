#!/usr/bin/env bash
# Inputs of the case study naive B cells against GM12878 (docs/case_study.md), all public:
#   GM12878 CHiCAGO calls (Cairns et al. 2016, GEO GSE81503)
#   hg19 HindIII promoter capture design (Javierre et al. 2016, OSF osf.io/u8tzp)
#   naive B cell calls and the DESeq2 table, as archived with EGA (Zenodo 10.5281/zenodo.23725671, record 03)
#   Roadmap Epigenomics fold change over input: E032 (B cells) and E116 (GM12878), four histone modifications
#   GENCODE v19 gene annotation
#   GM12878 SPRITE clusters (Quinodoz et al. 2018, GEO GSE114242; 7.2 GB)
# About 14 GB in total.
#
# Usage: download.sh <data dir>
set -euo pipefail
data=${1:?usage: download.sh <data dir>}
mkdir -p "$data/roadmap" "$data/sprite"
cd "$data"
get() { [[ -s "$2" ]] || curl -L --fail --retry 3 -o "$2" "$1"; }
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81503/suppl/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz \
    GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz
get https://osf.io/download/e594p/ human_PCHiC_hg19_HindIII_design.tar.gz
tar -xzf human_PCHiC_hg19_HindIII_design.tar.gz Human_hg19/Digest_Human_HindIII.rmap \
    Human_hg19/Digest_Human_HindIII_baits_e75_ID.baitmap
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
    gencode.v19.annotation.gtf.gz
get https://zenodo.org/records/23725671/files/03_processed_inputs.tar.gz 03_processed_inputs.tar.gz
tar -xzf 03_processed_inputs.tar.gz 03_processed_inputs/naiveB_javierre2016/naiveB_nB_washU.txt.gz \
    03_processed_inputs/expression/deseq2_GM12878_vs_naiveB.tsv
roadmap=https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange
for epigenome in E032 E116; do
    for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do
        get "$roadmap/$epigenome-$mark.fc.signal.bigwig" "roadmap/$epigenome-$mark.fc.signal.bigwig"
    done
done
get https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114242/suppl/GSE114242_human_combined_clusters.tar.gz \
    sprite/GSE114242_human_combined_clusters.tar.gz
tar -xzf sprite/GSE114242_human_combined_clusters.tar.gz -C sprite human.combined.mapq-ge30.clusters
echo "inputs in $data"
