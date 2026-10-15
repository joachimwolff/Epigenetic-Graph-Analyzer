#!/usr/bin/env bash
# Download every public dataset used in the EGA manuscript analyses.
#
# Usage: scripts/download_data.sh <target_dir>
#
# Layout written below <target_dir>:
#   jung2019/     promoter capture Hi-C (Jung et al. 2019, GSE86189), hg19, raw counts
#   mifsud2015/   GM12878 promoter capture Hi-C (Mifsud et al. 2015), hg19, significant interactions
#   sprite/       GM12878 DNA SPRITE clusters (Quinodoz et al. 2018, GSE114242), hg19
#   gam_h1/       H1-hESC GAM segregation tables (4DN), GRCh38
#   roadmap/      Roadmap Epigenomics consolidated fold-change bigWigs + RNA-seq RPKM, hg19
#   encode_gm12878/ ENCODE Broad histone / CTCF signal bigWigs for GM12878, hg19
#   annotation/   GENCODE v19 annotation, hg19 sequence and chromosome sizes, hg19->hg38 chain
#   fragments/    HindIII restriction fragments of hg19
set -euo pipefail

target=${1:?usage: download_data.sh <target_dir>}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$target"
cd "$target"

fetch() {
    # fetch <url> <output>; resumes partial downloads and skips finished ones
    local url=$1 out=$2
    if [[ -s "$out" && ! -e "$out.part" ]]; then
        return 0
    fi
    touch "$out.part"
    curl -sSfL -C - -o "$out" "$url"
    rm -f "$out.part"
}

# --- promoter capture Hi-C, Jung et al. 2019 ---------------------------------
mkdir -p jung2019
geo_jung=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86189/suppl
for sample in H1 GM ME npc TB msc; do
    for kind in po pp; do
        fetch "$geo_jung/GSE86189_${sample}.${kind}.all.txt.bz2" "jung2019/GSE86189_${sample}.${kind}.all.txt.bz2"
    done
done

# --- GM12878 promoter capture Hi-C, Mifsud et al. 2015 (ArrayExpress E-MTAB-2323)
mkdir -p mifsud2015
fetch https://www.ebi.ac.uk/biostudies/files/E-MTAB-2323/TS5_GM12878_promoter-other_significant_interactions.txt \
    mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt

# --- GM12878 SPRITE, Quinodoz et al. 2018 -------------------------------------
mkdir -p sprite
fetch https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114242/suppl/GSE114242_human_combined_clusters.tar.gz \
    sprite/GSE114242_human_combined_clusters.tar.gz

# --- H1-hESC GAM segregation tables, 4DN --------------------------------------
mkdir -p gam_h1
for entry in 4DNFI9FFNSMA:rep1_30kb 4DNFI6Z8MRTX:rep2_30kb \
             4DNFI8VQTXSR:rep1_40kb 4DNFIZAA4EX7:rep2_40kb \
             4DNFIULN6A82:rep1_50kb 4DNFIKYW3KLO:rep2_50kb; do
    acc=${entry%%:*}
    name=${entry##*:}
    url=$(curl -sSf -H 'Accept: application/json' \
        "https://data.4dnucleome.org/files-processed/${acc}/?format=json" |
        python3 -c 'import json,sys; print(json.load(sys.stdin)["open_data_url"])')
    fetch "$url" "gam_h1/H1_GAM_${name}.${acc}.csv"
done

# --- Roadmap Epigenomics: H1 (E003) and H1-derived lineages -------------------
# E004 mesendoderm, E005 trophoblast, E006 mesenchymal stem cells, E007 neural progenitors
mkdir -p roadmap
roadmap_fc=https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange
for eid in E003 E004 E005 E006 E007; do
    for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do
        fetch "$roadmap_fc/${eid}-${mark}.fc.signal.bigwig" "roadmap/${eid}-${mark}.fc.signal.bigwig"
    done
done
fetch https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.RPKM.pc.gz roadmap/57epigenomes.RPKM.pc.gz
fetch https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/Ensembl_v65.Gencode_v10.ENSG.gene_info roadmap/Ensembl_v65.Gencode_v10.ENSG.gene_info

# --- ENCODE GM12878 (Broad histone + CTCF), hg19 ------------------------------
mkdir -p encode_gm12878
encode_broad=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone
for track in H3k27ac H3k4me1 H3k4me3 H3k27me3 Ctcf; do
    fetch "$encode_broad/wgEncodeBroadHistoneGm12878${track}StdSig.bigWig" "encode_gm12878/GM12878_${track}.bigWig"
done

# --- annotation ---------------------------------------------------------------
mkdir -p annotation
fetch https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz annotation/gencode.v19.annotation.gtf.gz
fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes annotation/hg19.chrom.sizes
fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes annotation/hg38.chrom.sizes
fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz annotation/hg19ToHg38.over.chain.gz
fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz annotation/hg19.fa.gz

# --- HindIII restriction fragments of hg19 (the fragments of the Mifsud et al. capture design) ----------
mkdir -p fragments
if [[ ! -s fragments/hg19_HindIII_fragments.bed ]]; then
    python "$here/digest_genome.py" annotation/hg19.fa.gz AAGCTT 1 fragments/hg19_HindIII_fragments.bed \
        '^chr([0-9]+|X)$'
fi

echo "all downloads complete"
