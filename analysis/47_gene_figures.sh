#!/usr/bin/env bash
# Figures of single genes explained only beyond the promoter and its direct contacts: one track panel per level
# (figure_gene_level_tracks.py), the mark change per step against the mirrored loci (figure_gene_report.py) and, for
# NOL11 in GM12878, the multi-way support of the candidate paths (figure_gene_path_multiway.py).
#   Fig. 4 and Fig. 5      NOL11, naive B cells against GM12878 (35_gm12878_naiveB.sh, 36_gm12878_multiway_paths.py)
#   Supplementary Fig. S5  FAM172A, hESC against hNEC (31_step_decomposition.sh)
#   Supplementary Fig. S6  Arhgap26, Md against PA2 (31_step_decomposition.sh)
#   Supplementary Fig. S7  BRI3BP, naive B cells against GM12878
# RNA-seq coverage: rnaseq_coverage.sh. Outputs <out>/<gene>_level_step{1,2,3}.pdf, <gene>_levels.pdf and for NOL11
# NOL11_paths_{matrix,groups}.pdf (manuscript: fig_<gene>_*.pdf, fig_nol11_multiway_{matrix,groups}.pdf).
#
# Usage: 47_gene_figures.sh <data_dir> <results_dir> <out_dir>
set -euo pipefail
data=${1:?usage: 47_gene_figures.sh <data_dir> <results_dir> <out_dir>}
results=${2:?}
out=${3:?}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
py() { PYTHONPATH="$(dirname "$here")" python "$@"; }
roadmap() { printf '%s=%s,%s ' "$1" "$data/roadmap/$2-$1.fc.signal.bigwig" "$data/roadmap/$3-$1.fc.signal.bigwig"; }
kessler() { printf '%s=%s,%s ' "$1" "$data/kessler/E10.5_MD_$1_R1.bigwig" "$data/kessler/E10.5_PA2_$1_R1.bigwig"; }
gtf19="$data/annotation/gencode.v19.annotation.gtf.gz"
gtfm="$data/genomes/gencode.vM25.annotation.gtf.gz"
gm="$results/gm12878_naiveB"
dec="$results/revision_bioinf/decomposed"
mw="$results/multiway"

# NOL11 (Fig. 4, Fig. 5)
node=chr17:65705762-65729696
py "$here/figure_gene_level_tracks.py" "$gm/naiveB_vs_GM12878" "$node" "$gtf19" "$out/nol11" \
    --rna "$gm/rnaseq_coverage/naiveB_RNA-seq.rpm.bw,$gm/rnaseq_coverage/GM12878_RNA-seq.rpm.bw" \
    --tracks $(roadmap H3K27ac E032 E116) $(roadmap H3K4me1 E032 E116) $(roadmap H3K4me3 E032 E116) \
    $(roadmap H3K27me3 E032 E116) --conditions naiveB GM12878 --gene NOL11 --maxStep 3
py "$here/figure_gene_report.py" "$gm/naiveB_vs_GM12878" "$node" H3K27ac \
    "$data/roadmap/E032-H3K27ac.fc.signal.bigwig" "$data/roadmap/E116-H3K27ac.fc.signal.bigwig" "$out/nol11" \
    --name NOL11 --conditions naiveB GM12878 --step 2
py "$here/figure_gene_path_multiway.py" "$gm/naiveB_vs_GM12878" "$gm/multiway_paths" "$node" "$out/NOL11" \
    --data SPRITE="$mw/sprite_GM12878_hg19_HindIII_fragments.npz" HiPore-C="$mw/GM12878_hiporec_hg19_HindIII_fragments.npz" \
    Pore-C="$mw/GM12878_porec_hg19_HindIII_fragments.npz" --name NOL11 --maxPaths 4

# FAM172A (Supplementary Fig. S5)
node=chr5:93432979-93448526
rna="$results/revision_bioinf/freire_rnaseq_coverage"
py "$here/figure_gene_level_tracks.py" "$dec/hESC_vs_hNEC" "$node" "$gtf19" "$out/fam172a" \
    --rna "$rna/hESC_RNA-seq.rpm.bw,$rna/hNEC_RNA-seq.rpm.bw" \
    --tracks $(roadmap H3K4me1 E008 E009) $(roadmap H3K4me3 E008 E009) $(roadmap H3K27me3 E008 E009) \
    --conditions hESC hNEC --gene FAM172A --maxStep 3
py "$here/figure_gene_report.py" "$dec/hESC_vs_hNEC" "$node" H3K4me3 \
    "$data/roadmap/E008-H3K4me3.fc.signal.bigwig" "$data/roadmap/E009-H3K4me3.fc.signal.bigwig" "$out/fam172a" \
    --name FAM172A --conditions hESC hNEC --step 2

# Arhgap26 (Supplementary Fig. S6)
node=chr18:38595961-38603414
rna="$results/revision_bioinf/minoux_rnaseq_coverage"
py "$here/figure_gene_level_tracks.py" "$dec/Md_vs_PA2" "$node" "$gtfm" "$out/arhgap26" \
    --rna "$rna/E10.5_Md_RNA-seq.rpm.bw,$rna/E10.5_PA2_RNA-seq.rpm.bw" \
    --tracks $(kessler H3K27ac) $(kessler H3K4me2) $(kessler H3K27me3) --conditions Md PA2 --gene Arhgap26 --maxStep 3
py "$here/figure_gene_report.py" "$dec/Md_vs_PA2" "$node" H3K27me3 \
    "$data/kessler/E10.5_MD_H3K27me3_R1.bigwig" "$data/kessler/E10.5_PA2_H3K27me3_R1.bigwig" "$out/arhgap26" \
    --name Arhgap26 --conditions Md PA2 --step 3

# BRI3BP (Supplementary Fig. S7)
node=chr12:125476355-125499056
py "$here/figure_gene_level_tracks.py" "$gm/naiveB_vs_GM12878" "$node" "$gtf19" "$out/bri3bp" \
    --rna "$gm/rnaseq_coverage/naiveB_RNA-seq.rpm.bw,$gm/rnaseq_coverage/GM12878_RNA-seq.rpm.bw" \
    --tracks $(roadmap H3K27ac E032 E116) $(roadmap H3K4me1 E032 E116) $(roadmap H3K4me3 E032 E116) \
    $(roadmap H3K27me3 E032 E116) --conditions naiveB GM12878 --gene BRI3BP --maxStep 3
py "$here/figure_gene_report.py" "$gm/naiveB_vs_GM12878" "$node" H3K27ac \
    "$data/roadmap/E032-H3K27ac.fc.signal.bigwig" "$data/roadmap/E116-H3K27ac.fc.signal.bigwig" "$out/bri3bp" \
    --name BRI3BP --conditions naiveB GM12878 --step 2
