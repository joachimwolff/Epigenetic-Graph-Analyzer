#!/usr/bin/env bash
# Differential mode between naive B cells and the B-lymphoblastoid line GM12878 (EBV-transformed B cells), with the
# settings of 31_step_decomposition.sh.
#   promoter capture Hi-C: naive B cells of Javierre et al. 2016 (peak matrix column nB, javierre_to_washu.py) and
#     GM12878 of Cairns et al. 2016 (GSE81503); both hg19, HindIII, the same capture design; CHiCAGO score >= 5
#   histone marks: Roadmap consolidated fold change over input, E032 (B cells from peripheral blood) and E116
#     (GM12878), H3K27ac, H3K4me1, H3K4me3, H3K27me3
#   expression: ENCODE RSEM gene quantifications, three naive B donors and two GM12878 replicates (encode_deseq2.py)
# The multi-way data of GM12878 (SPRITE, Pore-C, HiPore-C on the HindIII fragments) test the explaining paths in
# 36_gm12878_multiway_paths.py.
set -euo pipefail
data=${1:?usage: 35_gm12878_naiveB.sh <data_dir> <out_dir> [gene permutations] [threads]}
out=${2:?}
gene_permutations=${3:-200000}
threads=${4:-6}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
nb="$data/naiveB"
[[ -s "$nb/javierre2016/naiveB_nB_washU.txt.gz" ]] || python "$here/javierre_to_washu.py" \
    "$nb/javierre2016/PCHiC_peak_matrix_cutoff5.txt.gz" nB "$nb/javierre2016/naiveB_nB_washU.txt.gz"
python "$here/encode_deseq2.py" "$nb/encode_rna/files.tsv" "$data/annotation/gencode.v19.annotation.gtf.gz" \
    naiveB GM12878 "$out/deseq2_GM12878_vs_naiveB.tsv"
roadmap() { for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
egaDifferential --verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 --chromatinThreshold 1 \
    --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --stepFeature decomposed \
    --interactionsA "$nb/javierre2016/naiveB_nB_washU.txt.gz" \
    --interactionsB "$data/chicago_calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz" \
    --format washu --baits "$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 5 --referenceScore 5 --names naiveB GM12878 \
    --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E032) --tracksB $(roadmap E116) \
    --markSign H3K27me3=-1 --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
    --expression "$out/deseq2_GM12878_vs_naiveB.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    -o "$out/naiveB_vs_GM12878" > "$out/naiveB_vs_GM12878.log" 2>&1
egaExplain -p "$out/naiveB_vs_GM12878" -o "$out/explain/naiveB_vs_GM12878" > "$out/explain_naiveB_vs_GM12878.log" 2>&1
