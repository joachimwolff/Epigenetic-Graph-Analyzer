#!/usr/bin/env bash
# Candidate trigger transcription factors (egaTriggers) of the genes explained by the differential runs of
# 16_hops_n.sh: motifs of JASPAR 2024 CORE vertebrates (non-redundant, 879 matrices) at the explaining nodes and
# promoter fragments, enrichment per gene group against GC- and length-matched background nodes, and per gene the
# motif-carrying factors whose own expression changes. The scanned sequence of a fragment is its 1 kb window with
# the highest active-mark signal (the active / enhancer marks of both conditions; Kasumi-1 and AML: DNase / DHS).
# Genomes are decompressed once to <results>/triggers/genomes (pyfaidx index next to them).
# Outputs: <results>/triggers/<comparison>_{motif_enrichment,gene_triggers}.tsv, figures, table_s31.tex.
set -euo pipefail
data=${1:?usage: 18_triggers.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> [hops_n dir] [threads]}
results=${2:?}
kasumi=${3:?}
assi=${4:?}
hops=${5:-$results/hops_n}
threads=${6:-8}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/triggers"
genomes="$out/genomes"
mkdir -p "$genomes"
[[ -s "$genomes/hg19.fa" ]] || zcat "$data/annotation/hg19.fa.gz" > "$genomes/hg19.fa"
[[ -s "$genomes/hg38.fa" ]] || zcat "$data/genomes/hg38.fa.gz" > "$genomes/hg38.fa"
[[ -s "$genomes/mm10.fa" ]] || zcat "$data/genomes/mm10.fa.gz" > "$genomes/mm10.fa"
motifs="$data/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
common=(--verbose --motifs "$motifs" --window 1000 --pvalue 1e-4 --alpha 0.05 --maxHops 3 --chromatinThreshold 1
        --roles declared --presence mean --backgroundRatio 3 --topPerGene 25 --threads "$threads" --seed 0 --cache)

# hESC / hNEC: H3K4me3 (active) and H3K4me1 (enhancer) of E008 and E009
egaTriggers "${common[@]}" -p "$hops/hESC_vs_hNEC" --genome "$genomes/hg19.fa" \
    --windowTracks "$data"/roadmap/E00{8,9}-H3K4me{1,3}.fc.signal.bigwig \
    --expression "$results/differential/deseq2_hNEC_vs_hESC.tsv" \
    -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1

# Md / PA2: H3K27ac and H3K4me2 of both tissues
k="$data/kessler"
egaTriggers "${common[@]}" -p "$hops/Md_vs_PA2" --genome "$genomes/mm10.fa" \
    --windowTracks "$k"/E10.5_{MD,PA2}_{H3K27ac,H3K4me2}_R1.bigwig \
    --expression "$results/kessler/deseq2_PA2_vs_Md.tsv" \
    -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1

# Kasumi-1 siMM / siRE: DNase day 2 and day 10 of both conditions
egaTriggers "${common[@]}" -p "$hops/siMM_vs_siRE" --genome "$genomes/hg38.fa" \
    --windowTracks "$kasumi"/dnase/bigwig/GSM343076{2,3,4,5}_DHS_kasumi1_*.bw \
    --expression "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" \
    -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1

# AML t(8;21) / FLT3-ITD: DHS of all patients; expression P as in 16_hops_n.sh (three patients per group)
egaTriggers "${common[@]}" -p "$hops/t821_vs_FLT3ITD" --genome "$genomes/hg38.fa" \
    --windowTracks "$assi"/bigwig/GSM28936{06,07,27,28,29,30}_*DHS.bw \
    --expression "$results/aml_assi/expression_FLT3ITD_vs_t821.tsv" --padjColumn pvalue --baseColumn mean_log2_B \
    -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1

python "$here/figure_triggers.py" "$out" "$out/figures"
python "$here/table_triggers.py" "$out" "$out/tables"
