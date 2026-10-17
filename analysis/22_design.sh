#!/usr/bin/env bash
# Experiment design (egaDesign) for the four differential runs of 16_hops_n.sh: every node of every explaining level
# of every explained gene is perturbed in silico (flatten: its mark change set to 0; silence: its marks set to 0;
# remove: node and edges deleted, for the ten strongest candidates of a gene) and the level statistic recomputed with
# the egaDifferential formulas; per-gene P before and after with 1000 permutations. Candidates ranked by the
# flatten effect and their evidence; BED windows of the three best nodes per gene from the 1 kb window with the
# highest signal of the explaining mark in the condition where it is higher.
# Validation (design_crispr_validation.py): single-condition ranking of the K562 graph (DNase) against the K562
# CRISPRi pairs of the ENCODE-rE2G training set (the crispr_benchmark.py pairs).
# Outputs: <results>/design/<comparison>_design*.{tsv,txt,bed,json}, <results>/design/validation/,
# <results>/design/figures/*.pdf, <results>/design/tables/table_s34.tex.
set -euo pipefail
data=${1:?usage: 22_design.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> [hops_n dir] [threads]}
results=${2:?}
kasumi=${3:?}
assi=${4:?}
hops=${5:-$results/hops_n}
threads=${6:-8}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/design"
mkdir -p "$out"
common=(--verbose --roles declared --presence mean --chromatinThreshold 1 --effect flatten --nRemove 10
        --permutations 1000 --strong 0.5 --top 10 --bedTop 3 --window 1000 --threads "$threads" --seed 0)

# hESC / hNEC: Roadmap E008 (hESC) and E009 (hNEC) fold enrichment tracks
roadmap() { for mark in H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }
egaDesign "${common[@]}" -p "$hops/hESC_vs_hNEC" --tss "$data/annotation/gencode.v19.annotation.gtf.gz" \
    --tracksA $(roadmap E008) --tracksB $(roadmap E009) -o "$out/hESC_vs_hNEC" > "$out/hESC_vs_hNEC.log" 2>&1

# Md / PA2
k="$data/kessler"
kessler() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$k/E10.5_$1_${mark}_R1.bigwig"; done; }
egaDesign "${common[@]}" -p "$hops/Md_vs_PA2" --tss "$k/mm10.ncbiRefSeq.gtf" \
    --tracksA $(kessler MD) --tracksB $(kessler PA2) -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1

# Kasumi-1 siMM / siRE
bw="$kasumi/dnase/bigwig"
egaDesign "${common[@]}" -p "$hops/siMM_vs_siRE" --tss "$data/annotation/gencode.v29.annotation.gtf.gz" \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1

# AML t(8;21) / FLT3-ITD
bw="$assi/bigwig"
join() { local IFS=,; echo "$*"; }
t821_dhs=$(join "$bw"/GSM289362[789]_t821*_DHS.bw "$bw"/GSM2893630_t821KIT-2_DHS.bw)
itd_dhs=$(join "$bw"/GSM289360[67]_ITD-*_DHS.bw)
egaDesign "${common[@]}" -p "$hops/t821_vs_FLT3ITD" --tss "$data/annotation/gencode.v29.annotation.gtf.gz" \
    --tracksA "DHS=$t821_dhs" --tracksB "DHS=$itd_dhs" -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1

# K562 CRISPRi validation: single-condition ranking on the K562 DMSO graph with ENCODE K562 DNase
python "$here/design_crispr_validation.py" "$results/k562/K562_dmso.graphml" \
    "$data/k562/crispr/github/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz" \
    "$data/accessibility/K562_DNase-seq.GRCh38.ENCFF413AHU.bigWig" "$out/validation/K562" "$threads" \
    > "$out/validation_K562.log" 2>&1

python "$here/figure_design.py" "$out" "$out/figures"
python "$here/table_design.py" "$out" "$out/tables"
