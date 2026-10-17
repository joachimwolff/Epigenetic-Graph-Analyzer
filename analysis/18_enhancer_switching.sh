#!/usr/bin/env bash
# Enhancer switching and promoter competition (egaSwitching) for the four comparisons of 16_hops_n.sh.
# Input per comparison: <hops_n>/<comparison>_union.graphml (and _elements.tsv for the element classes). If the
# egaDifferential run in <hops_n> is incomplete (union graph missing or older than the element table), the copy in
# <hops_n>/before_elements is used and the log says so.
# Output: <results>/enhancer_switching/<comparison>_switch_events.tsv, _switch_summary.tsv, _competition.tsv,
# _switch_hubs.tsv, _switch_biology.tsv, figures/ (arc plots of the top 3 events per condition, null distributions),
# table_s29.tex (analysis/table_switching.py).
# Known biology is checked, not forced: anterior Hoxa / Hoxb genes in PA2 (Md is Hox-negative), neural and
# pluripotency genes in hNEC, myeloid genes after RUNX1-ETO knockdown (siRE) and in FLT3-ITD vs t(8;21) AML.
set -euo pipefail
results=${1:?usage: 18_enhancer_switching.sh <results_dir> [permutations] [parallel comparisons]}
permutations=${2:-10000}
parallel=${3:-4}
here=$(cd "$(dirname "$0")" && pwd)
hops="$results/hops_n"
out="$results/enhancer_switching"
mkdir -p "$out/figures"
# numerical libraries: at most 4 threads per comparison (4 comparisons in parallel = 16)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4

input_of() {
    local c=$1
    if [[ -f "$hops/${c}_elements.tsv" && -f "$hops/${c}_union.graphml" && "$hops/${c}_union.graphml" -nt "$hops/${c}_elements.tsv" ]]; then
        echo "$hops/$c"
    else
        echo "$hops/before_elements/$c"
    fi
}

hoxa_anterior="Hoxa1,Hoxa2,Hoxa3,Hotairm1"
hoxa_posterior="Hoxa4,Hoxa5,Hoxa6,Hoxa7,Hoxa9,Hoxa10,Hoxa11,Hoxa13"
hoxb_anterior="Hoxb1,Hoxb2,Hoxb3"
neural="PAX6,SOX1,SOX3,NES,OTX2,LHX2,FOXG1,ZIC1,ZIC2,ZIC3,NR2F1,NR2F2,HES5,SIX3,POU3F2,DACH1,LMO3,ZBTB16,FABP7,MSI1,EMX2,SOX21,DLK1"
pluripotency="POU5F1,NANOG,LIN28A,DPPA4,PRDM14,ZFP42,SALL4,DNMT3B,L1TD1,TDGF1,LEFTY1,LEFTY2"
myeloid="CEBPA,CEBPB,CEBPE,SPI1,CSF1R,CD14,LYZ,MPO,ELANE,ITGAM,CSF3R,IRF8,KLF4,FCGR1A,S100A8,S100A9,CTSG,AZU1,PRTN3"
stem="CD34,KIT,LMO2,MEIS1,ERG"

declare -A expression level sets
expression[hESC_vs_hNEC]="$results/differential/deseq2_hNEC_vs_hESC.tsv"
expression[Md_vs_PA2]="$results/kessler/deseq2_PA2_vs_Md.tsv"
expression[siMM_vs_siRE]="$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv"
expression[t821_vs_FLT3ITD]="$results/aml_assi/expression_FLT3ITD_vs_t821.tsv"
sets[hESC_vs_hNEC]="neural=$neural:up pluripotency=$pluripotency:down"
sets[Md_vs_PA2]="Hoxa_anterior=$hoxa_anterior:up Hoxa_posterior=$hoxa_posterior:any Hoxb_anterior=$hoxb_anterior:up"
sets[siMM_vs_siRE]="myeloid=$myeloid:up stem=$stem:down"
sets[t821_vs_FLT3ITD]="myeloid=$myeloid:any stem=$stem:any"

run_one() {
    local c=$1 input
    input=$(input_of "$c")
    echo "input: $input" > "$out/$c.log"
    # shellcheck disable=SC2086
    python -m epigeneticgraphanalyzer.switching --verbose -p "$input" -o "$out/$c" \
        --expression "${expression[$c]}" --permutations "$permutations" --seed 0 \
        --geneSets ${sets[$c]} >> "$out/$c.log" 2>&1
}

pids=()
for c in hESC_vs_hNEC Md_vs_PA2 siMM_vs_siRE t821_vs_FLT3ITD; do
    run_one "$c" &
    pids+=($!)
    while (( $(jobs -rp | wc -l) >= parallel )); do wait -n; done
done
for pid in "${pids[@]}"; do wait "$pid"; done
python "$here/table_switching.py" "$out" "$out"
