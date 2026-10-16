#!/usr/bin/env bash
# Chicdiff on Kasumi-1 siMM vs siRE (two replicates each) and the overlap of its calls with the gained / lost edges
# of 08_aml_kasumi1.sh. Chicdiff 0.6 (github.com/RegulatoryGenomicsGroup/chicdiff, commit afd26d9) with Chicago 1.38
# and DESeq2 in the conda env 'chicdiff' (R 4.5). Inputs per replicate: chinput (bam2chicago) and the CHiCAGO
# background and scores of HiCExplorer v4 (kasumi_call_chicago.sh), converted by kasumi_chicdiff.R prepare.
# Steps run one after another, each only with >= 30 GB of available memory.
set -euo pipefail
kasumi=${1:?usage: kasumi_chicdiff.sh <kasumi1 dir> <output dir> <egaDifferential edges.tsv.gz>}
out=${2:?}
edges=${3:?}
here=$(cd "$(dirname "$0")" && pwd)
rscript=${RSCRIPT:-$HOME/miniconda3/envs/chicdiff/bin/Rscript}
python=${PYTHON:-python}
mkdir -p "$out"
wait_memory() {
    while (( $(free -g | awk '/^Mem:/{print $7}') < 30 )); do
        echo "$(date -Is) waiting for memory" >&2; sleep 300
    done
}
samples=(siMM_rep1 siMM_rep2 siRE_rep1 siRE_rep2)
for s in "${samples[@]}"; do
    [[ -s "$out/$s.Rds" ]] && continue
    wait_memory
    $rscript "$here/kasumi_chicdiff.R" prepare "$kasumi/chicago/${s}_scores.txt" "$kasumi/chicago/${s}_background.txt" \
        "$out/$s.Rds.tmp" > "$out/prepare_$s.log" 2>&1
    mv "$out/$s.Rds.tmp" "$out/$s.Rds"
done
design=$kasumi/design/kasumi1_hg38_promoters
if [[ ! -s "$out/peakmatrix.txt" ]]; then
    peaks=()
    for s in "${samples[@]}"; do peaks+=("$s=$kasumi/chicago/${s}_score1.txt"); done
    $rscript "$here/kasumi_chicdiff.R" peaks $design.rmap $design.baitmap "$out/peakmatrix.txt.tmp" "${peaks[@]}" \
        > "$out/peaks.log" 2>&1
    mv "$out/peakmatrix.txt.tmp" "$out/peakmatrix.txt"
fi
# design directory with exactly one .rmap and one .baitmap (the 4-column baitmap copy has another suffix)
mkdir -p "$out/design"
cp $design.rmap $design.baitmap "$out/design/"   # copies: the results share does not support symbolic links
spec() {   # condition:rep=chinput,Rds;...
    local c=$1 text="$1:"
    for r in rep1 rep2; do text+="${c}_$r=$kasumi/chinput/${c}_$r/${c}_$r.chinput,$out/${c}_$r.Rds;"; done
    printf '%s' "${text%;}"
}
if [[ ! -s "$out/siMM_vs_siRE_results.tsv.gz" ]]; then
    wait_memory
    (cd "$out" && /usr/bin/time -v $rscript "$here/kasumi_chicdiff.R" run "$out/design" "$out/siMM_vs_siRE" \
        "$out/peakmatrix.txt" "$(spec siMM)" "$(spec siRE)" > "$out/chicdiff.log" 2>&1)
fi
$python "$here/differential_robustness.py" chicdiff --edges "$edges" --chicdiff "$out/siMM_vs_siRE_results.tsv.gz" \
    --padj 0.05 -o "$out/siMM_vs_siRE" > "$out/overlap.log" 2>&1
