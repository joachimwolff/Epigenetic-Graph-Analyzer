#!/usr/bin/env bash
# Minimal example: the Md / PA2 comparison (mouse E10.5 mandibular arch versus second pharyngeal arch) restricted to
# the Hoxa locus, chr6:49.5-54.5 Mb (mm10). Inputs in example/data/ (prepare_example.py cuts them out of the public
# data of Kessler et al. 2023 and Minoux et al. 2017). Settings of the manuscript run (analysis/31_step_decomposition.sh:
# decomposed step features, up to three steps, 200,000 per-gene permutations, seed 0), with the PA2 domain calls of
# hicFindTADs for the reports (analysis/supporting_steps.sh reports).
#   1. egaCreateGraph  Md and PA2 graphs (replicate rule 'intersect', as egaDifferential)
#   2. egaDifferential Md versus PA2
#   3. egaExplain      integrated tables of the differential run
#   4. egaReport       main report and gene reports of Hoxa3 and Hoxa2 (HTML)
# Run time about one minute on one core.
#
# Usage: example/run_example.sh [output dir=example/output] [threads=2]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/output}
threads=${2:-2}
data="$here/data"
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
# ties in the community transitions follow Python's string hashing; the manuscript runs pin it (run_all.sh)
export PYTHONHASHSEED=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=${PYTHON:-python}
mkdir -p "$out/graphs" "$out/explain" "$out/reports"
md=("$data"/GSM6505208_*.txt.gz "$data"/GSM6505209_*.txt.gz)
pa2=("$data"/GSM6505210_*.txt.gz "$data"/GSM6505211_*.txt.gz)
gtf="$data/mm10.ncbiRefSeq.chr6_hoxa.gtf.gz"

echo "1. egaCreateGraph"
for x in Md PA2; do
    if [[ $x == Md ]]; then files=("${md[@]}"); else files=("${pa2[@]}"); fi
    "$py" -m epigeneticgraphanalyzer.egaCreateGraph -i "${files[@]}" -f washu --replicates intersect \
        --tss "$gtf" --geneTypes --maxDistance 2000000 --weight score --communityMethod leiden --resolution 1.0 \
        --seed 42 -o "$out/graphs/$x" > "$out/graphs/$x.log" 2>&1
done

echo "2. egaDifferential"
marks() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$data/E10.5_$1_${mark}_R1.bigwig"; done; }
fold() { for mark in H3K27ac H3K4me2 H3K27me3; do printf '%s=%s ' "$mark" "$data/$1_${mark}.fold.bw"; done; }
# shellcheck disable=SC2046
"$py" -m epigeneticgraphanalyzer.egaDifferential --verbose --maxHops 3 --genePermutations 200000 \
    --rnaThreshold 1 --chromatinThreshold 1 --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --stepFeature decomposed \
    --interactionsA "${md[@]}" --interactionsB "${pa2[@]}" \
    --format washu --replicates intersect --names Md PA2 --tss "$gtf" --geneTypes \
    --tracksA $(marks MD) --tracksB $(marks PA2) --markSign H3K27me3=-1 \
    --markRoles H3K27ac=active H3K4me2=active H3K27me3=repressive \
    --stateTracksA $(fold Md) --stateTracksB $(fold PA2) \
    --expression "$data/deseq2_PA2_vs_Md.chr6_hoxa.tsv" --log2FCColumn log2FoldChange \
    --padjColumn padj --alpha 0.05 --geneKey genes \
    -o "$out/Md_vs_PA2" > "$out/Md_vs_PA2.log" 2>&1

echo "3. egaExplain"
"$py" -m epigeneticgraphanalyzer.egaExplain -p "$out/Md_vs_PA2" -o "$out/explain/Md_vs_PA2" \
    > "$out/explain_Md_vs_PA2.log" 2>&1

echo "4. egaReport"
tads=(--tads "$data/PA2_25kb" --tadRegions "$data/kessler2023_HIRE_mm10.bed" --tadMinStrength 0.5
      --annotation "$gtf")
"$py" -W ignore -m epigeneticgraphanalyzer.egaReport -p "$out/Md_vs_PA2" "${tads[@]}" --pdf none \
    -o "$out/reports" > "$out/reports/Md_vs_PA2.log" 2>&1
"$py" -W ignore -m epigeneticgraphanalyzer.egaReport -p "$out/Md_vs_PA2" --gene Hoxa3 Hoxa2 "${tads[@]}" \
    --pdf none -o "$out/reports" > "$out/reports/Md_vs_PA2_genes.log" 2>&1

echo "5. key results"
"$py" "$here/summarize_example.py" "$out" "$out/summary"
echo "done: $out"
