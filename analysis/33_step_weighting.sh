#!/usr/bin/env bash
# Weighting of the step feature (egaDifferential --stepWeighting): the runs of 31_step_decomposition.sh (mark term,
# decomposed) for hESC / hNEC and Md / PA2 with the weightings 'unweighted' (every node at k steps once) and
# 'normalised' (path products divided by their sum: a weighted mean of the mark signal over the nodes at k steps),
# in <out>/<weighting>/, the variants run side by side with <threads> each. step_weighting_summary.py compares them
# with the product weighting of <product_dir> (the outputs of 31_step_decomposition.sh) and writes
# <out>/step_weighting_*.tsv. The package of this checkout is used (PYTHONPATH).
set -euo pipefail
data=${1:?usage: 33_step_weighting.sh <data_dir> <results_dir> <product_dir> <out_dir> [gene permutations] [threads per variant]}
results=${2:?}
product=${3:?}
out=${4:?}
gene_permutations=${5:-200000}
threads=${6:-4}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$out"
pids=()
for weighting in unweighted normalised; do
    bash "$here/31_step_decomposition.sh" "$data" "$results" "$out/$weighting" "$gene_permutations" "$threads" \
        "$weighting" > "$out/run_$weighting.log" 2>&1 &
    pids+=($!)
done
status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
python "$here/step_weighting_summary.py" "$out" "$product" unweighted normalised > "$out/summary.log" 2>&1
exit $status
