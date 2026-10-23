#!/usr/bin/env bash
# Domain boundaries from capture contacts, CTCF / cohesin binding and CTCF motif orientation, evaluated against Hi-C.
#   boundaries_hic.py        diamond insulation of the 4DN mcool files (GM12878, K562, H1; 10 / 25 / 50 kb); needs
#                            cooler (HICEXPLORER_PYTHON)
#   boundaries_features.py   CTCF summits, cohesin at the summit, CTCF motif orientation (JASPAR MA0139.2), capture
#                            insulation, convergent cohesin-CTCF loops among the graph edges
#   boundaries_evaluate.py   per-bin evidence (capture dip, cohesin-bound CTCF with its orientation pattern, loop
#                            anchors) and rule-based calls with fixed thresholds (no fitted model): dip, CTCF, CTCF
#                            pattern, dip OR / AND CTCF (pattern); precision / recall / F1 / AUPRC against Hi-C,
#                            threshold grid, bootstrap of the change when the dip is added to CTCF, domain agreement
#   figure_boundaries.py     one PDF per panel; table_boundaries.py: table_s35.tex
# Inputs: <data_dir> = raw directory (motifs, annotation, domains), <results_dir> with graphs/ and
# compartments_domains/hic_reference/ (23_hic_reference.sh). Outputs: <results>/boundaries/.
set -euo pipefail
mcool=${1:?usage: 27_boundaries.sh <mcool_dir> <data_dir> <results_dir> [threads] [bootstrap]}
data=${2:?}
results=${3:?}
threads=${4:-12}
n_boot=${5:-500}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=${PYTHON:-python}
hpy=${HICEXPLORER_PYTHON:-python}
out="$results/boundaries"
mkdir -p "$out"/{hic,cache,tables,models,calls,figures,logs}

"$hpy" "$here/boundaries_hic.py" "$mcool" "$out/hic" "$threads"
"$py" "$here/boundaries_features.py" "$data" "$results" "$threads"
"$py" "$here/boundaries_evaluate.py" "$results" "$data" "$threads" "$n_boot"
region=$("$py" "$here/figure_boundaries.py" select "$results")
"$hpy" "$here/boundaries_hic.py" "$mcool" "$out/hic" "$threads" \
    --regions "GM12878:25000:$region"
"$py" "$here/figure_boundaries.py" plot "$results"
"$py" "$here/table_boundaries.py" "$out" "$out/tables"
