#!/usr/bin/env bash
# Baselines of the expression models (expression_baselines.py): on the informative genes of the decomposed runs
# (31_step_decomposition.sh, 32_step_decomposition_aml.sh) and with the promoter model as base, the distance-weighted
# mark signal of the linear neighbourhood (window, 1/d, union-graph nodes and bigWig bins, no edges), the promoter
# degree, contact_direct, the direct and the graph features; graph features on top of the window and the reverse;
# concordant genes among the differential genes of the cascade. Outputs in <out>/ only. The package of this checkout
# is used (PYTHONPATH).
set -euo pipefail
decomposed=${1:?usage: 34_expression_baselines.sh <decomposed_dir> <out_dir> [threads]}
out=${2:?}
threads=${3:-8}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$out"
python "$here/expression_baselines.py" "$decomposed" "$out" --comparisons hESC_vs_hNEC Md_vs_PA2 siMM_vs_siRE \
    --threads "$threads" --binSize 5000 > "$out/expression_baselines.log" 2>&1
