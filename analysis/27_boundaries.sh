#!/usr/bin/env bash
# Domain boundaries from capture contacts, CTCF / cohesin binding and CTCF motif orientation, evaluated against Hi-C.
#   boundaries_hic.py        diamond insulation of the 4DN mcool files (GM12878, K562, H1; 10 / 25 / 50 kb) and of the
#                            Rubin et al. 2017 day matrices (50 kb, hg19); needs cooler (HICEXPLORER_PYTHON)
#   boundaries_features.py   CTCF summits, cohesin at the summit, CTCF motif orientation (JASPAR MA0139.2), capture
#                            insulation, convergent cohesin-CTCF loops among the graph edges, per-bin features
#   boundaries_evaluate.py   logistic boundary callers (capture, CTCF + cohesin, + orientation, combined), cross-cell
#                            and chromosome-split designs, precision / recall / F1 / AUPRC, bootstrap of the gain of
#                            the contacts over CTCF, domain agreement
#   boundaries_timecourse.py keratinocyte days 0 / 3 / 6: per-day boundaries, changes against Hi-C, expression
#   figure_boundaries.py     one PDF per panel; table_boundaries.py: table_s35.tex
# Inputs: <data_dir> = raw directory (motifs, annotation, domains, rubin2017_keratinocyte), <results_dir> with graphs/,
# timecourse_rubin2017/ and compartments_domains/hic_reference/ (23_hic_reference.sh). Outputs: <results>/boundaries/.
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

"$hpy" "$here/boundaries_hic.py" "$mcool" "$data/rubin2017_keratinocyte" "$out/hic" "$threads"
"$py" "$here/boundaries_features.py" "$data" "$results" "$threads"
"$py" "$here/boundaries_evaluate.py" "$results" "$data" "$threads" "$n_boot"
"$py" "$here/boundaries_timecourse.py" "$results" 200
region=$("$py" "$here/figure_boundaries.py" select "$results")
"$hpy" "$here/boundaries_hic.py" "$mcool" "$data/rubin2017_keratinocyte" "$out/hic" "$threads" \
    --regions "GM12878:25000:$region"
"$py" "$here/figure_boundaries.py" plot "$results"
"$py" "$here/table_boundaries.py" "$out" "$out/tables"
