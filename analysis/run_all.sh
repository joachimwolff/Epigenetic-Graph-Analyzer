#!/usr/bin/env bash
# Complete analysis of the manuscript: data download, calling, graphs, validation, differential mode, figures.
set -euo pipefail
data=${1:?usage: run_all.sh <data_dir> <results_dir>}
results=${2:?}
here=$(cd "$(dirname "$0")" && pwd)
"$here/../scripts/download_data.sh" "$data"
"$here/01_call_interactions.sh" "$data" "$results" 2
"$here/02_build_graphs.sh" "$data" "$results"
"$here/03_validate.sh" "$data" "$results"
"$here/04_differential.sh" "$data" "$results"
