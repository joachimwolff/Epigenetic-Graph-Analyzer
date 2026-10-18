#!/usr/bin/env bash
# Reviewer requests on the graph side of the multi-way support (Bioinformatics revision); outputs under
# <results>/revision_bioinf/<analysis>/, every run skipped if its output exists.
#   B1 degree-aware null: step curves with controls also matched on the graph degree class of both ends
#      (hop_distance.py --degreeClasses), and degree-preserving rewired graphs that keep chromosome, bait / other-end
#      kind and length bin of every edge (rewired_null.py; log10 length bins of 0.1, 0.3 and 1.0)
# Usage: 30_revision_graph.sh <data_dir> <results_dir> [stage] [parallel jobs]
set -euo pipefail
data=${1:?usage: 30_revision_graph.sh <data_dir> <results_dir> [stage] [jobs]}
results=${2:?}
stage=${3:-all}
jobs=${4:-4}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=${PYTHON:-python}
out="$results/revision_bioinf"
mkdir -p "$out"
sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
hiporec="$results/multiway/GM12878_hiporec_hg19_HindIII_fragments.npz"
gam="$results/multiway/gam_H1_hg38_30kb.npz"
chain="$data/annotation/hg19ToHg38.over.chain.gz"
gm="$results/graphs/GM_chicago.graphml"
hesc="$results/graphs/hESC.graphml"
h1="$results/graphs/H1.graphml"

run() {   # <output> <log> <command...>: run unless the output exists
    local output=$1 log=$2
    shift 2
    [[ -s "$output" ]] && return 0
    "$@" > "$log" 2>&1 || echo "FAILED: $output (see $log)" >&2
}
wait_slot() { while (( $(jobs -rp | wc -l) >= jobs )); do sleep 5; done; }

if [[ "$stage" == all || "$stage" == B1 ]]; then
    o="$out/B1_degree_null"
    mkdir -p "$o"
    pairings=("GM_chicago_sprite $gm $sprite -" "GM_chicago_hiporec $gm $hiporec -" "hESC_gam $hesc $gam $chain"
              "H1_gam $h1 $gam $chain")
    for p in "${pairings[@]}"; do
        read -r name graph multiway ch <<< "$p"
        chain_option=()
        [[ "$ch" != - ]] && chain_option=(--chain "$ch")
        for variant in "coverage:--coverageClasses 5" "coverage_degree3:--coverageClasses 5 --degreeClasses 3" \
                       "degree5_only:--coverageClasses 0 --degreeClasses 5"; do
            label=${variant%%:*}
            read -r -a options <<< "${variant#*:}"
            wait_slot
            run "$o/${name}_hops_${label}.tsv" "$o/${name}_hops_${label}.log" "$py" "$here/hop_distance.py" \
                "$graph" "$multiway" "$o/${name}_hops_${label}.tsv" "${chain_option[@]}" "${options[@]}" \
                --pairTable "$o/${name}_hops_${label}_pairs.tsv.gz" &
        done
    done
    wait
    for p in "${pairings[@]:0:3}"; do
        read -r name graph multiway ch <<< "$p"
        chain_option=()
        [[ "$ch" != - ]] && chain_option=(--chain "$ch")
        for width in 0.1 0.3 1.0; do
            wait_slot
            run "$o/${name}_rewired_w${width}_summary.tsv" "$o/${name}_rewired_w${width}.log" "$py" \
                "$here/rewired_null.py" "$graph" "$multiway" "$o/${name}_rewired_w${width}" "${chain_option[@]}" \
                --rewirings 10 --logBinWidth "$width" --seed 0 &
        done
    done
    wait
fi
