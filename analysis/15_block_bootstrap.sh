#!/usr/bin/env bash
# Chromosome block bootstrap (1,000 resamples of whole chromosomes, percentile 95% intervals) of results whose items
# share nodes and chromosomes (analysis/block_bootstrap_results.py):
#   * promoter hubs with the promoter and >= 2 partners in one GM12878 SPRITE cluster against shifted hubs
#     (per-hub table of analysis/hub_multiway.py), per hub size class and overall
#   * Delta R^2 of the direct and indirect levels of the nested expression models of the four egaDifferential
#     comparisons (gene feature tables of 04, 05, 07, 08), refitted in every resample
#   * ChAseR assortativity of the GM12878 pair networks at 1 to 5 steps (inputs of analysis/chaser_inputs.py) and
#     its difference to distance-matched randomisations of the step-1 network, written as edge lists by
#     analysis/chaser_random_networks.R (conda environment `chaser`)
set -euo pipefail
results=${1:?usage: 15_block_bootstrap.sh <results_dir> [resamples=1000]}
n_boot=${2:-1000}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/bootstrap"
mkdir -p "$out"
boot() { python "$here/block_bootstrap_results.py" "$@" --bootstrap "$n_boot" --seed 0; }

boot hubs "$results/hubs/GM_chicago_sprite_hubs.tsv.gz" "$out/GM_chicago_sprite" > "$out/GM_chicago_sprite_hubs.log" 2>&1

for comparison in differential/hESC_vs_hNEC kessler/Md_vs_PA2 aml_kasumi1/siMM_vs_siRE aml_assi/t821_vs_FLT3ITD; do
    name=$(basename "$comparison")
    boot nested "$results/${comparison}_gene_features.tsv" "$results/${comparison}_nested_models.tsv" "$out/$name" \
        > "$out/${name}_nested.log" 2>&1
done

chaser="$results/competitors/chaser"
[[ -s "$out/chaser_random/random_step1_chas.tsv" ]] || conda run -n chaser Rscript "$here/chaser_random_networks.R" \
    "$chaser" "$out/chaser_random" 10 1 > "$out/chaser_random.log" 2>&1
boot chaser "$chaser" "$chaser/chaser_assortativity_GM12878.tsv" "$out/chaser_random" "$out/GM12878" \
    > "$out/GM12878_chaser.log" 2>&1
