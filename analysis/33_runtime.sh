#!/usr/bin/env bash
# Wall time and peak memory (GNU time -v, maximum resident set size) of the largest runs of the manuscript, every run
# pinned to <cores> cores (taskset) and run one after the other:
#   1. egaCreateGraph of the largest graph, GM12878 strong core of Mifsud et al. 2015 (>= 20 reads; 332,078 edges,
#      Supplementary Table S1), settings of 02_build_graphs.sh
#   2. egaDifferential of the largest comparison, Kasumi-1 siMM / siRE (union of 183,165 and 220,874 edges),
#      settings of 32_step_decomposition_aml.sh: --stepFeature decomposed, 200,000 per-gene permutations
#   3. egaValidateMultiway of the GM12878 graph (CHiCAGO calls of Cairns et al. 2016) against GM12878 SPRITE, and of
#      the strong-core graph against SPRITE, settings of 03_validate.sh
# Inputs are read from <results_dir> (graphs, SPRITE matrix, DESeq2 table, HiPore-C step support of a previous
# run_all.sh) and the Kasumi-1 calls; outputs, logs and runtime.tsv go to <out_dir> only.
#
# Usage: 33_runtime.sh <data_dir> <results_dir> <kasumi1 dir> <out_dir> [cores=8]
set -euo pipefail
usage="usage: 33_runtime.sh <data_dir> <results_dir> <kasumi1 dir> <out_dir> [cores=8]"
data=${1:?$usage}
results=${2:?$usage}
kasumi=${3:?$usage}
out=${4:?$usage}
cores=${5:-8}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONHASHSEED=${PYTHONHASHSEED:-2}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4} \
       MKL_NUM_THREADS=${MKL_NUM_THREADS:-4} NUMEXPR_MAX_THREADS=${NUMEXPR_MAX_THREADS:-4}
py=${PYTHON:-python}
mkdir -p "$out"
table="$out/runtime.tsv"
printf 'run\tcommand\tcores\twall_seconds\tpeak_rss_gb\tuser_seconds\tsystem_seconds\texit_status\n' > "$table"

measure() {   # <name> <command...>
    local name=$1
    shift
    local status=0
    /usr/bin/time -v -o "$out/$name.time" taskset -c "0-$((cores - 1))" "$@" > "$out/$name.log" 2>&1 || status=$?
    awk -F': ' -v name="$name" -v cmd="$(basename "$3")" -v cores="$cores" -v status="$status" '
        /Elapsed \(wall clock\) time/ { n = split($2, t, ":"); wall = (n == 3) ? t[1] * 3600 + t[2] * 60 + t[3] : t[1] * 60 + t[2] }
        /Maximum resident set size/ { rss = $2 / 1048576 }
        /User time/ { user = $2 }
        /System time/ { sys = $2 }
        END { printf "%s\t%s\t%s\t%.1f\t%.2f\t%.1f\t%.1f\t%s\n", name, cmd, cores, wall, rss, user, sys, status }' \
        "$out/$name.time" >> "$table"
    return 0
}

# 1. largest graph
baits="$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap"
measure createGraph_GM_mifsud_c20 "$py" -m epigeneticgraphanalyzer.egaCreateGraph \
    -i "$data/mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt" -f mifsud --minCount 20 \
    --baits "$baits" --tss "$data/annotation/gencode.v19.annotation.gtf.gz" \
    --maxDistance 2000000 --weight score --communityMethod leiden --resolution 1.0 --seed 42 \
    -o "$out/GM_mifsud_c20_r1"

# 2. largest differential comparison
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
measure differential_Kasumi1_siMM_vs_siRE "$py" -m epigeneticgraphanalyzer.egaDifferential \
    --verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1 \
    --permutations 200 --bootstrap 1000 --seed 0 --threads "$cores" \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --stepFeature decomposed \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 5 --names siMM siRE \
    --markRoles DNase_d2=active DNase_d10=active --stateScale genomeMean \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    --expression "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    -o "$out/siMM_vs_siRE"

# 3. multi-way validation against SPRITE
sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
for graph in GM_chicago GM_mifsud_c20_r1; do
    measure "validate_${graph}_sprite" "$py" -m epigeneticgraphanalyzer.egaValidateMultiway --verbose \
        -g "$results/graphs/$graph.graphml" -m "$sprite" --confirmPairs \
        --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 --seed 0 --threads "$cores" \
        -o "$out/${graph}_sprite"
done
cat "$table"
