#!/usr/bin/env bash
# Calibration of the explanations beyond step 1 between naive B cells and GM12878, with the settings of
# 35_gm12878_naiveB.sh:
#   unchanged: genes whose expression does not change receive the real fold changes (pseudo_de.py); how often does
#     EGA offer them an explanation, level by level?
#   rewired: the interactions of both conditions are rewired together, keeping degrees, bait / other-end kinds,
#     chromosomes, edge lengths and the scores of both conditions (rewire_inputs.py); how often are the real
#     differentially expressed genes explained when the graph connects other fragments at the same distances?
# Each null is run with --seeds seeds; 43_calibration_summary.py compares the explanation levels with the real run.
#
# Usage: 42_calibration.sh <data_dir> <real run dir> <out_dir> [seeds] [gene permutations] [threads] [parallel]
#            [nulls: unchanged rewired] [rewiring length bin width]
set -euo pipefail
data=${1:?usage: 42_calibration.sh <data_dir> <real run dir> <out_dir> [seeds] [gene permutations] [threads] [parallel]}
real=${2:?}
out=${3:?}
seeds=${4:-3}
gene_permutations=${5:-200000}
threads=${6:-12}
parallel=${7:-2}
nulls=${8:-unchanged rewired}
bin_width=${9:-0.1}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
nb="$data/naiveB"
washu_a="$nb/javierre2016/naiveB_nB_washU.txt.gz"
washu_b="$data/chicago_calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz"
baitmap="$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap"
roadmap() { for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }

run() {  # <name> <interactions A> <interactions B> <expression>
    egaDifferential --verbose --maxHops 3 --genePermutations "$gene_permutations" --rnaThreshold 1 \
        --chromatinThreshold 1 --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
        --stepFeature decomposed \
        --interactionsA "$2" --interactionsB "$3" --format washu --baits "$baitmap" \
        --scoreThreshold 5 --referenceScore 5 --names naiveB GM12878 \
        --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --tracksA $(roadmap E032) --tracksB $(roadmap E116) \
        --markSign H3K27me3=-1 --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
        --expression "$4" --log2FCColumn log2FoldChange --padjColumn padj \
        --alpha 0.05 --geneKey genes -o "$out/$1/naiveB_vs_GM12878" > "$out/$1.log" 2>&1
}

jobs=()
for seed in $(seq 0 $((seeds - 1))); do
    if [[ " $nulls " == *" unchanged "* ]]; then
        mkdir -p "$out/unchanged_$seed"
        python "$here/pseudo_de.py" "$real/deseq2_GM12878_vs_naiveB.tsv" "$out/unchanged_$seed/expression.tsv" \
            --seed "$seed" > "$out/unchanged_$seed/expression.log"
        jobs+=("unchanged_$seed $washu_a $washu_b $out/unchanged_$seed/expression.tsv")
    fi
    if [[ " $nulls " == *" rewired "* ]]; then
        mkdir -p "$out/rewired_$seed"
        python "$here/rewire_inputs.py" "$washu_a" "$washu_b" "$baitmap" "$out/rewired_$seed/naiveB.washU.txt.gz" \
            "$out/rewired_$seed/GM12878.washU.txt.gz" --seed "$seed" --logBinWidth "$bin_width" \
            > "$out/rewired_$seed/rewiring.log" 2>&1
        jobs+=("rewired_$seed $out/rewired_$seed/naiveB.washU.txt.gz $out/rewired_$seed/GM12878.washU.txt.gz $real/deseq2_GM12878_vs_naiveB.tsv")
    fi
done
for job in "${jobs[@]}"; do
    while (($(jobs -rp | wc -l) >= parallel)); do wait -n; done
    run $job &
done
wait
[[ -f "$here/43_calibration_summary.py" ]] && python "$here/43_calibration_summary.py" "$real/naiveB_vs_GM12878" "$out" > "$out/summary.txt"
