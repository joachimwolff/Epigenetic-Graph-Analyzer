#!/usr/bin/env bash
# CRISPRi evaluation for the Bioinformatics revision.
# Pairs: crispr_benchmark.py on the K562 DMSO graph against the two CRISPRi pair sets of the ENCODE-rE2G benchmark
# (EPCrisprBenchmark, Gschwind et al.): the K562 training pairs and the K562 part of the held-out pairs, with ABC,
# ENCODE-rE2G and power-law contact scores of the published K562 ENCODE-rE2G predictions (ENCFF970QAX, ENCSR627ANP).
# B2 (crispr_cv.py): leakage-aware out-of-fold evaluation (leave-one-chromosome-out and GroupKFold by gene) of the
# nested logistic models, gene-cluster bootstrap CIs, coefficient tables and likelihood-ratio tests.
# B3 (crispr_head_to_head.py): head-to-head with the published predictions on identical pairs.
# Outputs: <out>/pairs/bench_{training,heldout}_*, <out>/B2_cv/, <out>/B3_head_to_head/.
set -euo pipefail
data=${1:?usage: 30_crispr_revision.sh <data_dir> <results_dir> <out_dir> [threads]}
results=${2:?}
out=${3:?}
threads=${4:-6}
here=$(cd "$(dirname "$0")" && pwd)
crispr="$data/k562/crispr"
graph="$results/k562/K562_dmso.graphml"
predictions="$crispr/encode/ENCFF970QAX.bed.gz"
mkdir -p "$out/pairs" "$out/B2_cv"

if [ ! -s "$out/pairs/bench_heldout_pairs.tsv.gz" ]; then
    python "$here/crispr_benchmark.py" "$graph" \
        "$crispr/github/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz" "$predictions" \
        "$out/pairs/bench_training" > "$out/pairs/bench_training.log" 2>&1
    python "$here/crispr_benchmark.py" "$graph" \
        "$crispr/github/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz" "$predictions" \
        "$out/pairs/bench_heldout" > "$out/pairs/bench_heldout.log" 2>&1
fi

python "$here/crispr_cv.py" "$out/pairs/bench_training_pairs.tsv.gz" "$out/pairs/bench_heldout_pairs.tsv.gz" \
    "$out/B2_cv" --threads "$threads" --bootstrap 1000 --seed 0 > "$out/B2_cv/crispr_cv.log" 2>&1

# B3 (crispr_head_to_head.py): distance, ABC, ENCODE-rE2G, power law, EGA direct / two-step / PPR and the
# out-of-sample EGA logistic models on identical pairs; coverage, intersection and missing = lowest score.
mkdir -p "$out/B3_head_to_head"
python "$here/crispr_head_to_head.py" "$out/pairs/bench_training_pairs.tsv.gz" \
    "$out/pairs/bench_heldout_pairs.tsv.gz" "$out/B2_cv" "$out/B3_head_to_head" --threads "$threads" \
    --bootstrap 1000 --seed 0 > "$out/B3_head_to_head/crispr_head_to_head.log" 2>&1
