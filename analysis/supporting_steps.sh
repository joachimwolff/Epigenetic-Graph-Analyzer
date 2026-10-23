#!/usr/bin/env bash
# Analyses whose outputs are read by later drivers, the supplementary tables and the figures, but that have no
# numbered driver of their own. One step per call (run_all.sh calls them in dependency order):
#   k562_graph            K562 DMSO and PON graphs (GSE250135, GRCh38; washU 'other end' column split into
#                         chrom,start,end), egaCreateGraph defaults, <results>/k562/          (09, 19, S1)
#   virtual_capture       graphs of the virtual capture calls (CHiCAGO on promoter bins of Hi-C / Micro-C, 5 kb, in
#                         <virtual capture dir>) and their agreement with the capture calls,
#                         <results>/virtual_capture/                                          (11, 19, S1, S9, S20)
#   hubs                  promoter hubs in GM12878 SPRITE and H1 GAM (hub_multiway.py), <results>/hubs/  (15, S4)
#   explain               egaExplain of the four differential runs of 04, 05, 07, 08, <results>/explain/ (S14)
#   competitors           ChAseR assortativity (chaser_inputs.py, chaser_assortativity.R, conda env `chaser`) and
#                         promoter hubs against HLM-Nbody polymer three-body predictions (HLM_NBODY: the sprite/
#                         directory of github.com/leiliu2015/HLM-Nbody, commit a49ea16; default
#                         <data_dir>/competitors/HLM-Nbody/sprite), <results>/competitors/     (15, S12, S13)
#   community_robustness  community_robustness.py, <results>/robustness_communities/          (S21)
#   crispr_benchmark      K562 CRISPRi benchmark of graph features (crispr_benchmark.py), training and held-out
#                         pairs, <results>/k562/bench_*                                         (S10)
#   sprite_resolution     SPRITE downsampling and pile-up on the Mifsud strong-core graph, <results>/diagnostics/
#   figures               figure_concept.py, figure_3d_support.py, figure_virtual_capture.py,
#                         figure_sprite_resolution.py, <results>/figures/
#   reports               egaReport main and gene reports of Md / PA2 and hESC / hNEC with the TAD calls of
#                         28_tad_crossing.sh, <results>/reports/
#
# Usage: supporting_steps.sh <step> <data_dir> <results_dir> <kasumi1 dir> <virtual capture dir> [threads]
set -euo pipefail
step=${1:?usage: supporting_steps.sh <step> <data_dir> <results_dir> <kasumi1 dir> <virtual capture dir> [threads]}
data=${2:?}
results=${3:?}
kasumi=${4:?}
vc=${5:?}
threads=${6:-8}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
py=${PYTHON:-python}
annotation="$data/annotation"

case $step in
k562_graph)
    out="$results/k562"
    mkdir -p "$out"
    for s in dmso:GSM7974410_k562_dmso pon:GSM7974411_k562_pon; do
        zcat "$data/k562/contacts/${s#*:}_washU.txt.gz" |
            awk -F'\t' -v OFS='\t' '{split($4,o,/[:,-]/); print $1","$2","$3, o[1]","o[2]","o[3], $5}' |
            gzip > "$out/K562_${s%%:*}_washu.txt.gz"
        "$py" -m epigeneticgraphanalyzer.egaCreateGraph -i "$out/K562_${s%%:*}_washu.txt.gz" -f washu --weight score \
            --maxDistance 2000000 --tss "$annotation/gencode.v29.annotation.gtf.gz" -o "$out/K562_${s%%:*}" \
            > "$out/K562_${s%%:*}.log" 2>&1
    done
    ;;
virtual_capture)
    out="$results/virtual_capture"
    g="$out/graphs"
    mkdir -p "$g"
    baits19="$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap"
    chain="$annotation/hg19ToHg38.over.chain.gz"
    gm_real="$data/chicago_calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz"
    h1_real="$data/chicago_calls/H1_jung2019_reprocessed/H1.ibed.gz"
    bench() { "$py" "$here/virtual_capture_benchmark.py" "$@"; }
    bench "$gm_real" washu "$baits19" "$chain" "$vc/GM12878/chicago_GM12878.ibed" "$vc/GM12878/GM12878.baitmap" 5000 \
        "$out/GM12878" > "$out/GM12878.log" 2>&1
    # two capture experiments of GM12878 on the same bins: CHiCAGO calls against the Mifsud strong core (>= 20 reads)
    awk -F'\t' 'NR==1 || $10>=20' "$data/mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt" \
        > "$out/mifsud_c20.txt"
    bench "$gm_real" washu "$baits19" "$chain" "$out/mifsud_c20.txt" "$vc/GM12878/GM12878.baitmap" 5000 \
        "$out/GM12878_cairns_vs_mifsud20" 5 mifsud "$baits19" "$chain" -1000 > "$out/GM12878_cairns_vs_mifsud20.log" 2>&1
    bench "$data/k562/contacts/GSM7974410_k562_dmso_washU.txt.gz" washu - - "$vc/K562/chicago_K562.ibed" \
        "$vc/K562/K562.baitmap" 5000 "$out/K562" > "$out/K562.log" 2>&1
    bench "$h1_real" ibed - "$chain" "$vc/H1/hicx_H1.ibed" "$vc/H1/H1.baitmap" 5000 "$out/H1_microc" \
        > "$out/H1_microc.log" 2>&1
    # matched number: the 39% highest-scoring cis calls of the Micro-C calls
    awk -F'\t' 'NR==1{print; next} $1==$5' "$vc/H1/hicx_H1.ibed" > "$out/H1_microc_cis.ibed"
    n=$(( $(wc -l < "$out/H1_microc_cis.ibed") - 1 ))
    (head -1 "$out/H1_microc_cis.ibed"; tail -n +2 "$out/H1_microc_cis.ibed" | sort -t$'\t' -k10,10gr | awk -v n=$((n * 39 / 100)) 'NR <= n') \
        > "$out/H1_microc_top.ibed"
    bench "$h1_real" ibed - "$chain" "$out/H1_microc_top.ibed" "$vc/H1/H1.baitmap" 5000 "$out/H1_microc_matched" \
        > "$out/H1_microc_matched.log" 2>&1
    rm -f "$out/H1_microc_cis.ibed"
    graph() {
        "$py" -m epigeneticgraphanalyzer.egaCreateGraph -i "$2" -f ibed --scoreThreshold 5 --maxDistance 2000000 \
            --weight score --tss "$annotation/gencode.v29.annotation.gtf.gz" -o "$g/$1" > "$g/$1.log" 2>&1
    }
    graph VC_H1_microc "$vc/H1/hicx_H1.ibed" &
    graph VC_GM12878 "$vc/GM12878/chicago_GM12878.ibed" &
    wait
    ;;
hubs)
    out="$results/hubs"
    mkdir -p "$out"
    "$py" "$here/hub_multiway.py" "$results/graphs/GM_chicago.graphml" \
        "$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz" "$out/GM_chicago_sprite" --threads "$threads" \
        > "$out/GM_chicago_sprite.log" 2>&1
    "$py" "$here/hub_multiway.py" "$results/graphs/hESC.graphml" "$results/multiway/gam_H1_hg38_30kb.npz" \
        "$out/hESC_gam" "$annotation/hg19ToHg38.over.chain.gz" --threads "$threads" > "$out/hESC_gam.log" 2>&1
    ;;
explain)
    out="$results/explain"
    mkdir -p "$out"
    for x in "differential/hESC_vs_hNEC hESC hNEC" "kessler/Md_vs_PA2 Md PA2" "aml_kasumi1/siMM_vs_siRE siMM siRE" \
             "aml_assi/t821_vs_FLT3ITD t821 FLT3ITD"; do
        set -- $x
        "$py" -m epigeneticgraphanalyzer.egaExplain --verbose -p "$results/$1" --names "$2" "$3" \
            -o "$out/$(basename "$1")" > "$out/$(basename "$1").log" 2>&1
    done
    ;;
competitors)
    chaser="$results/competitors/chaser"
    mkdir -p "$chaser" "$results/competitors/polymer"
    e="$data/encode_gm12878"
    "$py" "$here/chaser_inputs.py" "$results/graphs/GM_chicago.graphml" "$chaser" H3K27ac="$e/GM12878_H3k27ac.bigWig" \
        H3K4me1="$e/GM12878_H3k4me1.bigWig" H3K4me3="$e/GM12878_H3k4me3.bigWig" H3K27me3="$e/GM12878_H3k27me3.bigWig" \
        CTCF="$e/GM12878_Ctcf.bigWig" > "$chaser/inputs.log" 2>&1
    (cd "$chaser" && "${CHASER_RSCRIPT:-$HOME/miniconda3/envs/chaser/bin/Rscript}" "$here/chaser_assortativity.R" . \
        chaser_assortativity_GM12878.tsv 10 > chaser.log 2>&1)
    "$py" "$here/polymer_comparison.py" "${HLM_NBODY:-$data/competitors/HLM-Nbody/sprite}" \
        "$results/graphs/GM_chicago.graphml" "$results/competitors/polymer/GM12878" \
        > "$results/competitors/polymer/GM12878.log" 2>&1
    ;;
community_robustness)
    mkdir -p "$results/robustness_communities"
    "$py" "$here/community_robustness.py" "$data" "$results" "$results/robustness_communities" \
        > "$results/robustness_communities/run.log" 2>&1
    ;;
crispr_benchmark)
    c="$data/k562/crispr"
    for set in training:training_K562 heldout:heldout_5_cell_types; do
        "$py" "$here/crispr_benchmark.py" "$results/k562/K562_dmso.graphml" \
            "$c/github/EPCrisprBenchmark_combined_data.${set#*:}.GRCh38.tsv.gz" "$c/encode/ENCFF970QAX.bed.gz" \
            "$results/k562/bench_${set%%:*}" > "$results/k562/bench_${set%%:*}.log" 2>&1
    done
    ;;
sprite_resolution)
    out="$results/diagnostics"
    mkdir -p "$out"
    graph="$results/graphs/GM_mifsud_c20_r1.graphml"
    sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
    "$py" "$here/sprite_downsampling.py" "$graph" "$sprite" "$out/GM_mifsud_spritefrag_downsampling.tsv" 100000 \
        > "$out/downsampling.log" 2>&1 &
    "$py" "$here/sprite_pileup.py" "$graph" "$sprite" "$out/GM_mifsud_spritefrag_cov" 10 5000 > "$out/pileup_cov.log" 2>&1 &
    wait
    ;;
figures)
    out="$results/figures"
    mkdir -p "$out"
    "$py" "$here/figure_concept.py" "$out/figure1_concept"
    "$py" "$here/figure_3d_support.py" "$results" "$out/figure2"
    "$py" "$here/figure_virtual_capture.py" "$results" "$out/figure5"
    "$py" "$here/figure_sprite_resolution.py" "$results/diagnostics/GM_mifsud_spritefrag_downsampling.tsv" \
        "$results/diagnostics/GM_mifsud_spritefrag_cov_pileup.npz" "$out/figureS1_sprite_resolution.pdf"
    ;;
reports)
    out="$results/reports"
    mkdir -p "$out"
    tc="$results/tad_crossing"
    md=(--tads "$tc/hicfindtads/PA2_25kb" --tadRegions "$tc/kessler2023_HIRE_mm10.bed" --tadMinStrength 0.5
        --annotation "$data/kessler/mm10.ncbiRefSeq.gtf")
    hs=(--tads "$results/compartments_domains/hic_reference/H1_50kb" --tadChain "$annotation/hg38ToHg19.over.chain.gz"
        --tadMinStrength 0.5 --annotation "$annotation/gencode.v19.annotation.gtf.gz")
    report() { "$py" -W ignore -m epigeneticgraphanalyzer.egaReport --verbose "$@" -o "$out"; }
    report -p "$results/hops_n/Md_vs_PA2" "${md[@]}" > "$out/Md_vs_PA2.log" 2>&1
    report -p "$results/hops_n/Md_vs_PA2" --gene Hoxa3 Hoxa2 "${md[@]}" > "$out/Md_vs_PA2_genes.log" 2>&1
    report -p "$results/hops_n/hESC_vs_hNEC" "${hs[@]}" > "$out/hESC_vs_hNEC.log" 2>&1
    report -p "$results/hops_n/hESC_vs_hNEC" --gene STAMBPL1 CALCOCO1 NR2F1 EPHA7 "${hs[@]}" \
        > "$out/hESC_vs_hNEC_genes.log" 2>&1
    ;;
*)
    echo "unknown step: $step" >&2
    exit 1
    ;;
esac
