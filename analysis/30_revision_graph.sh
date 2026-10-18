#!/usr/bin/env bash
# Reviewer requests on the graph side of the multi-way support (Bioinformatics revision); outputs under
# <results>/revision_bioinf/<analysis>/, every run skipped if its output exists.
#   B1 degree-aware null: step curves with controls also matched on the graph degree class of both ends
#      (hop_distance.py --degreeClasses), and degree-preserving rewired graphs that keep chromosome, bait / other-end
#      kind and length bin of every edge (rewired_null.py; log10 length bins of 0.1, 0.3 and 1.0)
#   B6 input comparability: size, edge lengths and degrees of every graph (graph_comparability.py) and the GM12878
#      step curves with the graph built at CHiCAGO score 5, 8 and 12
#   B7 graph sensitivity: distance limit of the pairs (hop_distance.py --maxDistance), number of shortest paths and
#      personalised PageRank of the indirect pairs (path_ensemble.py)
#   B9 ascertainment: nodes and pairs lost in the mapping to the multi-way data (ascertainment.py), step curves on the
#      top four coverage quintiles of the graph loci (hop_distance.py --minCoverageQuantile 0.2)
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

if [[ "$stage" == all || "$stage" == B6 ]]; then
    # input comparability: size and shape of every graph, and the GM12878 graph at CHiCAGO score 5, 8 and 12 (the
    # published WashU calls of GSE81503 hold score >= 5 only, so lower thresholds cannot be built)
    o="$out/B6_comparability"
    mkdir -p "$o"
    r=$results
    jung="Jung et al. 2019 published po/pp interaction tables (GSE86189); differential driver not in this repository"
    {
        printf 'name\tpath\tassembly\tcalling\treplicate_rule\tsource\n'
        printf 'GM_chicago\t%s\thg19\tCHiCAGO score >= 5, edges <= 2 Mb\tmerged replicates (published)\tCairns et al. 2016, GSE81503\n' "$r/graphs/GM_chicago.graphml"
        printf 'GM_chicago_all\t%s\thg19\tCHiCAGO score >= 5, no distance limit\tmerged replicates (published)\tCairns et al. 2016, GSE81503\n' "$r/graphs_nolimit/GM_chicago_all.graphml"
        printf 'GM_mifsud_c20_r1\t%s\thg19\tsignificant (Mifsud TS5), >= 20 reads, edges <= 2 Mb\tas published\tMifsud et al. 2015\n' "$r/graphs/GM_mifsud_c20_r1.graphml"
        printf 'hESC\t%s\thg19\tCHiCAGO score >= 12, edges <= 2 Mb\tmerged replicates (published)\tFreire-Pritchett et al. 2017, GSE86821\n' "$r/graphs/hESC.graphml"
        printf 'hNPC\t%s\thg19\tCHiCAGO score >= 12, edges <= 2 Mb\tmerged replicates (published)\tFreire-Pritchett et al. 2017, GSE86821\n' "$r/graphs/hNPC.graphml"
        printf 'H1\t%s\thg19\tHiCUP + CHiCAGO, score >= 5, edges <= 2 Mb\tseven libraries pooled (01_h1_chicago.sh)\tJung et al. 2019 raw reads, reprocessed\n' "$r/graphs/H1.graphml"
        printf 'K562\t%s\tGRCh38\tCHiCAGO score >= 5\tsingle published call set\tGSE250135 (GSM7974410, DMSO)\n' "$r/graphs/K562.graphml"
        printf 'K562_dmso\t%s\tGRCh38\tCHiCAGO score >= 5\tsingle published call set\tGSE250135 (GSM7974410)\n' "$r/k562/K562_dmso.graphml"
        printf 'K562_pon\t%s\tGRCh38\tCHiCAGO score >= 5\tsingle published call set\tGSE250135 (GSM7974411)\n' "$r/k562/K562_pon.graphml"
        printf 'VC_GM12878\t%s\tGRCh38\tvirtual capture, CHiCAGO (HiCExplorer) 5 kb, score >= 5\tone Hi-C map\tRao et al. 2014 Hi-C\n' "$r/virtual_capture/graphs/VC_GM12878.graphml"
        printf 'VC_H1_microc\t%s\tGRCh38\tvirtual capture, CHiCAGO (HiCExplorer) 5 kb, score >= 5\tone Micro-C map\tKrietenstein et al. 2020 Micro-C\n' "$r/virtual_capture/graphs/VC_H1_microc.graphml"
        printf 'hESC_vs_hNEC\t%s\thg19\tCHiCAGO score >= 12 per condition\tmerged replicates (published)\tFreire-Pritchett et al. 2017 (04_differential.sh)\n' "$r/differential/hESC_vs_hNEC_union.graphml"
        printf 'Md_vs_PA2\t%s\tmm10\tCHiCAGO score >= 5\tintersect: score >= 5 in both replicates\tKessler et al. 2023, GSE211902 (05_kessler_md_pa2.sh)\n' "$r/kessler/Md_vs_PA2_union.graphml"
        printf 'siMM_vs_siRE\t%s\tGRCh38\tCHiCAGO (own calls) score >= 5\tall: present in both replicates, score >= 5 in one\tPtasinska et al., own processing (08_aml_kasumi1.sh)\n' "$r/aml_kasumi1/siMM_vs_siRE_union.graphml"
        printf 't821_vs_FLT3ITD\t%s\tGRCh38\tHOMER FDR <= 0.01\tone patient per subtype\tAssi et al. 2019, GSE108316 (07_aml_assi.sh)\n' "$r/aml_assi/t821_vs_FLT3ITD_union.graphml"
        for s in ME msc npc TB; do
            for p in 1e-2 1e-3; do
                printf 'H1_vs_%s_p%s\t%s\thg19\tp-value tier %s\tas published\t%s\n' "$s" "$p" "$r/differential/H1_vs_${s}_p${p}_union.graphml" "$p" "$jung"
            done
        done
        for d in D0 D3 D6; do
            printf 'keratinocyte_%s\t%s\thg19\tcount-derived score >= 5 (rubin_chic_to_ibed.py)\tall: both donors, score >= 5 in one\tRubin et al. 2017, GSE84660 (22_timecourse.sh)\n' "$d" "$r/timecourse_rubin2017/graphs/$d.graphml"
        done
    } > "$o/manifest.tsv"
    run "$o/graph_comparability.tsv" "$o/graph_comparability.log" "$py" "$here/graph_comparability.py" \
        "$o/manifest.tsv" "$o/graph_comparability.tsv" &
    calls="$data/chicago_calls"
    for t in 5 8 12; do
        prefix="$o/GM_chicago_score$t"
        [[ -s "$prefix.graphml" ]] || "$py" -m epigeneticgraphanalyzer.egaCreateGraph \
            -i "$calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz" -f washu \
            --scoreThreshold "$t" --baits "$calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap" \
            --tss "$data/annotation/gencode.v19.annotation.gtf.gz" --maxDistance 2000000 --weight score \
            --communityMethod leiden --resolution 1.0 --seed 42 -o "$prefix" > "$prefix.log" 2>&1
        for m in "sprite:$sprite" "hiporec:$hiporec"; do
            wait_slot
            run "${prefix}_${m%%:*}_hops.tsv" "${prefix}_${m%%:*}_hops.log" "$py" "$here/hop_distance.py" \
                "$prefix.graphml" "${m#*:}" "${prefix}_${m%%:*}_hops.tsv" &
        done
    done
    wait
fi

if [[ "$stage" == all || "$stage" == B7 ]]; then
    # graph sensitivity: (a) distance limit of the tested pairs (2 Mb, 5 Mb, all of one chromosome) on the graph with
    # edges <= 2 Mb and on the graph without an edge length limit; (b) number of shortest paths and personalised
    # PageRank of the tested pairs at two and three steps (pair tables of stage B1)
    o="$out/B7_sensitivity"
    mkdir -p "$o"
    for g in "GM_chicago:$gm" "GM_chicago_all:$results/graphs_nolimit/GM_chicago_all.graphml"; do
        for d in 2000000 5000000 0; do
            name="${g%%:*}_sprite_max$d"
            wait_slot
            run "$o/${name}_hops.tsv" "$o/${name}_hops.log" "$py" "$here/hop_distance.py" "${g#*:}" "$sprite" \
                "$o/${name}_hops.tsv" --maxDistance "$d" &
        done
    done
    for m in sprite hiporec; do
        wait_slot
        run "$o/GM_chicago_${m}_paths_strata.tsv" "$o/GM_chicago_${m}_paths.log" "$py" "$here/path_ensemble.py" \
            "$gm" "$out/B1_degree_null/GM_chicago_${m}_hops_coverage_pairs.tsv.gz" "$o/GM_chicago_${m}_paths" &
    done
    wait
fi

if [[ "$stage" == all || "$stage" == B9 ]]; then
    # ascertainment: nodes and pairs lost on the way to each multi-way data set, and the SPRITE / HiPore-C step curves
    # with the graph loci of the lowest coverage quintile removed from tests and controls
    o="$out/B9_ascertainment"
    mkdir -p "$o"
    porec="$results/multiway/GM12878_porec_hg19_HindIII_fragments.npz"
    k562_hiporec="$results/multiway/K562_hiporec_hg19_HindIII_fragments.npz"
    to_hg19="$data/annotation/hg38ToHg19.over.chain.gz"
    for p in "GM_chicago_sprite $gm $sprite -" "hESC_gam $hesc $gam $chain" "H1_gam $h1 $gam $chain" \
             "GM_chicago_hiporec $gm $hiporec -" "GM_chicago_porec $gm $porec -" \
             "K562_dmso_hiporec $results/k562/K562_dmso.graphml $k562_hiporec $to_hg19"; do
        read -r name graph multiway ch <<< "$p"
        chain_option=()
        [[ "$ch" != - ]] && chain_option=(--chain "$ch")
        wait_slot
        run "$o/${name}_pairs.tsv" "$o/${name}.log" "$py" "$here/ascertainment.py" "$graph" "$multiway" "$o/$name" \
            "${chain_option[@]}" &
    done
    for m in "sprite:$sprite" "hiporec:$hiporec"; do
        wait_slot
        run "$o/GM_chicago_${m%%:*}_hops_top4quintiles.tsv" "$o/GM_chicago_${m%%:*}_hops_top4quintiles.log" "$py" \
            "$here/hop_distance.py" "$gm" "${m#*:}" "$o/GM_chicago_${m%%:*}_hops_top4quintiles.tsv" \
            --minCoverageQuantile 0.2 &
    done
    wait
fi
