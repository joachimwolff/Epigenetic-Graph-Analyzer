#!/usr/bin/env bash
# Robustness of the gained / lost edges of the differential graph and of its communities (reviewer request), with
# analysis/differential_robustness.py:
#   threshold  edge statuses, coherent communities and community transitions recomputed at other calling thresholds
#              Kasumi-1 siMM vs siRE  CHiCAGO exports at score >= 1 per replicate, rule 'all' (present in both
#                                     replicates, score >= t in one), t = 5 (as run by 08_aml_kasumi1.sh), 3, 8, 10
#              Md vs PA2              published WashU calls (score >= 5) per replicate; reference = present in both
#                                     replicates (05_kessler_md_pa2.sh); score >= 5, 8, 10 in both or in either
#              plus, for the gained / lost edges of the reference, the largest score of the pair in the replicates of
#              the condition where it is absent
#   depth      Kasumi-1: binomial thinning of the CHiCAGO input counts (chinput) per replicate, (a) the deeper
#              condition siRE to the read count of the siMM replicate of the same number, (b) all four replicates
#              to 50%; CHiCAGO recalled as in kasumi_call_chicago.sh (HiCExplorer v4 chicChicago*, design unchanged,
#              export score >= 1) and the statuses recomputed with the rule of the run
#   chicdiff   Kasumi-1: Chicdiff (Cairns et al. 2019) on the per-replicate chinput files, overlap of its calls
#              (padj < 0.05) with the gained / lost edges (see the chicdiff step below)
# Large steps run one after another, each only with >= 30 GB of available memory, at most 8 threads.
set -euo pipefail
results=${1:?usage: 14_differential_robustness.sh <results_dir> <data_dir> <kasumi1 dir> [threshold|depth|chicdiff|all]}
data=${2:?}
kasumi=${3:?}
step=${4:-all}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$here/..${PYTHONPATH:+:$PYTHONPATH}" PYTHONHASHSEED=0
python=${PYTHON:-python}
threads=8
out="$results/robustness_differential"
mkdir -p "$out"
calls="$kasumi/chicago"
gtf="$data/annotation/gencode.v29.annotation.gtf.gz"
k="$data/kessler"
md="$k/GSM6505208_Md_CapHiC_Merged_All-Tr_2ndBioRep_Chicago_output_washU_text.txt.gz,$k/GSM6505209_Md_CapHiC_3rdBioRep_Chicago_output_washU_text.txt.gz"
pa2="$k/GSM6505210_PA2_CapHiC_2ndBioRep_Chicago_output_washU_text.txt.gz,$k/GSM6505211_PA2_E10.5wt_CapHiC_3rdBioRep_Tr1_Chicago_output_washU_text.txt.gz"
mm="$calls/siMM_rep1_score1.ibed.gz,$calls/siMM_rep2_score1.ibed.gz"
re="$calls/siRE_rep1_score1.ibed.gz,$calls/siRE_rep2_score1.ibed.gz"

wait_memory() {   # block until at least 30 GB are available
    while (( $(free -g | awk '/^Mem:/{print $7}') < 30 )); do
        echo "$(date -Is) waiting for memory: $(free -g | awk '/^Mem:/{print $7}') GB available" >&2
        sleep 300
    done
}

kasumi_statuses() {   # <output prefix> <variant> [<variant> ...]
    local prefix=$1; shift
    local variants=()
    for v in "$@"; do variants+=(--variant "$v"); done
    wait_memory
    $python "$here/differential_robustness.py" statuses "${variants[@]}" --format ibed --names siMM siRE \
        --rawA ${mm//,/ } --rawB ${re//,/ } --scoreBins 1 3 5 --referenceScore 5 --tss "$gtf" \
        --expression "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" --geneKey genes \
        --storedPrefix "$results/aml_kasumi1/siMM_vs_siRE" -o "$prefix" > "$prefix.log" 2>&1
}

if [[ $step == threshold || $step == all ]]; then
    wait_memory
    $python "$here/differential_robustness.py" statuses --format washu --names Md PA2 \
        --variant "both_5:intersect:none:$md:$pa2" \
        --variant "both_8:intersect:8:$md:$pa2" --variant "both_10:intersect:10:$md:$pa2" \
        --variant "either_5:union:5:$md:$pa2" --variant "either_8:union:8:$md:$pa2" \
        --variant "either_10:union:10:$md:$pa2" \
        --rawA ${md//,/ } --rawB ${pa2//,/ } --scoreBins 5 --referenceScore 5 --tss "$k/mm10.ncbiRefSeq.gtf" \
        --geneTypes --expression "$results/kessler/deseq2_PA2_vs_Md.tsv" --geneKey genes \
        --storedPrefix "$results/kessler/Md_vs_PA2" -o "$out/Md_vs_PA2_threshold" > "$out/Md_vs_PA2_threshold.log" 2>&1
    kasumi_statuses "$out/siMM_vs_siRE_threshold" "all_5:all:5:$mm:$re" "all_3:all:3:$mm:$re" \
        "all_8:all:8:$mm:$re" "all_10:all:10:$mm:$re"
fi

if [[ $step == depth || $step == all ]]; then
    design=$kasumi/design/kasumi1_hg38_promoters
    hicx=${HICX:-$HOME/src/HiCExplorer-v4/cpp/build/tools}
    depth="$out/depth"
    mkdir -p "$depth"
    reads() { awk -F'\t' 'NR>2{n+=$3} END{print n}' "$kasumi/chinput/$1/$1.chinput"; }
    [[ -s "$depth/read_counts.tsv" ]] || for s in siMM_rep1 siMM_rep2 siRE_rep1 siRE_rep2; do
        printf '%s\t%s\n' "$s" "$(reads $s)"; done > "$depth/read_counts.tsv"
    count() { awk -v s="$1" '$1==s{print $2}' "$depth/read_counts.tsv"; }
    recall() {   # <setting> <sample> <fraction>: thinned chinput -> CHiCAGO -> ibed at score >= 1
        local setting=$1 name=$2 fraction=$3 dir="$depth/$1"
        local prefix="$dir/${name}"
        mkdir -p "$dir"
        [[ -s "${prefix}_score1.ibed.gz" ]] && return 0
        wait_memory
        $python "$here/differential_robustness.py" downsample "$kasumi/chinput/$name/$name.chinput" "$fraction" \
            "$prefix.chinput" --seed 0 > "$prefix.downsample.log" 2>&1
        wait_memory
        "$hicx/chicChicagoBackgroundModel" --rmap $design.rmap --baitmap $design.baitmap --chinput "$prefix.chinput" \
            --nperbin $design.npb --nbaitsperbin $design.nbpb --proxOE $design.poe --threads $threads \
            -o "${prefix}_background.txt" > "${prefix}_background.log" 2>&1
        "$hicx/chicChicagoScores" --chinput "$prefix.chinput" --backgroundModel "${prefix}_background.txt" \
            --baitmap $design.baitmap --rmap $design.rmap --threads $threads -o "${prefix}_scores.txt" \
            > "${prefix}_scores.log" 2>&1
        "$hicx/chicChicagoSignificantInteractions" --scores "${prefix}_scores.txt" --rmap $design.rmap \
            --baitmap $design.baitmap --scoreThreshold 1 --threads $threads -o "${prefix}_score1.txt" \
            > "${prefix}_significant.log" 2>&1
        awk -F'\t' -v OFS='\t' 'NR==1{print "bait_chr","bait_start","bait_end","bait_name","otherEnd_chr","otherEnd_start","otherEnd_end","otherEnd_name","N_reads","score"; next} {print $1,$2,$3,$4,$5,$6,$7,".",$8,$10}' \
            "${prefix}_score1.txt" | gzip > "${prefix}_score1.ibed.gz.tmp"
        mv "${prefix}_score1.ibed.gz.tmp" "${prefix}_score1.ibed.gz"
        rm -f "$prefix.chinput" "${prefix}_scores.txt" "${prefix}_score1.txt"
    }
    # (a) siRE thinned to the read count of the siMM replicate of the same number (siMM unchanged)
    for r in 1 2; do
        recall matched "siRE_rep$r" "$(awk -v a="$(count siMM_rep$r)" -v b="$(count siRE_rep$r)" 'BEGIN{printf "%.6f", a/b}')"
    done
    # (b) all four replicates thinned to 50%
    for s in siMM_rep1 siMM_rep2 siRE_rep1 siRE_rep2; do recall half "$s" 0.5; done
    dm="$depth/matched"
    dh="$depth/half"
    kasumi_statuses "$out/siMM_vs_siRE_depth" "all_5:all:5:$mm:$re" \
        "matched_all_5:all:5:$mm:$dm/siRE_rep1_score1.ibed.gz,$dm/siRE_rep2_score1.ibed.gz" \
        "half_all_5:all:5:$dh/siMM_rep1_score1.ibed.gz,$dh/siMM_rep2_score1.ibed.gz:$dh/siRE_rep1_score1.ibed.gz,$dh/siRE_rep2_score1.ibed.gz"
fi

if [[ $step == chicdiff || $step == all ]]; then
    PYTHON=$python bash "$here/kasumi_chicdiff.sh" "$kasumi" "$out/chicdiff" "$results/aml_kasumi1/siMM_vs_siRE_edges.tsv.gz"
fi
