#!/usr/bin/env bash
# Evaluation of the interaction graphs with Pore-C and HiPore-C concatemers, as for SPRITE (03_validate.sh):
#   * GM12878 Pore-C, Deshpande et al. 2022 (GSE149117), GRCh38 fragment alignments of 21 libraries
#   * GM12878 and K562 HiPore-C, Zhong et al. 2023 (GSE202539), hg38 read alignments per flow cell
# Every read becomes one multi-way group (analysis/porec_to_clusters.py; fragments with mapping quality >= 1, lifted
# to hg19), assigned to the hg19 HindIII fragments of the capture design; graphs on GRCh38 (K562) are lifted to hg19.
# Pore-C is ligation-based: a read records several fragments ligated within one crosslinked complex.
set -euo pipefail
data=${1:?usage: 09_porec.sh <data_dir> <results_dir>}
results=${2:?}
here=$(cd "$(dirname "$0")" && pwd)
porec="$data/porec"
chain_to_hg19="$data/annotation/hg38ToHg19.over.chain.gz"
fragments="$data/fragments/hg19_HindIII_fragments.bed"
out="$results/porec"
mkdir -p "$out" "$results/multiway"

convert() {   # <name> <input files...>: per-read clusters, hg19
    local name=$1
    shift
    local clusters="$out/${name}_clusters.txt.gz"
    [[ -s "$clusters" ]] && return 0
    # one process per input file, ${JOBS:-8} at a time
    printf '%s\0' "$@" | xargs -0 -P "${JOBS:-8}" -I{} sh -c \
        'python "$1/porec_to_clusters.py" "$2" "$3.part.$(basename "$2").gz" $4 --chain "$5" --sep , 2>> "$6"' \
        _ "$here" {} "$clusters" "$READ_COLUMNS" "$chain_to_hg19" "$out/${name}_convert.log"
    cat "$clusters".part.*.gz > "$clusters" && rm "$clusters".part.*.gz
}

build() {   # <name>: incidence matrix on the HindIII fragments
    local name=$1
    local npz="$results/multiway/${name}_hg19_HindIII_fragments.npz"
    [[ -s "$npz" ]] || egaBuildMultiway --type sprite --input "$out/${name}_clusters.txt.gz" \
        --chromSizes "$data/annotation/hg19.chrom.sizes" --fragments "$fragments" -o "$npz" \
        > "$out/${name}_build.log" 2>&1
}

evaluate() {   # <graph prefix> <multiway name> <run name> [chain for the graph]
    local graph=$1 name=$2 run=$3
    shift 3
    local npz="$results/multiway/${name}_hg19_HindIII_fragments.npz"
    local chain=${1:-}
    [[ -s "$out/${run}_summary.tsv" ]] || egaValidateMultiway --verbose -g "$graph.graphml" -m "$npz" ${chain:+--chain "$chain"} --confirmPairs \
        --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 --seed 0 \
        -o "$out/$run" > "$out/$run.log" 2>&1
    [[ -s "$out/${run}_hops.tsv" ]] || python "$here/hop_distance.py" "$graph.graphml" "$npz" "$out/${run}_hops.tsv" ${chain:+"$chain"} \
        > "$out/${run}_hops.log" 2>&1
    [[ -s "$out/${run}_hubs_hub_summary.tsv" ]] || python "$here/hub_multiway.py" "$graph.graphml" "$npz" "$out/${run}_hubs" ${chain:+"$chain"} \
        > "$out/${run}_hubs.log" 2>&1
}

READ_COLUMNS="--read read_name --chrom chrom --start start --end end --mapq MapQual --minMapq 1"
convert GM12878_hiporec "$porec"/GM12878_hiporec_zhong2023/*_reads_alignment.csv.gz
convert K562_hiporec "$porec"/K562_hiporec_zhong2023/*_reads_alignment.csv.gz
# pore-c-py fragment tables: alignments that passed the pipeline filter only
READ_COLUMNS="--read read_name --chrom chrom --start start --end end --mapq mapping_quality --minMapq 1 --keep pass_filter=True"
convert GM12878_porec "$porec"/GM12878_porec_deshpande2022/*fragment_alignments.csv.gz
for name in GM12878_hiporec K562_hiporec GM12878_porec; do
    build "$name"
done
evaluate "$results/graphs/GM_chicago" GM12878_porec GM_chicago_porec
evaluate "$results/graphs/GM_chicago" GM12878_hiporec GM_chicago_hiporec
evaluate "$results/k562/K562_dmso" K562_hiporec K562_hiporec "$chain_to_hg19"
for run in GM12878_porec GM12878_hiporec; do
    [[ -s "$out/GM_chicago_${run}_community_pairs_summary.tsv" ]] || python "$here/sprite_community_pairs.py" "$results/multiway/${run}_hg19_HindIII_fragments.npz" \
        "$results/graphs/GM_chicago_nodes.tsv" "$out/GM_chicago_${run}" > "$out/GM_chicago_${run}_pairs.log" 2>&1
done
