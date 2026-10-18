#!/usr/bin/env bash
# Validation with ligation-free multi-way data.
#   GM12878 graphs (CHiCAGO calls, Mifsud strong core) against GM12878 SPRITE (hg19), SPRITE reads assigned to the
#   HindIII fragments of the capture design
#   embryonic stem cell graphs (H1 from analysis/01_h1_chicago.sh when present, H9 hESC of Freire-Pritchett et al.)
#   against H1 GAM (hg38, 30 kb; graph coordinates lifted from hg19)
#   controls: the GM12878 graph against H1 GAM, the stem cell graphs against GM12878 SPRITE, and the
#   hESC-derived neuroectoderm graph against H1 GAM
#   every run also tests each direct and indirect pair individually (--confirmPairs)
set -euo pipefail
data=${1:?usage: 03_validate.sh <data_dir> <results_dir> [threads=24]}
results=${2:?}
threads=${3:-24}
chain="$data/annotation/hg19ToHg38.over.chain.gz"
sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
gam="$results/multiway/gam_H1_hg38_30kb.npz"
mkdir -p "$results/multiway" "$results/validation" "$data/fragments"
[[ -s "$data/fragments/hg19_HindIII_fragments.bed" ]] || \
    python "$(dirname "$0")/../scripts/digest_genome.py" "$data/annotation/hg19.fa.gz" AAGCTT 1 \
        "$data/fragments/hg19_HindIII_fragments.bed" '^chr([0-9]+|X)$'
[[ -s "$sprite" ]] || egaBuildMultiway --type sprite --input "$data/sprite/human.combined.mapq-ge30.clusters.gz" \
    --chromSizes "$data/annotation/hg19.chrom.sizes" --fragments "$data/fragments/hg19_HindIII_fragments.bed" \
    -o "$sprite"
[[ -s "$gam" ]] || egaBuildMultiway --type gam \
    --input "$data"/gam_h1/H1_GAM_rep1_30kb.*.csv "$data"/gam_h1/H1_GAM_rep2_30kb.*.csv \
    --chromSizes "$data/annotation/hg38.chrom.sizes" -o "$gam"
run() {   # <graph> <multiway> <name> <worker processes> [options]
    local graph=$1 multiway=$2 name=$3 workers=$4
    shift 4
    [[ -s "$results/graphs/$graph.graphml" ]] || return 0
    egaValidateMultiway --verbose -g "$results/graphs/$graph.graphml" -m "$multiway" "$@" --confirmPairs \
        --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 --seed 0 --threads "$workers" \
        -o "$results/validation/$name" > "$results/validation/$name.log" 2>&1
}
# GAM runs are light and run four at a time; SPRITE runs hold the full cluster matrix and run two at a time; the
# worker processes of each run share the matrix (results do not depend on their number)
quarter=$(( threads / 4 > 0 ? threads / 4 : 1 ))
half=$(( threads / 2 > 0 ? threads / 2 : 1 ))
run H1 "$gam" H1_gam "$quarter" --chain "$chain" &
run hESC "$gam" hESC_gam "$quarter" --chain "$chain" &
run hNPC "$gam" hNEC_gam_crossstate "$quarter" --chain "$chain" &
run GM_chicago "$gam" GM_chicago_gam_crosscell "$quarter" --chain "$chain" &
wait
run GM_chicago "$sprite" GM_chicago_sprite "$half" &
run GM_mifsud_c20_r1 "$sprite" GM_mifsud_sprite "$half" &
wait
run hESC "$sprite" hESC_sprite_crosscell "$half" &
run H1 "$sprite" H1_sprite_crosscell "$half" &
wait
