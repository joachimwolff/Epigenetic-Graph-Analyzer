#!/usr/bin/env bash
# Time course of keratinocyte differentiation (Rubin et al. 2017, GSE84662, hg19; days 0, 3 and 6, two donors).
#   1. capture Hi-C contacts (GSE84660 b2b + b2g, CHiCAGO FDR < 0.01 at any day; b2g restricted to other ends on an
#      H3K27ac peak at any day) as one ibed per donor and day with a count-derived score (rubin_chic_to_ibed.py),
#      graphs per day with the replicate rule 'all' (both donors, score >= 5 in at least one)
#   2. H3K27ac and SMC1A bigWigs (bedtools coverage) divided by their genome-wide mean (scale_bigwig.py)
#   3. PyDESeq2 ~ donor + day: day 3 vs 0, day 6 vs 3, day 6 vs 0 (rubin_deseq2.py)
#   4. egaDifferential for the three pairs (maxHops 3; H3K27ac active, SMC1A architectural; learned roles too)
#   5. egaTriggers on day 3 vs day 0 (JASPAR 2024 CORE vertebrates, hg19)
#   6. egaTimecourse: trajectories, ordering of contact / mark / expression changes against a within-gene permutation
#      null, H3K27ac at the explaining nodes of candidate-trigger genes in the siRNA knockdowns (KLF4, ZNF750, EHF) at
#      day 3, UMI-4C viewpoints; figures (one PDF per panel) and table_s33.tex
# Outputs: <results>/timecourse_rubin2017/ only. Compartments / TADs of the Hi-C matrices are not part of this script.
set -euo pipefail
data=${1:?usage: 22_timecourse.sh <data_dir> <results_dir> [threads] [permutations]}
results=${2:?}
threads=${3:-12}
permutations=${4:-200}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
py=${PYTHON:-python}
raw="$data/rubin2017_keratinocyte"
out="$results/timecourse_rubin2017"
mkdir -p "$out"/{chic,bigwig,expression,graphs,differential,triggers}

# 1-3: inputs
[[ -s "$out/chic/D0_B1.ibed.gz" ]] || "$py" "$here/rubin_chic_to_ibed.py" "$raw/GSE84660_CHiC_contacts_b2b.bed.gz" \
    "$raw/GSE84660_CHiC_contacts_b2g.bed.gz" "$out/chic" > "$out/chic/conversion.log" 2>&1
for f in "$raw"/GSE84657_RAW/GSM22472{41,42,43,44,45,46,50,51,52,53,54,55}_*.bigWig; do
    b=$(basename "$f" .bigWig)
    [[ -s "$out/bigwig/$b.scaled.bw" ]] || { "$py" "$here/scale_bigwig.py" "$f" "$out/bigwig/$b.scaled.bw.tmp" \
        && mv "$out/bigwig/$b.scaled.bw.tmp" "$out/bigwig/$b.scaled.bw"; }
done
[[ -s "$out/expression/deseq2_D6_vs_D0.tsv" ]] || "$py" "$here/rubin_deseq2.py" \
    "$raw/GSE84661_GV19_TSS_RNAseq_RSEM_counts.bed.gz" "$out/expression" > "$out/expression/deseq2.log" 2>&1

tss="$data/annotation/gencode.v19.annotation.gtf.gz"
bw="$out/bigwig"
declare -A h3k27ac=([D0]=GSM2247241 [D3]=GSM2247242 [D6]=GSM2247243)
declare -A h3k27ac2=([D0]=GSM2247244 [D3]=GSM2247245 [D6]=GSM2247246)
declare -A smc1a=([D0]=GSM2247250 [D3]=GSM2247251 [D6]=GSM2247252)
declare -A smc1a2=([D0]=GSM2247253 [D3]=GSM2247254 [D6]=GSM2247255)
day_number() { echo "${1#D}"; }
tracks() {
    local d=$1 n
    n=$(day_number "$d")
    printf 'H3K27ac=%s,%s SMC1A=%s,%s' \
        "$bw/${h3k27ac[$d]}_H3K27ac_ChIPseq_Donor1_Day$n.scaled.bw" "$bw/${h3k27ac2[$d]}_H3K27ac_ChIPseq_Donor2_Day$n.scaled.bw" \
        "$bw/${smc1a[$d]}_SMC1A_ChIPseq_Donor1_Day$n.scaled.bw" "$bw/${smc1a2[$d]}_SMC1A_ChIPseq_Donor2_Day$n.scaled.bw"
}

# graphs per day
for d in D0 D3 D6; do
    [[ -s "$out/graphs/$d.graphml" ]] || "$py" -m epigeneticgraphanalyzer.egaCreateGraph --verbose \
        --interactions "$out/chic/${d}_B1.ibed.gz" "$out/chic/${d}_B2.ibed.gz" --format ibed \
        --baits "$out/chic/baits.bed" --replicates all --scoreThreshold 5 --tss "$tss" --weight score \
        --communityMethod leiden --resolution 1.0 --seed 42 -o "$out/graphs/$d" > "$out/graphs/$d.log" 2>&1
done

# 4: pairwise differential runs, three in parallel
per=$(( threads / 3 > 0 ? threads / 3 : 1 ))
run_pair() {
    local a=$1 b=$2
    [[ -s "$out/differential/${a}_vs_${b}_explanations_all_levels.tsv" ]] && return 0
    # shellcheck disable=SC2046
    "$py" -m epigeneticgraphanalyzer.egaDifferential --verbose --maxHops 3 --genePermutations 200000 \
        --rnaThreshold 1 --chromatinThreshold 1 --permutations "$permutations" --bootstrap 1000 --seed 0 \
        --threads "$per" --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 \
        --stateWindow 1000 \
        --interactionsA "$out/chic/${a}_B1.ibed.gz" "$out/chic/${a}_B2.ibed.gz" \
        --interactionsB "$out/chic/${b}_B1.ibed.gz" "$out/chic/${b}_B2.ibed.gz" \
        --format ibed --baits "$out/chic/baits.bed" --replicates all --scoreThreshold 5 --referenceScore 5 \
        --names "$a" "$b" --tss "$tss" --tracksA $(tracks "$a") --tracksB $(tracks "$b") \
        --markRoles H3K27ac=active SMC1A=architectural \
        --expression "$out/expression/deseq2_${b}_vs_${a}.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
        --alpha 0.05 --geneKey genes -o "$out/differential/${a}_vs_${b}" > "$out/differential/${a}_vs_${b}.log" 2>&1
}
run_pair D0 D3 &
p1=$!
run_pair D3 D6 &
p2=$!
run_pair D0 D6 &
p3=$!
wait $p1 $p2 $p3

# 5: candidate trigger factors of day 3 vs day 0
genome="$out/triggers/hg19.fa"
[[ -s "$genome" ]] || zcat "$data/annotation/hg19.fa.gz" > "$genome"
[[ -s "$out/triggers/D0_vs_D3_gene_triggers.tsv" ]] || "$py" -m epigeneticgraphanalyzer.egaTriggers --verbose \
    -p "$out/differential/D0_vs_D3" --genome "$genome" \
    --motifs "$data/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    --windowTracks "$bw"/GSM22472{41,42,44,45}_H3K27ac_*.scaled.bw \
    --expression "$out/expression/deseq2_D3_vs_D0.tsv" --window 1000 --pvalue 1e-4 --alpha 0.05 --maxHops 3 \
    --chromatinThreshold 1 --roles declared --presence mean --backgroundRatio 3 --topPerGene 25 \
    --threads "$threads" --seed 0 --cache -o "$out/triggers/D0_vs_D3" > "$out/triggers/D0_vs_D3.log" 2>&1

[[ ${SKIP_TIMECOURSE:-0} == 1 ]] && exit 0

# 6: time-course analysis
kd="$raw/GSE84662_RAW"
h3k27ac_day() { local n; n=$(day_number "$1"); printf '%s=%s,%s' "$1" \
    "$bw/${h3k27ac[$1]}_H3K27ac_ChIPseq_Donor1_Day$n.scaled.bw" "$bw/${h3k27ac2[$1]}_H3K27ac_ChIPseq_Donor2_Day$n.scaled.bw"; }
"$py" -m epigeneticgraphanalyzer.egaTimecourse --verbose --days D0 D3 D6 \
    --differential "$out/differential" --expressionDir "$out/expression" \
    --interactions D0="$out/chic/D0_B1.ibed.gz,$out/chic/D0_B2.ibed.gz" D3="$out/chic/D3_B1.ibed.gz,$out/chic/D3_B2.ibed.gz" \
                   D6="$out/chic/D6_B1.ibed.gz,$out/chic/D6_B2.ibed.gz" \
    --format ibed --baits "$out/chic/baits.bed" --replicates all --scoreThreshold 5 --referenceScore 5 \
    --markTracks "$(h3k27ac_day D0)" "$(h3k27ac_day D3)" "$(h3k27ac_day D6)" --orderMark H3K27ac \
    --maxHops 3 --rnaThreshold 1 --alpha 0.05 --chromatinThreshold 1 --contactThreshold 1 \
    --sensitivityThresholds 0.5 1.5 \
    --triggers "$out/triggers/D0_vs_D3" \
    --motifs "$data/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
    --knockdownControl "$kd/GSM2496159_H3K27ac_ChIPseq_KD_CTRi_Rep1.bw" "$kd/GSM2496160_H3K27ac_ChIPseq_KD_CTRi_Rep2.bw" \
    --knockdown "KLF4=$kd/GSM2496161_H3K27ac_ChIPseq_KD_KLF4i_Rep1.bw,$kd/GSM2496162_H3K27ac_ChIPseq_KD_KLF4i_Rep2.bw" \
                "ZNF750=$kd/GSM2496163_H3K27ac_ChIPseq_KD_ZNF750i_Rep1.bw,$kd/GSM2496164_H3K27ac_ChIPseq_KD_ZNF750i_Rep2.bw" \
                "EHF=$kd/GSM2496165_H3K27ac_ChIPseq_KD_EHFi_Rep1.bw,$kd/GSM2496166_H3K27ac_ChIPseq_KD_EHFi_Rep2.bw" \
    --umi4c "$kd"/GSM24961{43,44,45,46,47,48,49,50,51,52}_UMI4C_*.bed.gz \
    --permutations 10000 --bootstrap 2000 --seed 0 --threads "$threads" \
    -o "$out/timecourse" > "$out/timecourse.log" 2>&1

# sensitivity: H3K27ac of donor 2 only in the ordering test (the day 0 H3K27ac track of donor 1 correlates with every
# other H3K27ac track at Spearman ~0.64 against ~0.9 among the others); explaining nodes as in the main run
h3k27ac_donor2() { local n; n=$(day_number "$1"); printf '%s=%s' "$1" \
    "$bw/${h3k27ac2[$1]}_H3K27ac_ChIPseq_Donor2_Day$n.scaled.bw"; }
"$py" -m epigeneticgraphanalyzer.egaTimecourse --verbose --days D0 D3 D6 \
    --differential "$out/differential" --expressionDir "$out/expression" \
    --interactions D0="$out/chic/D0_B1.ibed.gz,$out/chic/D0_B2.ibed.gz" D3="$out/chic/D3_B1.ibed.gz,$out/chic/D3_B2.ibed.gz" \
                   D6="$out/chic/D6_B1.ibed.gz,$out/chic/D6_B2.ibed.gz" \
    --format ibed --baits "$out/chic/baits.bed" --replicates all --scoreThreshold 5 --referenceScore 5 \
    --markTracks "$(h3k27ac_donor2 D0)" "$(h3k27ac_donor2 D3)" "$(h3k27ac_donor2 D6)" --orderMark H3K27ac \
    --maxHops 3 --rnaThreshold 1 --alpha 0.05 --chromatinThreshold 1 --contactThreshold 1 \
    --sensitivityThresholds 0.5 1.5 --permutations 10000 --bootstrap 2000 --seed 0 --threads "$threads" \
    -o "$out/timecourse_donor2" > "$out/timecourse_donor2.log" 2>&1

# figures (one PDF per panel) and table_s33.tex
"$py" "$here/figure_timecourse.py" "$out" "$out/figures" > "$out/figures.log" 2>&1
"$py" "$here/table_timecourse.py" "$out" "$out/tables" > "$out/tables.log" 2>&1
