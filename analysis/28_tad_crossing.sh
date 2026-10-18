#!/usr/bin/env bash
# TAD boundaries crossed by gained, lost and shared contacts (egaTadCrossing), with real Hi-C TAD calls.
#   1. hicFindTADs on the balanced 25 kb coolers of Kessler et al. 2023 (GSE211901, Md and PA2 E10.5; converted by
#      the kessler_tad_crossing step, chromosome names without 'chr'): --minDepth 75 kb, --maxDepth 250 kb (3-10
#      bins), --step 25 kb, --thresholdComparisons 0.05, --delta 0.01, --correctForMultipleTesting fdr
#   2. egaTadCrossing on the Md / PA2 union graph (hops_n/Md_vs_PA2) with the PA2 and with the Md TADs, every
#      boundary ('all') and the boundaries at or above the median strength ('s50'); HIRE1 / HIRE2 of Kessler et al.
#      2023 name the partner regions; hESC / hNEC (hops_n/hESC_vs_hNEC, hg19) with the H1 hicFindTADs calls
#      (compartments_domains/hic_reference, GRCh38, lifted with hg38ToHg19)
#   3. egaBoundaries (rule-based, without Hi-C) on the Md and PA2 condition graphs with the Kessler CTCF tracks and
#      the CTCF motif (JASPAR MA0139); evaluated against the hicFindTADs boundaries of the same condition
#   4. summary tables, locus figure (one PDF per panel), table_s36.tex, egaExplain / egaReport / egaExplorer --tads
# Inputs: <data_dir> = raw directory (kessler/, genomes/, motifs/, annotation/, roadmap not needed), <results_dir>
# with hops_n/, kessler_tad_crossing/ and compartments_domains/hic_reference/. Outputs: <results>/tad_crossing/.
set -euo pipefail
data=${1:?usage: 28_tad_crossing.sh <data_dir> <results_dir> [threads]}
results=${2:?}
threads=${3:-12}
here=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py=${PYTHON:-python}
hbin=${HICEXPLORER_BIN:-}
out="$results/tad_crossing"
k="$data/kessler/kessler2023_hic_ctcf"
mkdir -p "$out"/{hicfindtads,Md_vs_PA2,hESC_vs_hNEC,egaboundaries,explain,reports,explorer,figures,tables}

# 1. hicFindTADs (Md and PA2 in parallel, half of the threads each)
half=$(( threads / 2 > 0 ? threads / 2 : 1 ))
for c in Md PA2; do
    if [ ! -s "$out/hicfindtads/${c}_25kb_boundaries.bed" ]; then
        "${hbin}hicFindTADs" -m "$results/kessler_tad_crossing/${c}_25kb_25000.cool" \
            --outPrefix "$out/hicfindtads/${c}_25kb" --minDepth 75000 --maxDepth 250000 --step 25000 \
            --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr \
            --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X \
            --numberOfProcessors "$half" > "$out/hicfindtads/${c}_hicFindTADs.log" 2>&1 &
    fi
done
wait

# HIRE1 / HIRE2 (Kessler et al. 2023, mm10)
printf 'chr6\t50789172\t50828639\tHIRE2\nchr6\t50913170\t51087888\tHIRE1\n' > "$out/kessler2023_HIRE_mm10.bed"

# 2. egaTadCrossing
for tads in PA2 Md; do
    for setting in all s50; do
        strength=0; [ "$setting" = s50 ] && strength=0.5
        "$py" -m epigeneticgraphanalyzer.egaTadCrossing --verbose -p "$results/hops_n/Md_vs_PA2" \
            --tads "$out/hicfindtads/${tads}_25kb" --tadRegions "$out/kessler2023_HIRE_mm10.bed" \
            --tadMinStrength "$strength" --gene Hoxa2 Hoxa3 Hoxa1 \
            -o "$out/Md_vs_PA2/${tads}_tads_${setting}" > "$out/Md_vs_PA2/${tads}_tads_${setting}.log" 2>&1
    done
done
for setting in all s50; do
    strength=0; [ "$setting" = s50 ] && strength=0.5
    "$py" -m epigeneticgraphanalyzer.egaTadCrossing --verbose -p "$results/hops_n/hESC_vs_hNEC" \
        --tads "$results/compartments_domains/hic_reference/H1_50kb" \
        --tadChain "$data/annotation/hg38ToHg19.over.chain.gz" --tadMinStrength "$strength" --gene EPHA7 \
        -o "$out/hESC_vs_hNEC/H1_tads_${setting}" > "$out/hESC_vs_hNEC/H1_tads_${setting}.log" 2>&1
done

# 3. egaBoundaries on the condition graphs (no Hi-C), evaluated against hicFindTADs
"$py" "$here/tad_crossing_inputs.py" "$results/hops_n/Md_vs_PA2" \
    "$k/GSM6505197_E10_Md_CTCF_normalized_by_autosomal_lib_1e6.scaled_0100.bw" \
    "$out/hicfindtads/Md_25kb_boundaries.bed" "$out/hicfindtads/PA2_25kb_boundaries.bed" "$out/egaboundaries/inputs"
declare -A ctcf=([Md]="$k/GSM6505197_E10_Md_CTCF_normalized_by_autosomal_lib_1e6.scaled_0100.bw"
                 [PA2]="$k/GSM6505196_E10_PA2_CTCF_normalized_by_autosomal_lib_1e6.scaled_0100.bw")
for c in Md PA2; do
    "$py" -m epigeneticgraphanalyzer.egaBoundaries --verbose --nodes "$out/egaboundaries/inputs/${c}_nodes.tsv" \
        --edges "$out/egaboundaries/inputs/${c}_edges.tsv" --chromSizes "$out/egaboundaries/inputs/chrom.sizes" \
        --binSize 25000 --ctcf "${ctcf[$c]}" --fasta "$data/genomes/mm10.fa.gz" \
        --motifs "$data/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt" \
        --reference "$out/egaboundaries/inputs/${c}_25kb_boundaries_chr.bed" \
        -o "$out/egaboundaries/${c}_ega" > "$out/egaboundaries/${c}_ega.log" 2>&1 &
done
wait

# 4. summaries, figure, table, displays
"$py" "$here/tad_crossing_summary.py" "$out" > "$out/summary.log" 2>&1
"$py" "$here/figure_tad_crossing.py" "$out" "$k" "$results/hops_n/Md_vs_PA2"
"$py" "$here/table_tad_crossing.py" "$out" "$out/tables"
md=(--tads "$out/hicfindtads/PA2_25kb" --tadRegions "$out/kessler2023_HIRE_mm10.bed" --tadMinStrength 0.5)
hs=(--tads "$results/compartments_domains/hic_reference/H1_50kb"
    --tadChain "$data/annotation/hg38ToHg19.over.chain.gz" --tadMinStrength 0.5)
"$py" -m epigeneticgraphanalyzer.egaExplain -p "$results/hops_n/Md_vs_PA2" "${md[@]}" -o "$out/explain/Md_vs_PA2"
"$py" -m epigeneticgraphanalyzer.egaExplain -p "$results/hops_n/hESC_vs_hNEC" "${hs[@]}" \
    -o "$out/explain/hESC_vs_hNEC"
"$py" -m epigeneticgraphanalyzer.egaReport -p "$results/hops_n/Md_vs_PA2" --gene Hoxa2 Hoxa3 "${md[@]}" \
    --annotation "$data/kessler/mm10.ncbiRefSeq.gtf" -o "$out/reports"
"$py" -m epigeneticgraphanalyzer.egaReport -p "$results/hops_n/hESC_vs_hNEC" --gene EPHA7 "${hs[@]}" \
    --annotation "$data/annotation/gencode.v19.annotation.gtf.gz" -o "$out/reports"
"$py" -m epigeneticgraphanalyzer.explorer -p "$results/hops_n/Md_vs_PA2" "${md[@]}" \
    --annotation "$data/kessler/mm10.ncbiRefSeq.gtf" -o "$out/explorer/"
