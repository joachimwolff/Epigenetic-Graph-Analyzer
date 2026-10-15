#!/usr/bin/env bash
# Differential mode on primary AML (Assi et al. 2019, GSE108316, hg38): t(8;21) against FLT3-ITD.
#   capture Hi-C: HOMER interaction calls of one patient per subtype (FDR <= 0.01), baits of the HindIII promoter
#   capture design (Mifsud et al. 2015) lifted to hg38
#   chromatin: DNase-seq of 4 t(8;21) and 2 FLT3-ITD patients (mean per subtype)
#   expression: RNA-seq coverage of 3 t(8;21) and 3 FLT3-ITD patients (analysis/assi_expression.py); with three
#   patients per subtype no gene passes an adjusted threshold, so differential genes are those with an unadjusted
#   Welch p < 0.05 and |log2FC| >= 1 (exploratory)
set -euo pipefail
assi=${1:?usage: 07_aml_assi.sh <assi2019 download dir> <data_dir> <results_dir> [permutations]}
data=${2:?}
results=${3:?}
permutations=${4:-200}
here=$(cd "$(dirname "$0")" && pwd)
gtf="$data/annotation/gencode.v29.annotation.gtf.gz"
sizes="$data/annotation/hg38.chrom.sizes"
baits="$data/chicago_calls/design_hg38_HindIII_lifted/baits_hg38_lifted.bed"
out="$results/aml_assi"
bw="$assi/bigwig"
mkdir -p "$out" "$bw"
to_bigwig() {
    local input=$1 output=$2
    [[ -s "$output" ]] && return 0
    # variableStep WIG (span 20) and bedGraph to sorted bedGraph, then bigWig
    zcat "$input" | awk -v OFS='\t' '
        /^track/ {next}
        /^variableStep/ {split($2, c, "="); chrom = c[2]; span = 1; for (i = 3; i <= NF; i++) if ($i ~ /^span=/) {split($i, s, "="); span = s[2]}; wig = 1; next}
        wig {print chrom, $1 - 1, $1 - 1 + span, $2; next}
        NF == 4 {print $1, $2, $3, $4}' | awk '$1 ~ /^chr([0-9]+|X|Y)$/' | LC_ALL=C sort -k1,1 -k2,2n -S 8G > "$output.bedGraph"
    bedGraphToBigWig "$output.bedGraph" "$sizes" "$output"
    rm "$output.bedGraph"
}
for f in "$assi"/*_DHS.wig.gz; do
    to_bigwig "$f" "$bw/$(basename "$f" .wig.gz).bw"
done
# RNA bedGraphs are stranded (minus strand as negative values): one bigWig per strand
for f in "$assi"/*.bedgraph.gz; do
    name=$(basename "$f" .bedgraph.gz)
    for strand in plus minus; do
        track="$bw/$name.$strand.bw"
        [[ -s "$track" ]] && continue
        zcat "$f" | awk -v OFS='\t' -v strand=$strand '!/^track/ && $1 ~ /^chr([0-9]+|X|Y)$/ {
            if (strand == "plus" && $4 > 0) print $1, $2, $3, $4;
            if (strand == "minus" && $4 < 0) print $1, $2, $3, -$4}' | LC_ALL=C sort -k1,1 -k2,2n -S 8G > "$track.bedGraph"
        bedGraphToBigWig "$track.bedGraph" "$sizes" "$track"
        rm "$track.bedGraph"
    done
done
join() { local IFS=,; echo "$*"; }
t821_dhs=$(join "$bw"/GSM289362[789]_t821*_DHS.bw "$bw"/GSM2893630_t821KIT-2_DHS.bw)
itd_dhs=$(join "$bw"/GSM289360[67]_ITD-*_DHS.bw)
stranded() { for f in "$@"; do printf '%s+%s\n' "${f%.plus.bw}.plus.bw" "${f%.plus.bw}.minus.bw"; done; }
# t821.ITDx2-5 carries both t(8;21) and FLT3-ITD and is left out of the contrast
t821_rna=$(join $(stranded "$bw"/GSM289406[78]_t821*.plus.bw "$bw"/GSM2894070_t821.KIT-2.plus.bw))
itd_rna=$(join $(stranded "$bw"/GSM289404[125]_ITD-*.plus.bw))
[[ -s "$out/expression_FLT3ITD_vs_t821.tsv" ]] || python "$here/assi_expression.py" "$gtf" \
    "$out/expression_FLT3ITD_vs_t821.tsv" "A=$t821_rna" "B=$itd_rna"
egaDifferential --verbose \
    --interactionsA "$assi/GSM2895201_t821-1_CHiC_interactions.txt.gz" \
    --interactionsB "$assi/GSM2895200_FLT3-ITD_CHiC_interactions.txt.gz" \
    --format homer --baits "$baits" --maxFdr 0.01 --referenceScore 2 --names t821 FLT3ITD \
    --tss "$gtf" --tracksA "DHS=$t821_dhs" --tracksB "DHS=$itd_dhs" \
    --expression "$out/expression_FLT3ITD_vs_t821.tsv" --log2FCColumn log2FoldChange --padjColumn pvalue \
    --alpha 0.05 --geneKey genes --rnaThreshold 1 --chromatinThreshold 1 --permutations "$permutations" \
    --bootstrap 1000 --seed 0 -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1
