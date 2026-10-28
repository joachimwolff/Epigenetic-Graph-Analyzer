#!/usr/bin/env bash
# Case study from A to Z: naive B cells against the B-lymphoblastoid line GM12878 (docs/case_study.md), with the
# settings of the manuscript (analysis/35_gm12878_naiveB.sh, analysis/03_validate.sh).
#   1. graph of GM12878 and its indirect contacts
#   2. SPRITE of GM12878 on the HindIII fragments, and the multi-way support of the graph by step
#   3. comparison of naive B cells and GM12878: explanation levels for every differentially expressed gene
#   4. integrated tables, main report, gene reports, explorer
#   5. figures: differential graph of NOL11, the region with its SPRITE clusters
#   6. calibration with unchanged genes (one run)
# Inputs from download.sh. Run time on 12 cores about one hour.
#
# Usage: run.sh <data dir> <output dir> [threads=12]
set -euo pipefail
data=${1:?usage: run.sh <data dir> <output dir> [threads]}
out=${2:?}
threads=${3:-12}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"/{graphs,multiway,validation,differential,explain,reports,explorer,figures,calibration}
design="$data/Human_hg19"
gm12878="$data/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz"
naiveb="$data/03_processed_inputs/naiveB_javierre2016/naiveB_nB_washU.txt.gz"
deseq2="$data/03_processed_inputs/expression/deseq2_GM12878_vs_naiveB.tsv"
gtf="$data/gencode.v19.annotation.gtf.gz"
sprite=$(ls "$data"/sprite/human.combined.mapq-ge30.clusters* | head -1)
roadmap() { for mark in H3K27ac H3K4me1 H3K4me3 H3K27me3; do printf '%s=%s ' "$mark" "$data/roadmap/$1-$mark.fc.signal.bigwig"; done; }

echo "1. graph of GM12878"
egaCreateGraph -i "$gm12878" -f washu --baits "$design/Digest_Human_HindIII_baits_e75_ID.baitmap" \
    --scoreThreshold 5 --tss "$gtf" --maxDistance 2000000 --weight score --communityMethod leiden \
    --resolution 1.0 --seed 42 -o "$out/graphs/GM12878" > "$out/graphs/GM12878.log" 2>&1
egaIndirectContacts -g "$out/graphs/GM12878.graphml" --maxDistance 2000000 --triplets -o "$out/graphs/GM12878" \
    > "$out/graphs/GM12878_indirect.log" 2>&1

echo "2. SPRITE and multi-way support"
awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3}' "$design/Digest_Human_HindIII.rmap" > "$out/multiway/hg19_HindIII_fragments.bed"
awk 'BEGIN{OFS="\t"} {if ($3 > end[$1]) end[$1] = $3} END{for (c in end) print c, end[c]}' \
    "$out/multiway/hg19_HindIII_fragments.bed" > "$out/multiway/hg19.chrom.sizes"
egaBuildMultiway --type sprite --input "$sprite" --chromSizes "$out/multiway/hg19.chrom.sizes" \
    --fragments "$out/multiway/hg19_HindIII_fragments.bed" -o "$out/multiway/GM12878_sprite.npz" \
    > "$out/multiway/build.log" 2>&1
egaValidateMultiway -g "$out/graphs/GM12878.graphml" -m "$out/multiway/GM12878_sprite.npz" --confirmPairs \
    --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 --seed 0 --threads "$threads" \
    -o "$out/validation/GM12878_sprite" > "$out/validation/GM12878_sprite.log" 2>&1

echo "3. naive B cells against GM12878"
differential() {  # <expression table> <output prefix>
    egaDifferential --verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1 \
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads" \
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
        --stepFeature decomposed \
        --interactionsA "$naiveb" --interactionsB "$gm12878" --format washu \
        --baits "$design/Digest_Human_HindIII_baits_e75_ID.baitmap" --scoreThreshold 5 --referenceScore 5 \
        --names naiveB GM12878 --tss "$gtf" --tracksA $(roadmap E032) --tracksB $(roadmap E116) \
        --markSign H3K27me3=-1 --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
        --expression "$1" --log2FCColumn log2FoldChange --padjColumn padj --alpha 0.05 --geneKey genes \
        -o "$2" > "$2.log" 2>&1
}
differential "$deseq2" "$out/differential/naiveB_vs_GM12878"

echo "4. tables, reports, explorer"
p="$out/differential/naiveB_vs_GM12878"
egaExplain -p "$p" -o "$out/explain/naiveB_vs_GM12878" > "$out/explain/explain.log" 2>&1
egaReport -p "$p" --annotation "$gtf" --labels "naive B cells" GM12878 --pdf none -o "$out/reports" \
    > "$out/reports/main.log" 2>&1
egaReport -p "$p" --gene NOL11 BRI3BP --annotation "$gtf" --labels "naive B cells" GM12878 --pdf none \
    -o "$out/reports" > "$out/reports/genes.log" 2>&1
egaExplorer -p "$p" --annotation "$gtf" --labels "naive B cells" GM12878 --differentialOnly -o "$out/explorer/" \
    > "$out/explorer/explorer.log" 2>&1

echo "5. figures"
egaPlotDifferential -g "${p}_union.graphml" --gene NOL11 --steps 2 --candidates "${p}_indirect_only_genes.tsv" \
    --mark H3K27ac -o "$out/figures/NOL11_differential.png" > "$out/figures/NOL11.log" 2>&1
egaPlotRegion -g "$out/graphs/GM12878.graphml" --region chr20:30000000-31000000 --noIndirect \
    --multiway "$out/multiway/GM12878_sprite.npz" -o "$out/figures/chr20_region" > "$out/figures/region.log" 2>&1

echo "6. calibration with unchanged genes"
python "$here/../../analysis/pseudo_de.py" "$deseq2" "$out/calibration/unchanged_expression.tsv" --seed 0 \
    > "$out/calibration/unchanged_expression.log"
differential "$out/calibration/unchanged_expression.tsv" "$out/calibration/naiveB_vs_GM12878"
echo "done: $out"
