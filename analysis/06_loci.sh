#!/usr/bin/env bash
# Example regions with several communities.
#   Figure 3 (GM12878) and Supplementary Figure S4 (H1-hESC): interactions as pyGenomeTracks arcs coloured by
#             community, the same region as a promoter-centred graph (egaPlotRegion), and the graph with the
#             multi-way co-occurrence of its interactions and single SPRITE clusters / GAM tubes.
#   Supplementary Figure S6: the SPRITE clusters / GAM tubes that contain loci of these communities and their
#             community purity; S5: a wider GM12878 region; S2, S3: single communities with their multi-way events.
# PGT is the pyGenomeTracks executable (default: pyGenomeTracks on PATH).
set -euo pipefail
data=${1:?usage: 06_loci.sh <data_dir> <results_dir> <figure_dir>}
results=${2:?}
figures=${3:?}
pgt=${PGT:-pyGenomeTracks}
chain="$data/annotation/hg19ToHg38.over.chain.gz"
sprite="$results/multiway/sprite_GM12878_hg19_HindIII_fragments.npz"
gam="$results/multiway/gam_H1_hg38_30kb.npz"
work="$results/regions"
mkdir -p "$figures" "$work"
# name | figure of arcs and graph | graph | multi-way data | region (hg19) | label | extra options
regions=(
    "GM|figure3|GM_chicago|$sprite|chr20:30000000-31000000|GM12878|--noIndirect --multiway $sprite"
    "H1|figureS4|H1|$gam|chr12:52550000-54620000|H1-hESC|--noIndirect --multiway $gam --chain $chain"
)
for entry in "${regions[@]}"; do
    IFS='|' read -r name figure graph multiway region label extra <<< "$entry"
    chrom=${region%%:*}
    span=${region#*:}
    zcat "$data/annotation/gencode.v19.annotation.gtf.gz" |
        awk -F'\t' -v c="$chrom" -v s="${span%-*}" -v e="${span#*-}" \
            '$1==c && $5>=s-300000 && $4<=e+300000 && ($3=="exon"||$3=="transcript"||$3=="gene") &&
             $9 ~ /gene_type "protein_coding"/' > "$work/${name}_genes.gtf"
    egaPlotRegion -g "$results/graphs/$graph.graphml" --region "$region" --genes "$work/${name}_genes.gtf" \
        --title "$label capture Hi-C" --graphSize 7.1 2.7 --geneFontSize 9 --arcHeight 5 --eventsHeight 7 $extra -o "$work/${name}_region"
    # plot from the first to the last fragment with a contact
    "$pgt" --tracks "$work/${name}_region_tracks.ini" --region "$(cat "$work/${name}_region_extent.txt")" \
        --width 30 --fontSize 9 \
        -o "$figures/${figure}_${name}_arcs.pdf"
    cp "$work/${name}_region_graph.pdf" "$figures/${figure}_${name}_graph.pdf"
    # promoter-centred graph with the multi-way evidence
    [[ -s "$work/${name}_region_graph_multiway.pdf" ]] && \
        cp "$work/${name}_region_graph_multiway.pdf" "$figures/${figure}_${name}_graph_multiway.pdf"
    chain_option=()
    [[ "$multiway" == "$gam" ]] && chain_option=(--chain "$chain")
    egaPlotLocus --style communities --panels evidence -g "$results/graphs/$graph.graphml" -m "$multiway" \
        "${chain_option[@]}" --region "$region" --title "$label" -o "$figures/figureS6_${name}_evidence.pdf"
done
# partition of a wider GM12878 region with seven communities
egaPlotLocus --style communities -g "$results/graphs/GM_chicago.graphml" -m "$sprite" \
    --region chr12:14370000-16850000 --title 'GM12878 SPRITE, chr12' -o "$figures/figureS5_GM_partition.pdf"
# community of a gene, looked up in the node table
community_of() {
    awk -F'\t' -v g="$2" 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
        $col["bait"] == "True" && ("," $col["genes"] ",") ~ ("," g ",") {print $col["community"]; exit}' "$1"
}
egaPlotLocus -g "$results/graphs/GM_chicago.graphml" -m "$sprite" \
    --highlightCommunity "$(community_of "$results/graphs/GM_chicago_nodes.tsv" MYOCD)" \
    --title 'GM12878 SPRITE, MYOCD community' -o "$figures/figureS2_locus_sprite.pdf"
egaPlotLocus -g "$results/graphs/H1.graphml" -m "$gam" --chain "$chain" \
    --highlightCommunity "$(community_of "$results/graphs/H1_nodes.tsv" HOXA3)" \
    --title 'H1-hESC GAM, HOXA community' -o "$figures/figureS3_locus_gam.pdf"
