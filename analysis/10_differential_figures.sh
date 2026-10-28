#!/usr/bin/env bash
# Differential graph figures, one PDF per panel (egaPlotDifferential --separatePanels):
#   <name>_arcs: gained, lost and shared interactions as arcs; <name>_A, _B: the graphs of the two conditions;
#   <name>_edges: change of the interactions; <name>_nodes_<mark>: change of expression and of each chromatin mark.
#   Figure 4: HOXB community, hESC to hNEC (chr17, hg19)
#   Figure 5: Hoxa3 neighbourhood, Md to PA2 (mouse, mm10)
#   Figure 6: MS4A community, Kasumi-1 siMM to siRE (chr11, hg38)
#   Supplementary Figure S2: FGF9 (hESC to hNEC), CEBPE and CCNA1 (Kasumi-1)
set -euo pipefail
results=${1:?usage: 10_differential_figures.sh <results_dir> <figure_dir>}
figures=${2:?}
mkdir -p "$figures"
plot() {
    local name=$1
    shift
    egaPlotDifferential "$@" --sideBySide --separatePanels -o "$figures/$name.pdf"
}
plot figure4_HOXB -g "$results/differential/hESC_vs_hNEC_union.graphml" --region chr17:46450000-46950000 \
    --candidates "$results/differential/hESC_vs_hNEC_indirect_only_genes.tsv"
plot figure5_Hoxa3 -g "$results/kessler/Md_vs_PA2_union.graphml" --gene Hoxa3 --steps 2 \
    --candidates "$results/kessler/Md_vs_PA2_indirect_only_genes.tsv"
plot figure6_MS4A -g "$results/aml_kasumi1/siMM_vs_siRE_union.graphml" --gene MS4A3 --steps 1 \
    --candidates "$results/aml_kasumi1/siMM_vs_siRE_indirect_only_genes.tsv"
plot figureS2a_FGF9 -g "$results/differential/hESC_vs_hNEC_union.graphml" --gene FGF9 --steps 2 \
    --candidates "$results/differential/hESC_vs_hNEC_indirect_only_genes.tsv"
plot figureS2b_CEBPE -g "$results/aml_kasumi1/siMM_vs_siRE_union.graphml" --gene CEBPE --steps 2
plot figureS2c_CCNA1 -g "$results/aml_kasumi1/siMM_vs_siRE_union.graphml" --gene CCNA1 --steps 2 \
    --candidates "$results/aml_kasumi1/siMM_vs_siRE_indirect_only_genes.tsv"
