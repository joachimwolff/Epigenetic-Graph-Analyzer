# Two conditions

`egaDifferential` compares two conditions A and B: it builds the union graph, attaches the marks of both conditions,
evaluates every differentially expressed gene at every explanation level and writes the tables of the reports. The
[tutorial](../tutorial.md) runs it on a small data set; this page explains the choices.

## Minimal call

```bash
egaDifferential --interactionsA A_rep1.ibed A_rep2.ibed --interactionsB B_rep1.ibed B_rep2.ibed \
    --format ibed --scoreThreshold 5 --names A B --tss genes.gtf.gz \
    --tracksA H3K27ac=A_H3K27ac.bw H3K27me3=A_H3K27me3.bw \
    --tracksB H3K27ac=B_H3K27ac.bw H3K27me3=B_H3K27me3.bw \
    --markRoles H3K27ac=active H3K27me3=repressive \
    --expression deseq2_B_vs_A.tsv --log2FCColumn log2FoldChange --padjColumn padj \
    --geneKey genes -o differential/A_vs_B
```

## Settings of the manuscript

```bash
egaDifferential ... --maxHops 3 --stepFeature decomposed --genePermutations 200000 \
    --rnaThreshold 1 --chromatinThreshold 1 --permutations 200 --bootstrap 1000 --seed 0 --threads 8 \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000 \
    --referenceScore 5
```

## Choices

| Question | Option |
|---|---|
| Which genes are differentially expressed? | `--log2FCColumn`, `--padjColumn`, `--alpha`, `--rnaThreshold` (log2); or `--expressionColumns A B` with `--minExpression` |
| When does a mark change count? | `--chromatinThreshold` (log2 change, default 1, two-fold) |
| What does a mark mean? | `--markRoles mark=active\|enhancer\|repressive\|architectural\|neutral`; `--markSign mark=-1` as a short form; `--roles learned` or `both` learns the roles per level |
| How far from the promoter? | `--maxHops` (default 3); the calibration shows that the third step rarely exceeds the rate of unchanged genes |
| Marks or contacts at the steps? | `--stepFeature decomposed` (recommended) separates the change of the marks from the change of the contacts; `combined`, `marks`, `contacts` |
| How are neighbours weighted? | `--stepWeighting product` (contact-weighted sum), `unweighted`, `normalised` (contact-weighted mean) |
| Chromatin states | `--stateThreshold` (fold enrichment for presence), `--stateTracksA/B` or `--statePeaksA/B`, `--stateScale genomeMean` for coverage tracks |
| Replicates | `--replicates all\|intersect\|union` |
| Threshold flags | give the full caller output (including pairs below the threshold) and `--supportFloor`; `--chicdiff` adds a Chicdiff table |
| Domain boundaries | `--tads` adds the boundary-crossing tables (as `egaTadCrossing`) |
| External elements | `--elements` (for example ENCODE cCREs or ChromHMM states), `--annotation` for element classes |
| Reproducibility | `--seed`; the results do not depend on `--threads` |

## Outputs

| File | Content |
|---|---|
| `<prefix>_explanation_levels.tsv` | per level: genes with a concordant observation (any, first), with *P* < 0.05 beyond the mirrored loci, with *q* < 0.05, threshold flags |
| `<prefix>_explanations_all_levels.tsv` | every promoter fragment at every level: result, strongest node, concordant marks, *P*, *q*, `beyond_mirror`, a text |
| `<prefix>_explanation_cascade.tsv` | first explaining level per gene |
| `<prefix>_indirect_only_genes.tsv` | genes explained only beyond step 1, with the explaining nodes |
| `<prefix>_edges.tsv.gz` | union edges: status (gained, lost, shared), strengths, threshold support |
| `<prefix>_union.graphml` | annotated union graph (input of `egaPlotDifferential`, `egaReport`, `egaExplorer`) |
| `<prefix>_mark_roles.tsv` | declared and learned roles with bootstrap intervals |
| `<prefix>_state_transitions.tsv`, `_bivalent_transitions.tsv` | chromatin state transitions of the promoters |
| `<prefix>_nested_models.tsv` | variance of expression changes explained by promoter, step 1, 2, 3 (population level) |
| `<prefix>_community_*.tsv` | communities of both conditions, their matching and coordinated expression |
| `<prefix>.log` | the run log with the command line, read by the reports |

## Next

* `egaExplain -p differential/A_vs_B -o explain/A_vs_B` joins the tables per gene, element and community.
* [Reports and plots](reports.md).
* [Calibration](calibration.md) gives the rate of explanations for unchanged genes.
