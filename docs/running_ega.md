# Running EGA from start to end

The commands below run a differential comparison of two conditions A and B from called promoter interactions,
signal tracks of histone marks and a table of differential expression, with the settings of the manuscript runs.
File names are placeholders. `REPRODUCE.md` maps every figure and table of the manuscript to its script, and
`example/run_example.sh` runs a complete comparison on the Hoxa locus in about one minute.

## Installation

```bash
mamba env create -f environment.yml && conda activate ega && pip install .
bash example/run_example.sh          # Md / PA2 comparison at the Hoxa locus
```

## 1. One graph per condition

Interactions called by CHiCAGO (WashU text; ibed, BEDPE and the table of Mifsud et al. 2015 are also read),
replicates intersected, communities for display.

```bash
for c in A B; do
  egaCreateGraph -i ${c}_rep1_washU.txt.gz ${c}_rep2_washU.txt.gz -f washu \
      --replicates intersect --scoreThreshold 5 --baits design.baitmap --tss genes.gtf.gz \
      --maxDistance 2000000 --weight score --communityMethod leiden --resolution 1.0 \
      --seed 42 -o graphs/$c
done
```

## 2. Indirect contacts

Pairs up to the maximal distance with their mediators, Adamic-Adar index, path score and triplets.

```bash
egaIndirectContacts -g graphs/B.graphml --maxDistance 2000000 --triplets -o graphs/B
```

## 3. Multi-way data

SPRITE clusters (and Pore-C concatemers in the same cluster format) on the restriction fragments, or GAM in bins,
and the evaluation of direct and indirect contacts against control pairs.

```bash
egaBuildMultiway --type sprite --input clusters.txt.gz --chromSizes genome.chrom.sizes \
    --fragments HindIII_fragments.bed -o multiway/B_sprite.npz
egaValidateMultiway -g graphs/B.graphml -m multiway/B_sprite.npz --confirmPairs \
    --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 \
    --seed 0 --threads 8 -o validation/B_sprite
```

## 4. Differential mode

Union graph, explanation levels (promoter, step 1 to step 3), mark roles, chromatin states, mirror control and
per-gene permutation test.

```bash
egaDifferential --maxHops 3 --stepFeature decomposed --genePermutations 200000 \
    --rnaThreshold 1 --chromatinThreshold 1 --permutations 200 --bootstrap 1000 \
    --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 \
    --stateWindow 1000 --seed 0 --threads 8 \
    --interactionsA A_rep1_washU.txt.gz A_rep2_washU.txt.gz \
    --interactionsB B_rep1_washU.txt.gz B_rep2_washU.txt.gz \
    --format washu --replicates intersect --scoreThreshold 5 --referenceScore 5 \
    --baits design.baitmap --names A B --tss genes.gtf.gz \
    --tracksA H3K27ac=A_H3K27ac.bw H3K4me1=A_H3K4me1.bw H3K4me3=A_H3K4me3.bw H3K27me3=A_H3K27me3.bw \
    --tracksB H3K27ac=B_H3K27ac.bw H3K4me1=B_H3K4me1.bw H3K4me3=B_H3K4me3.bw H3K27me3=B_H3K27me3.bw \
    --markSign H3K27me3=-1 \
    --markRoles H3K27ac=active H3K4me3=active H3K4me1=enhancer H3K27me3=repressive \
    --expression deseq2_B_vs_A.tsv --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes -o differential/A_vs_B
```

## 5. Integrated tables

```bash
egaExplain -p differential/A_vs_B -o explain/A_vs_B
```

## 6. Contacts across domain boundaries (optional)

With domains called on Hi-C of the same cells, for example by hicFindTADs.

```bash
egaTadCrossing -p differential/A_vs_B --tads B_domains --tadMinStrength 0.5 -o tads/A_vs_B
```

## 7. Reports

The main report, gene reports (HTML and PDF) and the interactive browser.

```bash
egaReport -p differential/A_vs_B --annotation genes.gtf.gz -o reports
egaReport -p differential/A_vs_B --gene GENE1 GENE2 --annotation genes.gtf.gz -o reports
egaExplorer -p differential/A_vs_B --annotation genes.gtf.gz -o explorer/
```
