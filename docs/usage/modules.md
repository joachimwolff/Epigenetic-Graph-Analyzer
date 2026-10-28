# Further modules

These modules work on the same graph and annotation. Each was evaluated against an independent reference; the key
number and the limitation are given with each.

## Contacts across domain boundaries

```bash
egaTadCrossing -p differential/A_vs_B --tads B_domains --tadMinStrength 0.5 --tadRegions named_regions.bed \
    --gene GENE1 -o tads/A_vs_B
```

For every edge, the boundaries it crosses; per promoter, the crossing gained and lost contacts and a sentence; a
genome-wide ranking of promoters (gained against lost among crossing edges) and of boundaries. On the mouse arches, it
ranked the published gain of Hoxa2 contacts into the HIRE1 region first among all promoters. Boundaries should come
from Hi-C of the same cells.

## Communities

`egaCreateGraph` and `egaDifferential` partition the graph into communities (Leiden or Louvain, `--resolution`).
`egaDifferential` matches the communities of both conditions and tests the coordinated expression change of each
community against random gene sets (`--communityPermutations`).

## Element classes

`egaDifferential` classifies the other ends as promoters (baited or not, from the transcription start sites of
`--annotation`, by default those of `--tss`), enhancers or other elements from the marks; `--elements` compares the classes with an external annotation such as ENCODE cCREs or
ChromHMM states.

## Compartments from marks (experimental)

```bash
egaCompartments --chromSizes hg19.chrom.sizes --binSize 100000 \
    -t H3K27ac=H3K27ac.bw H3K27me3=H3K27me3.bw --markRoles H3K27ac=active H3K27me3=repressive \
    --reference hic_eigenvector.bedGraph -o compartments/GM12878
```

A/B compartments from chromatin tracks. Against Hi-C compartments, the area under the curve was 0.85 to 0.91, about
as high as that of GC content.

## Boundaries from capture insulation and CTCF (experimental)

```bash
egaBoundaries --nodes graphs/B_nodes.tsv --edges graphs/B_edges.tsv --chromSizes mm10.chrom.sizes \
    --ctcf CTCF.bw --cohesin RAD21.bw --fasta mm10.fa.gz --motifs JASPAR2024_CORE.jaspar \
    --rule dip_or_ctcf_pattern --reference hic_boundaries.bed -o boundaries/B
```

Capture-insulation minima reached a precision of 0.64 to 0.69 at a recall of 0.24 to 0.33 against Hi-C boundaries.

## eQTL along graph paths

```bash
egaEqtl --leads leads.tsv --significant significant.tsv --nodes graphs/GM12878_nodes.tsv \
    --edges graphs/GM12878_edges.tsv --rmap HindIII.rmap --baits design.baitmap --gtf genes.gtf \
    --mirror -o eqtl/GM12878
```

eQTL variant-gene pairs at one, two and three steps against distance-matched pairs without an eQTL.

## Structural variants and copy number

`egaDifferential --structuralVariants` and `--copyNumber` add breakpoints and copy number changes of an external
caller (for example from whole-genome sequencing) as further levels. This option is not part of the evaluation of the
manuscript.
