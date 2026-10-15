# Epigenetic Graph Analyzer (EGA)

EGA represents promoter capture Hi-C (pcHi-C) interactions as a graph whose nodes are restriction
fragments and whose edges are significant interactions. Baited promoters and other ends share one node
space, so promoter-promoter interactions connect the neighbourhoods of different genes. On this graph EGA

* infers **indirect contacts**: pairs of fragments without a significant interaction that share a
  neighbour. Pairs of two non-baited fragments (`oe-oe`) cannot be measured by promoter capture at all;
* partitions the graph into **communities** of co-interacting promoters and regulatory elements;
* **validates** direct contacts, indirect contacts, two-hop triplets and communities with ligation-free
  multi-way contact data (SPRITE, GAM) against distance-matched controls;
* in **differential mode** compares two conditions: interaction changes, contact-weighted chromatin
  changes at the promoter, at direct and at indirect (two-hop) neighbours, and their relation to gene
  expression changes.

## Installation

```bash
pip install .            # core
pip install .[leiden]    # Leiden community detection (python-igraph, leidenalg), recommended
```

Python >= 3.9 with numpy, scipy, pandas, networkx, pyBigWig, matplotlib and pyliftover.

## Command line tools

| Tool | Purpose |
|---|---|
| `egaCallInteractions` | Significant interactions from raw pcHi-C counts (distance decay, bait and other-end factors, negative binomial test) |
| `egaCreateGraph` | Graph, gene annotation of baits, Leiden/Louvain communities (GraphML + TSV tables) |
| `egaIndirectContacts` | Indirect contacts with pair classes and scores, direct contacts, two-hop triplets |
| `egaAnnotateGraph` | Mean bigWig signal per node, expression values per bait |
| `egaBuildMultiway` | Bin SPRITE clusters or GAM segregation tables into a sparse incidence matrix |
| `egaValidateMultiway` | Validation of direct/indirect contacts, triplets and communities against multi-way data |
| `egaDifferential` | Differential interactions and promoter/direct/indirect chromatin changes versus expression changes |
| `egaPlotNeighbourhood` | Direct and indirect neighbourhood of a gene on the genomic axis |

Every tool prints its options with `--help`.

## Example

```bash
# 1. significant interactions (Jung et al. 2019 raw count tables)
egaCallInteractions -f jung -i H1.po.all.txt.bz2 H1.pp.all.txt.bz2 \
    --minDistance 15000 --pvalueThreshold 0.01 --minCount 6 -o H1.tsv.gz --summary H1.json

# 2. graph, communities and indirect contacts
egaCreateGraph -i H1.tsv.gz --tss gencode.v19.annotation.gtf.gz -o graphs/H1
egaIndirectContacts -g graphs/H1.graphml -o graphs/H1 --triplets

# 3. validation with GAM (hg38) of a graph built on hg19
egaBuildMultiway --type gam --input rep1_30kb.csv rep2_30kb.csv --chromSizes hg38.chrom.sizes -o gam.npz
egaValidateMultiway -g graphs/H1.graphml -m gam.npz --chain hg19ToHg38.over.chain.gz -o validation/H1_gam

# 4. differential mode
egaDifferential --interactionsA H1.tsv.gz --interactionsB NPC.tsv.gz --names H1 NPC \
    --tss gencode.v19.annotation.gtf.gz \
    --tracksA H3K27ac=E003-H3K27ac.bw H3K27me3=E003-H3K27me3.bw \
    --tracksB H3K27ac=E007-H3K27ac.bw H3K27me3=E007-H3K27me3.bw --markSign H3K27me3=-1 \
    --expression 57epigenomes.RPKM.pc.gz --expressionColumns E003 E007 -o differential/H1_vs_NPC
```

The analyses of the manuscript are scripted in `analysis/` (`scripts/download_data.sh` fetches all
public input data).

## Input formats

* Interactions: Jung et al. 2019 raw count tables (`jung`), Mifsud et al. 2015 tables (`mifsud`),
  CHiCAGO `ibed`, or a BEDPE-like table (`bedpe`: bait chrom/start/end, other chrom/start/end, score[, count]).
* Genes: GTF/GFF (transcript TSS) or BED of TSS.
* SPRITE: cluster files (`barcode<TAB>chrom:position<TAB>...`). GAM: segregation tables
  (`chrom start stop tube1 tube2 ...`, one table per replicate).

## Testing

```bash
pip install .[test,leiden]
pytest
```

The tests run every module and every command line tool on synthetic data with planted interactions,
hubs and expression effects.

## License

MIT
