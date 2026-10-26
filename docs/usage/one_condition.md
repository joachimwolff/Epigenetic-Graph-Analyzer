# One condition: graph, indirect contacts, annotation

For a single cell type, EGA builds the graph, lists the indirect contacts and attaches chromatin and expression to
the nodes. The commands below use the Md data of the tutorial; run them from `example/`.

## Build the graph

```bash
egaCreateGraph -i data/GSM6505208_*.txt.gz data/GSM6505209_*.txt.gz -f washu --replicates intersect \
    --tss data/mm10.ncbiRefSeq.chr6_hoxa.gtf.gz --geneTypes --maxDistance 2000000 \
    --weight score --communityMethod leiden --resolution 1.0 --seed 42 -o output/one/Md
```

Outputs: `Md.graphml` (the graph), `Md_nodes.tsv` (fragments with bait status, genes and community),
`Md_edges.tsv` (interactions with score and distance) and the communities. `--maxDistance` drops longer interactions;
`--communityMethod louvain` works without the `leiden` extra.

## Indirect contacts

```bash
egaIndirectContacts -g output/one/Md.graphml --maxDistance 2000000 --triplets -o output/one/Md
```

`Md_indirect.tsv` lists every pair at two steps with its class (`oe-oe`: two other ends, never measurable by promoter
capture; `bait-oe` and `bait-bait`: measurable but not called), the number of mediators, the Adamic-Adar index and a
path score; `Md_triplets.tsv` the two-step paths. `--maxMediatorDegree` excludes mediators with very many partners.

## Chromatin and expression on the nodes

```bash
egaAnnotateGraph -g output/one/Md.graphml \
    -t H3K27ac=data/E10.5_MD_H3K27ac_R1.bigwig H3K27me3=data/E10.5_MD_H3K27me3_R1.bigwig \
    -o output/one/Md_annotated
```

Every node receives the exact mean of each track over its fragment (`--flank` extends the fragments). An expression
table (`--expression`, `--expressionColumns`, `--geneKey`) attaches expression values to the baits.

## Plots

```bash
# neighbourhood of a gene: direct and indirect contacts on a genomic axis, nodes coloured by a mark
egaPlotNeighbourhood -g output/one/Md_annotated.graphml --gene Hoxa2 --steps 2 --colorBy H3K27ac \
    -o output/one/Hoxa2_neighbourhood.pdf

# a region as pyGenomeTracks input (arcs coloured by community) and as a promoter-centred graph
egaPlotRegion -g output/one/Md.graphml --region chr6:51000000-52500000 -o output/one/hoxa_region
```

`egaPlotRegion` writes `<prefix>_tracks.ini` for `pyGenomeTracks --tracks <prefix>_tracks.ini --region ...` and
`<prefix>_graph.pdf`. With `--multiway`, it adds the SPRITE clusters or GAM tubes that hold several loci of the
region (see [Multi-way data](multiway.md)).
