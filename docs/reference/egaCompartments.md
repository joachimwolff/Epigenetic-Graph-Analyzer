# egaCompartments

[Experimental] The evaluation of this module gave weak or null results; its output is exploratory. A/B compartment estimates per genomic bin from chromatin tracks (first principal component oriented by the active marks; active minus repressive role score) and community / edge statistics of a graph against TADs.

```text
egaCompartments [-h] [--verbose] [--version] --chromSizes
                             CHROMSIZES [--binSize BINSIZE]
                             [--tracks [TRACKS ...]]
                             [--markRoles [MARKROLES ...]] [--smooth SMOOTH]
                             [--genes GENES] [--reference REFERENCE]
                             [--nodes NODES] [--edges EDGES]
                             [--domains DOMAINS] [--shifts SHIFTS]
                             [--window WINDOW] [--seed SEED] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--chromSizes` `CHROMSIZES` (required) |  | Chromosome sizes (UCSC format) |
| `--binSize` `BINSIZE` | `100000` | Bin size in bp for the chromatin estimates |
| `--tracks`, `-t` `TRACKS ...` | `[]` | name=file.bw[,replicate.bw]; replicates are averaged |
| `--markRoles` `MARKROLES ...` | `[]` | name=role; roles active, enhancer, transcription (positive), repressive, heterochromatin (negative), other (PC1 only) |
| `--smooth` `SMOOTH` |  | Running mean over 2 * smooth + 1 bins |
| `--genes` `GENES` |  | BED of genes (or TSS); the gene density orients the PC1 when no active role is declared |
| `--reference` `REFERENCE` |  | Hi-C eigenvector (bedGraph, positive = A) to evaluate against |
| `--nodes` `NODES` |  | Node table of a graph (egaCreateGraph *_nodes.tsv) for the domain test |
| `--edges` `EDGES` |  | Edge table of the same graph (*_edges.tsv) |
| `--domains` `DOMAINS` |  | TADs (BED or Rao et al. domain list) for the domain test |
| `--shifts` `SHIFTS` | `1000` | Circular shifts of the domain test |
| `--window` `WINDOW` | `50000` | Distance in bp to count a boundary as matched |
| `--seed` `SEED` |  |  |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>_bins.tsv, <prefix>_compartments.bedGraph, <prefix>_evaluation.tsv, <prefix>_loadings.tsv and <prefix>_domains.tsv as applicable |
