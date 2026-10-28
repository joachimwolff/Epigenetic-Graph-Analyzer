# Input formats

## Called interactions

`egaCreateGraph` and `egaDifferential` read called interactions with `--format`. Only interactions on one
chromosome are used. Every reader returns the same table (bait, other end, distance, count, score).

| `--format` | Source | Layout |
|---|---|---|
| `ibed` | CHiCAGO `exportResults(..., outputFormat="interBed")` | header line; `bait_chr bait_start bait_end bait_name otherEnd_chr otherEnd_start otherEnd_end otherEnd_name N_reads score` |
| `washu` | CHiCAGO WashU text or WashU browser track | `chr1,831895,848168 TAB chr1,850619,874081 TAB 14.16` (bait first, score), or the long-range track `chr start end chr:start-end,score` |
| `bedpe` | any caller | `bait_chrom bait_start bait_end other_chrom other_start other_end score [count]` |
| `mifsud` | Table S5 of Mifsud et al. 2015 | as published |
| `homer` | HOMER `analyzeHiC` interaction tables | windows are snapped to a fixed grid; give the baits with `--baits` |
| `cool` | `.cool` / `.mcool::resolutions/<bin>` | every pixel touching a baited bin becomes a pair; `--baits` required, `--coolField` selects the pixel column |

**Scores and thresholds.** `--scoreThreshold` keeps pairs at or above a score (CHiCAGO: 5). `--referenceScore` is
the score that maps to edge strength 1. For the threshold flags of gained and lost interactions, the input of each
condition should contain the pairs below the calling threshold as well (for example the full CHiCAGO output); with
called pairs only, the flags read *not testable*.

**Replicates.** Give one file per replicate. `--replicates all` keeps pairs present in every replicate and reaching
the threshold in at least one, `intersect` pairs reaching it in every replicate, `union` pairs reaching it in any.

**Baits.** `--baits` takes the CHiCAGO `.baitmap` or a BED of baited fragments. Chromosome names with or without
`chr` are matched.

## Gene annotation

`--tss` takes a GTF (GENCODE, Ensembl, RefSeq; gzipped or not) or a BED of transcription start sites with a name
column. `--geneTypes` restricts the genes (default protein coding, lincRNA, antisense; without values, all types,
needed for RefSeq GTFs without a gene type).

## Chromatin marks

* **Signal tracks** (`--tracksA`, `--tracksB`, `egaAnnotateGraph --tracks`): `mark=file.bw`, with replicates as
  `mark=rep1.bw,rep2.bw` (averaged). The signal of a fragment is the exact mean of the track over the fragment. Any
  track works: ChIP-seq fold change over input, CUT&Tag, ATAC-seq, DNase-seq, DNA methylation.
* **States**: fold enrichment over input (`--stateScale fold`, the default; for example the Roadmap
  `.fc.signal.bigwig`), coverage without input (`--stateScale genomeMean`), separate state tracks
  (`--stateTracksA/B`) or peak calls (`--statePeaksA/B`: MACS2/MACS3 narrowPeak or broadPeak, SEACR, any BED).
* **Roles**: `--markRoles mark=role` or a TSV with the columns `mark` and `role`.

## Expression

A table with one row per gene:

* a log2 fold change of B against A (`--log2FCColumn`, for example DESeq2 `log2FoldChange`) and optionally an
  adjusted *P* (`--padjColumn`, with `--alpha`), or
* two expression columns of A and B (`--expressionColumns`, for example TPM), with `--minExpression`.

`--geneKey genes` matches the table by gene name, `gene_ids` by gene identifier (versions are ignored).

## Multi-way data

`egaBuildMultiway` writes a sparse matrix of loci by groups (`.npz`):

| `--type` | Input |
|---|---|
| `sprite` | SPRITE clusters, one line per cluster: `cluster_id TAB chrom:pos TAB chrom:pos ...`. With `--fragments` (a BED of the restriction fragments of the capture design) reads are assigned to fragments, else to bins of `--binSize` |
| `sprite` | Pore-C and HiPore-C concatemers in the same layout; `analysis/porec_to_clusters.py` converts per-fragment tables or pairs files with read ids (optionally lifted with a UCSC chain) |
| `gam` | GAM segregation tables with a header (`chrom start stop` and one 0/1 column per tube), one per replicate, as distributed by the 4D Nucleome Data Portal |

`--chromSizes` gives the chromosome sizes of the assembly of the multi-way data. If the graph and the multi-way data
use different assemblies, `egaValidateMultiway --chain` lifts the graph.

## Domains and boundaries

`--tads` (`egaDifferential`, `egaExplain`, `egaReport`, `egaExplorer`, `egaTadCrossing`) reads hicFindTADs output
(prefix or `_domains.bed`, `_boundaries.bed`, `_score.bedgraph`), Juicer Arrowhead or 4DN domain lists, cooltools
insulation tables or a BED of boundaries or domains; `--tadFormat` overrides the detection, `--tadChain` lifts them.
