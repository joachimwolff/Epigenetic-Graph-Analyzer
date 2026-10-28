# egaEqtl

Test eQTL variant-gene pairs for graph paths of 1, 2 and 3 steps against distance-matched non-eQTL pairs.

```text
egaEqtl [-h] [--verbose] [--version] --leads LEADS
                             --significant SIGNIFICANT --nodes NODES --edges
                             EDGES [--rmap RMAP] [--liftRmap LIFTRMAP]
                             [--binSize BINSIZE] [--chromSizes CHROMSIZES]
                             [--baits BAITS] [--baitGtf BAITGTF]
                             [--baitGeneTypes BAITGENETYPES] --gtf GTF
                             [--geneTypes GENETYPES]
                             [--maxDistance MAXDISTANCE] [--controls CONTROLS]
                             [--window WINDOW] [--lengthWindow LENGTHWINDOW]
                             [--noBaitMatching] [--bootstrap BOOTSTRAP]
                             [--modelBootstrap MODELBOOTSTRAP] [--mirror]
                             [--label LABEL] [--seed SEED] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--leads` `LEADS` (required) |  | Tested eQTL pairs (chrom, pos, gene_id) |
| `--significant` `SIGNIFICANT` (required) |  | All significant pairs of the data set (chrom, pos, gene_id); also the source of control fragments |
| `--nodes` `NODES` (required) |  | <graph>_nodes.tsv of egaCreateGraph |
| `--edges` `EDGES` (required) |  | <graph>_edges.tsv of egaCreateGraph |
| `--rmap` `RMAP` |  | Restriction fragments (rmap or BED, 1-based closed) |
| `--liftRmap` `LIFTRMAP` |  | Chain file to lift the --rmap fragments to the graph assembly |
| `--binSize` `BINSIZE` |  | Fixed bins as fragments (virtual capture) instead of --rmap |
| `--chromSizes` `CHROMSIZES` |  | Chromosome sizes for --binSize |
| `--baits` `BAITS` |  | Baited fragments (baitmap or BED) in the graph assembly |
| `--baitGtf` `BAITGTF` |  | Baits are the fragments holding a TSS of --baitGeneTypes in this GTF |
| `--baitGeneTypes` `BAITGENETYPES` | `protein_coding` |  |
| `--gtf` `GTF` (required) |  | Gene annotation (GTF) of the graph assembly for the promoters |
| `--geneTypes` `GENETYPES` | `protein_coding,lincRNA,antisense` |  |
| `--maxDistance` `MAXDISTANCE` | `1000000` |  |
| `--controls` `CONTROLS` | `5` | Controls per eQTL pair and control type |
| `--window` `WINDOW` | `0.1` | Distance matching window, log10 units |
| `--lengthWindow` `LENGTHWINDOW` | `0.15` | Fragment length matching window of control fragments, log10 units (negative: off) |
| `--noBaitMatching` |  | Do not match control fragments on the bait status of the eQTL fragment |
| `--bootstrap` `BOOTSTRAP` | `1000` | Block-bootstrap resamples for fractions, odds |
| `--modelBootstrap` `MODELBOOTSTRAP` | `200` | Block-bootstrap resamples for the models |
| `--mirror` |  | Also evaluate the mirrored control graph |
| `--label` `LABEL` | `` | Label written into every output table |
| `--seed` `SEED` |  |  |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  |  |
