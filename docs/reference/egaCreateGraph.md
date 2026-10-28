# egaCreateGraph

Build the promoter capture Hi-C interaction graph, annotate baited promoters with genes and partition the graph into communities.

```text
egaCreateGraph [-h] [--verbose] [--version] --interactions
                             INTERACTIONS [INTERACTIONS ...]
                             [--format {ibed,washu,bedpe,mifsud,homer,cool}]
                             [--baits BAITS] [--coolField COOLFIELD]
                             [--maxFdr MAXFDR]
                             [--scoreThreshold SCORETHRESHOLD]
                             [--minCount MINCOUNT]
                             [--replicates {all,intersect,union}]
                             [--maxDistance MAXDISTANCE] [--tss TSS]
                             [--geneTypes [GENETYPES ...]]
                             [--weight {score,count,unit}]
                             [--communityMethod {leiden,louvain}]
                             [--resolution RESOLUTION] [--seed SEED]
                             --outputPrefix OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--interactions`, `-i` `INTERACTIONS ...` (required) |  | Called interactions: CHiCAGO ibed or WashU text, BEDPE (bait first, score), the Mifsud et al. table, or cool matrices read out at the baited bins |
| `--format`, `-f` `{ibed, washu, bedpe, mifsud, homer, cool}` | `ibed` |  |
| `--baits` `BAITS` |  | Baited fragments of the capture design (BED or CHiCAGO .baitmap); required for --format cool, for other formats it sets the bait status and orientation of each pair |
| `--coolField` `COOLFIELD` | `count` | Pixel column of the cool matrices used as interaction score (e.g. count or a score) |
| `--maxFdr` `MAXFDR` |  | Keep pairs with an FDR at most this value (inputs with an FDR column, e.g. HOMER) |
| `--scoreThreshold` `SCORETHRESHOLD` |  | Keep pairs with score >= threshold (e.g. CHiCAGO 5, or a minimum pixel value) |
| `--minCount` `MINCOUNT` |  | Minimum read count of an edge (if available) |
| `--replicates` `{all, intersect, union}` |  | Treat the input files as replicates: 'all' keeps pairs present in every replicate and reaching --scoreThreshold in at least one (others may be below it), 'intersect' pairs reaching it in every replicate, 'union' in any; default: pool all files |
| `--maxDistance` `MAXDISTANCE` |  | Maximum bait-other distance of an edge (bp) |
| `--tss` `TSS` |  | Gene annotation: GTF(.gz) or BED of TSS/promoters |
| `--geneTypes` `GENETYPES ...` | `['protein_coding', 'lincRNA', 'antisense']` | Gene types taken from a GTF |
| `--weight` `{score, count, unit}` | `score` | Edge weight |
| `--communityMethod` `{leiden, louvain}` | `leiden` |  |
| `--resolution` `RESOLUTION` | `1.0` | Community resolution parameter |
| `--seed` `SEED` | `42` |  |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>.graphml, _nodes.tsv, _edges.tsv, _communities.tsv |
