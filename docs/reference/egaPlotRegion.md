# egaPlotRegion

Write pyGenomeTracks input (interaction arcs coloured by community, indirect contacts) for a region and draw the same region as a graph with the promoters in the centre.

```text
egaPlotRegion [-h] [--verbose] [--version] --graph GRAPH
                             --region REGION [--genes GENES]
                             [--maxCommunities MAXCOMMUNITIES] [--noIndirect]
                             [--title TITLE] [--geneFontSize GENEFONTSIZE]
                             [--arcHeight ARCHEIGHT]
                             [--genesHeight GENESHEIGHT] [--multiway MULTIWAY]
                             [--chain CHAIN] [--minLoci MINLOCI]
                             [--maxEvents MAXEVENTS]
                             [--eventsHeight EVENTSHEIGHT]
                             [--graphSize WIDTH HEIGHT] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  |  |
| `--region` `REGION` (required) |  | chrom:start-end |
| `--genes` `GENES` |  | Gene annotation for the gene track (BED12 or GTF, uncompressed) |
| `--maxCommunities` `MAXCOMMUNITIES` | `8` |  |
| `--noIndirect` |  | Do not draw indirect contacts |
| `--title` `TITLE` |  |  |
| `--geneFontSize` `GENEFONTSIZE` | `9` | Font size of gene names in the gene track |
| `--arcHeight` `ARCHEIGHT` | `7` | Height of the interaction track (cm) |
| `--genesHeight` `GENESHEIGHT` | `4` | Height of the gene track (cm) |
| `--multiway` `MULTIWAY` |  | SPRITE / GAM .npz on the assembly of the graph; adds a track with the clusters or tubes that contain several loci of the region |
| `--chain` `CHAIN` |  | Chain file from the graph assembly to the assembly of --multiway; the multi-way track of the pyGenomeTracks output is omitted in this case |
| `--minLoci` `MINLOCI` | `3` | Minimum region loci per cluster / tube |
| `--maxEvents` `MAXEVENTS` | `60` | Clusters / tubes drawn |
| `--eventsHeight` `EVENTSHEIGHT` | `6` | Height of the multi-way track (cm) |
| `--graphSize` `WIDTH HEIGHT` | `[7.0, 5.0]` | Size of the graph figure in inches |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>_tracks.ini, <prefix>_links_*.tsv and <prefix>_graph.pdf |
