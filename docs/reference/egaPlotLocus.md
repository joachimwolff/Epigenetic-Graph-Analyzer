# egaPlotLocus

Plot a region of the interaction graph (arcs of direct and indirect contacts) above the multi-way contact map of the same region with the graph contacts marked on it.

```text
egaPlotLocus [-h] [--verbose] [--version] --graph GRAPH
                             --multiway MULTIWAY [--region REGION]
                             [--chain CHAIN] [--metric {auto,cooc_oe,npmi}]
                             [--highlightCommunity HIGHLIGHTCOMMUNITY]
                             [--style {events,communities,map}]
                             [--panels {all,evidence}] [--minLoci MINLOCI]
                             [--maxRows MAXROWS] [--shifts SHIFTS]
                             [--title TITLE] --output OUTPUT
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  |  |
| `--multiway`, `-m` `MULTIWAY` (required) |  | .npz of egaBuildMultiway |
| `--region` `REGION` |  | chrom:start-end in the assembly of the graph (map style) |
| `--chain` `CHAIN` |  | Chain file from the graph assembly to the multi-way assembly |
| `--metric` `{auto, cooc_oe, npmi}` | `auto` | Map values: log2 observed/expected co-occurrence or NPMI ('auto': SPRITE -> cooc_oe, GAM -> npmi) |
| `--highlightCommunity` `HIGHLIGHTCOMMUNITY` |  | Community id drawn in colour |
| `--style` `{events, communities, map}` | `events` | 'events': the SPRITE clusters / GAM tubes that contain several loci of the community, 'map': pairwise co-occurrence map of the region |
| `--panels` `{all, evidence}` | `all` | communities: all panels or only the multi-way events and their purity |
| `--minLoci` `MINLOCI` | `3` | events: minimum community loci per cluster / tube |
| `--maxRows` `MAXROWS` | `60` | events: clusters / tubes drawn |
| `--shifts` `SHIFTS` | `50` | events: shifted control sets |
| `--title` `TITLE` |  |  |
| `--output`, `-o` `OUTPUT` (required) |  |  |
