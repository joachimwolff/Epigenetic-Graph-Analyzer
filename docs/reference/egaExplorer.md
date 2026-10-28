# egaExplorer

Self-contained offline HTML explorer of a differential run of egaDifferential: gene search with the text of egaReport, interactive arc view and neighbourhood graph, and browse views of the differential genes and communities.

```text
egaExplorer [-h] [--verbose] [--version] --prefix PREFIX
                             [--explainPrefix EXPLAINPREFIX]
                             [--names NAMES NAMES] [--labels LABELS LABELS]
                             [--annotation ANNOTATION] [--rechip RECHIP]
                             [--differentialOnly] [--maxNodes MAXNODES]
                             [--maxLog2 MAXLOG2] [--rnaThreshold RNATHRESHOLD]
                             [--chromatinThreshold CHROMATINTHRESHOLD]
                             [--log [LOG ...]] [--name NAME] [--tads TADS]
                             [--tadFormat {auto,hicexplorer,arrowhead,cooltools,bed,bed-boundaries,bed-domains,ega}]
                             [--tadWindow TADWINDOW] [--tadChain TADCHAIN]
                             [--tadRegions TADREGIONS]
                             [--tadMinStrength TADMINSTRENGTH] --output OUTPUT
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--prefix`, `-p` `PREFIX` (required) |  | Output prefix of egaDifferential |
| `--explainPrefix` `EXPLAINPREFIX` |  | Output prefix of egaExplain (default: <prefix dir>/explain/<prefix name> if present) |
| `--names` `NAMES` |  | Condition names A B (default: from the run) |
| `--labels` `LABELS` |  | Readable names of the conditions A B |
| `--annotation` `ANNOTATION` |  | Gene annotation (GTF, GTF.gz or BED with strand) for upstream / downstream |
| `--rechip` `RECHIP` |  | Re-ChIP signal per bait (default as egaReport) |
| `--differentialOnly` |  | Search index of the differential genes only (default: every baited gene of the graph) |
| `--maxNodes` `MAXNODES` | `80` | Default largest neighbourhood drawn (changeable in the page) |
| `--maxLog2` `MAXLOG2` | `2.0` | Colour scale limit of the log2 changes of nodes |
| `--rnaThreshold` `RNATHRESHOLD` |  | As egaReport |
| `--chromatinThreshold` `CHROMATINTHRESHOLD` |  | As egaReport |
| `--log` `LOG ...` |  | Run logs (default: <prefix>.log) |
| `--name` `NAME` |  | Name of the comparison (default: the prefix name) |
| `--tads` `TADS` |  | TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, _tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools insulation table, BED of boundaries or domains (optional score column), or egaBoundaries output |
| `--tadFormat` `{auto, hicexplorer, arrowhead, cooltools, bed, bed-boundaries, bed-domains, ega}` | `auto` | Format of --tads |
| `--tadWindow` `TADWINDOW` |  | cooltools: insulation window (default: the largest in the table) |
| `--tadChain` `TADCHAIN` |  | Chain file from the assembly of --tads to that of the graph (UCSC chromosome names) |
| `--tadRegions` `TADREGIONS` |  | BED of named regions (chrom start end name) used to name the partner regions |
| `--tadMinStrength` `TADMINSTRENGTH` |  | Only boundaries at or above this strength percentile (0-1) count as crossed |
| `--output`, `-o` `OUTPUT` (required) |  | Output HTML file, or a directory (then <dir>/<name>_explorer.html) |
