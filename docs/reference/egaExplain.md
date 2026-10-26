# egaExplain

egaExplain: an integrated report of a differential analysis, linking interaction changes, chromatin changes and expression changes.

```text
egaExplain [-h] [--verbose] [--version] --prefix PREFIX
                             [--names NAMES NAMES]
                             [--repressive [REPRESSIVE ...]] [--fdr FDR]
                             [--pseudo PSEUDO] [--tads TADS]
                             [--tadFormat {auto,hicexplorer,arrowhead,cooltools,bed,bed-boundaries,bed-domains,ega}]
                             [--tadWindow TADWINDOW] [--tadChain TADCHAIN]
                             [--tadRegions TADREGIONS]
                             [--tadMinStrength TADMINSTRENGTH] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--prefix`, `-p` `PREFIX` (required) |  | Output prefix of egaDifferential |
| `--names` `NAMES` |  | Condition names A B (default: from the graph) |
| `--repressive` `REPRESSIVE ...` | `['H3K27me3', 'H3K9me3']` | Repressive marks, expected to change opposite to activity |
| `--fdr` `FDR` | `0.1` | Communities with a coordinated change at this FDR |
| `--pseudo` `PSEUDO` | `1.0` | Pseudo count of the log2 mark change |
| `--tads` `TADS` |  | TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, _tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools insulation table, BED of boundaries or domains (optional score column), or egaBoundaries output |
| `--tadFormat` `{auto, hicexplorer, arrowhead, cooltools, bed, bed-boundaries, bed-domains, ega}` | `auto` | Format of --tads |
| `--tadWindow` `TADWINDOW` |  | cooltools: insulation window (default: the largest in the table) |
| `--tadChain` `TADCHAIN` |  | Chain file from the assembly of --tads to that of the graph (UCSC chromosome names) |
| `--tadRegions` `TADREGIONS` |  | BED of named regions (chrom start end name) used to name the partner regions |
| `--tadMinStrength` `TADMINSTRENGTH` |  | Only boundaries at or above this strength percentile (0-1) count as crossed |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  |  |
