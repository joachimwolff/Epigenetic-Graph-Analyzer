# egaTadCrossing

Per edge of a differential run (egaDifferential): the TAD boundaries it crosses, the strongest of them and same- vs cross-domain, per condition and by status (gained / lost / shared); per promoter the crossing contacts, the boundary they cross, their partner region and a sentence; genome-wide ranking of promoters (Fisher: gained / lost among crossing vs non-crossing edges, BH) and of boundaries (gained / lost crossing edges vs the expectation from edges of the same length).

```text
egaTadCrossing [-h] [--verbose] [--version] --prefix PREFIX
                             --tads TADS
                             [--tadFormat {auto,hicexplorer,arrowhead,cooltools,bed,bed-boundaries,bed-domains,ega}]
                             [--tadWindow TADWINDOW] [--tadChain TADCHAIN]
                             [--tadRegions TADREGIONS]
                             [--tadMinStrength TADMINSTRENGTH]
                             [--names NAMES NAMES] [--gene [GENE ...]]
                             --outputPrefix OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--prefix`, `-p` `PREFIX` (required) |  | Output prefix of egaDifferential |
| `--tads` `TADS` (required) |  | TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, _tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools insulation table, BED of boundaries or domains (optional score column), or egaBoundaries output |
| `--tadFormat` `{auto, hicexplorer, arrowhead, cooltools, bed, bed-boundaries, bed-domains, ega}` | `auto` | Format of --tads |
| `--tadWindow` `TADWINDOW` |  | cooltools: insulation window (default: the largest in the table) |
| `--tadChain` `TADCHAIN` |  | Chain file from the assembly of --tads to that of the graph (UCSC chromosome names) |
| `--tadRegions` `TADREGIONS` |  | BED of named regions (chrom start end name) used to name the partner regions |
| `--tadMinStrength` `TADMINSTRENGTH` |  | Only boundaries at or above this strength percentile (0-1) count as crossed |
| `--names` `NAMES` |  | Condition names (default: from the run) |
| `--gene` `GENE ...` | `[]` | Print the sentence of these genes |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>_edges.tsv.gz, _summary.tsv, _promoters.tsv, _boundaries.tsv, _genes.tsv |
