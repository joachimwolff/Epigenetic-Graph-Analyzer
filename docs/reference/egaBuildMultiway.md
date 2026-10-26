# egaBuildMultiway

Bin SPRITE clusters or GAM segregation tables into a sparse bin x group incidence matrix (.npz) used by egaValidateMultiway.

```text
egaBuildMultiway [-h] [--verbose] [--version] --type {sprite,gam}
                             --input INPUT [INPUT ...] --chromSizes CHROMSIZES
                             [--binSize BINSIZE] [--fragments FRAGMENTS]
                             [--minClusterSize MINCLUSTERSIZE]
                             [--maxClusterSize MAXCLUSTERSIZE] --output OUTPUT
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--type` `{sprite, gam}` (required) |  |  |
| `--input` `INPUT ...` (required) |  | SPRITE cluster file, or one GAM segregation table per replicate |
| `--chromSizes` `CHROMSIZES` (required) |  | Chromosome sizes of the data assembly |
| `--binSize` `BINSIZE` | `25000` | SPRITE bin size (GAM uses its windows) |
| `--fragments` `FRAGMENTS` |  | SPRITE: BED of intervals (e.g. restriction fragments of the capture Hi-C design); reads are assigned to the interval containing them instead of fixed bins |
| `--minClusterSize` `MINCLUSTERSIZE` | `2` | SPRITE: minimum reads per cluster |
| `--maxClusterSize` `MAXCLUSTERSIZE` | `1000` | SPRITE: maximum reads per cluster |
| `--output`, `-o` `OUTPUT` (required) |  | Output .npz |
