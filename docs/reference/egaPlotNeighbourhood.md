# egaPlotNeighbourhood

Plot the neighbourhood of a gene: nodes on a genomic axis, arcs for direct contacts, dashed arcs for indirect contacts of the gene promoter.

```text
egaPlotNeighbourhood [-h] [--verbose] [--version] --graph GRAPH --gene
                             GENE [--steps {1,2}] [--colorBy COLORBY] --output
                             OUTPUT
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  |  |
| `--gene` `GENE` (required) |  | Gene name as annotated on the bait node |
| `--steps` `{1, 2}` | `2` |  |
| `--colorBy` `COLORBY` |  | Numeric node attribute used for node colour (e.g. H3K27ac) |
| `--output`, `-o` `OUTPUT` (required) |  | Image file (pdf, png, svg) |
