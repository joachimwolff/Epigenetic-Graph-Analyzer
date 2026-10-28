# egaPlotDifferential

Plot the differential interaction graph of egaDifferential for a region or the neighbourhood of a gene: gained, lost and shared interactions as arcs and as a promoter-centred graph with chromatin and expression changes.

```text
egaPlotDifferential [-h] [--verbose] [--version] --graph GRAPH
                             (--region REGION | --gene GENE) [--steps STEPS]
                             [--names NAMES NAMES] [--mark MARK [MARK ...]]
                             [--candidates CANDIDATES] [--sideBySide]
                             [--separatePanels] [--maxLog2 MAXLOG2]
                             [--maxDistance MAXDISTANCE] [--title TITLE]
                             --outFileName OUTFILENAME [--dpi DPI]
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  | <prefix>_union.graphml of egaDifferential |
| `--region` `REGION` |  | chrom:start-end |
| `--gene` `GENE` |  | Gene name; draws its neighbourhood up to --steps |
| `--steps` `STEPS` | `2` | Neighbourhood size around --gene (graph steps); to show the nodes that explain a candidate, at least its explaining level (explaining_hops of --candidates; a warning is logged otherwise) |
| `--names` `NAMES` |  | Condition names A B (default: from the graph) |
| `--mark` `MARK ...` |  | Chromatin marks whose log2 change colours the other ends, one panel each (default: all marks of the graph; the combined figure uses the first) |
| `--candidates` `CANDIDATES` |  | <prefix>_indirect_only_genes.tsv of egaDifferential; outlines the candidates and the nodes at their explaining level (2, 3, ... hops) that explain them |
| `--sideBySide` |  | Also draw condition A and condition B |
| `--separatePanels` |  | Write every panel as its own file (<out>_arcs, _A, _B, _edges, _nodes with the extension of --outFileName), for figures assembled in a manuscript |
| `--maxLog2` `MAXLOG2` | `2.0` | Colour scale limit of log2 changes |
| `--maxDistance` `MAXDISTANCE` | `2000000` | Interactions longer than this (bp) are not drawn (default: 2 Mb, as in the graphs) |
| `--title` `TITLE` |  |  |
| `--outFileName`, `-o` `OUTFILENAME` (required) |  | Figure file; the format follows the extension: pdf, svg, eps, ps (vector) or png, jpg, tif (raster) |
| `--dpi` `DPI` | `300` | Resolution of raster formats (default: 300) |
