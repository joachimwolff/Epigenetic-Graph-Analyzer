# egaAnnotateGraph

Annotate graph nodes with the mean bigWig signal over each fragment and baits with gene expression values.

```text
egaAnnotateGraph [-h] [--verbose] [--version] --graph GRAPH
                             [--tracks [TRACKS ...]] [--flank FLANK]
                             [--expression EXPRESSION]
                             [--expressionColumns [EXPRESSIONCOLUMNS ...]]
                             [--geneKey {gene_ids,genes}] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  | GraphML input |
| `--tracks`, `-t` `TRACKS ...` | `[]` | name=file.bw[,replicate.bw ...]; replicates are averaged |
| `--flank` `FLANK` |  | Extend fragments by this many bp on both sides |
| `--expression` `EXPRESSION` |  | Expression table (gene identifier column + value columns) |
| `--expressionColumns` `EXPRESSIONCOLUMNS ...` |  | Value columns of the expression table to attach (default: all numeric) |
| `--geneKey` `{gene_ids, genes}` | `gene_ids` | Node attribute used to match the expression table |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>.graphml and _nodes.tsv |
