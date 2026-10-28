# egaIndirectContacts

Infer indirect contacts (pairs sharing a neighbour without a direct significant interaction) and list the two-hop paths (triplets) of a graph.

```text
egaIndirectContacts [-h] [--verbose] [--version] --graph GRAPH
                             [--maxDistance MAXDISTANCE]
                             [--minDistance MINDISTANCE]
                             [--maxMediatorDegree MAXMEDIATORDEGREE]
                             [--triplets] --outputPrefix OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  | GraphML of egaCreateGraph |
| `--maxDistance` `MAXDISTANCE` | `2000000` | Maximum distance of a pair (bp) |
| `--minDistance` `MINDISTANCE` |  | Minimum distance of a pair (bp) |
| `--maxMediatorDegree` `MAXMEDIATORDEGREE` |  | Ignore mediators with a larger degree (limits hub explosion) |
| `--triplets` |  | Also write the open two-hop triplets |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>_indirect.tsv, _direct.tsv and optionally _triplets.tsv |
