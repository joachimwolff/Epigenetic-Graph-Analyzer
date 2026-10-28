# egaValidateMultiway

Validate graph-derived direct and indirect contacts, two-hop triplets and communities with ligation-free multi-way contact data against distance-matched controls.

```text
egaValidateMultiway [-h] [--verbose] [--version] --graph GRAPH
                             --multiway MULTIWAY [--chain CHAIN]
                             [--maxDistance MAXDISTANCE]
                             [--maxMediatorDegree MAXMEDIATORDEGREE]
                             [--controls CONTROLS]
                             [--compartments COMPARTMENTS] [--domains DOMAINS]
                             [--coverageClasses COVERAGECLASSES]
                             [--confirmPairs] [--maxPairs MAXPAIRS]
                             [--maxTriplets MAXTRIPLETS]
                             [--communityMinMembers COMMUNITYMINMEMBERS]
                             [--skipPairs] [--skipTriplets]
                             [--skipCommunities] [--seed SEED]
                             [--threads THREADS] --outputPrefix OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--graph`, `-g` `GRAPH` (required) |  | GraphML of egaCreateGraph |
| `--multiway`, `-m` `MULTIWAY` (required) |  | .npz of egaBuildMultiway |
| `--chain` `CHAIN` |  | UCSC chain file to lift graph coordinates to the multi-way assembly |
| `--maxDistance` `MAXDISTANCE` | `2000000` |  |
| `--maxMediatorDegree` `MAXMEDIATORDEGREE` |  |  |
| `--controls` `CONTROLS` | `10` | Matched controls per tested pair |
| `--compartments` `COMPARTMENTS` |  | A/B compartments in the multi-way assembly (BED with names A*/B* or eigenvector bedGraph); controls match the compartments of both ends |
| `--domains` `DOMAINS` |  | Topological domains in the multi-way assembly (BED); controls match whether both ends lie in the same domain |
| `--coverageClasses` `COVERAGECLASSES` |  | Also match pair controls on the coverage class (quantiles over graph loci) of both ends; 0 matches on offset and locus pools only |
| `--confirmPairs` |  | Test every direct and indirect pair individually against background pairs of the same distance and coverage (writes <prefix>_confirm_*.tsv) |
| `--maxPairs` `MAXPAIRS` | `500000` | Random subset size of direct and of indirect pairs (all pairs if fewer) |
| `--maxTriplets` `MAXTRIPLETS` | `500000` | Random subset size of triplets |
| `--communityMinMembers` `COMMUNITYMINMEMBERS` | `3` | A group supports a community if it holds at least this many of its bins |
| `--skipPairs` |  |  |
| `--skipTriplets` |  |  |
| `--skipCommunities` |  |  |
| `--seed` `SEED` |  |  |
| `--threads` `THREADS` |  | Worker processes for the pair, triplet and community queries (default: available CPUs, at most 24); the results do not depend on the number |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  |  |
