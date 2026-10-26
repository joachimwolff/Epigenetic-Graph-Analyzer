# egaBoundaries

[Experimental] The evaluation of this module gave weak or null results; its output is exploratory. Per-bin boundary evidence: capture insulation dip (graph edges spanning each bin against a distance-matched mirror expectation), cohesin-bound CTCF with motif orientation (reverse sites on the left, forward on the right), convergent CTCF loops among the edges; boundaries by a rule with fixed thresholds (no fitted weights).

```text
egaBoundaries [-h] [--verbose] [--version] --nodes NODES
                             --edges EDGES --chromSizes CHROMSIZES
                             [--chromosomes [CHROMOSOMES ...]]
                             [--binSize BINSIZE] [--minExpected MINEXPECTED]
                             [--dipProminence DIPPROMINENCE] [--dipMax DIPMAX]
                             [--ctcfMin CTCFMIN] [--window WINDOW]
                             [--rule {dip,ctcf,ctcf_pattern,loop_anchor,dip_or_ctcf,dip_and_ctcf,dip_or_ctcf_pattern,dip_and_ctcf_pattern}]
                             [--ctcf [CTCF ...]] [--cohesin [COHESIN ...]]
                             [--fasta FASTA] [--motifs MOTIFS]
                             [--motifId MOTIFID] [--motifPvalue MOTIFPVALUE]
                             [--minEnrichment MINENRICHMENT]
                             [--reference REFERENCE] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--nodes` `NODES` (required) |  | Node table of a graph (egaCreateGraph *_nodes.tsv) |
| `--edges` `EDGES` (required) |  | Edge table of the same graph (*_edges.tsv) |
| `--chromSizes` `CHROMSIZES` (required) |  | Chromosome sizes (UCSC format) |
| `--chromosomes` `CHROMOSOMES ...` |  | Restrict to these chromosomes |
| `--binSize` `BINSIZE` | `25000` |  |
| `--minExpected` `MINEXPECTED` | `5.0` | Bins with fewer expected spanning edges are not covered (no insulation) |
| `--dipProminence` `DIPPROMINENCE` | `0.3` | Dip: local minimum of the capture insulation with at least this prominence (log2) |
| `--dipMax` `DIPMAX` |  | Dip: capture insulation (log2 o/e) at most this |
| `--ctcfMin` `CTCFMIN` | `10.0` | CTCF evidence: cohesin-bound sites with CTCF enrichment >= this |
| `--window` `WINDOW` | `1` | CTCF and loop evidence within +- this many bins |
| `--rule` `{dip, ctcf, ctcf_pattern, loop_anchor, dip_or_ctcf, dip_and_ctcf, dip_or_ctcf_pattern, dip_and_ctcf_pattern}` |  | Rule for the boundary calls (default: dip_and_ctcf_pattern with CTCF, dip without); the calls of every rule are columns of <prefix>_bins.tsv |
| `--ctcf` `CTCF ...` | `[]` | CTCF ChIP bigWig(s), averaged |
| `--cohesin` `COHESIN ...` | `[]` | Cohesin (RAD21 / SMC1A / SMC3) bigWig(s) |
| `--fasta` `FASTA` |  | Genome FASTA (plain or gzipped) for the motif orientation |
| `--motifs` `MOTIFS` |  | JASPAR-format motif file |
| `--motifId` `MOTIFID` | `MA0139` | CTCF motif ID (prefix match) |
| `--motifPvalue` `MOTIFPVALUE` | `0.0001` |  |
| `--minEnrichment` `MINENRICHMENT` | `10.0` | CTCF summits: minimum signal over the genome-wide mean |
| `--reference` `REFERENCE` |  | BED of reference (Hi-C) boundaries; boundary = interval midpoint |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Writes <prefix>_bins.tsv, <prefix>_capture_insulation.bedGraph, <prefix>_boundaries.bed, <prefix>_domains.bed and, with CTCF, <prefix>_ctcf_sites.bed and <prefix>_loops.tsv; with --reference <prefix>_evaluation.tsv |
