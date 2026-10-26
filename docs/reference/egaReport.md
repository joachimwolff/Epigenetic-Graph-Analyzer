# egaReport

Self-contained HTML and PDF report of a differential run of egaDifferential (and egaExplain): summary, explanation cascade, major findings with figures, hubs and communities, methods; or reports of single genes (--gene) or of one community (--community).

```text
egaReport [-h] [--verbose] [--version] --prefix PREFIX
                             [--explainPrefix EXPLAINPREFIX]
                             [--names NAMES NAMES] [--gene GENE [GENE ...]]
                             [--community COMMUNITY] [--communitySide {A,B}]
                             [--labels LABELS LABELS]
                             [--annotation ANNOTATION] [--rechip RECHIP]
                             [--top TOP] [--maxCommunities MAXCOMMUNITIES]
                             [--maxNodes MAXNODES] [--overlay]
                             [--overlayTracksA [OVERLAYTRACKSA ...]]
                             [--overlayTracksB [OVERLAYTRACKSB ...]]
                             [--rnaThreshold RNATHRESHOLD]
                             [--chromatinThreshold CHROMATINTHRESHOLD]
                             [--maxLog2 MAXLOG2] [--log [LOG ...]]
                             [--name NAME]
                             [--pdf {auto,chromium,matplotlib,none}]
                             [--chromium CHROMIUM] [--tads TADS]
                             [--tadFormat {auto,hicexplorer,arrowhead,cooltools,bed,bed-boundaries,bed-domains,ega}]
                             [--tadWindow TADWINDOW] [--tadChain TADCHAIN]
                             [--tadRegions TADREGIONS]
                             [--tadMinStrength TADMINSTRENGTH] --outputDir
                             OUTPUTDIR
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--prefix`, `-p` `PREFIX` (required) |  | Output prefix of egaDifferential |
| `--explainPrefix` `EXPLAINPREFIX` |  | Output prefix of egaExplain (default: <prefix dir>/explain/<prefix name> if present) |
| `--names` `NAMES` |  | Condition names A B (default: read from <prefix>_state_marks.tsv, the edges or the graph) |
| `--gene` `GENE ...` |  | Gene report(s) for these gene names |
| `--community` `COMMUNITY` |  | Community report for this community ID |
| `--communitySide` `{A, B}` | `B` | Condition whose community IDs --community refers to |
| `--labels` `LABELS` |  | Readable names of the conditions A B for the title (default: known names of common data sets, else the condition names) |
| `--annotation` `ANNOTATION` |  | Gene annotation (GTF, GTF.gz or BED with strand) for upstream / downstream of the report gene; without it only the distance is given |
| `--rechip` `RECHIP` |  | Re-ChIP signal per bait (default: <prefix dir>/rechip/rechip_baits.tsv.gz if its columns match the condition names) |
| `--top` `TOP` | `2` | Genes shown per cascade level in the gene cards of the main report |
| `--maxCommunities` `MAXCOMMUNITIES` | `3` | Coordinated communities drawn in the main report |
| `--maxNodes` `MAXNODES` | `80` | Largest neighbourhood drawn; larger neighbourhoods keep the promoter, its strongest partners and the paths to the explaining nodes |
| `--overlay` |  | Gene reports: add the plain overlay view, one figure of two panels for the region of the gene: the inputs as a track overlay (called interactions per condition, mark bigWig tracks of both conditions, expression change of the promoters) next to the elements of the report (ranked levels with provenance and contact support, distance control, 3D support, boundary crossings); every panel is also written as a PDF (<name>_gene_<gene>_overlay_inputs.pdf, _overlay_ega.pdf) |
| `--overlayTracksA` `OVERLAYTRACKSA ...` |  | mark=file.bw[,rep.bw] of condition A for the overlay (default: the --tracksA of egaDifferential from <prefix>_support.json; without bigWig files the fragment means of the graph are drawn) |
| `--overlayTracksB` `OVERLAYTRACKSB ...` |  | as --overlayTracksA for condition B |
| `--rnaThreshold` `RNATHRESHOLD` |  | \|log2FC\| of a differential gene (default: from the run log, else 1, as egaDifferential) |
| `--chromatinThreshold` `CHROMATINTHRESHOLD` |  | \|log2FC\| of an explaining mark change (default: from the run log, else 1) |
| `--maxLog2` `MAXLOG2` | `2.0` | Colour scale limit of the log2 changes of nodes |
| `--log` `LOG ...` |  | Run logs or scripts with the egaDifferential command line (default: <prefix>.log) |
| `--name` `NAME` |  | Name of the comparison (default: the prefix name) |
| `--pdf` `{auto, chromium, matplotlib, none}` | `auto` | PDF engine: Chromium headless (auto: if available), matplotlib, or no PDF |
| `--chromium` `CHROMIUM` |  | Chromium executable (default: snap or PATH) |
| `--tads` `TADS` |  | TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, _tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools insulation table, BED of boundaries or domains (optional score column), or egaBoundaries output |
| `--tadFormat` `{auto, hicexplorer, arrowhead, cooltools, bed, bed-boundaries, bed-domains, ega}` | `auto` | Format of --tads |
| `--tadWindow` `TADWINDOW` |  | cooltools: insulation window (default: the largest in the table) |
| `--tadChain` `TADCHAIN` |  | Chain file from the assembly of --tads to that of the graph (UCSC chromosome names) |
| `--tadRegions` `TADREGIONS` |  | BED of named regions (chrom start end name) used to name the partner regions |
| `--tadMinStrength` `TADMINSTRENGTH` |  | Only boundaries at or above this strength percentile (0-1) count as crossed |
| `--outputDir`, `-o` `OUTPUTDIR` (required) |  | Output directory |
