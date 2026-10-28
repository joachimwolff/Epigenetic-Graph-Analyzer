# Epigenetic Graph Analyzer (EGA)

EGA represents promoter capture Hi-C (pcHi-C) interactions as a graph whose nodes are restriction
fragments and whose edges are significant interactions. Baited promoters and other ends share one node
space, so promoter-promoter interactions connect the neighbourhoods of different genes. On this graph EGA

* infers **indirect contacts**: pairs of fragments without a significant interaction that share a
  neighbour. Pairs of two non-baited fragments (`oe-oe`) cannot be measured by promoter capture at all;
* partitions the graph into **communities** of co-interacting promoters and regulatory elements;
* **validates** direct contacts, indirect contacts, two-hop triplets and communities with ligation-free
  multi-way contact data (SPRITE, GAM) against distance-matched controls;
* in **differential mode** compares two conditions: interaction changes, contact-weighted chromatin
  changes at the promoter, at direct and at indirect (two-hop) neighbours, and their relation to gene
  expression changes.

## Installation

```bash
pip install .            # core
pip install .[leiden]    # Leiden community detection (python-igraph, leidenalg), recommended
```

Python >= 3.9 with numpy, scipy, pandas, networkx, pyBigWig, matplotlib and pyliftover. The exact environment of
the manuscript analyses is `environment.yml` (`mamba env create -f environment.yml`); the external tools of the
analysis drivers are in `envs/`. `REPRODUCE.md` maps every figure and table of the manuscript to its script, and
`example/run_example.sh` runs a small Md / PA2 comparison at the Hoxa locus (about one minute).

## Quick start

```bash
mamba env create -f environment.yml && conda activate ega && pip install .
bash example/run_example.sh          # Md / PA2 comparison at the Hoxa locus, about one minute
```

The outputs are compared with `example/expected/` by `pytest -m slow tests/test_example.py`.

## Command line tools

| Tool | Purpose |
|---|---|
| `egaCreateGraph` | Graph from called interactions (CHiCAGO, HOMER, other callers) or cool matrices; gene annotation of baits; Leiden communities |
| `egaIndirectContacts` | Indirect contacts with pair classes and scores, direct contacts, two-step triplets |
| `egaAnnotateGraph` | Mean bigWig signal per node, expression values per bait |
| `egaBuildMultiway` | SPRITE clusters, GAM segregation tables or concatemers as a sparse locus-by-group matrix |
| `egaValidateMultiway` | Direct and indirect contacts, triplets and hubs against multi-way data with matched controls |
| `egaEqtl` | eQTL enrichment of direct edges and paths against matched non-eQTL pairs |
| `egaDifferential` | Two conditions: union graph, chromatin states, mark roles, contact-weighted mark changes |
| `egaExplain` | Gene-centred integrative reports: all levels from the promoter to step n |
| `egaReport`, `egaExplorer` | Reports as HTML or PDF, interactive browser page |
| `egaTadCrossing` | Gained and lost contacts across domain boundaries (HiCExplorer, Arrowhead, cooltools, BED) |
| `egaPlotNeighbourhood`, `egaPlotLocus`, `egaPlotRegion`, `egaPlotDifferential` | Plots |

Experimental modules (their evaluation gave weak or null results; the output is exploratory and the `--help` text says so):

| Tool | Purpose |
|---|---|
| `egaCompartments` | A/B compartments from chromatin tracks |
| `egaBoundaries` | Rule-based domain boundaries from capture contacts and CTCF |

Every tool prints its options with `--help`.

## Example

```bash
# 1. interactions called by CHiCAGO (ibed; analysis/01_h1_chicago.sh for the H1-hESC reads of Jung et al. 2019)
# 2. graph, communities and indirect contacts
egaCreateGraph -i H1.ibed.gz -f ibed --scoreThreshold 5 --tss gencode.v19.annotation.gtf.gz -o graphs/H1
egaIndirectContacts -g graphs/H1.graphml -o graphs/H1 --triplets

# 3. validation with GAM (hg38) of a graph built on hg19
egaBuildMultiway --type gam --input rep1_30kb.csv rep2_30kb.csv --chromSizes hg38.chrom.sizes -o gam.npz
egaValidateMultiway -g graphs/H1.graphml -m gam.npz --chain hg19ToHg38.over.chain.gz -o validation/H1_gam

# 4. differential mode
egaDifferential --interactionsA H1.ibed.gz --interactionsB NPC.ibed.gz --format ibed --scoreThreshold 5 --names H1 NPC \
    --tss gencode.v19.annotation.gtf.gz \
    --tracksA H3K27ac=E003-H3K27ac.bw H3K27me3=E003-H3K27me3.bw \
    --tracksB H3K27ac=E007-H3K27ac.bw H3K27me3=E007-H3K27me3.bw --markSign H3K27me3=-1 \
    --expression 57epigenomes.RPKM.pc.gz --expressionColumns E003 E007 -o differential/H1_vs_NPC
```

### Mark presence in the chromatin states

The chromatin states of `egaDifferential` (bivalent, active and repressed promoters, poised and active
enhancers, element classes, learned mark combinations) need the presence of every mark on every fragment.
By default a mark is present when its signal over the fragment reaches `--stateThreshold` (fold enrichment
over input, default 2; statistics: fragment mean and maximum 1 kb window, `--stateWindow`; separate tracks
with `--stateTracksA` / `--stateTracksB`).

Peak calls of any caller (MACS2 / MACS3 narrowPeak or broadPeak, SEACR, generic BED; the columns chrom,
start, end are read) can define the presence instead, per mark and condition:

```bash
egaDifferential ... \
    --statePeaksA H3K4me3=H1_H3K4me3_peaks.narrowPeak H3K27me3=H1_H3K27me3_peaks.broadPeak \
    --statePeaksB H3K4me3=NPC_H3K4me3_peaks.narrowPeak H3K27me3=NPC_H3K27me3_peaks.broadPeak \
    --statePeakMinOverlap 1
```

A mark with peaks is present on a fragment when at least one peak overlaps the fragment by at least
`--statePeakMinOverlap` bp (default 1); several files of one mark (`mark=a.bed,b.bed`) are pooled. A
mark given with peaks needs peaks in both conditions. Marks without peaks keep the fold-enrichment rule,
so both rules can be mixed (e.g. peaks for H3K27me3, fold enrichment for H3K4me3). Peaks define presence
only: a state transition still needs the signal of every crossing mark to change by `--stateMinChange`
in the direction of the crossing, the state summaries show the signal values, and the mark changes of
the explanation levels (promoter, steps 1 .. `--maxHops`) always come from the signal tracks. In
`<prefix>_state_marks.tsv` a peak-defined mark has the scale `peaks` and the same fraction of present
fragments at every threshold.

The analyses of the manuscript are scripted in `analysis/` (`scripts/download_data.sh` fetches all
public input data).

## Input formats

* Interactions: Jung et al. 2019 raw count tables (`jung`), Mifsud et al. 2015 tables (`mifsud`),
  CHiCAGO `ibed`, or a BEDPE-like table (`bedpe`: bait chrom/start/end, other chrom/start/end, score[, count]).
* Genes: GTF/GFF (transcript TSS) or BED of TSS.
* SPRITE: cluster files (`barcode<TAB>chrom:position<TAB>...`). GAM: segregation tables
  (`chrom start stop tube1 tube2 ...`, one table per replicate).

## Testing

```bash
pip install .[test,leiden]
pytest
```

The tests run every module and every command line tool on synthetic data with planted interactions,
hubs and expression effects.

## License

MIT
