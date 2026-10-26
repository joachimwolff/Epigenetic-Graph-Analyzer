# Epigenetic Graph Analyzer

The Epigenetic Graph Analyzer (EGA) represents promoter capture Hi-C interactions as a graph of restriction
fragments, annotates the graph with chromatin marks and gene expression, measures with multi-way contact data
(SPRITE, GAM, Pore-C, HiPore-C) how far proximity derived from the graph extends, and writes a gene-centred report
for every differentially expressed gene of a comparison of two conditions: from the promoter over its direct
contacts to elements two and three steps away.

The reports list observations and their support. Which candidate explanation to follow, for example in a
perturbation experiment, is the decision of the analyst.

* Source code: <https://github.com/joachimwolff/Epigenetic-Graph-Analyzer> (MIT licence)
* Every command prints its options with `--help`; the [command reference](reference/index.md) lists them all.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
tutorial
concepts
inputs
```

```{toctree}
:maxdepth: 2
:caption: Use cases

usage/one_condition
usage/differential
usage/multiway
usage/reports
usage/calibration
usage/virtual_capture
usage/modules
running_ega
```

```{toctree}
:maxdepth: 2
:caption: Reference

reference/index
```
