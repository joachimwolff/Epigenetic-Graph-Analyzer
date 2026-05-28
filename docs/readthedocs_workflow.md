# Epigenetic Graph Analyzer: Concept and Workflow

## Overview
Epigenetic Graph Analyzer is a graph-based framework for integrating capture Hi-C interactions with RNA-seq and ChIP-seq differential signals.

The central idea, consistent with the manuscript draft and data club presentation, is:

1. Represent capture Hi-C contacts as a graph where bait and interacting elements are nodes and interaction strengths are edge weights.
2. Build a differential graph between two biological conditions.
3. Annotate graph nodes and edges with multi-omics changes.
4. Analyze direct and indirect interaction paths to explain gene regulation.

This documentation separates the workflow into:

1. Differential multi-omics analysis workflow.
2. Standalone SPRITE validation workflow for indirect mappings.

## Scientific Motivation
Capture Hi-C gives pairwise contact evidence. Biological regulation, however, often involves multi-way chromatin hubs.

The framework addresses this by combining:

1. Differential edge behavior in interaction graphs.
2. Differential RNA expression.
3. Differential histone or protein occupancy at genomic loci.
4. Path-based interpretation (bait -> OE -> OE...) for indirect influence.

The standalone SPRITE validation then checks whether graph-inferred indirect triplets are supported by independent 3D contact observations.

## Data Inputs
The project supports use with data assets such as:

1. Capture Hi-C derived graphs for paired conditions.
2. Genome size files for interval/bin operations.
3. Promoter annotation tables.
4. RNA-seq differential results (log2 fold change).
5. ChIP-seq differential peak summaries (log2 fold change and significance).
6. SPRITE cluster files for orthogonal 3D validation.

## Differential Workflow

### 1) Build condition-specific interaction graphs
Generate graph structures from capture Hi-C style interaction tables.

```{mermaid}
flowchart LR
A[Capture Hi-C interactions] --> B[egaPrepareInteractions]
B --> C[egaCreateGraph condition A]
B --> D[egaCreateGraph condition B]
C --> E[Graph A]
D --> F[Graph B]
```

### 2) Merge and compare graph communities
Align communities across conditions and compute interaction fold changes.

```{mermaid}
flowchart LR
E[Graph A] --> G[egaMergeGraphs]
F[Graph B] --> G
G --> H[Merged differential community graphs]
```

### 3) Add RNA and chromatin annotations
Annotate nodes with RNA-seq and bigWig-derived signal summaries, then compute relation-to-local-neighborhood metrics.

```{mermaid}
flowchart LR
H[Merged graph] --> I[egaComputeRLM]
I --> J[Node-level and neighborhood-level tables]
J --> K[egaClusterBaits]
K --> L[Regulatory pattern clusters]
```

### 4) Interpret direct and indirect mechanisms
Interpretation is performed in path-length space:

1. Step 1: direct bait-neighbor effects.
2. Step 2: one intermediate node.
3. Step 3: extended local hubs.

This supports hypothesis classes highlighted in the manuscript/presentation:

1. Expression change with interaction and chromatin change.
2. Expression change with interaction change but weak chromatin shift.
3. Expression change explained through indirect graph paths.

## Standalone SPRITE Validation Workflow
This branch validates whether graph-inferred indirect triplets are present in real multi-way 3D data.

### Core question
If A interacts with B and B interacts with C in the capture Hi-C graph, does SPRITE report co-occurrence of A, B, and C in the same cluster?

### Procedure

```{mermaid}
flowchart LR
A[GraphML with node positions] --> B[Map nodes to genomic bins]
C[SPRITE clusters] --> D[Count triplet support]
B --> E[Graph-derived indirect triplets]
D --> F[Triplet evidence]
E --> G[egaValidateWithSPRITE]
F --> G
G --> H[Supported vs unsupported indirect triplets]
```

### Modes
`egaValidateWithSPRITE` provides two counting modes:

1. `targeted` (default): low-memory mode that counts only graph-derived triplets.
2. `full`: counts all SPRITE triplets (higher memory footprint).

### Outputs
The validation step yields:

1. Detailed table of each indirect graph triplet with SPRITE support count.
2. Summary table with total triplets, supported fraction, and thresholds used.

## Recommended End-to-End Execution

1. Standardize interaction input: `egaPrepareInteractions`.
2. Build graphs per condition: `egaCreateGraph`.
3. Merge and derive differential communities: `egaMergeGraphs`.
4. Annotate and compute local regulation metrics: `egaComputeRLM`, optionally `egaMergeRLM`.
5. Cluster and visualize regulatory patterns: `egaClusterBaits`, `egaPlotlyGraph`.
6. Validate indirect triplets with orthogonal 3D data: `egaValidateWithSPRITE`.

## Practical Notes

1. For large SPRITE files, use `targeted` mode and optionally set `--maxClusterBins`.
2. Ensure graph nodes carry valid genomic `position` attributes before SPRITE validation.
3. Keep RNA-seq and ChIP-seq contrasts matched to the same condition pair used for graph comparison.

## Scope Boundaries

1. Differential branch explains condition-dependent rewiring and regulatory consequences.
2. Standalone SPRITE branch is validation-focused and can run independently once a graph is available.
3. SPRITE support strengthens confidence in indirect mappings but does not by itself establish causality.
