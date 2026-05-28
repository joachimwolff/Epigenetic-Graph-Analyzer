# epigenetic_graph_analyzer

Epigenetic Graph Analyzer integrates capture Hi-C interactions with RNA-seq and epigenetic signals and supports graph-based downstream analysis.

## Installation

```bash
pip install -r requirements.txt
pip install .
```

## CLI Tools

Core tools:

- `egaCreateGraph`
- `egaMergeGraphs`
- `egaClusterBaits`
- `egaPlotlyGraph`
- `egaComputeRLM`
- `egaMergeRLM`
- `egaValidateWithSPRITE`

Utility tools:

- `egaPrepareInteractions`
- `egaFetchPromoters`
- `egaReconstructCutSites`
- `egaTADStats`
- `egaCheckDownloads`

## Examples

Validate indirect graph links with SPRITE (low-memory mode by default):

```bash
egaValidateWithSPRITE \
	--sprite /path/to/sprite.clusters.gz \
	--graph /path/to/graph.graphml \
	--genomeSize /path/to/chrom.sizes \
	--outputPrefix /path/to/output/sprite_validation \
	--binSize 10000 \
	--minSpriteCount 1
```

Compute RLM table from graph + bigWig signals:

```bash
egaComputeRLM \
	--graph /path/to/graph.graphml \
	--conditions condA condB \
	--bigwigFiles condA_rep1.bw,condA_rep2.bw condB_rep1.bw,condB_rep2.bw \
	--deseqCounts /path/to/deseq.tsv \
	--steps 2 \
	--outputFolder /path/to/output
```

## Testing

```bash
python -m unittest discover -s tests -p "test_*.py"
```

## License

This project is licensed under the MIT License.