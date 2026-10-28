# Command reference

Every command prints the same information with `--help`. These pages are generated from the argument parsers by `docs/generate_reference.py`.

## Graphs

* [egaCreateGraph](egaCreateGraph.md): Build the promoter capture Hi-C interaction graph, annotate baited promoters with genes and partition the graph into communities
* [egaIndirectContacts](egaIndirectContacts.md): Infer indirect contacts (pairs sharing a neighbour without a direct significant interaction) and list the two-hop paths (triplets) of a graph
* [egaAnnotateGraph](egaAnnotateGraph.md): Annotate graph nodes with the mean bigWig signal over each fragment and baits with gene expression values

## Multi-way data

* [egaBuildMultiway](egaBuildMultiway.md): Bin SPRITE clusters or GAM segregation tables into a sparse bin x group incidence matrix (.npz) used by egaValidateMultiway
* [egaValidateMultiway](egaValidateMultiway.md): Validate graph-derived direct and indirect contacts, two-hop triplets and communities with ligation-free multi-way contact data against distance-matched controls

## Two conditions and reports

* [egaDifferential](egaDifferential.md): Compare the interaction graphs of two conditions: gained, lost and shared interactions (optionally tested by Chicdiff), contact-weighted chromatin changes at promoter, direct and indirect (two-hop) level and up to --maxHops graph steps, and their relation to expression changes
* [egaExplain](egaExplain.md): egaExplain: an integrated report of a differential analysis, linking interaction changes, chromatin changes and expression changes
* [egaReport](egaReport.md): Self-contained HTML and PDF report of a differential run of egaDifferential (and egaExplain): summary, explanation cascade, major findings with figures, hubs and communities, methods; or reports of single genes (--gene) or of one community (--community)
* [egaExplorer](egaExplorer.md): Self-contained offline HTML explorer of a differential run of egaDifferential: gene search with the text of egaReport, interactive arc view and neighbourhood graph, and browse views of the differential genes and communities

## Plots

* [egaPlotNeighbourhood](egaPlotNeighbourhood.md): Plot the neighbourhood of a gene: nodes on a genomic axis, arcs for direct contacts, dashed arcs for indirect contacts of the gene promoter
* [egaPlotLocus](egaPlotLocus.md): Plot a region of the interaction graph (arcs of direct and indirect contacts) above the multi-way contact map of the same region with the graph contacts marked on it
* [egaPlotRegion](egaPlotRegion.md): Write pyGenomeTracks input (interaction arcs coloured by community, indirect contacts) for a region and draw the same region as a graph with the promoters in the centre
* [egaPlotDifferential](egaPlotDifferential.md): Plot the differential interaction graph of egaDifferential for a region or the neighbourhood of a gene: gained, lost and shared interactions as arcs and as a promoter-centred graph with chromatin and expression changes

## Further modules

* [egaCompartments](egaCompartments.md): A/B compartment estimates per genomic bin from chromatin tracks (first principal component oriented by the active marks; active minus repressive role score) and community / edge statistics of a graph against TADs (experimental)
* [egaBoundaries](egaBoundaries.md): Per-bin boundary evidence: capture insulation dip (graph edges spanning each bin against a distance-matched mirror expectation), cohesin-bound CTCF with motif orientation (reverse sites on the left, forward on the right), convergent CTCF loops among the edges; boundaries by a rule with fixed thresholds (no fitted weights) (experimental)
* [egaTadCrossing](egaTadCrossing.md): Per edge of a differential run (egaDifferential): the TAD boundaries it crosses, the strongest of them and same- vs cross-domain, per condition and by status (gained / lost / shared); per promoter the crossing contacts, the boundary they cross, their partner region and a sentence; genome-wide ranking of promoters (Fisher: gained / lost among crossing vs non-crossing edges, BH) and of boundaries (gained / lost crossing edges vs the expectation from edges of the same length)
* [egaEqtl](egaEqtl.md): Test eQTL variant-gene pairs for graph paths of 1, 2 and 3 steps against distance-matched non-eQTL pairs
