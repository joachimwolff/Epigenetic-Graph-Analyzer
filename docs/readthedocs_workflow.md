# Concept and workflow

## Graph model

Promoter capture Hi-C measures contacts between a baited promoter fragment and any other fragment. EGA
turns the significant interactions into an undirected graph:

* **nodes** are restriction fragments, keyed by `chrom:start-end`; a fragment is a *bait* node if it is
  captured, otherwise an *other-end* (OE) node. A baited fragment that is the other end of a
  promoter-promoter interaction is the same node as the bait, so paths can pass through promoters;
* **edges** are significant interactions, weighted by log2(observed/expected) or the input score.

Because nodes are keyed by coordinates, graphs of different conditions built from the same fragment
set share their nodes and can be compared edge by edge.

## Indirect contacts

Two nodes on the same chromosome without a significant interaction but with a common neighbour (the
mediator) form an indirect contact. Pairs are classified by the capture design:

| Class | Meaning |
|---|---|
| `oe-oe` | neither end baited: the contact cannot be measured by promoter capture |
| `bait-oe` | measurable, not called significant |
| `bait-bait` | both ends baited, measurable, not called significant |

Each pair carries the number of mediators, the Adamic-Adar index and a path score (sum over mediators of
the weaker edge). Two-hop paths u - v - w are reported as triplets.

```{mermaid}
flowchart LR
A[raw pcHi-C counts] --> B[egaCallInteractions]
B --> C[egaCreateGraph]
C --> D[egaIndirectContacts]
C --> E[egaAnnotateGraph]
F[SPRITE clusters / GAM tables] --> G[egaBuildMultiway]
C --> H[egaValidateMultiway]
G --> H
B --> I[egaDifferential]
J[bigWig marks, expression] --> I
```

## Interaction calling

`egaCallInteractions` models the expected count of a bait-other pair as
`f(distance) * s(bait) * t(other end)`. The distance decay is estimated in log-spaced distance bins,
bait factors by empirical Bayes shrinkage, and the other-end factor as a power law of the fragment
coverage (or per fragment, or not at all). Counts are tested against a negative binomial with a
maximum likelihood dispersion. Candidate interactions are excluded from the background fit in two
refit rounds. Pairs are called by BH-FDR or, as in Jung et al. (2019), by a nominal p-value and a
minimum count; in the latter case the expected number of false calls under the fitted background is
reported.

## Validation with multi-way data

SPRITE clusters and GAM tubes record sets of loci that were together in one crosslinked complex or one
nuclear slice, without proximity ligation. `egaBuildMultiway` bins them into a sparse incidence matrix
(SPRITE clusters weighted 2/n). `egaValidateMultiway` then tests

1. **pairs** (direct edges and indirect contacts): co-occurrence, co-occurrence relative to the
   genome-wide expectation at the same distance, and for GAM the NPMI, each compared with control pairs
   at the identical bin offset drawn from the same pools of graph loci (bait or OE) and not connected
   within two hops;
2. **triplets**: the number of groups holding all three loci of a two-hop path, compared with control
   triplets around the same mediator;
3. **communities**: groups holding at least k bins of a community, compared with the same bin set
   circularly shifted along the chromosome.

Graph coordinates can be lifted to the assembly of the multi-way data with `--chain`.

## Differential mode

`egaDifferential` takes the complete call tables of two conditions, chromatin marks (bigWig) and
expression for both:

1. pairs significant in either condition are tested for a change of enrichment with a conditional
   binomial test;
2. on the union graph each bait gene receives, per mark, the log2 change of the mark at the promoter,
   of the contact-weighted mark sum over direct neighbours, and over indirect (two-hop) nodes;
3. expression changes are related to these features by nested linear models (promoter; + direct;
   + indirect) and by a hierarchical classification (promoter > direct > indirect > unexplained);
4. the contribution of the indirect level is tested by permuting mark values over non-promoter nodes.

Communities of both conditions are matched by Jaccard index.
