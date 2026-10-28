# Validation with multi-way data

Multi-way data show whether the indirect contacts of a graph are spatial proximity in single nuclei. EGA uses SPRITE,
GAM, Pore-C and HiPore-C.

## 1. Multi-way data as a matrix

```bash
# SPRITE clusters, reads assigned to the restriction fragments of the capture design
egaBuildMultiway --type sprite --input GM12878.clusters.gz --chromSizes hg19.chrom.sizes \
    --fragments HindIII_fragments.bed -o multiway/GM12878_sprite.npz

# Pore-C / HiPore-C: convert the concatemers to the cluster layout first
python analysis/porec_to_clusters.py porec_fragments.parquet porec_clusters.txt.gz \
    --read read_name --chrom chrom --start start --end end --minMapq 1 --mapq mapping_quality \
    --chain hg38ToHg19.over.chain.gz
egaBuildMultiway --type sprite --input porec_clusters.txt.gz --chromSizes hg19.chrom.sizes \
    --fragments HindIII_fragments.bed -o multiway/GM12878_porec.npz

# GAM: one segregation table per replicate, in the windows of the table
egaBuildMultiway --type gam --input H1_gam_rep1.txt.gz H1_gam_rep2.txt.gz --chromSizes hg38.chrom.sizes \
    -o multiway/H1_gam.npz
```

SPRITE clusters with 2 to 1,000 reads are kept (`--minClusterSize`, `--maxClusterSize`).

## 2. Pairs, paths and communities against controls

```bash
egaValidateMultiway -g graphs/GM12878.graphml -m multiway/GM12878_sprite.npz \
    --confirmPairs --coverageClasses 5 --controls 10 --maxPairs 300000 --maxTriplets 100000 \
    --seed 0 --threads 8 -o validation/GM12878_sprite
```

* Every direct and indirect pair is compared with `--controls` pairs at the identical offset, of the same node type,
  in the same coverage classes (`--coverageClasses`) and not connected within two steps.
* `--compartments` and `--domains` add the matching of A/B compartments and domain membership.
* Two-step paths (triplets) are compared with control triplets around the same mediator, communities with shifted
  copies.
* `--chain` lifts the graph to the assembly of the multi-way data.

The outputs give, per step, the probability that a graph pair co-occurs more often than its controls, and the
absolute share of pairs found together.

## 3. Plots

```bash
# a community or region of the graph above its multi-way events
egaPlotLocus -g graphs/GM12878.graphml -m multiway/GM12878_sprite.npz --region chr20:30000000-31000000 \
    --style events -o figures/locus.pdf
# arcs with the SPRITE clusters below
egaPlotRegion -g graphs/GM12878.graphml --region chr20:30000000-31000000 \
    --multiway multiway/GM12878_sprite.npz -o figures/region
```

## Paths of single genes

The explaining paths of a gene (promoter, mediator, element) can be tested in the multi-way data of the same cell
type; `analysis/36_gm12878_multiway_paths.py` counts the groups holding all three fragments against controls shifted
by a whole number of restriction fragments.
