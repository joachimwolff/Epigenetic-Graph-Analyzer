# Reproducing the manuscript

This file maps every figure and table of the manuscript and its supplement to the script that produces it, the
driver that calls the script and the output file. Paths of outputs are relative to the results directory
(`<results>`); `analysis/` is the folder of this repository.

## Software

* Package and analysis scripts: `environment.yml` (exact versions of the environment that produced the results;
  Python 3.11.17). `mamba env create -f environment.yml`, then `pip install --no-deps -e .`
* External tools, each in its own environment (`envs/*.yml`, all solve with `mamba env create --dry-run`):

| Tool | Version | Environment | Used by |
|---|---|---|---|
| HiCUP | 0.9.2 | `envs/capchic.yml` | `01_h1_chicago.sh`, `kasumi_map_runs.sh`, `26_structural_variants.sh` |
| bowtie2 | 2.5.5 | `envs/capchic.yml` | HiCUP |
| samtools | 1.24 | `envs/capchic.yml` | HiCUP, `capture_sv_*.py` (`SAMTOOLS`) |
| R / CHiCAGO | 4.5.3 / 1.38.0 | `envs/capchic.yml` | `run_chicago.R` (H1-hESC calls) |
| chicagoTools `bam2chicago.sh` | not versioned (CHiCAGO repository) | | `01_h1_chicago.sh`, `kasumi_call_chicago.sh` |
| HiCExplorer CHiCAGO tools (C++) | 4.0.0, local build of the HiCExplorer 4 C++ branch | not packaged | `kasumi_call_chicago.sh` (Kasumi-1 calls) |
| HiCExplorer (hicFindTADs, hicPCA) | 3.7.5 | `envs/hicexplorer.yml` | `23_hic_reference.sh`, `27_boundaries.sh`, `28_tad_crossing.sh` |
| cooler / hic2cool | 0.10.2 / 1.0.1 | `envs/hicexplorer.yml` | `hic_eigenvectors.py`, `boundaries_hic.py`, Kessler et al. 2023 .hic conversion |
| pyGenomeTracks | 3.9 | `envs/hicdiff.yml` | `06_loci.sh` (`PGT`) |
| UCSC bedGraphToBigWig, liftOver | 482 | `envs/ucsc_kent.yml` | `07_aml_assi.sh`, `08_aml_kasumi1.sh` |
| salmon | 2.8.0 | `envs/salmon.yml` | `kasumi_rna.sh` |
| PyDESeq2 | 0.5.4 | `environment.yml` | `hesc_deseq2.py`, `kessler_deseq2.py`, `kasumi_deseq2.py`, `rubin_deseq2.py` |
| MOODS-python | 1.9.4.1 | `environment.yml` | `egaTriggers` (`20_triggers.sh`) |
| chromosight / mustache-hic | 1.6.3 / 1.3.3 | `envs/loopcallers.yml` | `12_loopcaller_benchmark.sh` |
| Chicdiff (R) | 0.6, installed from source | `envs/chicdiff.yml` | `kasumi_chicdiff.sh` |
| ChAseR (R package `chaser`) | 0.0.0.9, installed from source | `envs/chaser.yml` | `supporting_steps.sh competitors`, `15_block_bootstrap.sh` |
| HLM-Nbody | commit a49ea16 | not packaged | `polymer_comparison.py` |

## Complete run

Downloads: `scripts/download_data.sh <data_dir>`; Kasumi-1 and Assi et al. 2019 raw data are mapped with
`analysis/kasumi_map_runs.sh`, `kasumi_call_chicago.sh`, `kasumi_rna.sh` and the HiCUP runs of `26_structural_variants.sh`.
The results of the manuscript were produced by

```bash
THREADS=24 analysis/run_all.sh /mnt/processed/wolffjoa/ega/raw /mnt/processed/wolffjoa/ega/results_final \
    ~/data/ega_aml/kasumi1 ~/data/ega_aml/assi2019 ~/data/ega_virtual_capture ~/data/ega_virtual_capture/mcool
```

(started at commit d548762 on 2026-10-18 and finished the same day, including resumed steps; step list and
durations in `<results>/run_all_steps.tsv`; resume with a start step as 7th argument). `run_all.sh` sets `PYTHONHASHSEED=2` for every step (ties in the community transitions of
`egaDifferential` follow the string hash order) and limits numerical libraries to 4 threads per process.
Raw-data processing (HiCUP, CHiCAGO, salmon, Pore-C conversion, .hic conversion, competitor calls) is reused from
`REUSE_FROM` and listed in `<results>/reused_intermediates.tsv`.

The analyses added for the Bioinformatics version run after `run_all.sh`, on its results directory:

```bash
analysis/31_step_decomposition.sh <data> <results> <results>/revision_bioinf/decomposed 200000 6
analysis/32_step_decomposition_aml.sh <results> <kasumi1 dir> <assi2019 dir> <results>/revision_bioinf/decomposed 200000 6
analysis/30_revision_graph.sh <data> <results>              # B1, B6, B7, B9 to B13 -> <results>/revision_bioinf/
analysis/30_crispr_revision.sh <data> <results> <results>/revision_bioinf   # pairs, B2, B3 (CRISPRi)
analysis/33_runtime.sh <data> <results> <kasumi1 dir> <results>/revision_bioinf/runtime 8   # run time and memory
```

## Random seeds

| Where | Seed |
|---|---|
| `run_all.sh` (all steps) | `PYTHONHASHSEED=2` |
| `egaCreateGraph` (Leiden communities, `02_build_graphs.sh` and others) | `--seed 42` |
| `egaDifferential` (community permutations, bootstrap, per-gene permutations, role bootstrap) | `--seed 0` |
| `egaValidateMultiway` (`03_validate.sh`, `09_porec.sh`) | `--seed 0` |
| `hop_distance.py` (pairs, controls, chromosome bootstrap) | seed = number of steps |
| `block_bootstrap_results.py` (`15_block_bootstrap.sh`), `differential_robustness.py` (`14`), `rewired_null.py` (`30_revision_graph.sh` B1) | `--seed 0` |
| `hub_multiway.py`, `sprite_downsampling.py`, `sprite_pileup.py` | `numpy.random.default_rng(0)` (pile-up controls seed 1) |
| `egaPlotRegion` layout (Fig. 2) | 0 |
| `figure_differential_graph_fa2.py` (ForceAtlas2, Fig. 4e) | `--seed 1` |

## Main text

| Item | Script | Driver / command | Output |
|---|---|---|---|
| Fig. 1a-c (concept) | `figure_concept.py` | `run_all.sh` step `figures` (`supporting_steps.sh figures`) | `figures/figure1_concept{a,b,c}.pdf` (manuscript: `figure1{a,b,c}.pdf`) |
| Fig. 1d (pair test schematic) | `figure_multiway_summary.py` | `python analysis/figure_multiway_summary.py <results> <prefix>`; not in a driver | `<prefix>_test.pdf` (manuscript: `figure1d_test.pdf`) |
| Fig. 2a (chr20:30-31 Mb arcs and SPRITE clusters) | `egaPlotRegion`, pyGenomeTracks | `06_loci.sh` (step `06_loci`) | `figures/figure3_GM_arcs.pdf` |
| Fig. 2b-g (graph, multi-way support, four clusters) | `egaPlotRegion` | `06_loci.sh` | `figures/figure3_GM_graph.pdf`, `figure3_GM_graph_multiway.pdf`, `figure3_GM_graph_multiway_example{1-4}.pdf` |
| Fig. 3a (re-ChIP by state) | `rechip_validation.py` | `16_hops_n.sh` (step `16_hops_n`): `rechip_validation.py <results> <results>/hops_n/rechip` | `hops_n/rechip/rechip_validation.pdf` (manuscript: `fig_rechip_validation.pdf`) |
| Fig. 3b, c (cascades hESC / hNEC, Md / PA2, mark term) | `figure_explanation_cascade.py` on the `31_step_decomposition.sh` runs | `python analysis/figure_explanation_cascade.py <results>/revision_bioinf/decomposed <figdir> hESC_hNEC=hESC_vs_hNEC Md_PA2=Md_vs_PA2` | `<figdir>/explanation_cascade_{hESC_hNEC,Md_PA2}.pdf` (manuscript: `fig_cascade_*.pdf`) |
| Fig. 4a, b (Hi-C of Md / PA2 with domains) | `figure_tad_crossing.py` | `28_tad_crossing.sh` (step `28_tad_crossing`) | `tad_crossing/figures/hic_{Md,PA2}.pdf` (manuscript: `fig_tad_hic_*.pdf`) |
| Fig. 4c, d (Hoxa2 arcs) | `figure_tad_crossing.py` | `28_tad_crossing.sh` | `tad_crossing/figures/arcs_{Md,PA2}.pdf` (manuscript: `fig_tad_arcs_*.pdf`) |
| Fig. 4e (differential graph, ForceAtlas2) | `figure_differential_graph_fa2.py` | `python analysis/figure_differential_graph_fa2.py <results>/hops_n/Md_vs_PA2_union.graphml chr6:52207514 fig_hoxa_locus_fa2.pdf --seed 1` | `fig_hoxa_locus_fa2.pdf` |
| Fig. 4f (Hoxa2 boundary graph) | `figure_promoter_boundary_graph.py` | `python analysis/figure_promoter_boundary_graph.py <results>/hops_n/Md_vs_PA2_union.graphml chr6:52162863 51175000 fig_hoxa2_boundary_graph.pdf --far-label "across the boundary"` | `fig_hoxa2_boundary_graph.pdf` |
| Table 1 (what the parts of EGA are checked against) | none, written by hand | | |
| Table 2 (co-occurrence by step, shared at two steps) | `hop_distance.py --rates` | steps 1-5 equal `porec/{GM_chicago_hiporec,GM_chicago_porec,K562_hiporec}_hops.tsv` (`09_porec.sh`) and `reviewer_controls/{GM_chicago_sprite,hESC_gam,H1_gam}_hops_coverage.tsv` (`11_reviewer_controls.sh`); the shared fractions (columns `together`, `together_control`; GAM `share`, `share_control`, step 2) from `python analysis/hop_distance.py <graph> <multiway.npz> <out.tsv> [--chain <chain>] --rates` for GM_chicago x SPRITE / HiPore-C / Pore-C, hESC and H1 x GAM (`--chain hg19ToHg38`), K562_dmso x K562 HiPore-C (`--chain hg38ToHg19`) | `<results>/absolute_rates/*.tsv` |

## Supplement

| Item | Script | Driver | Output |
|---|---|---|---|
| Fig. S1a-c (multi-way support by step, hubs) | `figure_multiway_summary.py` | `python analysis/figure_multiway_summary.py <results> <prefix>`; inputs from `09_porec.sh`, `11_reviewer_controls.sh`, step `hubs` | `<prefix>_{a,b,c}.pdf` (manuscript: `figure2_support_{a,b,c}.pdf`) |
| Fig. S2 (SPRITE at fragment resolution) | `sprite_downsampling.py`, `sprite_pileup.py`, `figure_sprite_resolution.py` | steps `sprite_resolution` and `figures` | `diagnostics/GM_mifsud_spritefrag_downsampling.tsv`, `diagnostics/GM_mifsud_spritefrag_cov_pileup.npz`, `figures/figureS1_sprite_resolution.pdf` |
| Fig. S3a, b (TAD-separation score, CTCF) | `figure_tad_crossing.py` | `28_tad_crossing.sh` | `tad_crossing/figures/{insulation,ctcf}.pdf` (manuscript: `fig_tad_insulation.pdf`, `fig_tad_ctcf.pdf`) |
| Table S1 (graphs) | `supplement_tables.py` | `02_build_graphs.sh`, step `k562_graph`, `04`/`05`/`08` (condition graphs); table: step `tables` | `graphs/*_{nodes,edges,communities}.tsv`, `tables/` |
| Table S2 (data sets) | none, written by hand | | |
| Table S3 (multi-way support by step with chromosome bootstrap intervals) | `hop_distance.py` | `11_reviewer_controls.sh` (coverage and compartment / domain matching), `09_porec.sh`; absolute rates as Table 2 | `reviewer_controls/*_hops_{coverage,domains}.tsv`, `porec/*_hops.tsv`, `absolute_rates/*.tsv` |
| Table S4 (by pair class) | `egaValidateMultiway`, `celltype_specific_edges.py`; `supplement_tables.py` | `03_validate.sh`, `09_porec.sh`, `11_reviewer_controls.sh`; step `tables` | `validation/*_summary.tsv`, `porec/*_summary.tsv`, `reviewer_controls/GM_vs_hESC_domains.tsv` |
| Table S5 (per-contact support) | `egaValidateMultiway --confirmPairs`; `supplement_tables.py` | `03_validate.sh`, `09_porec.sh`; step `tables` | `validation/*_confirm_summary.tsv`, `porec/*_confirm_summary.tsv` |
| Table S6 (hubs) | `hub_multiway.py`, `egaValidateMultiway` (concatemers), `block_bootstrap_results.py` | step `hubs`, `09_porec.sh`, `15_block_bootstrap.sh` | `hubs/*_hub_summary.tsv`, `porec/*_hubs_hub_summary.tsv`, `bootstrap/` |
| Table S7 (expression and chromatin models) | `egaDifferential --stepFeature decomposed` | `31_step_decomposition.sh`, `32_step_decomposition_aml.sh` | `revision_bioinf/decomposed/*_nested_models.tsv`, `*_permutation.json`, `*_gene_classes.tsv`, `*_indirect_gene_test.tsv` |
| Table S8a (state transitions) | `egaDifferential`; `supplement_tables.py` | `16_hops_n.sh`; step `tables` | `hops_n/*_bivalent_transitions.tsv` |
| Table S8b (re-ChIP) | `rechip_validation.py` | `16_hops_n.sh` | `hops_n/rechip/rechip_validation.tsv` |
| Table S9 (mark roles) | `blind_roles.py`, `egaDifferential --roles both`, `step_decomposition_summary.py` | `21_blind_roles.sh`, `31_step_decomposition.sh`; the held-out rows come from runs without a driver in this repository | `hops_n/blind/blind_roles.tsv`, `revision_bioinf/decomposed/*_mark_roles.tsv`, `revision_bioinf/heldout_roles/heldout_roles_agreement*.tsv` |
| Table S10 (report levels) | `egaDifferential --stepFeature decomposed`, `egaExplain` | `31_step_decomposition.sh`, `32_step_decomposition_aml.sh` | `revision_bioinf/decomposed/*_explanation_levels.tsv`, `*_explanation_cascade.tsv` |
| Table S11 (robustness to threshold and depth) | `differential_robustness.py`, `support_summary.py`, `subthreshold_contacts.py`; `supplement_tables.py` | `14_differential_robustness.sh`, `29_threshold_support.sh`, `13_candidate_robustness.sh`; step `tables` | `robustness_differential/*_{variants,status_transfer}.tsv`, `flags_demo/support_summary.tsv`, `robustness_candidates/*.json` |
| Table S12 (boundaries crossed) | `egaTadCrossing`, `tad_crossing_summary.py`, `table_tad_crossing.py` | `28_tad_crossing.sh`; step `tables` | `tad_crossing/tables/table_s36.tex` |
| Table S13 (further modules) | none, written by hand; numbers from `18`, `19`, `20`, `22`, `24`, `26`, `27` and step `virtual_capture` | | |
| Table S14 (features of related methods) | none, written by hand (literature) | | |
| Table S15 (sensitivity at two steps) | `hop_distance.py`, `rewired_null.py`, `graph_comparability.py`, `path_ensemble.py`, `ascertainment.py`, `dependence.py`, `balance.py`, `gam_windows.py`, `liftover_report.py` | `30_revision_graph.sh` (B1, B6, B7, B9 to B13) | `revision_bioinf/B*/`, summary `revision_bioinf/SUMMARY_graph.md` |

## Minimal example

`example/run_example.sh` runs `egaCreateGraph`, `egaDifferential`, `egaExplain` and `egaReport` on the Md / PA2
comparison cut to the Hoxa locus (chr6:49.5-54.5 Mb, mm10; 2.2 MB of public data in `example/data/`, cut by
`example/prepare_example.py`) with the settings of `31_step_decomposition.sh`, in about 45 s on two cores. The key
outputs are compared with `example/expected/` by `tests/test_example.py` (`pytest -m slow tests/test_example.py`).
As in the genome-wide run, Hoxa3 is explained at step 2 by the H3K27ac gain at the Hoxa2 promoter; P values differ
from the genome-wide run because the per-gene null draws from the window only.

## Run time and memory

`analysis/33_runtime.sh` (GNU time; every run pinned to 8 cores with taskset and run alone; workstation with 32
cores and 125 GB of memory, 2026-10-19). Peak memory is the maximum resident set size of the largest single process;
the worker processes share the input by fork. Outputs in
`/mnt/processed/wolffjoa/ega/results/revision_bioinf/runtime/` (`runtime.tsv`); the differential and validation
outputs are identical to those of the manuscript runs.

| Run | Input | Wall time | Peak memory |
|---|---|---|---|
| `egaCreateGraph` with Leiden communities | largest graph: GM12878 strong core of Mifsud et al. 2015, 207,957 nodes, 332,078 edges | 32 s | 2.9 GB |
| `egaDifferential --stepFeature decomposed`, 200,000 per-gene permutations, 200 community permutations, 1,000 bootstrap samples, `--threads 8` | largest comparison: Kasumi-1 siMM / siRE, two replicates each, union graph of 147,525 nodes and 286,493 edges | 8.1 min | 4.2 GB |
| `egaValidateMultiway` against GM12878 SPRITE, `--threads 8` | GM12878 graph (CHiCAGO, 88,124 edges) | 3.2 min | 4.1 GB |
| `egaValidateMultiway` against GM12878 SPRITE, `--threads 8` | GM12878 strong core (332,078 edges) | 12.9 min | 8.7 GB |
