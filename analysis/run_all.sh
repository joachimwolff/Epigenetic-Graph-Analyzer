#!/usr/bin/env bash
# Complete analysis of the manuscript in dependency order, from the interaction calls to the supplementary tables,
# figures, reports and explorers, all written to one results directory.
#
# Usage: run_all.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> <virtual capture dir> <mcool dir>
#                   [start step] [stop step]
#   data_dir             downloads (scripts/download_data.sh), e.g. /mnt/processed/wolffjoa/ega/raw
#   results_dir          output directory (created; a new directory gives a complete run)
#   kasumi1 dir          Kasumi-1 HiCUP BAMs, chinput, CHiCAGO calls, DNase bigWigs, salmon quantification
#                        (kasumi_map_runs.sh, kasumi_call_chicago.sh, kasumi_rna.sh)
#   assi2019 dir         Assi et al. 2019 downloads, bigWigs and the HiCUP runs of the capture Hi-C read subsets
#   virtual capture dir  CHiCAGO calls on promoter bins of the 4DN Hi-C / Micro-C matrices (virtual_capture.py)
#   mcool dir            4DN mcool files (4DNFIIV4M7TF GM12878, 4DNFI18UHVRO K562, 4DNFI9GMP2J8 H1)
#   start step           first step to run (name or number from the list below; default: the first), to resume a
#                        failed run; stop step: last step to run (default: the last)
#
# Raw-data processing is not repeated when its output exists (REUSE_FROM, default the directory of the previous
# run, /mnt/processed/wolffjoa/ega/results): the step 'reuse' copies these intermediates into <results_dir> and
# lists them in <results_dir>/reused_intermediates.tsv:
#   H1 HiCUP / bam2chicago / CHiCAGO calls (01_h1_chicago.sh; <data_dir>/chicago_calls/H1_jung2019_reprocessed),
#   Kasumi-1 HiCUP, CHiCAGO calls and salmon (kasumi1 dir), Assi capture Hi-C read subsets mapped with HiCUP
#   (assi2019 dir), virtual capture CHiCAGO calls (virtual capture dir), Pore-C / HiPore-C read conversion
#   (porec/*_clusters.txt.gz), the Kessler et al. 2023 .hic files converted to 25 kb coolers
#   (kessler_tad_crossing/*.cool), the thinned-depth CHiCAGO recalls of Kasumi-1
#   (robustness_differential/depth), and the competitor calls (Mustache and Chromosight loops of
#   12_loopcaller_benchmark.sh, Chicdiff of kasumi_chicdiff.sh). Every analysis from the graphs onward is recomputed.
#
# Environment: EGA_PREFIX (conda env with the package dependencies, default ~/miniconda3/envs/ega),
# HICEXPLORER_PREFIX (default ~/miniconda3/envs/hicexplorer), THREADS (largest number of threads of one step,
# default 24; worker processes of the multi-way validation), DIFFERENTIAL_PERMUTATIONS (community permutations of
# egaDifferential, default 200; PERMUTATIONS is read by 13_candidate_robustness.sh), H1_WORK / BOWTIE2_HG19 /
# HICUP_DIGEST_HG19 / CHICAGO_TOOLS (inputs of 01_h1_chicago.sh when its calls are missing), SAMTOOLS, PGT
# (pyGenomeTracks), HLM_NBODY.
# Every step is logged with start, end and duration in <results_dir>/run_all_steps.tsv; the output of the drivers
# goes to the log files of each step and to standard output.
set -euo pipefail
usage="usage: run_all.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> <virtual capture dir> <mcool dir> [start step] [stop step]"
data=$(realpath "${1:?$usage}")
results=${2:?$usage}
kasumi=$(realpath "${3:?$usage}")
assi=$(realpath "${4:?$usage}")
vc=$(realpath "${5:?$usage}")
mcool=$(realpath "${6:?$usage}")
start=${7:-}
stop=${8:-}
mkdir -p "$results"
results=$(realpath "$results")
here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
threads=${THREADS:-24}
permutations=${DIFFERENTIAL_PERMUTATIONS:-200}
reuse_from=${REUSE_FROM:-/mnt/processed/wolffjoa/ega/results}
ega=${EGA_PREFIX:-$HOME/miniconda3/envs/ega}
hicexplorer=${HICEXPLORER_PREFIX:-$HOME/miniconda3/envs/hicexplorer}

# the package of this checkout, the tools of the ega environment, UCSC tools for bigWig conversion
export PATH="$ega/bin:$PATH:$HOME/miniconda3/envs/__ucsc_kent/bin"
export PYTHONPATH="$repo${PYTHONPATH:+:$PYTHONPATH}"
export PYTHON="$ega/bin/python"
# Python's per-process string hashing orders sets of node names; ties in the community transitions of egaDifferential
# follow that order (see community_robustness.py, where 2 reproduces the stored results), so the run pins it
export PYTHONHASHSEED=${PYTHONHASHSEED:-2}
export HICEXPLORER_BIN="$hicexplorer/bin/" HICEXPLORER_PYTHON="$hicexplorer/bin/python"
export EGA_RAW="$data" EGA_RESULTS="$results" EGA_HOPS="$results/hops_n" EGA_MCOOL="$mcool"
export SAMTOOLS=${SAMTOOLS:-$HOME/miniconda3/envs/capchic/bin/samtools}
export PGT=${PGT:-$HOME/miniconda3/envs/hicdiff/bin/pyGenomeTracks}
# numerical libraries: at most 4 threads per process unless a driver sets its own limit
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4} \
       MKL_NUM_THREADS=${MKL_NUM_THREADS:-4} NUMEXPR_MAX_THREADS=${NUMEXPR_MAX_THREADS:-4}
figures="$results/figures"
support() { bash "$here/supporting_steps.sh" "$1" "$data" "$results" "$kasumi" "$vc" "${2:-8}"; }

# ---------------------------------------------------------------------------------------------------------- steps
step_reuse() {
    local list="$results/reused_intermediates.tsv"
    printf 'intermediate\tsource\treason\n' > "$list"
    copy() {   # <relative path glob> <reason>
        local pattern=$1 reason=$2 f rel
        for f in $reuse_from/$pattern; do
            [[ -e "$f" ]] || { echo "missing intermediate: $f" >&2; return 1; }
            rel=${f#"$reuse_from"/}
            mkdir -p "$results/$(dirname "$rel")"
            [[ -e "$results/$rel" ]] || cp -r "$f" "$results/$rel"
            printf '%s\t%s\t%s\n' "$rel" "$f" "$reason" >> "$list"
        done
    }
    copy 'porec/*_clusters.txt.gz' 'Pore-C / HiPore-C read conversion (porec_to_clusters.py, unchanged since 2a43283)'
    copy 'porec/*_convert.log' 'log of the Pore-C read conversion'
    copy 'kessler_tad_crossing/*_25kb_25000.cool' 'Kessler et al. 2023 .hic converted to 25 kb coolers (hic2cool, balanced)'
    copy 'kessler_tad_crossing/*_convert.log' 'log of the .hic conversion'
    copy 'kessler_tad_crossing/*_balance.log' 'log of the cooler balancing'
    copy 'loopcallers/*_mustache.tsv' 'competitor calls (Mustache 1.3.3 loops on the mcool files)'
    copy 'loopcallers/*_mustache_chr*' 'competitor calls per chromosome and logs'
    copy 'loopcallers/*_chromosight.tsv' 'competitor calls (Chromosight 1.6.3 loops on the mcool files)'
    copy 'loopcallers/*_chromosight.log' 'competitor log'
    copy 'robustness_differential/depth/read_counts.tsv' 'read counts of the Kasumi-1 chinput files'
    copy 'robustness_differential/depth/matched' 'CHiCAGO recalls of the thinned Kasumi-1 chinput (seed 0)'
    copy 'robustness_differential/depth/half' 'CHiCAGO recalls of the thinned Kasumi-1 chinput (seed 0)'
    copy 'robustness_differential/chicdiff/*.Rds' 'competitor calls (Chicdiff inputs and results)'
    copy 'robustness_differential/chicdiff/peakmatrix.txt' 'competitor input (Chicdiff peak matrix)'
    copy 'robustness_differential/chicdiff/siMM_vs_siRE_results.tsv.gz' 'competitor calls (Chicdiff)'
    copy 'robustness_differential/chicdiff/design' 'Chicdiff design files'
    copy 'robustness_differential/chicdiff/*.log' 'Chicdiff logs'
    for f in "$data/chicago_calls/H1_jung2019_reprocessed/H1.ibed.gz" "$kasumi"/chicago/si*_score1.ibed.gz \
             "$kasumi"/hicup/SRR*/*.hicup.bam "$assi"/chic_subset/hicup/SRR*/*.hicup.bam \
             "$vc"/GM12878/chicago_GM12878.ibed "$vc"/K562/chicago_K562.ibed "$vc"/H1/hicx_H1.ibed \
             "$kasumi"/rna_salmon/SRR*/quant.sf; do
        [[ -s "$f" ]] || { echo "missing raw-processing output: $f" >&2; return 1; }
        printf '%s\t%s\t%s\n' "(input, not copied)" "$f" 'raw-data processing outside the results directory' >> "$list"
    done
}
step_01_h1_chicago() {
    local calls="$data/chicago_calls/H1_jung2019_reprocessed/H1.ibed.gz"
    if [[ -s "$calls" ]]; then
        echo "H1 CHiCAGO calls from the raw reads present, reused: $calls"
        return 0
    fi
    local work=${H1_WORK:-$HOME/data/ega_h1_pchic}
    PATH="$HOME/miniconda3/envs/capchic/bin:$PATH" bash "$here/01_h1_chicago.sh" "$work" \
        "${BOWTIE2_HG19:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/hicup/bowtie2_hg19/hg19}" \
        "${HICUP_DIGEST_HG19:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/hicup/Digest_hg19_HindIII_None_15-12-58_17-09-2026.txt}" \
        "${CHICAGO_TOOLS:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/chicagoTools}" "$threads"
    mkdir -p "$(dirname "$calls")"
    gzip -c "$work/chicago/H1.ibed" > "$calls"
}
step_02_graphs() { bash "$here/02_build_graphs.sh" "$data" "$results"; }
step_k562_graph() { support k562_graph; }
step_03_validate() { bash "$here/03_validate.sh" "$data" "$results" "$threads"; }
step_virtual_capture() { support virtual_capture; }
step_kasumi_rna() {
    # PyDESeq2 on the salmon quantification of kasumi_rna.sh (salmon itself is raw-data processing, reused)
    mkdir -p "$results/aml_kasumi1"
    python "$here/kasumi_deseq2.py" "$kasumi/rna/samples.tsv" "$kasumi/rna_salmon" \
        "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" \
        "${GENCODE_V29_TRANSCRIPTS:-$HOME/data/genomes/hg38/gencode_v29/gencode.v29.transcripts.fa.gz}" \
        > "$results/aml_kasumi1/rna.log" 2>&1
}
step_04_differential() { bash "$here/04_differential.sh" "$data" "$results" "$permutations"; }
step_05_kessler() { bash "$here/05_kessler_md_pa2.sh" "$data" "$results"; }
step_07_aml_assi() { bash "$here/07_aml_assi.sh" "$assi" "$data" "$results" "$permutations"; }
step_08_aml_kasumi1() { bash "$here/08_aml_kasumi1.sh" "$kasumi" "$data" "$results" "$permutations"; }
step_06_loci() { bash "$here/06_loci.sh" "$data" "$results" "$figures"; }
step_09_porec() { JOBS=8 bash "$here/09_porec.sh" "$data" "$results" "$threads"; }
step_hubs() { support hubs "$threads"; }
step_10_differential_figures() { bash "$here/10_differential_figures.sh" "$results" "$figures"; }
step_11_reviewer_controls() { bash "$here/11_reviewer_controls.sh" "$data" "$results" "$threads"; }
step_12_loopcallers() { bash "$here/12_loopcaller_benchmark.sh" "$data" "$vc" "$results"; }
step_13_candidate_robustness() { bash "$here/13_candidate_robustness.sh" "$results" "$data" "$kasumi" all; }
step_14_differential_robustness() { bash "$here/14_differential_robustness.sh" "$results" "$data" "$kasumi" all; }
step_explain() { support explain; }
step_competitors() { support competitors; }
step_15_block_bootstrap() { bash "$here/15_block_bootstrap.sh" "$results" 1000; }
step_17_minoux_fold() { bash "$here/17_minoux_fold.sh" "$data" "$results"; }
step_16_hops_n() { bash "$here/16_hops_n.sh" "$data" "$results" "$kasumi" "$assi" "$permutations" "$threads"; }
step_18_switching() { bash "$here/18_enhancer_switching.sh" "$results" 10000 4; }
step_19_eqtl() { bash "$here/19_eqtl.sh" 8 "$data" "$results"; }
step_20_triggers() { bash "$here/20_triggers.sh" "$data" "$results" "$kasumi" "$assi" "$results/hops_n" 8; }
step_21_blind_roles() { bash "$here/21_blind_roles.sh" "$data" "$results" "$permutations" "$threads"; }
step_22_timecourse() { bash "$here/22_timecourse.sh" "$data" "$results" 12 "$permutations"; }
step_23_hic_reference() { bash "$here/23_hic_reference.sh" "$mcool" "$data" "$results" 8; }
step_24_compartments_domains() { bash "$here/24_compartments_domains.sh" 1000; }
step_26_structural_variants() { bash "$here/26_structural_variants.sh" "$data" "$results" "$kasumi" "$assi" 12 "$SAMTOOLS"; }
step_27_boundaries() { bash "$here/27_boundaries.sh" "$mcool" "$data" "$results" 12 500; }
step_28_tad_crossing() { bash "$here/28_tad_crossing.sh" "$data" "$results" 12; }
# threshold sensitivity / provenance demo (Kasumi-1, Md / PA2): needs the raw score files, 16 (settings), 17 (fold
# tracks), 23/24 (K562 TADs) and 28 (PA2 TADs)
step_29_threshold_support() { bash "$here/29_threshold_support.sh" "$data" "$results" "$kasumi" all 16; }
step_community_robustness() { support community_robustness; }
step_crispr_benchmark() { support crispr_benchmark; }
step_kasumi_sv_cnv() { support kasumi_sv_cnv; }
step_sprite_resolution() { support sprite_resolution "$threads"; }
step_tables() {
    local out="$results/tables" t
    mkdir -p "$out"
    python "$here/supplement_tables.py" "$results" "$out" > "$out/supplement_tables.log" 2>&1
    python "$here/table_switching.py" "$results/enhancer_switching" "$results/enhancer_switching"
    python "$here/table_eqtl.py" "$results/eqtl_validation" "$results/eqtl_validation/tables"
    python "$here/table_triggers.py" "$results/triggers" "$results/triggers/tables"
    python "$here/table_timecourse.py" "$results/timecourse_rubin2017" "$results/timecourse_rubin2017/tables"
    python "$here/table_sv.py" "$results" "$results/sv/tables"
    python "$here/table_boundaries.py" "$results/boundaries" "$results/boundaries/tables"
    python "$here/table_tad_crossing.py" "$results/tad_crossing" "$results/tad_crossing/tables"
    for t in "$results"/enhancer_switching/table_s*.tex "$results"/{eqtl_validation,triggers,timecourse_rubin2017,sv,boundaries,tad_crossing}/tables/table_s*.tex; do
        [[ -e "$t" ]] && cp "$t" "$out/"
    done
    ls "$out"
}
step_figures() { support figures; }
step_reports() { support reports; }

steps=(reuse 01_h1_chicago 02_graphs k562_graph 03_validate virtual_capture kasumi_rna 04_differential 05_kessler
       07_aml_assi 08_aml_kasumi1 06_loci 09_porec hubs 10_differential_figures 11_reviewer_controls 12_loopcallers
       13_candidate_robustness 14_differential_robustness explain competitors 15_block_bootstrap 17_minoux_fold
       16_hops_n 18_switching 19_eqtl 20_triggers 21_blind_roles 22_timecourse 23_hic_reference
       24_compartments_domains 26_structural_variants 27_boundaries 28_tad_crossing 29_threshold_support
       community_robustness
       crispr_benchmark kasumi_sv_cnv sprite_resolution tables figures reports)

index_of() {   # step name or its number prefix -> position in the list
    local i
    for i in "${!steps[@]}"; do
        [[ "${steps[$i]}" == "$1" || "${steps[$i]}" == "$1"_* ]] && { echo "$i"; return 0; }
    done
    echo "unknown step: $1 (steps: ${steps[*]})" >&2
    return 1
}
first=0
last=$(( ${#steps[@]} - 1 ))
[[ -n "$start" ]] && first=$(index_of "$start")
[[ -n "$stop" ]] && last=$(index_of "$stop")

timing="$results/run_all_steps.tsv"
[[ -s "$timing" ]] || printf 'step\tstart\tend\tseconds\tstatus\n' > "$timing"
echo "run_all.sh $(date -Is) commit $(git -C "$repo" rev-parse --short HEAD 2>/dev/null || echo unknown)" \
    "$(git -C "$repo" status --porcelain 2>/dev/null | grep -q . && echo '(uncommitted changes)')"
echo "data $data, results $results, kasumi1 $kasumi, assi2019 $assi, virtual capture $vc, mcool $mcool"
echo "steps ${steps[$first]} .. ${steps[$last]}, threads $threads, permutations $permutations"
for (( i = first; i <= last; i++ )); do
    name=${steps[$i]}
    t0=$(date +%s)
    echo "=== $(date -Is) step $name"
    status=0
    # a background subshell keeps errexit active inside the step function
    ( "step_$name" ) &
    wait $! || status=$?
    t1=$(date +%s)
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$(date -d @"$t0" -Is)" "$(date -d @"$t1" -Is)" $((t1 - t0)) "$status" >> "$timing"
    echo "=== $(date -Is) step $name finished with status $status after $((t1 - t0)) s"
    if (( status != 0 )); then
        echo "run_all.sh: step $name failed; resume with start step $name" >&2
        exit "$status"
    fi
done
echo "RUN_ALL_DONE $(date -Is)"
