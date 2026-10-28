#!/usr/bin/env bash
# eQTL validation of graph paths: eQTL variant-gene pairs against distance-matched non-eQTL pairs at 1, 2 and 3
# steps, on the GM12878 graphs (matched cell type for LCL eQTLs), the K562 graph (mismatch control) and the
# virtual-capture GM12878 graph, each also mirrored. Four eQTL data sets: GEUVADIS LCL SuSiE (PIP >= 0.5), GTEx v8
# LCL and whole blood (lead per gene), eQTLGen whole blood (lead per gene).
# Usage: bash 19_eqtl.sh [jobs=8] [data_dir] [results_dir]
# (data_dir and results_dir default to EGA_RAW / EGA_RESULTS, else /mnt/processed/wolffjoa/ega/{raw,results})
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
export PYTHONPATH="$(dirname "$here")${PYTHONPATH:+:$PYTHONPATH}"
python="${PYTHON:-$HOME/miniconda3/envs/ega/bin/python}"
jobs="${1:-8}"
root=/mnt/processed/wolffjoa/ega
raw="${2:-${EGA_RAW:-$root/raw}}"
res="${3:-${EGA_RESULTS:-$root/results}}"
out="$res/eqtl_validation"
inputs="$out/inputs"
mkdir -p "$out/runs" "$inputs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ ! -s "$inputs/prepare_summary.tsv" ]; then
    "$python" "$here/eqtl_prepare.py" "$raw/eqtl" "$raw/annotation" "$inputs"
fi

design="$raw/chicago_calls/design_hg19_HindIII"
gtf19="$raw/annotation/gencode.v19.annotation.gtf.gz"
gtf38="$raw/annotation/gencode.v29.annotation.gtf.gz"
declare -A graph_args graph_build
graph_args[GM_chicago]="--nodes $res/graphs/GM_chicago_nodes.tsv --edges $res/graphs/GM_chicago_edges.tsv --rmap $design/Digest_Human_HindIII.rmap --baits $design/Digest_Human_HindIII_baits_e75_ID.baitmap --gtf $gtf19"
graph_build[GM_chicago]=hg19
graph_args[GM_mifsud]="--nodes $res/graphs/GM_mifsud_c20_r1_nodes.tsv --edges $res/graphs/GM_mifsud_c20_r1_edges.tsv --rmap $design/Digest_Human_HindIII.rmap --baits $design/Digest_Human_HindIII_baits_e75_ID.baitmap --gtf $gtf19"
graph_build[GM_mifsud]=hg19
graph_args[K562]="--nodes $res/k562/K562_dmso_nodes.tsv --edges $res/k562/K562_dmso_edges.tsv --rmap $design/Digest_Human_HindIII.rmap --liftRmap $raw/annotation/hg19ToHg38.over.chain.gz --baits $raw/chicago_calls/design_hg38_HindIII_lifted/baits_hg38_lifted.bed --gtf $gtf38"
graph_build[K562]=hg38
graph_args[VC_GM12878]="--nodes $res/virtual_capture/graphs/VC_GM12878_nodes.tsv --edges $res/virtual_capture/graphs/VC_GM12878_edges.tsv --binSize 5000 --chromSizes $raw/annotation/hg38.chrom.sizes --baitGtf $gtf38 --gtf $gtf38"
graph_build[VC_GM12878]=hg38

commands="$out/runs/commands.txt"
: > "$commands"
for graph in GM_chicago GM_mifsud K562 VC_GM12878; do
    build="${graph_build[$graph]}"
    for data in geuvadis_susie gtex_lcl gtex_wb eqtlgen; do
        prefix="$out/runs/${graph}__${data}"
        [ -s "${prefix}_models.tsv" ] && continue
        echo "$python -m epigeneticgraphanalyzer.egaEqtl --verbose ${graph_args[$graph]} --leads $inputs/$data.$build.leads.tsv.gz --significant $inputs/$data.$build.significant.tsv.gz --mirror --label ${graph}__${data} -o $prefix > $prefix.log 2>&1" >> "$commands"
    done
done
xargs -P "$jobs" -I{} bash -c '{}' < "$commands"
"$python" "$here/table_eqtl.py" "$out" "$out/tables"
"$python" "$here/figure_eqtl.py" "$out" "$out/figures"
echo EQTL_RUNS_DONE
