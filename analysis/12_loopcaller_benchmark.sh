#!/usr/bin/env bash
# Virtual capture (CHiCAGO on promoter bins of Hi-C / Micro-C) against established loop callers on the same matrices:
# Mustache and Chromosight at 5 kb, promoter-anchored loops (analysis/loops_to_ibed.py), benchmarked against the
# promoter capture calls of the same cell type exactly as the virtual calls (analysis/virtual_capture_benchmark.py).
# Loop callers run in the conda environment `loopcallers` (chromosight 1.6.3, mustache-hic 1.3.3 from PyPI).
set -euo pipefail
data=${1:?usage: 12_loopcaller_benchmark.sh <data_dir> <virtual capture dir> <results_dir>}
vc=${2:?}
results=${3:?}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/loopcallers"
mkdir -p "$out"
res=5000
chroms=$(seq -f 'chr%g' 1 22; echo chrX)
# name | mcool | real calls | real format | real baits | chain real -> matrix assembly | virtual calls
sets=(
    "GM12878|$vc/mcool/4DNFIIV4M7TF.mcool|$data/chicago_calls/GSE81503_GM12878_CHiCAGO/GSE81503_GM12878_PCHiC_merge_final_washU_text.txt.gz|washu|$data/chicago_calls/design_hg19_HindIII/Digest_Human_HindIII_baits_e75_ID.baitmap|$data/annotation/hg19ToHg38.over.chain.gz|$vc/GM12878/chicago_GM12878.ibed"
    "K562|$vc/mcool/4DNFI18UHVRO.mcool|$data/k562/contacts/GSM7974410_k562_dmso_washU.txt.gz|washu|-|-|$vc/K562/chicago_K562.ibed"
    "H1|$vc/mcool/4DNFI9GMP2J8.mcool|$data/chicago_calls/H1_jung2019_reprocessed/H1.ibed.gz|ibed|-|$data/annotation/hg19ToHg38.over.chain.gz|$vc/H1/hicx_H1.ibed"
)
for entry in "${sets[@]}"; do
    IFS='|' read -r name mcool real format baits chain virtual <<< "$entry"
    baitmap="$vc/$name/$name.baitmap"
    # Mustache, per chromosome
    if [[ ! -s "$out/${name}_mustache.tsv" ]]; then
        for chrom in $chroms; do
            conda run -n loopcallers mustache -f "$mcool" -ch "$chrom" -r $res -p 4 \
                -o "$out/${name}_mustache_${chrom}.tsv" > "$out/${name}_mustache_${chrom}.log" 2>&1 || true
        done
        awk 'FNR > 1 || NR == 1' "$out/${name}"_mustache_chr*.tsv > "$out/${name}_mustache.tsv"
    fi
    # Chromosight, genome-wide (cis loops)
    [[ -s "$out/${name}_chromosight.tsv" ]] || conda run -n loopcallers chromosight detect --pattern loops \
        --threads 4 --no-plotting "$mcool::resolutions/$res" "$out/${name}_chromosight" > "$out/${name}_chromosight.log" 2>&1
    for tool in mustache chromosight; do
        python "$here/loops_to_ibed.py" "$out/${name}_${tool}.tsv" "$tool" "$baitmap" $res "$out/${name}_${tool}.ibed" \
            > "$out/${name}_${tool}_ibed.log" 2>&1
        python "$here/virtual_capture_benchmark.py" "$real" "$format" "$baits" "$chain" "$out/${name}_${tool}.ibed" \
            "$baitmap" $res "$out/${name}_${tool}" 5 ibed - - 0 > "$out/${name}_${tool}_benchmark.log" 2>&1
        # the same number of the highest-scoring virtual calls (cis, within 2 Mb)
        python "$here/top_calls.py" "$virtual" "$out/${name}_${tool}.ibed" "$out/${name}_virtual_top_${tool}.ibed"
        python "$here/virtual_capture_benchmark.py" "$real" "$format" "$baits" "$chain" \
            "$out/${name}_virtual_top_${tool}.ibed" "$baitmap" $res "$out/${name}_virtual_top_${tool}" 5 \
            > "$out/${name}_virtual_top_${tool}_benchmark.log" 2>&1
    done
done
