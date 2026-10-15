#!/usr/bin/env bash
# RNA-seq of Kasumi-1 four days after siRNA knockdown of RUNX1-ETO (siRE) or mismatch control (siMM), three
# replicates each (Ptasinska et al. 2014, GSE54478): salmon quantification against GENCODE v29 transcripts, gene
# level counts and PyDESeq2 (analysis/kasumi_deseq2.py).
#
# Usage: kasumi_rna.sh <rna dir with samples.tsv and <run>_{1,2}.fastq.gz> <gencode.v29.transcripts.fa.gz>
#                      <salmon binary> <output dir> [threads]
set -euo pipefail
rna=${1:?usage: kasumi_rna.sh <rna dir> <transcripts.fa.gz> <salmon> <output dir> [threads]}
transcripts=${2:?}
salmon=${3:?}
out=${4:?}
threads=${5:-16}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
index="$out/salmon_index_gencode_v29"
[[ -s "$index/versionInfo.json" ]] || "$salmon" index --gencode -t "$transcripts" -i "$index" -p "$threads" \
    > "$out/salmon_index.log" 2>&1
tail -n +2 "$rna/samples.tsv" | while IFS=$'\t' read -r run gsm title condition replicate; do
    [[ -s "$out/$run/quant.sf" ]] || "$salmon" quant -i "$index" -l A -1 "$rna/${run}_1.fastq.gz" \
        -2 "$rna/${run}_2.fastq.gz" -p "$threads" --validateMappings -o "$out/$run" > "$out/$run.log" 2>&1
done
python "$here/kasumi_deseq2.py" "$rna/samples.tsv" "$out" "$out/deseq2_siRE_vs_siMM.tsv" "$transcripts"
