#!/usr/bin/env bash
# HiCUP (HindIII, hg38) for the given Kasumi-1 capture Hi-C runs, in parallel
export PATH=$HOME/miniconda3/envs/capchic/bin:$PATH:$HOME/miniconda3/bin
work=$HOME/data/ega_aml/kasumi1
index=$HOME/data/genomes/hg38/bowtie2/hg38
digest=$(ls $HOME/data/genomes/hg38/hicup_digest/Digest_hg38_HindIII*.txt* | head -1)
threads=${THREADS:-8}
for run in "$@"; do
    out=$work/hicup/$run
    ls $out/*.hicup.bam > /dev/null 2>&1 && continue
    mkdir -p $out
    hicup --bowtie2 $(command -v bowtie2) --index $index --digest $digest --format Sanger --threads $threads --zip \
        --outdir $out $work/chic/${run}_1.fastq.gz $work/chic/${run}_2.fastq.gz > $out/hicup.log 2>&1 &
done
wait
echo MAPPING_DONE "$@"
