#!/usr/bin/env bash
# Kasumi-1 capture Hi-C per replicate: bam2chicago, then the CHiCAGO pipeline of HiCExplorer v4 (C++), export score >= 1
set -euo pipefail
export PATH=$HOME/miniconda3/envs/capchic/bin:$PATH:$HOME/miniconda3/bin
work=$HOME/data/ega_aml/kasumi1
design=$work/design/kasumi1_hg38_promoters
tools=/mnt/processed/wolffjoa/capturehic-to-hic-raw/chicagoTools
hicx=$HOME/src/HiCExplorer-v4/cpp/build/tools
declare -A sample=([SRR7521620]=siMM_rep1 [SRR7521621]=siMM_rep2 [SRR7521622]=siRE_rep1 [SRR7521623]=siRE_rep2)
mkdir -p $work/chinput $work/chicago
for run in "${!sample[@]}"; do
    name=${sample[$run]}
    (
    [[ -s $work/chinput/$name/$name.chinput ]] || (cd $work/chinput && bash $tools/bam2chicago.sh \
        $work/hicup/$run/${run}_1_2.hicup.bam $design.baitmap $design.rmap $name > $name.log 2>&1)
    cd $work/chicago
    $hicx/chicChicagoBackgroundModel --rmap $design.rmap --baitmap $design.baitmap --chinput $work/chinput/$name/$name.chinput \
        --nperbin $design.npb --nbaitsperbin $design.nbpb --proxOE $design.poe --threads 6 -o ${name}_background.txt > ${name}_background.log 2>&1
    $hicx/chicChicagoScores --chinput $work/chinput/$name/$name.chinput --backgroundModel ${name}_background.txt \
        --baitmap $design.baitmap --rmap $design.rmap --threads 6 -o ${name}_scores.txt > ${name}_scores.log 2>&1
    $hicx/chicChicagoSignificantInteractions --scores ${name}_scores.txt --rmap $design.rmap --baitmap $design.baitmap \
        --scoreThreshold 1 --threads 6 -o ${name}_score1.txt > ${name}_significant.log 2>&1
    awk -F'\t' -v OFS='\t' 'NR==1{print "bait_chr","bait_start","bait_end","bait_name","otherEnd_chr","otherEnd_start","otherEnd_end","otherEnd_name","N_reads","score"; next} {print $1,$2,$3,$4,$5,$6,$7,".",$8,$10}' ${name}_score1.txt | gzip > ${name}_score1.ibed.gz
    echo "done $name"
    ) &
done
wait
echo KASUMI_CHICAGO_DONE
