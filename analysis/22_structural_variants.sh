#!/usr/bin/env bash
# Structural variants and copy number as explanation levels (egaDifferential --structuralVariants, --svContacts,
# --copyNumber) for the two AML comparisons, written to <results>/sv/ (the outputs of 16_hops_n.sh are not touched).
#   Kasumi-1 siMM / siRE: capture_sv_calls.py on the HiCUP BAMs of both replicates per condition (8_aml_kasumi1.sh,
#     kasumi1/map_runs.sh); the t(8;21) is expected in both conditions.
#   AML t(8;21) / FLT3-ITD (Assi et al. 2019): only HOMER calls (cis) are deposited as processed data, so the first
#     30 M read pairs of each capture Hi-C run (SRR6392078 t(8;21) patient, SRR6392077 FLT3-ITD patient; ENA) are
#     mapped with HiCUP (HindIII, hg38) here and passed to capture_sv_calls.py with the lifted promoter baits.
# Both differential runs use the settings of 16_hops_n.sh plus the SV inputs; then egaExplain and the supplementary
# table (table_sv.py -> table_s32.tex).
set -euo pipefail
data=${1:?usage: 22_structural_variants.sh <data_dir> <results_dir> <kasumi1 dir> <assi2019 dir> [threads] [samtools]}
results=${2:?}
kasumi=${3:?}
assi=${4:?}
threads=${5:-12}
samtools=${6:-samtools}
here=$(cd "$(dirname "$0")" && pwd)
out="$results/sv"
sizes="$data/annotation/hg38.chrom.sizes"
mkdir -p "$out"

# 1. calls per comparison
h="$kasumi/hicup"
[[ -s "$out/kasumi1_sv.bedpe" ]] || python "$here/capture_sv_calls.py" \
    --bam "siMM=$h/SRR7521620/SRR7521620_1_2.hicup.bam,$h/SRR7521621/SRR7521621_1_2.hicup.bam" \
    --bam "siRE=$h/SRR7521622/SRR7521622_1_2.hicup.bam,$h/SRR7521623/SRR7521623_1_2.hicup.bam" \
    --baits "$kasumi/design/kasumi1_hg38_promoters.baitmap" --chromSizes "$sizes" --samtools "$samtools" \
    --threads "$threads" -o "$out/kasumi1" > "$out/kasumi1.log" 2>&1
sub="$assi/chic_subset"
# bowtie2 index and HiCUP digest of hg38 as for Kasumi-1 (kasumi_map_runs.sh)
genomes=${GENOMES:-$HOME/data/genomes}
index=$genomes/hg38/bowtie2/hg38
digest=$(ls "$genomes"/hg38/hicup_digest/Digest_hg38_HindIII*.txt* | head -1)
for run in SRR6392078 SRR6392077; do
    if [[ ! -s "$sub/${run}_1_30M.fastq.gz" ]]; then
        for mate in 1 2; do
            path=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp" |
                   tail -1 | cut -f2 | tr ';' '\n' | grep "_$mate.fastq.gz")
            curl -s --retry 5 "https://$path" | zcat | head -n 120000000 | gzip -1 > "$sub/${run}_${mate}_30M.fastq.gz" || true
        done
    fi
    ls "$sub/hicup/$run"/*.hicup.bam > /dev/null 2>&1 && continue
    mkdir -p "$sub/hicup/$run"
    hicup --bowtie2 "$(command -v bowtie2)" --index "$index" --digest "$digest" --format Sanger \
        --threads $((threads / 2)) --zip --outdir "$sub/hicup/$run" "$sub/${run}_1_30M.fastq.gz" \
        "$sub/${run}_2_30M.fastq.gz" > "$sub/hicup/$run/hicup.log" 2>&1 &
done
wait
[[ -s "$out/aml_sv.bedpe" ]] || python "$here/capture_sv_calls.py" \
    --bam "t821=$(ls "$sub"/hicup/SRR6392078/*.hicup.bam)" --bam "FLT3ITD=$(ls "$sub"/hicup/SRR6392077/*.hicup.bam)" \
    --baits "$data/chicago_calls/design_hg38_HindIII_lifted/baits_hg38_lifted.bed" --chromSizes "$sizes" \
    --samtools "$samtools" --threads "$threads" -o "$out/aml" > "$out/aml.log" 2>&1

# 2. differential runs with the structural variant levels (settings of 16_hops_n.sh)
common=(--verbose --maxHops 3 --genePermutations 200000 --rnaThreshold 1 --chromatinThreshold 1
        --permutations 200 --bootstrap 1000 --seed 0 --threads "$threads"
        --roles both --roleBootstrap 1000 --stateThreshold 2 --stateThresholds 1.5 2 3 --stateWindow 1000
        --svWindow 1)
calls="$kasumi/chicago"
bw="$kasumi/dnase/bigwig"
egaDifferential "${common[@]}" \
    --interactionsA "$calls/siMM_rep1_score1.ibed.gz" "$calls/siMM_rep2_score1.ibed.gz" \
    --interactionsB "$calls/siRE_rep1_score1.ibed.gz" "$calls/siRE_rep2_score1.ibed.gz" \
    --format ibed --replicates all --scoreThreshold 5 --names siMM siRE \
    --markRoles DNase_d2=active DNase_d10=active --stateScale genomeMean \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" \
    --tracksA "DNase_d2=$bw/GSM3430764_DHS_kasumi1_siMM_D2.bw" "DNase_d10=$bw/GSM3430762_DHS_kasumi1_siMM_D10.bw" \
    --tracksB "DNase_d2=$bw/GSM3430765_DHS_kasumi1_siRE_D2.bw" "DNase_d10=$bw/GSM3430763_DHS_kasumi1_siRE_D10.bw" \
    --expression "$results/aml_kasumi1/deseq2_siRE_vs_siMM.tsv" --log2FCColumn log2FoldChange --padjColumn padj \
    --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    --structuralVariants "$out/kasumi1_sv.bedpe" --svContacts "$out/kasumi1_sv_contacts.tsv" \
    --copyNumber "siMM=$out/kasumi1_siMM_copy_number.bedGraph" "siRE=$out/kasumi1_siRE_copy_number.bedGraph" \
    -o "$out/siMM_vs_siRE" > "$out/siMM_vs_siRE.log" 2>&1

bw="$assi/bigwig"
join() { local IFS=,; echo "$*"; }
t821_dhs=$(join "$bw"/GSM289362[789]_t821*_DHS.bw "$bw"/GSM2893630_t821KIT-2_DHS.bw)
itd_dhs=$(join "$bw"/GSM289360[67]_ITD-*_DHS.bw)
egaDifferential "${common[@]}" \
    --interactionsA "$assi/GSM2895201_t821-1_CHiC_interactions.txt.gz" \
    --interactionsB "$assi/GSM2895200_FLT3-ITD_CHiC_interactions.txt.gz" \
    --format homer --baits "$data/chicago_calls/design_hg38_HindIII_lifted/baits_hg38_lifted.bed" --maxFdr 0.01 \
    --referenceScore 2 --names t821 FLT3ITD --markRoles DHS=active --stateScale genomeMean \
    --tss "$data/annotation/gencode.v29.annotation.gtf.gz" --tracksA "DHS=$t821_dhs" --tracksB "DHS=$itd_dhs" \
    --expression "$results/aml_assi/expression_FLT3ITD_vs_t821.tsv" --log2FCColumn log2FoldChange \
    --padjColumn pvalue --alpha 0.05 --geneKey genes \
    --proximitySupport "$results/porec/K562_hiporec_hops.tsv" \
    --structuralVariants "$out/aml_sv.bedpe" --svContacts "$out/aml_sv_contacts.tsv" \
    --copyNumber "t821=$out/aml_t821_copy_number.bedGraph" "FLT3ITD=$out/aml_FLT3ITD_copy_number.bedGraph" \
    -o "$out/t821_vs_FLT3ITD" > "$out/t821_vs_FLT3ITD.log" 2>&1

# 3. per-gene reports and the supplementary table
for comparison in siMM_vs_siRE t821_vs_FLT3ITD; do
    egaExplain -p "$out/$comparison" -o "$out/explain/$comparison" > "$out/explain_$comparison.log" 2>&1
done
python "$here/table_sv.py" "$results" "$results/sv/tables"
