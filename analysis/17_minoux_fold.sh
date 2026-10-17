#!/usr/bin/env bash
# Fold enrichment over input of the Minoux et al. 2017 ChIP-seq tracks (GEO GSE89435) of the mandibular (Md) and
# second (PA2) pharyngeal arch at E10.5: H3K4me2, H3K27me3, H3K27ac and the sequential ChIP H3K27me3 then H3K4me2
# (re-ChIP, a direct measurement of bivalency). The three input libraries are averaged (minoux_fold_enrichment.py).
# Output: <results>/minoux_fold/{Md,PA2}_{H3K4me2,H3K27me3,H3K27ac,reChIP}.fold.bw, used for the chromatin states of
# the Md / PA2 comparison (16_hops_n.sh --stateTracksA/B) and by rechip_validation.py.
set -euo pipefail
data=${1:?usage: 17_minoux_fold.sh <data_dir> <results_dir>}
results=${2:?}
here=$(cd "$(dirname "$0")" && pwd)
geo="$data/kessler/minoux_gse89435"
reference="$data/kessler/E10.5_MD_H3K27ac_R1.bigwig"   # chromosome sizes (mm10)
out="$results/minoux_fold"
mkdir -p "$out"
inputs=("$geo"/GSM2371702_*.wig.gz "$geo"/GSM2371706_*.wig.gz "$geo"/GSM2371725_*.wig.gz)
# condition mark GSM
while read -r condition mark gsm; do
    python "$here/minoux_fold_enrichment.py" "$reference" "$out/${condition}_${mark}.fold.bw" "$geo"/${gsm}_*.wig.gz \
        --input "${inputs[@]}"
done <<'EOF'
Md H3K4me2 GSM2371710
PA2 H3K4me2 GSM2371711
Md H3K27me3 GSM2371715
PA2 H3K27me3 GSM2371716
Md H3K27ac GSM2371719
PA2 H3K27ac GSM2371720
Md reChIP GSM2371723
PA2 reChIP GSM2371724
EOF
