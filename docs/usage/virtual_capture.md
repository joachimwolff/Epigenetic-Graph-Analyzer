# Virtual promoter capture from Hi-C or Micro-C

Without promoter capture Hi-C, a Hi-C or Micro-C map can serve as the input of EGA. Two routes exist.

## Called interactions (recommended)

The bins of one resolution act as restriction fragments, and the bins holding a transcription start site are the
baits. `analysis/virtual_capture.py` writes the CHiCAGO design and input from an `.mcool`; CHiCAGO then calls the
interactions as for capture data (environment `envs/capchic.yml`):

```bash
python analysis/virtual_capture.py hic.mcool::resolutions/5000 genes.gtf.gz virtual/ GM12878 protein_coding
# design files: makeDesignFiles_py3.py of chicagoTools on virtual/GM12878.rmap and .baitmap
# CHiCAGO: runChicago.R with the design and virtual/GM12878.chinput, output ibed
egaCreateGraph -i virtual/GM12878.ibed -f ibed --scoreThreshold 5 --tss genes.gtf.gz -o graphs/VC_GM12878
```

## Pixels of the matrix

`--format cool` reads the pixels of a `.cool` or `.mcool` that touch a baited bin directly, without calling:

```bash
egaCreateGraph -i hic.mcool::resolutions/5000 -f cool --baits promoter_bins.bed --coolField count \
    --tss genes.gtf.gz -o graphs/GM12878_pixels
```

With `--scoreThreshold` on a score matrix (for example a CHiCAGO score matrix stored as cool), only pixels at or above
the threshold are kept.

## How well it works

In GM12878, 19.3% of the virtual calls were among the capture calls, with a recall of 12.6% of the capture calls, and
the two-step support in SPRITE was 0.517 against 0.55 for the capture graph. Virtual capture is a usable input where no
capture data exist, with fewer and weaker indirect contacts.
