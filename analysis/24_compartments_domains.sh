#!/usr/bin/env bash
# Can EGA estimate A/B compartments and domains? Evaluation against Hi-C (after 23_hic_reference.sh).
#   hic_eigenvectors.py        Hi-C eigenvectors phased by GC (needs cooler: HICEXPLORER_PYTHON, default python)
#   compartments_chromatin.py  compartments from chromatin tracks (PC1, role score, baselines) vs Hi-C eigenvectors
#   domains_graph.py           community / hub / edge / cascade-explanation statistics against TADs
#   compartment_switch.py      condition-specific compartment switches and expression concordance (prototype)
#   figure_compartments.py     one PDF per panel
# Results: <results>/compartments_domains/{tables,switch,figures}. Inputs and paths: compartments_common.py
# (EGA_RAW, EGA_RESULTS, EGA_HOPS override the defaults).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
shifts=${1:-1000}
"${HICEXPLORER_PYTHON:-python}" "$here/hic_eigenvectors.py" 100000 250000
python "$here/compartments_chromatin.py" 100000
python "$here/domains_graph.py" "$shifts" 50000
python "$here/compartment_switch.py" 0.5 0.05
python "$here/figure_compartments.py"
