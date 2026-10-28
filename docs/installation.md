# Installation

EGA is a Python package (Python 3.9 or newer). The conda environment of the repository fixes the versions with which
the results of the manuscript were produced.

## With conda / mamba (recommended)

```bash
git clone https://github.com/joachimwolff/Epigenetic-Graph-Analyzer.git
cd Epigenetic-Graph-Analyzer
mamba env create -f environment.yml
conda activate ega
pip install .
```

`pip install -e .` instead of `pip install .` installs the package in development mode.

## With pip only

```bash
pip install "epigeneticgraphanalyzer[leiden,cool,eqtl,compare,test] @ git+https://github.com/joachimwolff/Epigenetic-Graph-Analyzer.git"
```

`pyBigWig` needs a C compiler and the curl and zlib headers when no binary wheel exists for the platform; the conda
route avoids this.

| Extra | Needed for |
|---|---|
| `leiden` | Leiden communities (`--communityMethod leiden`; Louvain works without it) |
| `cool` | `.cool` / `.mcool` input (`egaCreateGraph --format cool`) |
| `eqtl` | regression models of `egaEqtl` |
| `compare` | partition comparison of the community robustness |
| `test` | the test suite |

## Check the installation

```bash
egaCreateGraph --version
python -m pytest -q tests          # from the repository root, about 2 minutes
bash example/run_example.sh        # the example of the tutorial, about 30 seconds
```

The example writes its key results to `example/output/summary/`; they have to be identical to
`example/expected/`:

```bash
diff -r example/expected example/output/summary && echo identical
```

## External tools

EGA starts from called interactions. Calling them from reads needs HiCUP and CHiCAGO, which have their own conda
environments (`envs/*.yml` in the repository). PDF output of `egaReport` uses a headless Chromium if available and
falls back to matplotlib.
