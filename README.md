# TreeFlow paper

Reproducibility repository for:

Christiaan Swanepoel, Mathieu Fourment, Xiang Ji, Hassan Nasif, Marc A Suchard, Frederick A Matsen IV, Alexei Drummond. ["TreeFlow: probabilistic programming and automatic differentiation for phylogenetics"](https://arxiv.org/abs/2211.05220). arXiv preprint arXiv:2211.05220 (2022).

## Repository structure

```
config/              Snakemake and model configuration files
data/                Sequence alignments (carnivores, H3N2)
manuscript/          LaTeX source and figures
scripts/             R plotting scripts
supplementary-data/  Pre-computed topologies, starting values, and BEAST XMLs
treeflow_pipeline/   Python package for pipeline logic and Snakemake rules
workflow/            Snakemake workflow files (data.smk, ms.smk, sim.smk)
```

## Requirements

### Python

- Python >= 3.9 (3.12 recommended)
- [`treeflow`](https://github.com/christiaanjs/treeflow) and its dependencies (TensorFlow, TensorFlow Probability)
- Python dependencies listed in `setup.py`

### External tools

- [BEAST 2](https://www.beast2.org/) >= 2.7, with the [Feast](https://github.com/tgvaughan/feast) package installed
- [RAxML](https://github.com/stamatak/standard-RAxML) (for topology inference)
- [LSD2](https://github.com/tothuhien/lsd2) (for tree rooting and dating)
- [BEAGLE](https://github.com/beagle-dev/beagle-lib) (optional, for accelerated BEAST runs)

### R

- R >= 4.0
- Install required CRAN and Bioconductor packages:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "purrr", "stringr",
  "ape", "reticulate", "gridExtra", "glue", "tibble",
  "broom", "magrittr"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("treeio")
```

- Install the local `treeflowbenchmarksr` package (from this workspace):

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_local("../treeflow-benchmarks/treeflowbenchmarksr")
```

### Benchmarks (optional)

To reproduce the benchmark comparison (Fig 7), you also need:
- [`treeflow-benchmarks`](https://github.com/christiaanjs/treeflow-benchmarks)
- [`phylojax`](https://github.com/christiaanjs/phylojax)
- [`bito`](https://github.com/phylovi/bito) (must be built from source; see its README for build instructions)
- The R package `treeflowbenchmarksr` (included in the treeflow-benchmarks repo)

### LaTeX

Install with `tlmgr install <PACKAGE>`
* `minted`
* `multirow`
* `lastpage`
* `standalone`

## Installation

Install this package and its dependencies:

```bash
pip install -e .                             # this package (treeflow-paper)
pip install -e /path/to/treeflow             # treeflow library
pip install -e /path/to/treeflow-benchmarks  # benchmarks (optional)
pip install arviz                            # additional dependency for results
```

Ensure `beast`, `raxmlHPC`, and `lsd` are available on your `PATH`.

## Reproducing results

The analysis is orchestrated with [Snakemake](https://snakemake.readthedocs.io/) (v9+). There are two main workflows:

### Data pipeline (`workflow/data.smk`)

Runs topology inference, BEAST MCMC, TreeFlow VI, and generates diagnostic plots for each dataset (carnivores, H3N2).

```bash
# Run everything (topology, BEAST, VI, plots)
snakemake -s workflow/data.smk --cores 1

# Run specific targets
snakemake -s workflow/data.smk out/carnivores/variational-samples.csv --cores 1
snakemake -s workflow/data.smk out/h3n2/beast.xml --cores 1
```

#### H3N2 multi-run VI fit

For the manuscript figure, the H3N2 VI fit is run 4 times (60,000 iterations each,
different seeds) rather than once, so a pooled posterior estimate with an inter-run
Monte Carlo error band can be shown -- the same treatment the carnivores base-model
figure gets from repeated runs of its example notebook. This is a separate target
from the single-run `out/h3n2/variational-samples.csv` above (which is still produced
and used for the quick local `marginals.png`/`traces.png`/`timing-data.csv`
diagnostics):

```bash
snakemake -s workflow/data.smk \
    out/h3n2/variational-multi-run/samples.csv \
    out/h3n2/variational-multi-run/tree-samples.nexus \
    out/h3n2/variational-multi-run/timing.csv \
    --cores 1
```

This produces, per seed, `variational-multi-run/{trace,samples,tree-samples}-run{1..4}.pickle|csv|nexus`
plus a `variational-benchmark.txt`-style timing file per run, and pools them into:
- `variational-multi-run/samples.csv` -- all runs' parameter samples, annotated with a `run` column
- `variational-multi-run/tree-samples.nexus` -- all runs' tree samples, concatenated
- `variational-multi-run/timing.csv` -- per-run seed/iteration-count/wall-clock-time

`workflow/ms.smk`'s `flu_marginals_plot`/`flu_tree_plot` rules consume these directly
(taking precedence over the generic single-run `data_marginals_plot`/`data_tree_plot`
rules for the H3N2 dataset).

### Manuscript pipeline (`workflow/ms.smk`)

Generates publication figures and compiles the manuscript. Depends on outputs from the data pipeline and the treeflow-benchmarks pipeline.

```bash
snakemake -s workflow/ms.smk --cores 1
```

To generate individual figures:

```bash
snakemake -s workflow/ms.smk manuscript/figures/carnivores-marginals.png --cores 1
snakemake -s workflow/ms.smk manuscript/figures/benchmark-log-scale-plot.png --cores 1
```

### BEAST runs

The BEAST MCMC analyses can be run directly using the provided BEAST 2.7 XMLs:

```bash
beast -seed 123 supplementary-data/carnivores/beast-2.7.xml
beast -seed 123 supplementary-data/h3n2/beast-2.7.xml
```

These write output to `out/{dataset}/beast.log` and `out/{dataset}/beast.trees`. Each runs a 30,000,000-iteration chain and may take several hours. Adding `-beagle` uses the BEAGLE library for acceleration if installed.

Alternatively, the data pipeline can generate BEAST XMLs from the model configuration and run them automatically via the `beast_run` rule.

### Benchmark pipeline

The benchmark comparison requires a separate pipeline in the `treeflow-benchmarks` repository:

```bash
cd /path/to/treeflow-benchmarks
snakemake --cores 1
```

This simulates phylogenetic data at multiple taxon counts and benchmarks likelihood computation across methods (TreeFlow, BEAGLE/bito, JAX). The output `plot-data.csv` is consumed by the manuscript pipeline for Fig 7.

### Carnivores example (Fig 6)

The branch-specific kappa plot (Fig 6a) requires running the `carnivores.ipynb` notebook in the `treeflow` repository's `examples/` directory, which produces `demo-out/carnivores-alt-trees.nexus`.

## Configuration

- `config/data-config.yaml` — seed, output directory, data URLs
- `config/data-models.yaml` — phylogenetic model specifications for each dataset
- `config/beast-config.yaml` — BEAST chain length and logging frequency
- `config/ms-config.yaml` — manuscript template and figure options

## Supplementary data

Pre-computed intermediate results are provided in `supplementary-data/` for each dataset:

- `beast-2.7.xml` — BEAST 2.7 analysis XML (ready to run)
- `topology.nwk` — inferred tree topology (Newick)
- `model.yaml` — model specification
