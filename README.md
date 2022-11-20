# TreeFlow writeups

## TreeFlow paper

Christiaan Swanepoel, Mathieu Fourment, Xiang Ji, Hassan Nasif, Marc A Suchard, Frederick A Matsen IV, Alexei Drummond. ["TreeFlow: probabilistic programming and automatic differentiation for phylogenetics"](https://arxiv.org/abs/2211.05220). arXiv preprint arXiv:2211.05220 (2022).

### Data

Sequence alignments used in the examples in the paper (carnivores and H3N2 datasets) are in the `data` directory.

Tree topologies, starting values, TreeFlow model definition files, and BEAST XMLs are in the `supplementary-data` directory.

### Benchmarks

See the [`treeflow-benchmarks` repository](https://github.com/christiaanjs/treeflow-benchmarks)

## Installation
See `setup/setup.sh` for installation (out of date). The `TREEFLOW_LIB` environment variable must be set to specify the location where dependencies are installed.
This script assumes you have:
* conda
* gcc >= 7.5
* JDK 8+
* ant 

## Requirements
* [`treeflow`](https://github.com/christiaanjs/treeflow)
* Python requirements in `setup.py`
* BEAST 2
  * [Feast](https://github.com/tgvaughan/feast)
  * [Experimenter](https://github.com/christiaanjs/beast-validation)
  * [bwest](https://github.com/4ment/bwest)
* RAxML
* [LSD 0.2](http://www.atgc-montpellier.fr/LSD/)