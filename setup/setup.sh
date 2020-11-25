#!/bin/bash
set -euo pipefail

SCRIPT=`realpath $0`
TREEFLOW_PAPER_DIR=`dirname $SCRIPT`

cd $TREEFLOW_LIB
git clone https://github.com/christiaanjs/treeflow.git
cd treeflow
bash setup.sh

cd $TREEFLOW_PAPER_DIR
pip install -e .
