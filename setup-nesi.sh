#!/bin/bash
set -euo pipefail

module load Miniconda3/4.8.3
module load GCC/9.2.0

mkdir -p ~/treeflow-lib
export TREEFLOW_LIB=~/treeflow-lib

SCRIPT=`realpath $0`
SETUP_DIR=`dirname $SCRIPT`

bash $SETUP_DIR/setup.sh
