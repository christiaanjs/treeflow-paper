#!/bin/bash
set -eo pipefail

conda init bash
source ~/.bashrc

mkdir -p ~/treeflow-lib/bin
export TREEFLOW_LIB=~/treeflow-lib

SCRIPT=`realpath $0`
SETUP_DIR=`dirname $SCRIPT`

bash $SETUP_DIR/load-modules-nesi.sh
bash $SETUP_DIR/setup.sh
bash $SETUP_DIR/setup-beast.sh
