#!/bin/bash

module load Miniconda3/4.8.3
module load GCC/9.2.0
module load Java/1.8.0_144
module load RAxML/8.2.12-gimkl-2020a
module load ant/1.10.1-Java-1.8.0_144
export TREEFLOW_LIB=~/treeflow-lib
export PATH=$TREEFLOW_LIB/bin:$PATH
export DISPLAY=""

conda activate libsbn