#!/bin/bash
set -euo pipefail

cd $TREEFLOW_LIB
git clone https://github.com/tothuhien/lsd-0.3beta.git
ln -s $TREEFLOW_LIB/lsd-0.3beta/bin/lsd_unix $TREEFLOW_LIB/bin/lsd