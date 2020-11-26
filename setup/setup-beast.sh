#!/bin/bash
set -eo pipefail

cd $TREEFLOW_LIB
wget https://github.com/CompEvol/beast2/releases/download/v2.6.3/BEAST.v2.6.3.Linux.tgz
tar -xvzf BEAST.v2.6.3.Linux.tgz
rm BEAST.v2.6.3.Linux.tgz
ln -s $TREEFLOW_LIB/beast/bin/* $TREEFLOW_LIB/bin
packagemanager -add feast
packagemanager -add BEASTLabs
packagemanager -add MASTER

# Build beast-validation
git clone https://github.com/CompEvol/beast2.git
cd beast2
ant dist_all_BEAST
cd ..
git clone https://github.com/BEAST2-Dev/BEASTLabs.git
cd BEASTLabs
ant dist_all_BEASTlabs
cd ..
git clone https://github.com/tgvaughan/feast.git
ant build-jar
cd ..
git clone https://github.com/tgvaughan/MASTER.git
cd MASTER
ant build-package
git clone https://github.com/christiaanjs/beast-validation.git
ant addon
# Add beast-validation package
EXPERIMENTER_DIR=~/.beast/2.6/Experimenter
EXPERIMENTER_ZIP=Experimenter.addon.v0.0.1.zip
mkdir -p $EXPERIMENTER_DIR
cp build/dist/$EXPERIMENTER_ZIP $EXPERIMENTER_DIR
cd $EXPERIMENTER_DIR
unzip $EXPERIMENTER_ZIP
rm $EXPERIMENTER_ZIP