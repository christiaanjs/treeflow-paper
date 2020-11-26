#!/bin/bash
set -eo pipefail

# Let us use conda inside this shell
source $(which conda | xargs dirname)/../etc/profile.d/conda.sh

source setup/load-modules-nesi.sh

snakemake \
    -j 300 \
    --cluster-config config/nesi-config.yaml \
    --cluster "sbatch --mem={cluster.mem} --job-name={rule} --time={cluster.time} --output=/dev/null --error=/dev/null" \
    -s workflow/sim.smk