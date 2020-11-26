#!/bin/bash
set -eo pipefail

# Let us use conda inside this shell
source $(which conda | xargs dirname)/../etc/profile.d/conda.sh

source setup/load-modules-nesi.sh
LOG_FILE="/nesi/nobackup/uoa02746/log/slurm-%j.out"
snakemake \
    -j 300 \
    --cluster-config config/nesi-config.yaml \
    --cluster "sbatch --mem={cluster.mem} --job-name={rule} --time={cluster.time} --output=$LOG_FILE --error=$LOG_FILE" \
    -s workflow/sim.smk