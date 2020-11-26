#!/bin/bash
set -euo pipefail

bash setup/load-modules-nesi.sh
snakemake -j 100 --cluster-config config/nesi-config.yaml --cluster "sbatch --mem={cluster.mem} --job-name={rule} --time={cluster.time} --output=/dev/null --error=/dev/null" -s workflow/sim.smk