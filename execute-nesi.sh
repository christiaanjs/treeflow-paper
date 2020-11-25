#!/bin/bash
set -euo pipefail

bash setup/load-modules-nesi.sh
snakemake -j 100 --cluster-config config/nesi-config.yaml --cluster "sbatch --mem={cluster.mem} --time={cluster.time}"