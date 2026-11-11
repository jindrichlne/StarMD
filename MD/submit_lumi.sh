#!/bin/bash
# Submit all 6 replicas (2 states x 3) as separate single-GCD jobs on LUMI small-g.
set -euo pipefail
cd "$(dirname "$0")"
for s in state1_lumen-open state2_cytosol-open; do
  for r in rep1 rep2 rep3; do
    sbatch --chdir="$s/$r" --job-name="${s%%_*}_${r}" run_lumi.sbatch
  done
done
