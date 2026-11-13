#!/bin/bash
# Submit all 6 replicas (2 states x 3) as separate single-GCD jobs on LUMI small-g.
#   SBATCH_ACCOUNT=project_<number> EBU_USER_PREFIX=/project/project_<number>/EasyBuild ./submit_lumi.sh
# (GMX_MODULE overrides the GROMACS module name; see run_lumi.sbatch.)
set -euo pipefail
: "${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT to your LUMI project, e.g. project_<number>}"
export SBATCH_ACCOUNT
cd "$(dirname "$0")"
for s in state1_lumen-open state2_cytosol-open; do
  for r in rep1 rep2 rep3; do
    sbatch --chdir="$s/$r" --job-name="${s%%_*}_${r}" run_lumi.sbatch
  done
done
