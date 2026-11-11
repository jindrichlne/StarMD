#!/bin/bash
# Run all 6 replicas (2 states x 3) one after another on this machine.
# Each replica: nvt -> npt1..npt4 -> 100 ns prod (see scripts/run_md.sh).
# On a cluster, submit one job per replica instead, e.g.
#   cd state1_lumen-open/rep1 && GMX=gmx_mpi MDRUN_FLAGS="-ntomp 16" ./run_md.sh
set -euo pipefail
cd "$(dirname "$0")"
for s in state1_lumen-open state2_cytosol-open; do
  for r in rep1 rep2 rep3; do
    echo "=== $s/$r"
    (cd "$s/$r" && ./run_md.sh)
  done
done
