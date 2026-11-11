#!/bin/bash
# Equilibration + 100 ns production for ONE replica. Run it from the replica directory
# (e.g. state1_lumen-open/rep1). Restartable: finished stages (<stage>.gro present) are
# skipped, an interrupted stage continues from its checkpoint.
#
#   stage  ens  dt    length  restraints BB/SC/lipid-z (kJ/mol/nm^2)
#   nvt    NVT  1 fs  250 ps  4000/2000/1000   (new velocities, replica seed)
#   npt1   NPT  1 fs  250 ps  2000/1000/400
#   npt2   NPT  2 fs    1 ns  1000/500/200
#   npt3   NPT  2 fs    1 ns   500/200/40
#   npt4   NPT  2 fs    5 ns    50/0/0
#   prod   NPT  2 fs  100 ns   none
#
# Environment:
#   GMX          gmx binary (default: gmx; e.g. "gmx_mpi" or "srun gmx_mpi")
#   MDRUN_FLAGS  extra mdrun flags, e.g. "-ntmpi 14" (local GROMACS 2026.3 has no OpenMP,
#                use thread-MPI ranks), "-ntomp 16 -pin on" on a cluster, or on a GPU node
#                "-ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu"
#   STAGES       subset/order of stages to run (default: all)
set -euo pipefail

GMX=${GMX:-gmx}
MDRUN_FLAGS=${MDRUN_FLAGS:-}
STAGES=${STAGES:-"nvt npt1 npt2 npt3 npt4 prod"}
TOP=../topol.top
NDX=../index.ndx
REF=../01_min/min.gro          # restraint reference = minimised structure

prev_gro=../01_min/min.gro
prev_cpt=""
for stage in nvt npt1 npt2 npt3 npt4 prod; do
  if [[ -f $stage.gro ]]; then         # stage already finished
    prev_gro=$stage.gro; prev_cpt=$stage.cpt; continue
  fi
  [[ " $STAGES " == *" $stage "* ]] || { echo "stopping before $stage (not in STAGES)"; exit 0; }
  if [[ ! -f $stage.tpr ]]; then
    tflag=()
    [[ -n $prev_cpt ]] && tflag=(-t "$prev_cpt")   # carry velocities + coupling state
    $GMX grompp -f $stage.mdp -c "$prev_gro" -r $REF ${tflag[@]+"${tflag[@]}"} -p $TOP -n $NDX \
         -o $stage.tpr -po ${stage}_out.mdp -maxwarn 0
  fi
  # shellcheck disable=SC2086
  $GMX mdrun -deffnm $stage -cpi $stage.cpt $MDRUN_FLAGS
  [[ -f $stage.gro ]] || { echo "$stage did not finish"; exit 1; }
  prev_gro=$stage.gro; prev_cpt=$stage.cpt
done
echo "replica $(basename "$PWD") done"
