#!/bin/bash
# Build one SLC35A2 state: caps -> pdb2gmx (amber14sb, all His = HIE) -> packmol-memgen
# (POPC:POPE:CHL 40:40:20, TIP3P) -> GROMACS topology (Lipid21 port) -> cavity water +
# neutralise/0.15 M NaCl -> restrained EM -> structure checks -> 3 replica directories
# with mdp files and nvt.tpr.
#
# usage: scripts/build_state.sh <state>      e.g. state1_lumen-open
# Run from the MD/ directory. Steps whose outputs exist are skipped.
set -euo pipefail

STATE=$1
MD=$(cd "$(dirname "$0")/.." && pwd)
S=$MD/$STATE
B=$S/00_build
GMX=${GMX:-gmx}
NTMPI=${NTMPI:-7}       # this GROMACS build has no OpenMP: parallelise with thread-MPI ranks
case $STATE in
  state1_*) SID=1 ;;
  state2_*) SID=2 ;;
  *) echo "unknown state $STATE"; exit 1 ;;
esac
mkdir -p "$B"

# 1. ACE/NME caps on the artificial termini (res 29 / 341), UniProt numbering kept
[[ -f $B/protein_capped.pdb ]] || python3 "$MD/scripts/cap_termini.py" "$MD/structures/SLC35A2_${STATE}.pdb" "$B/protein_capped.pdb"

# 2. pdb2gmx: amber14sb, TIP3P, H rebuilt, standard charge states (PROPKA: no shifts at pH 7), HIE x4, +13 e
if [[ ! -f $B/protein.gro ]]; then
  (cd "$B" && printf "1\n1\n1\n1\n" | $GMX pdb2gmx -f protein_capped.pdb -o protein.gro -p protein.top \
      -i posre_pdb2gmx.itp -ff amber14sb -water tip3p -ignh -his > pdb2gmx.log 2>&1)
  (cd "$B" && $GMX editconf -f protein.gro -o protein_H.pdb > /dev/null 2>&1 \
      && sed -e 's/HISE/HIE /' protein_H.pdb | grep -E "^(ATOM|TER|END)" > protein_H_memgen.pdb)
fi

# 3. membrane + solvent
packed() {  # did the LAST packmol-memgen run in this directory finish? (its log is appended)
  [[ -f $B/packmol-memgen.log ]] || return 1
  awk '/Execution line/{ok=0} /Packing process finished/{ok=1} END{exit !ok}' "$B/packmol-memgen.log"
}
if ! packed; then
  "$MD/scripts/run_packmol_memgen.sh" "$S"
fi

# 4. GROMACS system + 5. cavity hydration and ions (neutralise + 0.15 M NaCl)
if [[ ! -f $S/topol.top ]]; then
  python3 "$MD/scripts/packmol_to_gromacs.py" "$S" ../lipid21.ff | tee "$B/convert.log"
  python3 "$MD/scripts/hydrate_and_ionize.py" "$S" 0.15 2026 | tee -a "$B/convert.log"
fi

# 6. energy minimisation (restrained protein heavy atoms + lipid head z; reference = packed input)
#    packmol leaves lipid tails < 1 A apart; a plain EM either overflows or forces cis->trans
#    flips / ring piercing. Ladder: soft-core lambda 0.9 (water/ions frozen) -> 0.99 -> normal.
mkdir -p "$S/01_min"
cp "$MD/mdp/min_sc1.mdp" "$MD/mdp/min_sc2.mdp" "$MD/mdp/min.mdp" "$S/01_min/"
if [[ ! -f $S/01_min/min.gro ]]; then
  prev=../system.gro
  for st in min_sc1 min_sc2 min; do
    (cd "$S/01_min" \
     && $GMX grompp -f $st.mdp -c $prev -r ../system.gro -p ../topol.top -n ../index.ndx \
          -o $st.tpr -po ${st}_out.mdp -maxwarn 0 > grompp_$st.log 2>&1 \
     && $GMX mdrun -deffnm $st -ntmpi "$NTMPI" > mdrun_$st.log 2>&1) \
      || { echo "$st failed"; exit 1; }
    grep -q "converged to Fmax" "$S/01_min/$st.log" || { echo "$st did not converge to Fmax < emtol"; exit 1; }
    prev=$st.gro
  done
fi
for f in min_sc1 min_sc2 min; do
  echo "$f: $(grep -E 'converged' "$S/01_min/$f.log" | tail -1)"
  grep -E "Potential Energy|Maximum force" "$S/01_min/$f.log" | tail -2
done

# 7. checks on the minimised structure
python3 "$MD/scripts/check_structure.py" "$S/topol.top" "$S/system.gro" "$S/01_min/min.gro" | tee "$S/01_min/check_structure.log"

# 8. replicas: own seeds, mdp copies, nvt.tpr ready to run
for r in 1 2 3; do
  R=$S/rep$r
  mkdir -p "$R"
  SEED=$((SID * 1000 + r))
  python3 "$MD/scripts/write_mdp.py" "$R" "$SEED"
  rm -f "$R"/min*.mdp
  cp "$MD/scripts/run_md.sh" "$R/"
  (cd "$R" && $GMX grompp -f nvt.mdp -c ../01_min/min.gro -r ../01_min/min.gro -p ../topol.top -n ../index.ndx \
      -o nvt.tpr -po nvt_out.mdp -maxwarn 0 > grompp_nvt.log 2>&1) || { echo "grompp nvt failed for rep$r"; exit 1; }
  echo "rep$r: seed $SEED, nvt.tpr ready"
done
