#!/bin/bash
# Pack POPC:POPE:CHL1 40:40:20 (both leaflets) + TIP3P + 0.15 M NaCl around a pre-oriented,
# already protonated (pdb2gmx) protein. Membrane normal = z, bilayer centre at z = 0 (README).
# --pbc: packmol packs against the periodic images, so the box is seamless.
# usage: run_packmol_memgen.sh <state_dir>
set -euo pipefail
cd "$1/00_build"
packmol-memgen --pdb protein_H_memgen.pdb --preoriented --notprotonate \
  --lipids POPC:POPE:CHL1 --ratio 4:4:2 \
  --dist 15 --dist_wat 17.5 \
  --salt --salt_c Na+ --salt_a Cl- --saltcon 0.15 \
  --pbc --nloop 30 --noprogress --overwrite \
  --output bilayer_protein.pdb --log packmol-memgen.log > packmol-memgen.stdout 2>&1
