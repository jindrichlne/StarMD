# SLC35A2 MD — GROMACS 2026 inputs

Two AlphaFold3 states of SLC35A2 (res 29–341, see `structures/README.md`) in a
POPC:POPE:CHL 40:40:20 bilayer. **3 replicas × 100 ns per state, 6 runs in total.**

| | |
|---|---|
| Protein | amber14sb.ff (ff14SB, built into GROMACS 2026), ACE/NME caps, UniProt numbering kept |
| Lipids | AMBER Lipid21, local GROMACS port `lipid21.ff/` (added through `lipid21_params.itp`) |
| Water / ions | TIP3P, Joung–Cheatham Na⁺/Cl⁻ (amber14sb.ff), neutralised + 0.15 M NaCl |
| Non-bonded | rvdw = rcoulomb = 1.0 nm (Amber convention), plain LJ cut-off (Potential-shift, forces = Amber truncation), PME, **DispCorr = no** |
| Ensemble | 310 K V-rescale (τ = 1 ps; SOLU / MEMB / SOLV), 1 bar C-rescale **semi-isotropic** (τ = 5 ps) |
| Integration | dt = 2 fs, h-bond constraints (LINCS), rigid TIP3P (SETTLE) |
| Output | whole system (`compressed-x-grps = System`) every 20 ps in prod (5000 frames / 100 ns) |

| System | box (nm) | atoms | POPC / POPE / CHL (lumen + cytosol leaflet) | TIP3P (cavity) | Na⁺ / Cl⁻ |
|---|---|---|---|---|---|
| state1 lumen-open | 10.546 × 10.546 × 9.472 | 106 650 | 159 / 159 / 79 (80+79, 80+79, 40+39) | 18 213 (201) | 49 / 62 |
| state2 cytosol-open | 10.436 × 10.436 × 9.225 | 101 548 | 155 / 155 / 77 (78+77, 78+77, 39+38) | 16 909 (166) | 46 / 59 |

Protein charge +13 → 0.15 M NaCl from the water count plus 13 extra Cl⁻.

## Layout

```
MD/
├── structures/              input AF3 models (+ README)
├── lipid21.ff/              Lipid21 GROMACS port (provided)
├── mdp/                     template mdp files (min, nvt, npt1-4, prod)
├── scripts/                 the whole build pipeline (see below)
├── run_all.sh               runs all 6 replicas one after another
└── state{1_lumen-open,2_cytosol-open}/
    ├── system.gro           packed, hydrated, ionised start structure
    ├── topol.top            topology (includes ../lipid21.ff and amber14sb.ff)
    ├── protein.itp, posre_protein.itp, index.ndx
    ├── 00_build/            intermediates (capped PDB, pdb2gmx, packmol-memgen logs)
    ├── 01_min/              energy minimisation ladder (DONE), min.gro, check_structure.log
    └── rep{1,2,3}/          nvt.mdp npt{1..4}.mdp prod.mdp, nvt.tpr (ready), run_md.sh
```

Keep the `MD/` tree together when copying to a cluster: `topol.top` includes `../lipid21.ff/`.

## Running

Per replica (restartable; finished stages are skipped, interrupted ones continue from the checkpoint):

```bash
cd state1_lumen-open/rep1
MDRUN_FLAGS="-ntmpi 14" ./run_md.sh                          # this Mac (no OpenMP build)
GMX=gmx_mpi MDRUN_FLAGS="-ntomp 16 -pin on" ./run_md.sh      # cluster CPU
MDRUN_FLAGS="-ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu" ./run_md.sh   # GPU node
```

Equilibration/production schedule (restraints in kJ mol⁻¹ nm⁻²; BB = backbone N/CA/C/O,
SC = other protein heavy atoms, lipid = z-only on P31 of POPC/POPE and O1 of CHL;
reference = minimised structure, `refcoord-scaling = com`):

| stage | ensemble | dt | length | BB | SC | lipid-z |
|---|---|---|---|---|---|---|
| min_sc1 (done) | steepest descent, soft-core λ = 0.9, water/ions frozen | – | Fmax < 1000 | 4000 | 2000 | 1000 |
| min_sc2 (done) | steepest descent, soft-core λ = 0.99 | – | Fmax < 1000 | 4000 | 2000 | 1000 |
| min (done) | steepest descent, full potential | – | Fmax < 1000 | 4000 | 2000 | 1000 |
| nvt | NVT, new velocities | 1 fs | 250 ps | 4000 | 2000 | 1000 |
| npt1 | NPT | 1 fs | 250 ps | 2000 | 1000 | 400 |
| npt2 | NPT | 2 fs | 1 ns | 1000 | 500 | 200 |
| npt3 | NPT | 2 fs | 1 ns | 500 | 200 | 40 |
| npt4 | NPT | 2 fs | 5 ns | 50 | 0 | 0 |
| prod | NPT | 2 fs | 100 ns | – | – | – |

Replicas share the minimised structure and differ in their velocity / thermostat seeds
(`gen-seed` = `ld-seed` = 1001, 1002, 1003 for state1 and 2001, 2002, 2003 for state2).

## How the systems were built (`scripts/build_state.sh <state>`)

1. `cap_termini.py` – ACE/NME heavy atoms on the truncation points (res 28 / 342).
2. `pdb2gmx` – amber14sb, TIP3P, H rebuilt. PROPKA3 (pH 7): no residue shifts its charge
   state (E75 4.8, D321 2.6–3.5, K148 7.5–8.2, K297 8.6–8.7, K78 ≥ 10.3, R314 12.4), so all
   Asp/Glu/Lys/Arg are charged in both states; His32/96/102/319 neutral (pKa < 6.4) and
   HIE in both states (pdb2gmx H-bond analysis also gives HIE). No disulfides. Net charge +13.
3. `run_packmol_memgen.sh` – packmol-memgen (AmberTools26), protein pre-oriented
   (membrane normal z, centre z = 0), POPC:POPE:CHL1 4:4:2 in both leaflets,
   ≥ 15 Å protein–box edge in xy, 17.5 Å water beyond the protein in z, packed with
   `--pbc` (packmol packs against the periodic images; the GROMACS box is the packmol box).
4. `packmol_to_gromacs.py` – maps the Amber-named packmol output onto the pdb2gmx protein
   and the Lipid21 .itp atom order, writes `system.gro`, `topol.top`, posres, index.
5. `hydrate_and_ionize.py` – packmol-memgen leaves the transporter cavity inside the bilayer
   empty; buried cavity water is added from `gmx solvate` (ray-casting burial test, ≥ 3 Å
   from lipids, so lipid-facing grooves stay dry). Ions are then set to "neutralise + 0.15 M
   NaCl" from the final water count.
6. Restrained minimisation ladder. packmol does not converge at this lipid density
   (~550 lipid–lipid atom pairs < 1 Å). A plain EM either overflows (state1, two H atoms
   0.1 Å apart) or resolves the overlaps by force (state2: an oleoyl bond flipped to trans,
   a POPE ester bond threaded through a cholesterol ring). Soft-core EM at λ = 0.9 (water and
   ions frozen, otherwise softened Na⁺/Cl⁻ pairs collapse) → λ = 0.99 → normal EM gives clean
   structures. Then `check_structure.py`: ring piercing (CHL rings, aromatic/Pro rings),
   chirality vs. the input, oleoyl cis double bond — all clean in both states; no
   intermolecular contact < 1.2 Å; protein Cα RMSD to the AF3 model 0.09 Å.
7. `write_mdp.py` – mdp files per replica with replica seeds; `nvt.tpr` pre-built.

## Checks done

- Lipid21 port vs AmberTools (`tleap` + `sander`, single POPC / POPE / CHL): bond, angle,
  dihedral, LJ-14 identical; Coulomb-14 within 0.003 % (Amber vs GROMACS Coulomb constant).
- amber14sb + `lipid21_params.itp` vs the standalone `lipid21.ff/forcefield.itp` on the whole
  lipid patch: all energy terms identical (no atom-type cross-talk).
- Smoke test of the whole chain (nvt → npt1–4 → prod, shortened) for both states: all grompp
  with `-maxwarn 0` (only the note on COM removal with position restraints), no LINCS/SETTLE
  warnings, T = 310 K, density ≈ 986 kg m⁻³.
- Speed on this Mac (M4 Pro, CPU only, `-ntmpi 7`, two runs side by side): ≈ 8 ns/day per
  replica at 2 fs, i.e. ~12–13 days per 100 ns replica. Use GPU/HPC for production.

## Notes

- `mdrun` writes residue numbers from the topology: the protein keeps UniProt numbering
  (ACE 28 … NME 342, so `resid 59 196` selections in `structures/README.md` work), lipid
  residues are numbered 1–3 within each lipid molecule (select lipids by `resname`).
- The membrane is symmetric (README leaves leaflet asymmetry unspecified).
