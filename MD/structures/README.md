# SLC35A2 — MD starting structures

Two AlphaFold3 conformations of the human UDP-galactose translocator SLC35A2 (UniProt P78381). They were mined from 290 models and are the two alternating-access states. Simulate each with 3 replicas: **6 MD runs** in total.

| File | State | Open to | Narrowest point of the open side | Models in state |
|---|---|---|---|---|
| `SLC35A2_state1_lumen-open.pdb` | state1 | Golgi lumen (cytosol sealed) | 3.6 Å radius | 265 of 290 |
| `SLC35A2_state2_cytosol-open.pdb` | state2 | cytosol (lumen sealed) | 2.0 Å radius | 25 of 290 |

The two models differ by a C-alpha distance-matrix RMSD (dRMSD) of 0.28 nm. Each file is the most typical model (medoid) of its state, taken from the better-confidence half of that state.

## What is in the files

- **Residues 29–341**, the same construct in both files: 313 aa, chain A, heavy atoms only, no hydrogens, ligands or waters. The disordered tails (1–28 and 342–396, pLDDT mostly < 60) are trimmed. Untrimmed copies are in `full_length/`.
- **UniProt numbering, not renumbered.** Keep it through system building (no `-renum` in `pdb2gmx`), so the switch pairs below match the trajectory.
- **Oriented for membrane building:** membrane normal along z, bilayer centre at z = 0, **cytosol at −z, Golgi lumen at +z**. The normal and centre come from a hydrophobic-slab fit (±15 Å). As a check, 81–88% of Trp/Tyr sit in the 8–18 Å interface belts. The cytosolic side was set by the positive-inside rule (20 vs 6 Lys+Arg) and the N- and C-termini, which are cytosolic in nucleotide-sugar transporters. If you re-orient with PPM in CHARMM-GUI, keep the termini on the cytosolic side.
- **B-factor column = pLDDT** (AF3 confidence, 0–100).
- **Lower-confidence internal loops** (pLDDT < 70): 137–142 and the Gly-rich 183–195 in both files; also 224–226 and 289–302 in state2. Expect these to move early in MD.

## System to build

- **Membrane: POPC:POPE:CHOL 40:40:20.** Per-leaflet asymmetry is not specified.
- **Termini:** 29 and 341 are artificial truncation points. Cap them with neutral groups (ACE/NME; in CHARMM-GUI, the ACE and CT3 patches) rather than charged termini.
- **Charge:** residues 29–341 carry **+13** at pH 7 (Lys/Arg +1, Asp/Glu −1, His neutral). Neutralise and add 0.15 M NaCl.
- **Protonation:** charged residues near the membrane centre (within ±10 Å) are E75, K78, K148 and K297, plus R314 and D321 in state1. They likely line the substrate pocket. Worth a pKa check (e.g. PROPKA) before choosing protonation states.
- MD settings and analysis follow `project.md` steps 4–5. Replace the water box of step 4b with the bilayer above.

## What to monitor

`switch_pairs_to_monitor.csv` lists the 25 Cα–Cα distances that best separate the two states. For each it gives the mean and spread per state across the AF3 models. The top ones:

| Pair | Lumen-open (nm) | Cytosol-open (nm) |
|---|---|---|
| P61–D196 | 3.25 | 1.23 |
| T59–D196 | 2.98 | 1.01 |
| D196–F316 | 3.38 | 1.56 |
| L60–D196 | 3.25 | 1.28 |

The cytosol-open state closes its luminal side by folding the loops around N135 (TM4–5) and D196 (TM6–7) onto the TM2 end (residues 56–62) and the 315–317 region. A replica that drifts to the other state's value has switched state.

For example, in GROMACS:

```bash
gmx distance -s md.tpr -f md.xtc -select 'name CA and resid 59 196' -oall T59_D196.xvg
```

## Provenance

`provenance/` holds the mining log, the MD selection, the per-model open/closed table (`run1_gates.csv`) and the overview figures. Everything was produced by `run_example.sh` in the parent folder, using AlphaFold3 at commit 28c2d29, on 2026-11-10. To regenerate this folder:

```bash
cd .. && ../.venv/bin/python scripts/prepare_md_inputs.py mining/run1 -o md_input --name SLC35A2
```
