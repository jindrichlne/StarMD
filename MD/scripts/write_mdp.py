#!/usr/bin/env python3
"""Write the GROMACS .mdp files (min, nvt, npt1-4, prod) for the SLC35A2 membrane systems.

Force fields: amber14sb (protein), TIP3P, Lipid21 (lipids).
Cut-offs: Amber convention, 1.0 nm for both LJ (plain cut-off; Potential-shift only
offsets the energy, forces equal Amber's truncation) and PME real space. No DispCorr.

Equilibration (CHARMM-GUI-like stepwise release of restraints):
  stage  ens  dt     length   BB    SC    lipid-z (P31 of PC/PE, O1 of CHL)
  min_sc1 EM soft-core lambda 0.9, water/ions frozen     4000  2000  1000
  min_sc2 EM soft-core lambda 0.99                        4000  2000  1000
  min    EM   -      -        4000  2000  1000
  nvt    NVT  1 fs   250 ps   4000  2000  1000
  npt1   NPT  1 fs   250 ps   2000  1000   400
  npt2   NPT  2 fs     1 ns   1000   500   200
  npt3   NPT  2 fs     1 ns    500   200    40
  npt4   NPT  2 fs     5 ns     50     0     0
  prod   NPT  2 fs   100 ns      -     -     -
(force constants in kJ mol^-1 nm^-2)

usage: write_mdp.py <outdir> [gen_seed]   (gen_seed is used for nvt gen-seed and ld-seed)
"""
import os
import sys

TEMP = 310  # K

HEADER = """; {title}
; SLC35A2 in POPC:POPE:CHL 40:40:20 | amber14sb + TIP3P + Lipid21 | GROMACS 2026
"""

NONBONDED = """
; ---- Non-bonded: Amber ff14SB / Lipid21 conventions ------------------------
; 1.0 nm LJ cut-off (Potential-shift only shifts the energy; forces are those of
; Amber's plain truncation), PME with 1.0 nm real-space cut-off.
; NO long-range dispersion correction.
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
vdw-modifier             = Potential-shift
rvdw                     = 1.0
DispCorr                 = no
fourierspacing           = 0.12
pme-order                = 4
ewald-rtol               = 1e-5
"""

CONSTRAINTS = """
; ---- Constraints (Amber SHAKE-equivalent: bonds to H; water rigid via SETTLE)
constraints              = h-bonds
constraint-algorithm     = LINCS
lincs-order              = 4
lincs-iter               = 1
"""

TCOUPL = """
; ---- Temperature coupling --------------------------------------------------
tcoupl                   = V-rescale
tc-grps                  = SOLU MEMB SOLV
tau-t                    = 1.0  1.0  1.0
ref-t                    = {T}  {T}  {T}
""".format(T=TEMP)

COMM = """
; ---- Centre-of-mass motion removal: protein+bilayer and solvent separately
comm-mode                = Linear
comm-grps                = SOLU_MEMB SOLV
nstcomm                  = 100
"""

PCOUPL = """
; ---- Pressure coupling: semi-isotropic (tension-free bilayer) --------------
pcoupl                   = C-rescale
pcoupltype               = semiisotropic
tau-p                    = 5.0
ref-p                    = 1.0     1.0
compressibility          = 4.5e-5  4.5e-5
{refcoord}"""

NO_PCOUPL = """
; ---- No pressure coupling (NVT) --------------------------------------------
pcoupl                   = no
"""


def restraint_define(bb, sc, lip):
    if bb == 0 and sc == 0 and lip == 0:
        return ""
    d = f"-DPOSRES -DPOSRES_FC_BB={bb:.1f} -DPOSRES_FC_SC={sc:.1f}"
    if lip > 0:
        d += f" -DPOSRES_LIPID -DPOSRES_FC_LIPID={lip:.1f}"
    return d


def output_block(nxtc, nener, nlog, xtc_grp="System"):
    return f"""
; ---- Output (whole system is written to the compressed trajectory) --------
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = {nxtc}
compressed-x-grps        = {xtc_grp}
compressed-x-precision   = 1000
nstenergy                = {nener}
nstcalcenergy            = 100
nstlog                   = {nlog}
"""


def write(outdir, name, text):
    with open(os.path.join(outdir, name), "w") as f:
        f.write(text.lstrip("\n"))


def main(outdir, seed):
    os.makedirs(outdir, exist_ok=True)

    # ---------------- minimisation
    write(outdir, "min.mdp", HEADER.format(title="Energy minimisation (steepest descent), restrained protein heavy atoms + lipid head z") + f"""
define                   = {restraint_define(4000, 2000, 1000)}
integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = 50000
nstenergy                = 500
nstlog                   = 500
""" + NONBONDED + """
; ---- No constraints during EM (water stays rigid via SETTLE)
constraints              = none
""")

    # ---------------- soft-core pre-minimisation (packmol leaves overlapping lipid tails)
    # sc1: lambda 0.9 (very soft, overlaps resolve without cis->trans flips / ring piercing);
    #      water + ions frozen, otherwise softened Na+/Cl- pairs can collapse onto each other
    # sc2: lambda 0.99 (stiff soft-core wall, everything free), then the normal min.mdp
    SOFTCORE = """
; ---- all non-bonded interactions at lambda = {lam} with soft-core (finite at r = 0)
free-energy              = yes
couple-moltype           = system
couple-lambda0           = none
couple-lambda1           = vdw-q
couple-intramol          = yes
init-lambda-state        = 0
vdw-lambdas              = {lam}
coul-lambdas             = {lam}
sc-alpha                 = 0.5
sc-power                 = 1
sc-sigma                 = 0.3
sc-coul                  = yes
nstdhdl                  = 0
"""
    for name, lam, freeze, title in [
            ("min_sc1.mdp", 0.9, True, "Soft-core pre-minimisation 1 (lambda 0.9, water/ions frozen): removes packmol lipid-tail overlaps"),
            ("min_sc2.mdp", 0.99, False, "Soft-core pre-minimisation 2 (lambda 0.99, all atoms free)")]:
        write(outdir, name, HEADER.format(title=title) + f"""
define                   = {restraint_define(4000, 2000, 1000)}
integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = 10000
nstenergy                = 500
nstlog                   = 500
""" + NONBONDED + """
constraints              = none
""" + SOFTCORE.format(lam=lam) + ("""
; ---- water and ions frozen (only protein/lipids overlap after packing)
freezegrps               = SOLV
freezedim                = Y Y Y
""" if freeze else ""))

    # ---------------- NVT
    nvt_steps = 250000  # 250 ps at 1 fs
    write(outdir, "nvt.mdp", HEADER.format(title="NVT equilibration, 250 ps, dt = 1 fs, strong restraints, new velocities") + f"""
define                   = {restraint_define(4000, 2000, 1000)}
integrator               = md
dt                       = 0.001
nsteps                   = {nvt_steps}
continuation             = no
gen-vel                  = yes
gen-temp                 = {TEMP}
gen-seed                 = {seed}
ld-seed                  = {seed}
""" + output_block(10000, 5000, 5000) + NONBONDED + CONSTRAINTS + TCOUPL + COMM + NO_PCOUPL)

    # ---------------- NPT stages
    npt = [  # name, dt (ps), nsteps, BB, SC, LIPID, description
        ("npt1", 0.001, 250000, 2000, 1000, 400, "250 ps, dt = 1 fs"),
        ("npt2", 0.002, 500000, 1000, 500, 200, "1 ns, dt = 2 fs"),
        ("npt3", 0.002, 500000, 500, 200, 40, "1 ns, dt = 2 fs"),
        ("npt4", 0.002, 2500000, 50, 0, 0, "5 ns, dt = 2 fs, weak backbone restraints only"),
    ]
    for name, dt, nsteps, bb, sc, lip, desc in npt:
        write(outdir, f"{name}.mdp", HEADER.format(title=f"NPT equilibration {name}: {desc}; restraints BB {bb} / SC {sc} / lipid-z {lip} kJ/mol/nm^2") + f"""
define                   = {restraint_define(bb, sc, lip)}
integrator               = md
dt                       = {dt}
nsteps                   = {nsteps}
continuation             = yes
gen-vel                  = no
ld-seed                  = {seed}
""" + output_block(25000 if dt == 0.002 else 10000, 5000, 5000) + NONBONDED + CONSTRAINTS + TCOUPL + COMM
              + PCOUPL.format(refcoord="refcoord-scaling         = com\n"))

    # ---------------- production
    write(outdir, "prod.mdp", HEADER.format(title="Production, 100 ns, dt = 2 fs, NPT 310 K / 1 bar semi-isotropic, no restraints") + f"""
integrator               = md
dt                       = 0.002
nsteps                   = 50000000
continuation             = yes
gen-vel                  = no
ld-seed                  = {seed}
""" + output_block(10000, 5000, 5000) + NONBONDED + CONSTRAINTS + TCOUPL + COMM
          + PCOUPL.format(refcoord=""))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else -1)
