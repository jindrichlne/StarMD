#!/usr/bin/env python3
"""Convert a packmol-memgen system (Amber naming) into GROMACS input for
amber14sb.ff (protein, from pdb2gmx) + Lipid21 port (lipid21.ff) + TIP3P/ions.

Writes into <state_dir>:
  system.gro           coordinates, ordered Protein | POPC | POPE | CHL | SOL | NA | CL,
                       UniProt residue numbers kept for the protein (ACE 28 ... NME 342)
  protein.itp          protein moleculetype extracted from the pdb2gmx topology
  posre_protein.itp    protein heavy-atom restraints (backbone / side chain force constants as macros)
  topol.top            full topology
  index.ndx            SOLU, MEMB, SOLV, SOLU_MEMB (+ System, Protein, Backbone, C-alpha)

usage: packmol_to_gromacs.py <state_dir> <lipid21.ff dir relative to state_dir>
       (expects <state_dir>/00_build/{bilayer_protein.pdb, protein.gro, protein.top})
"""
import os
import re
import sys
from collections import Counter

import numpy as np

LIPIDS = {  # molecule -> (residues in itp order, head residue name)
    "POPC": ["PA", "PC", "OL"],
    "POPE": ["PA", "PE", "OL"],
    "CHL": ["CHL"],
}
LIPID_RESTRAINT_ATOM = {"POPC": ("PC", "P31"), "POPE": ("PE", "P31"), "CHL": ("CHL", "O1")}
BACKBONE = {"N", "CA", "C", "O"}


# ------------------------------------------------------------------ readers
def read_pdb(fn):
    """Return list of residues: dict(resn, resi, atoms=[(name, xyz)]) and the CRYST1 box (A)."""
    residues, box, cur, key = [], None, None, None
    for line in open(fn):
        if line.startswith("CRYST1"):
            box = np.array([float(line[6:15]), float(line[15:24]), float(line[24:33])])
        elif line.startswith(("ATOM", "HETATM")):
            name, resn = line[12:16].strip(), line[17:21].strip()
            resi = line[22:27]
            k = (resn, resi, line[21])
            if k != key or cur is None:
                cur = dict(resn=resn, resi=resi, atoms=[])
                residues.append(cur)
                key = k
            cur["atoms"].append((name, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
        elif line.startswith("TER"):
            key = None  # force a new residue after TER
    return residues, box


def read_gro(fn):
    lines = open(fn).read().splitlines()
    n = int(lines[1])
    atoms = []
    for l in lines[2:2 + n]:
        atoms.append(dict(resi=int(l[0:5]), resn=l[5:10].strip(), name=l[10:15].strip(),
                          xyz=np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])])))
    return atoms


def itp_atoms(fn):
    """[(resnr, resname, atomname)] of the first [ atoms ] section."""
    out, sec = [], None
    for line in open(fn):
        s = line.split(";")[0].strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            if sec != "atoms" and out:
                break
            continue
        if sec == "atoms" and s:
            f = s.split()
            out.append((int(f[2]), f[3], f[4]))
    return out


# ------------------------------------------------------------------ main
def main(state_dir, ffdir_rel):
    build = os.path.join(state_dir, "00_build")
    residues, box = read_pdb(os.path.join(build, "bilayer_protein.pdb"))
    centre = None
    inp = open(os.path.join(build, "packmol.inp")).read()
    six = r"\s+".join([r"(\S+)"] * 6)
    pbc = re.search(r"^\s*pbc\s+" + six, inp, re.M)
    if pbc:
        # packed with packmol-memgen --pbc: the periodic box is exactly the packmol pbc box
        v = np.array([float(x) for x in pbc.groups()])
        lo, hi = v[:3], v[3:]
        box, centre = hi - lo, (hi + lo) / 2.0
        print(f"periodic box from packmol 'pbc' {lo} .. {hi} A")
    elif box is None:
        # no PBC packing and no CRYST1: box = union of the 'inside box' constraints + tolerance
        tol = float(re.search(r"^tolerance\s+(\S+)", inp, re.M).group(1))
        cons = np.array([[float(x) for x in m] for m in re.findall(r"inside box\s+" + six, inp)])
        lo, hi = cons[:, :3].min(0), cons[:, 3:].max(0)
        box, centre = hi - lo + tol, (hi + lo) / 2.0
        print(f"no CRYST1: box from packmol constraints {lo} .. {hi} + tolerance {tol} A")
    prot_gro = read_gro(os.path.join(build, "protein.gro"))
    nprot = len(prot_gro)

    # ---- split packmol residues into protein / lipids / water / ions
    flat = [(r, a) for r in residues for a in r["atoms"]]
    prot_pm = flat[:nprot]
    alias = {("ILE", "CD1"): "CD"}  # packmol-memgen renames GROMACS Ile CD to Amber CD1
    for (r, (name, _)), g in zip(prot_pm, prot_gro):
        if alias.get((r["resn"], name), name) != g["name"]:
            sys.exit(f"protein atom mismatch: packmol {r['resn']}{r['resi']}:{name} vs gro {g['resn']}{g['resi']}:{g['name']}")
    prot_xyz = np.array([a[1] for _, a in prot_pm])
    rest = residues
    # residues after the protein
    count = 0
    for i, r in enumerate(residues):
        count += len(r["atoms"])
        if count == nprot:
            rest = residues[i + 1:]
            break
    else:
        sys.exit("protein atom count does not end on a residue boundary")

    mols = {k: [] for k in ["POPC", "POPE", "CHL", "SOL", "NA", "CL"]}
    i = 0
    while i < len(rest):
        r = rest[i]
        rn = r["resn"]
        if rn == "PA":
            head, tail = rest[i + 1], rest[i + 2]
            mol = {"PC": "POPC", "PE": "POPE"}[head["resn"]]
            assert tail["resn"] == "OL", (head["resn"], tail["resn"])
            mols[mol].append([r, head, tail])
            i += 3
        elif rn == "CHL":
            mols["CHL"].append([r])
            i += 1
        elif rn in ("WAT", "HOH", "SOL"):
            d = dict(r["atoms"])
            mols["SOL"].append([d["O"], d["H1"], d["H2"]])
            i += 1
        elif rn == "Na+":
            mols["NA"].append([r["atoms"][0][1]])
            i += 1
        elif rn == "Cl-":
            mols["CL"].append([r["atoms"][0][1]])
            i += 1
        else:
            sys.exit(f"unknown residue {rn}")

    # ---- build ordered coordinate list (Angstrom), with gro naming
    ff = os.path.normpath(os.path.join(state_dir, ffdir_rel))
    out = []  # (resi, resn, name, xyz)
    for g, xyz in zip(prot_gro, prot_xyz):
        out.append((g["resi"], g["resn"], g["name"], xyz))
    resi = prot_gro[-1]["resi"]
    for mol in ["POPC", "POPE", "CHL"]:
        tmpl = itp_atoms(os.path.join(ff, f"{mol}.itp"))
        for m in mols[mol]:
            lookup = {}
            for res in m:
                for name, xyz in res["atoms"]:
                    lookup[(res["resn"], name)] = xyz
            if len(lookup) != len(tmpl):
                sys.exit(f"{mol}: {len(lookup)} atoms in packmol vs {len(tmpl)} in itp")
            base = resi
            for rnr, rname, aname in tmpl:
                out.append((base + rnr, rname, aname, lookup[(rname, aname)]))
            resi = base + max(t[0] for t in tmpl)
    for xyzs in mols["SOL"]:
        resi += 1
        for name, xyz in zip(["OW", "HW1", "HW2"], xyzs):
            out.append((resi, "SOL", name, xyz))
    for ion, nm in [("NA", "NA"), ("CL", "CL")]:
        for xyzs in mols[ion]:
            resi += 1
            out.append((resi, ion, nm, xyzs[0]))

    # ---- shift into the GROMACS box [0, L] (packmol-memgen box is centred on the origin in xy)
    allxyz = np.array([o[3] for o in out])
    if centre is None:
        centre = (allxyz.min(0) + allxyz.max(0)) / 2.0
    shift = box / 2.0 - centre
    with open(os.path.join(state_dir, "system.gro"), "w") as f:
        f.write("SLC35A2 in POPC:POPE:CHL 40:40:20, TIP3P, 0.15 M NaCl (packmol-memgen -> GROMACS)\n")
        f.write(f"{len(out)}\n")
        for n, (ri, rn, an, xyz) in enumerate(out, 1):
            x = (xyz + shift) / 10.0
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (ri % 100000, rn[:5], an[:5], n % 100000, *x))
        f.write("%10.5f%10.5f%10.5f\n" % tuple(box / 10.0))

    # ---- protein.itp from the pdb2gmx topology (moleculetype ... up to the posre include)
    top = open(os.path.join(build, "protein.top")).read()
    start = top.index("[ moleculetype ]")
    end = top.index("; Include Position restraint file")
    mt = top[start:end].replace("Protein_chain_A", "Protein        ", 1)
    with open(os.path.join(state_dir, "protein.itp"), "w") as f:
        f.write("; SLC35A2 residues 29-341 (UniProt numbering), ACE/NME caps, amber14sb (pdb2gmx, all His = HIE)\n")
        f.write(mt.rstrip() + "\n")
    patoms = itp_atoms(os.path.join(state_dir, "protein.itp"))
    with open(os.path.join(state_dir, "posre_protein.itp"), "w") as f:
        f.write("; Protein heavy-atom position restraints; force constants set through mdp 'define':\n")
        f.write(";   -DPOSRES -DPOSRES_FC_BB=<kJ/mol/nm2> -DPOSRES_FC_SC=<kJ/mol/nm2>\n")
        f.write("; backbone = N CA C O of amino acids; side chain = all other heavy atoms (incl. ACE/NME caps)\n")
        f.write("[ position_restraints ]\n;  ai  funct  fcx           fcy           fcz\n")
        for n, (rnr, rn, an) in enumerate(patoms, 1):
            if an.startswith("H"):
                continue
            fc = "POSRES_FC_BB" if (an in BACKBONE and rn not in ("ACE", "NME")) else "POSRES_FC_SC"
            f.write(f"{n:6d}    1  {fc:12s}  {fc:12s}  {fc:12s}\n")

    # ---- topology
    counts = [("Protein", 1)] + [(m, len(mols[m])) for m in ["POPC", "POPE", "CHL", "SOL", "NA", "CL"] if mols[m]]
    lipid_blocks = ""
    for mol in ["POPC", "POPE", "CHL"]:
        rn, an = LIPID_RESTRAINT_ATOM[mol]
        idx = next(n for n, t in enumerate(itp_atoms(os.path.join(ff, f"{mol}.itp")), 1) if t[1] == rn and t[2] == an)
        lipid_blocks += f"""
#include "{ffdir_rel}/{mol}.itp"
#ifdef POSRES_LIPID
; z-only restraint on {an} ({rn}) to keep the leaflets flat during equilibration
[ position_restraints ]
;  ai  funct  fcx  fcy  fcz
{idx:6d}    1     0    0   POSRES_FC_LIPID
#endif
"""
    with open(os.path.join(state_dir, "topol.top"), "w") as f:
        f.write(f"""; SLC35A2 (res 29-341, ACE/NME caps) in POPC:POPE:CHL 40:40:20 bilayer, TIP3P, 0.15 M NaCl
; protein: amber14sb.ff (GROMACS 2026 built-in ff14SB port)
; lipids:  AMBER Lipid21 GROMACS port ({ffdir_rel}), added with lipid21_params.itp
; water/ions: TIP3P + Joung-Cheatham ions from amber14sb.ff

; default restraint force constants (overridden by -D... in the mdp 'define')
#ifndef POSRES_FC_BB
#define POSRES_FC_BB 1000.0
#endif
#ifndef POSRES_FC_SC
#define POSRES_FC_SC 1000.0
#endif
#ifndef POSRES_FC_LIPID
#define POSRES_FC_LIPID 1000.0
#endif

#include "amber14sb.ff/forcefield.itp"
#include "{ffdir_rel}/lipid21_params.itp"

#include "protein.itp"
#ifdef POSRES
#include "posre_protein.itp"
#endif
{lipid_blocks}
#include "amber14sb.ff/tip3p.itp"
#include "amber14sb.ff/ions_tip3p.itp"

[ system ]
SLC35A2 in POPC:POPE:CHL 40:40:20 bilayer, TIP3P, 0.15 M NaCl

[ molecules ]
""")
        for m, n in counts:
            f.write(f"{m:10s} {n:8d}\n")

    # ---- index
    groups = {"System": [], "Protein": [], "Backbone": [], "C-alpha": [], "SOLU": [], "MEMB": [], "SOLV": [], "SOLU_MEMB": []}
    for n, (ri, rn, an, _) in enumerate(out, 1):
        groups["System"].append(n)
        if n <= nprot:
            groups["Protein"].append(n)
            groups["SOLU"].append(n)
            groups["SOLU_MEMB"].append(n)
            if an in ("N", "CA", "C") and rn not in ("ACE", "NME"):
                groups["Backbone"].append(n)
            if an == "CA" and rn not in ("ACE", "NME"):
                groups["C-alpha"].append(n)
        elif rn in ("PA", "PC", "PE", "OL", "CHL"):
            groups["MEMB"].append(n)
            groups["SOLU_MEMB"].append(n)
        else:
            groups["SOLV"].append(n)
    with open(os.path.join(state_dir, "index.ndx"), "w") as f:
        for g, idx in groups.items():
            f.write(f"[ {g} ]\n")
            for k in range(0, len(idx), 15):
                f.write(" ".join(f"{x:6d}" for x in idx[k:k + 15]) + "\n")

    # ---- report
    up = lambda m: sum(1 for x in mols[m] if np.mean([a[1][2] for r in x for a in r["atoms"]]) > 0)
    print(f"box (nm): {box[0]/10:.3f} {box[1]/10:.3f} {box[2]/10:.3f}   atoms: {len(out)}")
    for m, n in counts:
        extra = f"   (upper/lumen {up(m)}, lower/cytosol {n - up(m)})" if m in LIPIDS else ""
        print(f"  {m:8s} {n:6d}{extra}")
    print("  groups: " + ", ".join(f"{g}={len(v)}" for g, v in groups.items()))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
