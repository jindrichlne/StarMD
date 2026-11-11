#!/usr/bin/env python3
"""Post-process the packmol-memgen -> GROMACS system.

1. Hydrate the protein cavity inside the membrane slab. packmol-memgen fills water
   only outside the bilayer (|z| > leaflet width), so the open vestibule of the
   transporter starts empty. Candidate waters come from `gmx solvate` (fills all voids);
   only those that are inside the membrane slab AND buried by the protein in the xy-plane
   (ray casting: protein heavy atoms block >= MIN_BLOCKED of N_RAYS horizontal rays) AND
   not close to lipid atoms are kept. Lipid-facing grooves are therefore not hydrated.
2. Set the salt to "neutralise + 0.15 M NaCl" (README): N_NaCl = round(0.15 M * N_wat / 55.51 M),
   Na+ = N_NaCl, Cl- = N_NaCl + Q_protein. Missing ions replace random bulk waters
   (> 6 A from the protein, > 5 A from other ions, outside the membrane), pairs added
   alternately to the upper and lower water compartment.

Rewrites <state_dir>/system.gro, the [ molecules ] counts of topol.top and index.ndx.
usage: hydrate_and_ionize.py <state_dir> [conc_M=0.15] [seed=2026]
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

N_RAYS, MIN_BLOCKED, RAY_LEN, RAY_RAD = 16, 14, 20.0, 2.5   # Angstrom
LIPID_CLEAR = 3.0          # min distance of a cavity water O to lipid atoms (A)
SLAB = 24.0                # half-thickness of the region considered "membrane" (A)


def read_gro(fn):
    L = open(fn).read().splitlines()
    n = int(L[1])
    atoms = [(int(l[0:5]), l[5:10].strip(), l[10:15].strip(),
              float(l[20:28]) * 10, float(l[28:36]) * 10, float(l[36:44]) * 10) for l in L[2:2 + n]]
    box = np.array([float(v) for v in L[2 + n].split()[:3]]) * 10
    return L[0], atoms, box


def write_gro(fn, title, atoms, box):
    with open(fn, "w") as f:
        f.write(title + "\n%d\n" % len(atoms))
        for k, (ri, rn, an, x, y, z) in enumerate(atoms, 1):
            f.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (ri % 100000, rn[:5], an[:5], k % 100000, x / 10, y / 10, z / 10))
        f.write("%10.5f%10.5f%10.5f\n" % tuple(box / 10))


def mic(d, box):
    return d - box * np.round(d / box)


def buried(points, prot, box):
    """number of blocked horizontal rays for each point"""
    ang = np.linspace(0, 2 * np.pi, N_RAYS, endpoint=False)
    dirs = np.c_[np.cos(ang), np.sin(ang), np.zeros(N_RAYS)]
    out = np.zeros(len(points), int)
    for k, p in enumerate(points):
        d = mic(prot - p, box)
        near = d[np.linalg.norm(d, axis=1) < RAY_LEN + RAY_RAD]
        if len(near) == 0:
            continue
        t = near @ dirs.T                                   # projection on each ray
        perp2 = (near ** 2).sum(1)[:, None] - t ** 2
        hit = (t > 0) & (t < RAY_LEN) & (perp2 < RAY_RAD ** 2)
        out[k] = hit.any(0).sum()
    return out


def protein_charge(itp):
    q, sec = 0.0, None
    for line in open(itp):
        s = line.split(";")[0].strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            continue
        if sec == "atoms" and s:
            q += float(s.split()[6])
    return int(round(q))


def main(state_dir, conc=0.15, seed=2026):
    rng = np.random.default_rng(seed)
    gro = os.path.join(state_dir, "system.gro")
    title, atoms, box = read_gro(gro)
    resn = np.array([a[1] for a in atoms])
    xyz = np.array([a[3:] for a in atoms])
    is_prot = np.zeros(len(atoms), bool)
    # protein = everything before the first lipid residue
    first_lip = np.where(np.isin(resn, ["PA", "PC", "PE", "OL", "CHL"]))[0][0]
    is_prot[:first_lip] = True
    is_lip = np.isin(resn, ["PA", "PC", "PE", "OL", "CHL"])
    zc = xyz[is_lip & np.array([a[2] in ("P31", "O1") for a in atoms]), 2].mean()  # bilayer centre
    prot_heavy = xyz[is_prot & ~np.array([a[2].startswith("H") for a in atoms])]

    # ---- 1. cavity water
    tmp = os.path.join(state_dir, "00_build", "solvate_tmp.gro")
    subprocess.run(["gmx", "solvate", "-cp", gro, "-cs", "spc216.gro", "-o", tmp],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, satoms, _ = read_gro(tmp)
    new = satoms[len(atoms):]
    assert all(a[1] == "SOL" for a in new)
    wat = [new[k:k + 3] for k in range(0, len(new), 3)]
    O = np.array([w[0][3:] for w in wat])
    inslab = np.abs(O[:, 2] - zc) < SLAB
    cand = np.where(inslab)[0]
    nb = buried(O[cand], prot_heavy, box)
    cand = cand[nb >= MIN_BLOCKED]
    lip = xyz[is_lip]
    keep = []
    for k in cand:
        if np.min(np.linalg.norm(mic(lip - O[k], box), axis=1)) > LIPID_CLEAR:
            keep.append(k)
    cav = [wat[k] for k in keep]
    print(f"bilayer centre z = {zc / 10:.3f} nm; gmx solvate candidates in slab: {inslab.sum()}, "
          f"buried: {len(cand)}, kept as cavity water: {len(cav)}")

    # ---- assemble blocks
    blocks = {k: [] for k in ["PROT_LIP", "SOL", "NA", "CL"]}
    i = 0
    while i < len(atoms):
        rn = atoms[i][1]
        if rn == "SOL":
            blocks["SOL"].append(atoms[i:i + 3]); i += 3
        elif rn in ("NA", "CL"):
            blocks[rn].append(atoms[i]); i += 1
        else:
            blocks["PROT_LIP"].append(atoms[i]); i += 1
    # cavity waters get TIP3P names
    for w in cav:
        blocks["SOL"].append([(0, "SOL", n, *a[3:]) for n, a in zip(["OW", "HW1", "HW2"], w)])

    # ---- 2. ions
    q = protein_charge(os.path.join(state_dir, "protein.itp"))
    na, cl = len(blocks["NA"]), len(blocks["CL"])
    while True:
        nw = len(blocks["SOL"])
        n_nacl = int(round(conc * nw / 55.51))
        add_na, add_cl = n_nacl - na, n_nacl + q - cl
        if add_na <= 0 and add_cl <= 0:
            break
        if add_na < 0 or add_cl < 0 or add_na != add_cl:
            sys.exit(f"unexpected ion balance: have Na {na} Cl {cl}, want Na {n_nacl} Cl {n_nacl + q}")
        # one pair at a time (target depends weakly on the number of waters left)
        ions = np.array([a[3:] for a in blocks["NA"] + blocks["CL"]], float).reshape(-1, 3)
        Ow = np.array([w[0][3:] for w in blocks["SOL"]])
        upper = (na + cl) % 4 < 2       # alternate compartments pairwise
        dz = mic(Ow[:, 2:3] - zc, box[2:3])[:, 0]
        ok = (np.abs(dz) > 28.0) & ((dz > 0) == upper)
        idx = rng.permutation(np.where(ok)[0])
        placed = []
        for k in idx:
            p = Ow[k]
            if np.min(np.linalg.norm(mic(prot_heavy - p, box), axis=1)) < 6.0:
                continue
            allions = np.vstack([ions] + [Ow[j][None] for j in placed])
            if len(allions) and np.min(np.linalg.norm(mic(allions - p, box), axis=1)) < 5.0:
                continue
            placed.append(k)
            if len(placed) == 2:
                break
        if len(placed) < 2:
            sys.exit("could not place ions")
        kna, kcl = placed
        blocks["NA"].append((0, "NA", "NA", *Ow[kna]))
        blocks["CL"].append((0, "CL", "CL", *Ow[kcl]))
        blocks["SOL"] = [w for j, w in enumerate(blocks["SOL"]) if j not in (kna, kcl)]
        na, cl = na + 1, cl + 1
    nw = len(blocks["SOL"])
    print(f"waters {nw} (+{len(cav)} cavity), Na+ {na}, Cl- {cl}, protein charge {q:+d}, "
          f"NaCl = {na} pairs -> {na / nw * 55.51:.3f} M (Cl- {cl / nw * 55.51:.3f} M)")

    # ---- renumber and write
    out = list(blocks["PROT_LIP"])
    resi = out[-1][0]
    for w in blocks["SOL"]:
        resi += 1
        out += [(resi, "SOL", n, *a[3:]) for n, a in zip(["OW", "HW1", "HW2"], w)]
    for ion in ("NA", "CL"):
        for a in blocks[ion]:
            resi += 1
            out.append((resi, ion, ion, *a[3:]))
    write_gro(gro, title, out, box)

    top = os.path.join(state_dir, "topol.top")
    lines = open(top).read().splitlines()
    k = lines.index("[ molecules ]")
    counts = {"SOL": nw, "NA": na, "CL": cl}
    for j in range(k + 1, len(lines)):
        f = lines[j].split()
        if f and f[0] in counts:
            lines[j] = f"{f[0]:10s} {counts[f[0]]:8d}"
    open(top, "w").write("\n".join(lines) + "\n")

    # ---- index
    groups = {"System": [], "Protein": [], "Backbone": [], "C-alpha": [], "SOLU": [], "MEMB": [], "SOLV": [], "SOLU_MEMB": []}
    for n, (ri, rn, an, *_) in enumerate(out, 1):
        groups["System"].append(n)
        if n <= first_lip:
            for g in ("Protein", "SOLU", "SOLU_MEMB"):
                groups[g].append(n)
            if rn not in ("ACE", "NME") and an in ("N", "CA", "C"):
                groups["Backbone"].append(n)
            if rn not in ("ACE", "NME") and an == "CA":
                groups["C-alpha"].append(n)
        elif rn in ("PA", "PC", "PE", "OL", "CHL"):
            groups["MEMB"].append(n); groups["SOLU_MEMB"].append(n)
        else:
            groups["SOLV"].append(n)
    with open(os.path.join(state_dir, "index.ndx"), "w") as f:
        for g, idx in groups.items():
            f.write(f"[ {g} ]\n")
            for j in range(0, len(idx), 15):
                f.write(" ".join(f"{x:6d}" for x in idx[j:j + 15]) + "\n")
    print(f"atoms {len(out)}; groups: " + ", ".join(f"{g}={len(v)}" for g, v in groups.items()))


if __name__ == "__main__":
    main(sys.argv[1], *(float(a) if i == 0 else int(a) for i, a in enumerate(sys.argv[2:])))
