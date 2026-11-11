#!/usr/bin/env python3
"""Sanity checks for a packed + minimised membrane system (packmol artefacts).

1. ring piercing: any bond threading a cholesterol ring or a protein aromatic/Pro ring
2. chirality: sign of the signed volume at every heavy-atom centre with >= 3 heavy
   neighbours, compared with a reference structure (the packed input, whose internal
   geometry comes from rigid templates / pdb2gmx)
3. oleoyl cis double bond (OL C18-C19=C110-C111)

usage: check_structure.py <topol.top> <reference.gro> <test.gro>
Exit status 1 if any problem is found.
"""
import os
import re
import sys

import numpy as np

RINGS = {
    "CHL": [["C1", "C2", "C3", "C4", "C5", "C10"], ["C5", "C6", "C7", "C8", "C9", "C10"],
            ["C8", "C9", "C11", "C12", "C13", "C14"], ["C13", "C14", "C15", "C16", "C17"]],
    "PHE": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "TYR": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"], ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]],
    "HISE": [["CG", "ND1", "CE1", "NE2", "CD2"]], "HIE": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "HID": [["CG", "ND1", "CE1", "NE2", "CD2"]], "HIP": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "PRO": [["N", "CA", "CB", "CG", "CD"]],
}


def parse_itp_molecules(path, seen=None, defines=None):
    """Very small topology reader: returns {molname: (atoms[(resnr,resn,name)], bonds[(i,j)])}
    and the [ molecules ] list. Follows #include and honours #ifdef/#ifndef/#else/#endif."""
    mols, order = {}, []
    defines = set() if defines is None else defines
    stack, cur, sec = [], None, None

    def active():
        return all(stack)

    def read(fn):
        nonlocal cur, sec
        d = os.path.dirname(fn)
        for raw in open(fn):
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                tok = line.split()
                if tok[0] == "#ifdef":
                    stack.append(tok[1] in defines)
                elif tok[0] == "#ifndef":
                    stack.append(tok[1] not in defines)
                elif tok[0] == "#else":
                    stack[-1] = not stack[-1]
                elif tok[0] == "#endif":
                    stack.pop()
                elif tok[0] == "#define" and active():
                    defines.add(tok[1])
                elif tok[0] == "#include" and active():
                    inc = tok[1].strip('"')
                    p = os.path.join(d, inc)
                    if not os.path.exists(p):
                        gmxdata = os.environ.get("GMXDATA") or os.path.join(os.path.dirname(os.path.dirname(
                            os.popen("which gmx").read().strip())), "share", "gromacs")
                        p = os.path.join(gmxdata, "top", inc)
                    read(p)
                continue
            if not active():
                continue
            if line.startswith("["):
                sec = line.strip("[] ").lower()
                continue
            f = line.split()
            if sec == "moleculetype":
                cur = f[0]
                mols[cur] = ([], [])
            elif sec == "atoms" and cur:
                mols[cur][0].append((int(f[2]), f[3], f[4]))
            elif sec in ("bonds", "settles") and cur:
                if sec == "settles":
                    i = int(f[0])
                    mols[cur][1].extend([(i - 1, i), (i - 1, i + 1)])
                else:
                    mols[cur][1].append((int(f[0]) - 1, int(f[1]) - 1))
            elif sec == "molecules":
                order.append((f[0], int(f[1])))
    read(path)
    return mols, order


def read_gro(fn):
    L = open(fn).read().splitlines()
    n = int(L[1])
    x = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in L[2:2 + n]])
    box = np.array([float(v) for v in L[2 + n].split()[:3]])
    return x * 10.0, box * 10.0  # Angstrom


def mic(d, box):
    return d - box * np.round(d / box)


def build_system(top):
    mols, order = parse_itp_molecules(top)
    names, resn, resid, molid, bonds = [], [], [], [], []
    off, rid, mid = 0, 0, 0
    for m, n in order:
        atoms, b = mols[m]
        for _ in range(n):
            last = None
            for (rnr, rn, an) in atoms:
                if rnr != last:
                    rid += 1
                    last = rnr
                names.append(an); resn.append(rn); resid.append(rid); molid.append(mid)
            bonds.extend((i + off, j + off) for i, j in b)
            off += len(atoms)
            mid += 1
    return (np.array(names), np.array(resn), np.array(resid), np.array(molid), np.array(bonds), order)


def point_in_polygon(pt, poly):
    inside, n = False, len(poly)
    for k in range(n):
        a, b = poly[k], poly[(k + 1) % n]
        if (a[1] > pt[1]) != (b[1] > pt[1]):
            xint = a[0] + (pt[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if pt[0] < xint:
                inside = not inside
    return inside


def check_rings(x, box, names, resn, resid, bonds):
    problems = []
    bond_mid = x[bonds[:, 0]] + 0.5 * mic(x[bonds[:, 1]] - x[bonds[:, 0]], box)
    ring_list = []
    # find rings residue by residue (vectorised by residue boundaries)
    starts = np.r_[0, np.where(np.diff(resid) != 0)[0] + 1]
    ends = np.r_[starts[1:], len(resid)]
    for s, e in zip(starts, ends):
        rn = resn[s]
        if rn not in RINGS:
            continue
        nm = {names[k]: k for k in range(s, e)}
        for ring in RINGS[rn]:
            ring_list.append(([nm[a] for a in ring], f"{rn}{resid[s]}"))
    for idx, label in ring_list:
        p = x[idx[0]] + mic(x[idx] - x[idx[0]], box)
        c = p.mean(0)
        u, sv, vt = np.linalg.svd(p - c)
        nrm, e1, e2 = vt[2], vt[0], vt[1]
        poly = np.c_[(p - c) @ e1, (p - c) @ e2]
        own = resid[idx[0]]
        near = np.where(np.linalg.norm(mic(bond_mid - c, box), axis=1) < 3.0)[0]
        for b in near:
            i, j = bonds[b]
            if resid[i] == own and resid[j] == own:
                continue
            if i in idx or j in idx:
                continue
            pi = c + mic(x[i] - c, box)
            pj = pi + mic(x[j] - x[i], box)
            di, dj = (pi - c) @ nrm, (pj - c) @ nrm
            if di * dj >= 0:
                continue
            q = pi + (pj - pi) * di / (di - dj)
            if point_in_polygon(np.array([(q - c) @ e1, (q - c) @ e2]), poly):
                problems.append(f"ring {label} pierced by bond {resn[i]}{resid[i]}:{names[i]}-{resn[j]}{resid[j]}:{names[j]}")
    return problems, len(ring_list)


def chiral_centres(names, bonds):
    nb = {}
    for i, j in bonds:
        nb.setdefault(i, []).append(j)
        nb.setdefault(j, []).append(i)
    heavy = lambda k: not names[k].startswith("H")
    out = []
    for a, lst in nb.items():
        if not heavy(a) or len(lst) != 4:
            continue  # only sp3 centres (4 bonded partners)
        h = sorted(k for k in lst if heavy(k))
        if len(h) >= 3:
            out.append((a, h[0], h[1], h[2]))
    return np.array(out)


def signed_vol(x, box, c):
    a = x[c[:, 0]]
    v1, v2, v3 = (mic(x[c[:, k]] - a, box) for k in (1, 2, 3))
    return np.einsum("ij,ij->i", v1, np.cross(v2, v3))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - np.einsum("ij,ij->i", b0, b1)[:, None] * b1
    w = b2 - np.einsum("ij,ij->i", b2, b1)[:, None] * b1
    xx = np.einsum("ij,ij->i", v, w)
    yy = np.einsum("ij,ij->i", np.cross(b1, v), w)
    return np.degrees(np.arctan2(yy, xx))


def main(top, ref, test):
    names, resn, resid, molid, bonds, order = build_system(top)
    xr, boxr = read_gro(ref)
    xt, boxt = read_gro(test)
    assert len(xr) == len(xt) == len(names), (len(xr), len(xt), len(names))
    problems = []

    rp, nring = check_rings(xt, boxt, names, resn, resid, bonds)
    problems += rp
    print(f"rings checked: {nring}; pierced: {len(rp)}")

    cc = chiral_centres(names, bonds)
    vr, vt = signed_vol(xr, boxr, cc), signed_vol(xt, boxt, cc)
    flip = np.where((np.sign(vr) != np.sign(vt)) & (np.abs(vr) > 0.5))[0]
    for k in flip:
        a = cc[k, 0]
        problems.append(f"chirality flipped at {resn[a]}{resid[a]}:{names[a]}")
    print(f"sp3 centres checked: {len(cc)}; inverted: {len(flip)}")

    sel = lambda nm: np.where((resn == "OL") & (names == nm))[0]
    p = [xt[sel(n)] for n in ("C18", "C19", "C110", "C111")]
    for k in (1, 2, 3):  # make each quadruplet whole across PBC
        p[k] = p[k - 1] + mic(p[k] - p[k - 1], boxt)
    d = dihedral(*p)
    trans = np.where(np.abs(d) > 90)[0]
    for k in trans:
        problems.append(f"oleoyl double bond not cis (dihedral {d[k]:.0f}) in OL{resid[sel('C19')[k]]}")
    print(f"oleoyl double bonds checked: {len(d)}; not cis: {len(trans)}")

    for p in problems:
        print("PROBLEM:", p)
    print("OK" if not problems else f"{len(problems)} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main(*sys.argv[1:4])
