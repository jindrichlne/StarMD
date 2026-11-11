#!/usr/bin/env python3
"""Add ACE (N-terminal) and NME (C-terminal) cap heavy atoms to a heavy-atom PDB.

The README asks for neutral caps on the artificial truncation points (res 29 and 341).
Cap atoms are built with ideal Amber/peptide internal coordinates (NeRF):
  ACE: C, O, CH3 attached to N of the first residue (omega = 180 deg); the free
       phi of the first residue is scanned and the most clash-free value is kept.
  NME: N in the C/CA/O plane (sp2, opposite the CA-O bisector), C(methyl) trans to CA.
Hydrogens are added later by pdb2gmx (amber14sb.ff ACE/NME .hdb entries).
UniProt numbering is kept: ACE gets resid first-1, NME gets resid last+1.

usage: cap_termini.py in.pdb out.pdb
"""
import sys
import numpy as np


def nerf(a, b, c, bond, angle, torsion):
    """Place D so that |CD|=bond, angle(B,C,D)=angle, dihedral(A,B,C,D)=torsion (deg)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(angle),
                   bond * np.sin(angle) * np.cos(torsion),
                   bond * np.sin(angle) * np.sin(torsion)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def fmt(serial, name, resn, chain, resi, xyz):
    nm = name if len(name) == 4 else " " + name.ljust(3)
    return ("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s\n"
            % (serial % 100000, nm, resn, chain, resi, xyz[0], xyz[1], xyz[2], name[0]))


def main(fin, fout):
    atoms = []
    for line in open(fin):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(dict(name=line[12:16].strip(), resn=line[17:20].strip(),
                              chain=line[21], resi=int(line[22:26]),
                              xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                              line=line))
    chain = atoms[0]["chain"]
    first, last = atoms[0]["resi"], atoms[-1]["resi"]
    get = lambda r, n: next(a["xyz"] for a in atoms if a["resi"] == r and a["name"] == n)
    allxyz = np.array([a["xyz"] for a in atoms])

    # ---- ACE: scan phi of the first residue and keep the least clashing cap
    N, CA, C = get(first, "N"), get(first, "CA"), get(first, "C")
    best = None
    for phi in range(-180, 180, 10):
        c0 = nerf(C, CA, N, 1.335, 121.9, phi)          # C(i-1): C-N-CA angle, phi = C(i-1)-N-CA-C
        ch3 = nerf(CA, N, c0, 1.522, 116.6, 180.0)      # CA(i-1)-C(i-1)-N-CA = omega = 180
        o = nerf(CA, N, c0, 1.229, 122.9, 0.0)          # O(i-1) cis to CA(i)
        cap = np.array([ch3, c0, o])
        far = np.linalg.norm(allxyz - CA, axis=1) > 2.6  # ignore the bonded first residue backbone
        dmin = np.min(np.linalg.norm(allxyz[far][:, None, :] - cap[None], axis=2))
        if best is None or dmin > best[0]:
            best = (dmin, phi, cap)
    dmin_ace, phi_ace, ace = best

    # ---- NME: N in plane, methyl trans to CA
    CAl, Cl, Ol = get(last, "CA"), get(last, "C"), get(last, "O")
    u = (CAl - Cl) / np.linalg.norm(CAl - Cl) + (Ol - Cl) / np.linalg.norm(Ol - Cl)
    nN = Cl - 1.335 * u / np.linalg.norm(u)
    nC = nerf(CAl, Cl, nN, 1.449, 121.9, 180.0)
    nme = np.array([nN, nC])
    dmin_nme = np.min(np.linalg.norm(allxyz[np.linalg.norm(allxyz - Cl, axis=1) > 1.6][:, None, :] - nme[None], axis=2))

    out, s = [], 1
    for name, xyz in zip(["CH3", "C", "O"], ace):
        out.append(fmt(s, name, "ACE", chain, first - 1, xyz)); s += 1
    for a in atoms:
        out.append(a["line"][:6] + "%5d" % (s % 100000) + a["line"][11:]); s += 1
    for name, xyz in zip(["N", "C"], nme):
        out.append(fmt(s, name, "NME", chain, last + 1, xyz)); s += 1
    with open(fout, "w") as f:
        f.writelines(out)
        f.write("TER\nEND\n")
    print(f"ACE: phi({first})={phi_ace} deg, min dist to non-bonded atoms {dmin_ace:.2f} A")
    print(f"NME: min dist to non-bonded atoms {dmin_nme:.2f} A")


if __name__ == "__main__":
    main(*sys.argv[1:3])
