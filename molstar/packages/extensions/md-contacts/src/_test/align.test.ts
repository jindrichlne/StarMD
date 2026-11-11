/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { describe, expect, it } from 'vitest';
import { alignResidues, superpose } from '../align.js';
import type { ResidueInfo, TrajectoryData } from '../analysis.js';

const residues = (seq: string, first = 1): ResidueInfo[] => {
  const three: Record<string, string> = {
    A: 'ALA',
    G: 'GLY',
    L: 'LEU',
    K: 'LYS',
    D: 'ASP',
    S: 'SER',
    T: 'THR',
    V: 'VAL',
    F: 'PHE',
    E: 'GLU',
  };
  return Array.from(seq, (c, i) => ({ seq: first + i, name: three[c], chain: 'A', atom: i, hasCA: true }));
};

function traj(res: ResidueInfo[], xyz: number[][]): TrajectoryData {
  return { nFrames: 1, nRes: res.length, times: Float32Array.of(0), xyz: Float32Array.from(xyz.flat()), residues: res };
}

describe('md-contacts alignment', () => {
  it('pairs identical residues as they are', () => {
    const a = residues('AGLKDSTV');
    const out = alignResidues(a, residues('AGLKDSTV'));
    expect(out.method).toBe('same residues');
    expect(out.pairs.length / 2).toBe(8);
    expect(out.identity).toBe(1);
  });

  it('aligns a renumbered copy with a deletion and a mutation by sequence', () => {
    const a = residues('AGLKDSTVFEAGLKDSTVFE');
    // residues 5-6 deleted, one mutation, numbering shifted
    const b = residues('AGLKTVFEAGLKASTVFE', 101);
    const out = alignResidues(a, b);
    expect(out.method).toBe('sequence alignment');
    expect(out.columns.filter((c) => c.b < 0).length).toBe(2);
    expect(out.pairs.length / 2).toBe(18);
    expect(out.identity).toBeCloseTo(17 / 18, 5);
  });

  it('superposes a rotated and shifted copy exactly, and reports what moved', () => {
    const res = residues('AGLKDSTV');
    const pts = [
      [0, 0, 0],
      [0.38, 0, 0],
      [0.6, 0.3, 0],
      [0.9, 0.5, 0.2],
      [1.2, 0.4, 0.5],
      [1.4, 0.7, 0.6],
      [1.7, 0.6, 0.9],
      [2.0, 0.9, 1.0],
    ];
    // 90 degrees about z, then a shift
    const moved = pts.map(([x, y, z]) => [-y + 3, x - 2, z + 1]);
    const a = traj(res, pts);
    const { pairs } = alignResidues(res, res);
    const same = superpose(a, 0, traj(res, moved), 0, pairs);
    // float32 coordinates
    expect(same.rmsd).toBeLessThan(1e-3);
    const m = same.transform;
    const [x, y, z] = moved[3];
    expect(m[0] * x + m[4] * y + m[8] * z + m[12]).toBeCloseTo(pts[3][0], 3);
    // one residue displaced: it stands out
    const bent = moved.map((p, i) => (i === 7 ? [p[0], p[1], p[2] + 1] : p));
    const fit = superpose(a, 0, traj(res, bent), 0, pairs);
    expect(fit.displacement.indexOf(Math.max(...fit.displacement))).toBe(7);
  });
});
