/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { AlignmentOp, globalAlign } from '@molstar/core/algorithms/sequence/needleman-wunsch';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { MinimizeRmsd } from '@molstar/core/math/linear-algebra/3d/minimize-rmsd';
import { type ResidueInfo, residueKey, type TrajectoryData } from './analysis.js';

/** Sequence alignment of two simulations' residues, and the superposition of two of their frames on it. */

// biome-ignore format: a table
const ONE: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  // force-field names of protonation states and bonded forms
  HID: 'H', HIE: 'H', HIP: 'H', HISD: 'H', HISE: 'H', HISH: 'H', HSD: 'H', HSE: 'H', HSP: 'H',
  CYX: 'C', CYM: 'C', ASH: 'D', GLH: 'E', LYN: 'K', MSE: 'M', SEC: 'U',
};

export function oneLetter(name: string) {
  return ONE[name] ?? 'X';
}

export interface SequenceAlignment {
  /** residue index in `a` and in `b` of every column; -1 is a gap */
  columns: { a: number; b: number }[];
  /** the aligned residue pairs: `pairs[2k]` in `a`, `pairs[2k + 1]` in `b` */
  pairs: Int32Array;
  /** share of the aligned pairs with the same residue type */
  identity: number;
  /** how the columns were found */
  method: 'same residues' | 'sequence alignment';
}

/**
 * Aligns the residues of two simulations: as they are when both have the same residues (chain and number), else by a
 * global sequence alignment (so a mutant, a homologue or a renumbered copy still pairs up).
 */
export function alignResidues(a: ResidueInfo[], b: ResidueInfo[]): SequenceAlignment {
  const columns: { a: number; b: number }[] = [];
  const same = a.length === b.length && a.every((r, i) => residueKey(r) === residueKey(b[i]));
  if (same) {
    for (let i = 0; i < a.length; i++) columns.push({ a: i, b: i });
  } else {
    const code = (r: ResidueInfo) => oneLetter(r.name).charCodeAt(0);
    const { ops } = globalAlign(a.map(code), b.map(code), {
      match: 5,
      mismatch: -4,
      gapOpen: 10,
      gapExtend: 1,
      leadingGaps: 'free',
      trailingGaps: 'free',
      gapTies: 'open',
    });
    let i = 0,
      j = 0;
    for (const op of ops) {
      if (op === AlignmentOp.Pair) columns.push({ a: i++, b: j++ });
      else if (op === AlignmentOp.OnlyA) columns.push({ a: i++, b: -1 });
      else columns.push({ a: -1, b: j++ });
    }
  }
  const pairs: number[] = [];
  let identical = 0;
  for (const c of columns) {
    if (c.a < 0 || c.b < 0) continue;
    pairs.push(c.a, c.b);
    if (oneLetter(a[c.a].name) === oneLetter(b[c.b].name)) identical++;
  }
  return {
    columns,
    pairs: Int32Array.from(pairs),
    identity: pairs.length ? identical / (pairs.length / 2) : 0,
    method: same ? 'same residues' : 'sequence alignment',
  };
}

export interface Superposition {
  /** nm, over the aligned residues after the fit */
  rmsd: number;
  /** moves `b`'s coordinates (nm) onto `a` */
  transform: Mat4;
  /** nm: how far each aligned pair is apart after the fit, in the order of `pairs` */
  displacement: Float32Array;
}

/** Superposes frame `fb` of `b` on frame `fa` of `a` over the aligned residue pairs (representative atoms). */
export function superpose(
  a: TrajectoryData,
  fa: number,
  b: TrajectoryData,
  fb: number,
  pairs: Int32Array,
): Superposition {
  const np = pairs.length / 2;
  const pa = MinimizeRmsd.Positions.empty(np),
    pb = MinimizeRmsd.Positions.empty(np);
  const oa = fa * a.nRes * 3,
    ob = fb * b.nRes * 3;
  for (let k = 0; k < np; k++) {
    const ia = oa + pairs[2 * k] * 3,
      ib = ob + pairs[2 * k + 1] * 3;
    (pa.x as Float64Array)[k] = a.xyz[ia];
    (pa.y as Float64Array)[k] = a.xyz[ia + 1];
    (pa.z as Float64Array)[k] = a.xyz[ia + 2];
    (pb.x as Float64Array)[k] = b.xyz[ib];
    (pb.y as Float64Array)[k] = b.xyz[ib + 1];
    (pb.z as Float64Array)[k] = b.xyz[ib + 2];
  }
  const displacement = new Float32Array(np);
  if (np < 3) return { rmsd: 0, transform: Mat4.identity(), displacement };
  const { bTransform: m, rmsd } = MinimizeRmsd.compute({ a: pa, b: pb });
  for (let k = 0; k < np; k++) {
    const x = pb.x[k],
      y = pb.y[k],
      z = pb.z[k];
    const dx = m[0] * x + m[4] * y + m[8] * z + m[12] - pa.x[k],
      dy = m[1] * x + m[5] * y + m[9] * z + m[13] - pa.y[k],
      dz = m[2] * x + m[6] * y + m[10] * z + m[14] - pa.z[k];
    displacement[k] = Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  return { rmsd, transform: Mat4.clone(m), displacement };
}
