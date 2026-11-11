/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { describe, expect, it } from 'vitest';
import { findConformers, type TrajectoryData } from '../analysis.js';
import { localSwitches, parseThreeDi, summarizeThreeDi } from '../threedi.js';

const identity = (n: number) => Int32Array.from({ length: n }, (_, i) => i);

describe('3Di trajectories', () => {
  it('reads one line per frame and picks the analysed residues', () => {
    const t = parseThreeDi('DVLQ\nDVPQ\nDALQ\n', 3, 4, 4, identity(4));
    expect(String.fromCharCode(t.letters[1 * 4 + 2])).toBe('P');
    // capped termini: Foldseek saw 6 residues, 4 are analysed
    const capped = parseThreeDi('XDVLQX\nXDVPQX\n', 2, 4, 6, [1, 2, 3, 4]);
    expect(String.fromCharCode(...capped.letters.slice(0, 4))).toBe('DVLQ');
    expect(() => parseThreeDi('DVLQ\n', 2, 4, 4, identity(4))).toThrow(/frames/);
    expect(() => parseThreeDi('DVL\nDVL\n', 2, 4, 6, identity(4))).toThrow(/letters/);
  });

  it('summarises columns: entropy, usual letter, change since the start', () => {
    const t = parseThreeDi('DVLQ\nDVPQ\nDVLQ\nDVPA\n', 4, 4, 4, identity(4));
    const s = summarizeThreeDi(t);
    expect(s.entropy[0]).toBe(0);
    expect(s.entropy[2]).toBeCloseTo(1, 5);
    expect(String.fromCharCode(s.consensus[3])).toBe('Q');
    expect(s.consensusShare[3]).toBeCloseTo(0.75, 5);
    expect(Array.from(s.changed)).toEqual([0, 0.25, 0, 0.5]);
  });

  it('finds the residues whose letter differs between two sets of frames', () => {
    const t = parseThreeDi('DVLQ\nDVLQ\nDVPQ\nDVPQ\n', 4, 4, 4, identity(4));
    const out = localSwitches(
      [0, 1].map((frame) => ({ t, frame })),
      [2, 3].map((frame) => ({ t, frame })),
      4,
    );
    expect(out.length).toBe(1);
    expect(out[0]).toMatchObject({ residue: 2, letterA: 'L', letterB: 'P' });
    expect(out[0].divergence).toBeCloseTo(1, 5);
  });
});

describe('conformer tree', () => {
  // two shapes of a 12-residue chain: straight, and bent at residue 6
  const n = 12,
    nFrames = 60;
  const xyz = new Float32Array(nFrames * n * 3);
  const letters = new Uint8Array(nFrames * n);
  for (let f = 0; f < nFrames; f++) {
    const bent = f >= 40;
    for (let i = 0; i < n; i++) {
      const o = (f * n + i) * 3;
      const d = 0.38 * (i - 6);
      xyz[o] = (i < 6 || !bent ? d : 0) + 0.004 * Math.sin(f * 7 + i);
      xyz[o + 1] = i >= 6 && bent ? d : 0;
      letters[f * n + i] = (bent && i >= 5 && i <= 7 ? 'P' : 'D').charCodeAt(0);
    }
  }
  const t: TrajectoryData = {
    nFrames,
    nRes: n,
    times: Float32Array.from({ length: nFrames }, (_, f) => f * 100),
    xyz,
    residues: Array.from({ length: n }, (_, i) => ({ seq: i + 1, name: 'ALA', chain: 'A', atom: i, hasCA: true })),
  };

  it('keeps the tree: a root that splits into the two shapes by their share of frames', async () => {
    const r = await findConformers([t], { framesPerTrajectory: 60 });
    expect(r.unit).toBe('nm');
    const root = r.tree[0];
    expect(root.members.length).toBe(60);
    const [a, b] = root.children!.map((id) => r.tree[id]);
    expect([a.members.length, b.members.length].sort((x, y) => x - y)).toEqual([20, 40]);
    expect(a.x0).toBe(0);
    expect(b.x1).toBeCloseTo(1, 6);
    expect(a.conformer).not.toBe(b.conformer);
    expect(root.conformer).toBe(-1);
    // the cut lies between the split of the two shapes and the splits inside them
    expect(r.cutHeight).toBeLessThan(root.height);
    expect(r.cutHeight).toBeGreaterThan(Math.max(a.height, b.height));
    for (const node of r.tree) expect(node.members).toContain(node.medoid);
  });

  it('clusters by 3Di letters when asked', async () => {
    const r = await findConformers([t], { framesPerTrajectory: 60, fingerprint: '3di' }, { letters: [letters] });
    expect(r.unit).toBe('letters');
    expect(r.suggested).toBe(2);
    expect(r.conformers.map((c) => c.size).sort((x, y) => x - y)).toEqual([20, 40]);
    // three of twelve residues change their letter
    expect(r.tree[0].height).toBeCloseTo(0.25, 5);
    await expect(findConformers([t], { fingerprint: '3di' })).rejects.toThrow(/3Di/);
  });
});
