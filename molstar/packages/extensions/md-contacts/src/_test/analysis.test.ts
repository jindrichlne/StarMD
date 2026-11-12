/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { describe, expect, it } from 'vitest';
import {
  evaluateStability,
  findConformers,
  frameDistances,
  mapStats,
  replicaAgreement,
  switchPairs,
  type TrajectoryData,
} from '../analysis.js';
import { parseXpm, parseXvg } from '../xpm.js';

/** A chain of `n` residues on a line, with residues from `hinge` on swung by `angle(frame)` and a little noise. */
function hinged(nFrames: number, n: number, hinge: number, angle: (f: number) => number, seed = 1): TrajectoryData {
  let s = seed;
  const rnd = () => {
    s = (s * 1664525 + 1013904223) % 4294967296;
    return s / 4294967296 - 0.5;
  };
  const xyz = new Float32Array(nFrames * n * 3);
  for (let f = 0; f < nFrames; f++) {
    const a = angle(f);
    for (let i = 0; i < n; i++) {
      const d = 0.38 * (i - hinge);
      const x = i < hinge ? d : d * Math.cos(a),
        y = i < hinge ? 0 : d * Math.sin(a);
      const o = (f * n + i) * 3;
      xyz[o] = x + rnd() * 0.02;
      xyz[o + 1] = y + rnd() * 0.02;
      xyz[o + 2] = rnd() * 0.02;
    }
  }
  return {
    nFrames,
    nRes: n,
    times: Float32Array.from({ length: nFrames }, (_, f) => f * 100),
    xyz,
    residues: Array.from({ length: n }, (_, i) => ({ seq: i + 1, name: 'ALA', chain: 'A', atom: i, hasCA: true })),
  };
}

/** A compact zig-zag fold of `n` residues; residues 20 to 22 are displaced by `shift(frame)` nm; uniform noise. */
function folded(nFrames: number, n: number, shift: (f: number) => number, noise: number, seed = 5): TrajectoryData {
  let s = seed;
  const rnd = () => {
    s = (s * 1664525 + 1013904223) % 4294967296;
    return (s / 4294967296 - 0.5) * noise;
  };
  const xyz = new Float32Array(nFrames * n * 3);
  for (let f = 0; f < nFrames; f++) {
    const d = shift(f);
    for (let i = 0; i < n; i++) {
      const row = Math.floor(i / 10),
        col = row % 2 ? 9 - (i % 10) : i % 10;
      const o = (f * n + i) * 3;
      xyz[o] = col * 0.38 + rnd();
      xyz[o + 1] = row * 0.5 + rnd();
      xyz[o + 2] = (i >= 20 && i <= 22 ? d : 0) + rnd();
    }
  }
  return {
    nFrames,
    nRes: n,
    times: Float32Array.from({ length: nFrames }, (_, f) => f * 100),
    xyz,
    residues: Array.from({ length: n }, (_, i) => ({ seq: i + 1, name: 'ALA', chain: 'A', atom: i, hasCA: true })),
  };
}

describe('md-contacts analysis', () => {
  it('computes symmetric frame distances', () => {
    const t = hinged(2, 10, 5, () => 0);
    const d = frameDistances(t, 0);
    expect(d[0 * 10 + 9]).toBeCloseTo(9 * 0.38, 1);
    expect(d[3 * 10 + 7]).toBe(d[7 * 10 + 3]);
    expect(d[4 * 10 + 4]).toBe(0);
  });

  it('judges a steady shape as held and an opening hinge as departed', () => {
    const steady = evaluateStability(hinged(200, 30, 15, () => 0.2));
    expect(steady.verdict).toBe('holds');
    const opening = evaluateStability(hinged(200, 30, 15, (f) => 0.2 + (f / 200) * 2.4));
    expect(opening.verdict).not.toBe('holds');
    expect(opening.drmsdLate).toBeGreaterThan(steady.drmsdLate + 0.3);
    expect(evaluateStability(hinged(5, 30, 15, () => 0)).verdict).toBe('short');
    // moves in the first tenth, then sits still: flat at the end, but not where it started
    const settled = evaluateStability(hinged(400, 30, 15, (f) => 0.2 + Math.min(1, f / 40) * 0.5));
    expect(settled.verdict).toBe('shifted');
    expect(Math.abs(settled.drift)).toBeLessThan(0.03);
  });

  it('tells agreeing replicas from disagreeing ones', () => {
    const a = mapStats(
      hinged(100, 30, 15, () => 0.2, 1),
      0.8,
    );
    const b = mapStats(
      hinged(100, 30, 15, () => 0.2, 2),
      0.8,
    );
    const c = mapStats(
      hinged(100, 30, 15, () => 1.6, 3),
      0.8,
    );
    expect(replicaAgreement([a, b]).agree).toBe(true);
    expect(replicaAgreement([a, c]).agree).toBe(false);
  });

  it('finds two conformers of a two-state hinge and the pairs across it', async () => {
    const t = hinged(300, 30, 15, (f) => (f < 150 ? 0.2 : 1.8));
    const result = await findConformers([t], { framesPerTrajectory: 100 });
    expect(result.suggested).toBe(2);
    expect(result.conformers.map((c) => c.size).sort()).toEqual([150, 150]);
    expect(result.assignment[0][10]).not.toBe(result.assignment[0][290]);
    const top = switchPairs(result.conformers[0].stats, result.conformers[1].stats, { top: 5 });
    // the best pairs span the hinge
    for (const p of top) expect(p.i < 15 && p.j >= 15).toBe(true);
    const one = await findConformers([hinged(200, 30, 15, () => 0.2)], { framesPerTrajectory: 80 });
    expect(one.suggested).toBe(1);
  });

  it('finds a local rearrangement under thermal noise', async () => {
    // three residues of a loop move 1 nm in the last two thirds; every atom jitters by up to 0.08 nm
    const t = folded(600, 60, (f) => (f < 200 ? 0 : 1), 0.16);
    const r = await findConformers([t], { framesPerTrajectory: 150 });
    expect(r.suggested).toBe(2);
    expect(r.conformers.map((c) => c.size)).toEqual([400, 200]);
    const top = switchPairs(r.conformers[0].stats, r.conformers[1].stats, { top: 5 });
    for (const p of top) expect((p.i >= 20 && p.i <= 22) || (p.j >= 20 && p.j <= 22)).toBe(true);
  });

  it('does not cut one continuous motion into conformers', async () => {
    // a steady drift
    const drift = await findConformers([hinged(400, 30, 15, (f) => 0.2 + f / 400)], { framesPerTrajectory: 100 });
    expect(drift.suggested).toBe(1);
    expect(drift.splits.some((s) => s.accepted)).toBe(false);
    // a slow wobble around one shape
    let x = 0,
      seed = 7;
    const wobble = Array.from({ length: 600 }, () => {
      seed = (seed * 1664525 + 1013904223) % 4294967296;
      x = 0.95 * x + (seed / 4294967296 - 0.5) * 0.1;
      return x;
    });
    const one = await findConformers([hinged(600, 30, 15, (f) => 0.6 + wobble[f])], { framesPerTrajectory: 150 });
    expect(one.suggested).toBe(1);
    // asking for a number still cuts it
    const forced = await findConformers([hinged(400, 30, 15, (f) => 0.2 + f / 400)], {
      framesPerTrajectory: 100,
      nClusters: 3,
    });
    expect(forced.conformers.length).toBe(3);
    expect(forced.suggested).toBe(1);
  });

  it('finds conformers at different scales and across replicas', async () => {
    // two close shapes and one far away: the far one must not hide the close pair
    const t = hinged(600, 30, 15, (f) => (f < 200 ? 0.2 : f < 400 ? 0.4 : 1.8));
    const r = await findConformers([t], { framesPerTrajectory: 150 });
    expect(r.suggested).toBe(3);
    expect(new Set([r.assignment[0][100], r.assignment[0][300], r.assignment[0][500]]).size).toBe(3);
    // replicas that sit in different shapes, one of them crossing over
    const reps = [
      hinged(300, 30, 15, () => 0.2, 1),
      hinged(300, 30, 15, () => 0.9, 2),
      hinged(300, 30, 15, (f) => (f < 100 ? 0.2 : 0.9), 3),
    ];
    const across = await findConformers(reps, { framesPerTrajectory: 100 });
    expect(across.suggested).toBe(2);
    expect(across.assignment[0][150]).toBe(across.assignment[2][50]);
    expect(across.assignment[1][150]).toBe(across.assignment[2][250]);
    // a state between two others, on one line with them: both gaps are found, although neither end is dense at the
    // centre of "the other two"
    const line = await findConformers([hinged(600, 30, 15, (f) => (f < 250 ? 0.2 : f < 400 ? 0.9 : 1.6))], {
      framesPerTrajectory: 150,
    });
    expect(line.suggested).toBe(3);
    expect(line.splits.filter((s) => s.accepted).length).toBe(3);
    // a steady motion and then a jump: the motion stays one conformer, the jump makes the second
    const jump = await findConformers([hinged(600, 30, 15, (f) => (f < 400 ? 0.2 + f / 500 : 2.2))], {
      framesPerTrajectory: 150,
    });
    expect(jump.suggested).toBe(2);
    expect(jump.assignment[0][10]).toBe(jump.assignment[0][390]);
    expect(jump.assignment[0][10]).not.toBe(jump.assignment[0][500]);
    // a coarser scale never finds more
    expect(line.nearby[0].count).toBeGreaterThanOrEqual(line.suggested);
    expect(line.nearby[1].count).toBeLessThanOrEqual(line.suggested);
    // a few stray frames are not a conformer
    const stray = await findConformers([hinged(400, 30, 15, (f) => (f === 200 || f === 201 ? 2.5 : 0.2))], {
      framesPerTrajectory: 400,
    });
    expect(stray.suggested).toBe(1);
  });
});

describe('GROMACS file readers', () => {
  it('reads an XPM matrix with its levels and axes, y flipped', () => {
    const xpm = `/* XPM */
/* title:   "Mean smallest distance" */
/* legend:  "Distance (nm)" */
/* x-label: "Residue Index" */
/* y-label: "Residue Index" */
/* type:    "Continuous" */
static char *gromacs_xpm[] = {
"2 2   3 1",
"A  c #FFFFFF " /* "0" */,
"B  c #808080 " /* "0.5" */,
"C  c #000000 " /* "1" */,
/* x-axis:  1 2 */
/* y-axis:  1 2 */
"BC",
"AB"
}`;
    const [m] = parseXpm(xpm);
    expect(m.title).toBe('Mean smallest distance');
    expect([m.width, m.height]).toEqual([2, 2]);
    // the last text row is y = 1 (row 0)
    expect(Array.from(m.values)).toEqual([0, 0.5, 0.5, 1]);
    expect(m.xAxis).toEqual([1, 2]);
    expect(m.max).toBe(1);
  });

  it('reads an XVG file', () => {
    const xvg = parseXvg('# c\n@ title "RMSD"\n@ xaxis label "Time (ps)"\n@ s0 legend "Backbone"\n0 0.1\n10 0.2\n');
    expect(xvg.title).toBe('RMSD');
    expect(xvg.legends).toEqual(['Backbone']);
    expect(Array.from(xvg.columns[1])).toEqual([0.1, 0.2]);
  });
});
