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
