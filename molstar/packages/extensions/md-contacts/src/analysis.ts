/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { MinimizeRmsd } from '@molstar/core/math/linear-algebra/3d/minimize-rmsd';

/**
 * The numeric core of the MD contacts extension: residue-residue distance maps of trajectory frames, their mean and
 * spread, stability measures per replica, agreement of replicas, and conformers (clusters of frames by distance-map
 * similarity) with the residue pairs that tell them apart. No plugin or DOM imports: it runs in a worker or in Node.
 *
 * Distances are in nm (the GROMACS unit), times in ps.
 */

/** One residue of the analysed chain(s), with the atom that stands for it (CA when it has one). */
export interface ResidueInfo {
  /** residue number as in the structure file (never renumbered) */
  seq: number;
  name: string;
  chain: string;
  /** index of the representative atom in the trajectory frames */
  atom: number;
  hasCA: boolean;
}

/** Representative-atom coordinates of every frame of one trajectory. */
export interface TrajectoryData {
  nFrames: number;
  nRes: number;
  /** ps */
  times: Float32Array;
  /** nm, `xyz[(f * nRes + r) * 3 + c]` */
  xyz: Float32Array;
  residues: ResidueInfo[];
}

export function residueKey(r: ResidueInfo) {
  return `${r.chain}:${r.seq}`;
}

export function residueLabel(r: ResidueInfo) {
  return `${r.name} ${r.seq}${r.chain ? ` ${r.chain}` : ''}`;
}

/** Distance matrix (n x n, symmetric) of one frame, written into `out`. */
export function frameDistances(t: TrajectoryData, frame: number, out?: Float32Array) {
  const n = t.nRes;
  out ??= new Float32Array(n * n);
  const xyz = t.xyz,
    o = frame * n * 3;
  for (let i = 0; i < n; i++) {
    const xi = xyz[o + i * 3],
      yi = xyz[o + i * 3 + 1],
      zi = xyz[o + i * 3 + 2];
    out[i * n + i] = 0;
    for (let j = i + 1; j < n; j++) {
      const dx = xyz[o + j * 3] - xi,
        dy = xyz[o + j * 3 + 1] - yi,
        dz = xyz[o + j * 3 + 2] - zi;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      out[i * n + j] = d;
      out[j * n + i] = d;
    }
  }
  return out;
}

export interface MapStats {
  n: number;
  /** frames accumulated */
  count: number;
  mean: Float32Array;
  std: Float32Array;
  /** fraction of frames with the pair closer than the contact cutoff */
  contact: Float32Array;
  cutoff: number;
}

/** Running sums of a distance map over frames; `finish` turns them into mean, std and contact frequency. */
export class MapAccumulator {
  private sum: Float64Array;
  private sumSq: Float64Array;
  private hits: Uint32Array;
  count = 0;

  constructor(
    readonly n: number,
    readonly cutoff: number,
  ) {
    this.sum = new Float64Array(n * n);
    this.sumSq = new Float64Array(n * n);
    this.hits = new Uint32Array(n * n);
  }

  /** `map`: residue indices of `t` for each of the `n` rows (identity when omitted) */
  add(t: TrajectoryData, frame: number, map?: ArrayLike<number>) {
    const { n, cutoff, sum, sumSq, hits } = this;
    const xyz = t.xyz,
      o = frame * t.nRes * 3;
    for (let i = 0; i < n; i++) {
      const a = o + (map ? map[i] : i) * 3;
      const xi = xyz[a],
        yi = xyz[a + 1],
        zi = xyz[a + 2];
      for (let j = i + 1; j < n; j++) {
        const b = o + (map ? map[j] : j) * 3;
        const dx = xyz[b] - xi,
          dy = xyz[b + 1] - yi,
          dz = xyz[b + 2] - zi;
        const d2 = dx * dx + dy * dy + dz * dz;
        const d = Math.sqrt(d2);
        const k = i * n + j;
        sum[k] += d;
        sumSq[k] += d2;
        if (d < cutoff) hits[k]++;
      }
    }
    this.count++;
  }

  finish(): MapStats {
    const { n, count, cutoff } = this;
    const mean = new Float32Array(n * n),
      std = new Float32Array(n * n),
      contact = new Float32Array(n * n);
    const c = Math.max(1, count);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const k = i * n + j,
          k2 = j * n + i;
        const m = this.sum[k] / c;
        const v = Math.max(0, this.sumSq[k] / c - m * m);
        mean[k] = mean[k2] = m;
        std[k] = std[k2] = Math.sqrt(v);
        contact[k] = contact[k2] = this.hits[k] / c;
      }
      contact[i * n + i] = 1;
    }
    return { n, count, mean, std, contact, cutoff };
  }
}

/** Mean, std and contact frequency of the distance map over the frames `[first, last]` of a trajectory, every `stride`. */
export function mapStats(
  t: TrajectoryData,
  cutoff: number,
  options?: { first?: number; last?: number; stride?: number; filter?: (frame: number) => boolean },
) {
  const acc = new MapAccumulator(t.nRes, cutoff);
  const first = options?.first ?? 0,
    last = Math.min(options?.last ?? t.nFrames - 1, t.nFrames - 1),
    stride = Math.max(1, options?.stride ?? 1);
  for (let f = first; f <= last; f += stride) {
    if (options?.filter && !options.filter(f)) continue;
    acc.add(t, f);
  }
  return acc.finish();
}

/** Pools already computed stats (same residues) weighted by their frame counts. */
export function poolStats(stats: MapStats[]): MapStats {
  const n = stats[0].n,
    cutoff = stats[0].cutoff;
  const mean = new Float32Array(n * n),
    std = new Float32Array(n * n),
    contact = new Float32Array(n * n);
  let count = 0;
  for (const s of stats) count += s.count;
  const c = Math.max(1, count);
  for (let k = 0; k < n * n; k++) {
    let m = 0,
      sq = 0,
      ct = 0;
    for (const s of stats) {
      m += s.mean[k] * s.count;
      sq += (s.std[k] * s.std[k] + s.mean[k] * s.mean[k]) * s.count;
      ct += s.contact[k] * s.count;
    }
    m /= c;
    mean[k] = m;
    std[k] = Math.sqrt(Math.max(0, sq / c - m * m));
    contact[k] = ct / c;
  }
  return { n, count, mean, std, contact, cutoff };
}

/** RMS difference of two distance maps over pairs at least `minSep` residues apart (dRMSD, nm). */
export function mapDrmsd(a: ArrayLike<number>, b: ArrayLike<number>, n: number, minSep = 3) {
  let s = 0,
    c = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + minSep; j < n; j++) {
      const d = a[i * n + j] - b[i * n + j];
      s += d * d;
      c++;
    }
  }
  return c ? Math.sqrt(s / c) : 0;
}

/**
 * dRMSD of every frame to a reference distance map. `map` gives the trajectory residue of each reference row, for a
 * reference with fewer residues than the trajectory (the residues two simulations share).
 */
export function drmsdSeries(
  t: TrajectoryData,
  ref: ArrayLike<number>,
  nRef: number,
  map?: ArrayLike<number>,
  minSep = 3,
) {
  const out = new Float32Array(t.nFrames);
  const xyz = t.xyz;
  for (let f = 0; f < t.nFrames; f++) {
    const o = f * t.nRes * 3;
    let s = 0,
      c = 0;
    for (let i = 0; i < nRef; i++) {
      const a = o + (map ? map[i] : i) * 3;
      const xi = xyz[a],
        yi = xyz[a + 1],
        zi = xyz[a + 2];
      for (let j = i + minSep; j < nRef; j++) {
        const b = o + (map ? map[j] : j) * 3;
        const dx = xyz[b] - xi,
          dy = xyz[b + 1] - yi,
          dz = xyz[b + 2] - zi;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz) - ref[i * nRef + j];
        s += d * d;
        c++;
      }
    }
    out[f] = c ? Math.sqrt(s / c) : 0;
  }
  return out;
}

/** Distance of one residue pair in every frame. */
export function pairSeries(t: TrajectoryData, i: number, j: number) {
  const out = new Float32Array(t.nFrames);
  const xyz = t.xyz,
    n = t.nRes;
  for (let f = 0; f < t.nFrames; f++) {
    const a = (f * n + i) * 3,
      b = (f * n + j) * 3;
    const dx = xyz[a] - xyz[b],
      dy = xyz[a + 1] - xyz[b + 1],
      dz = xyz[a + 2] - xyz[b + 2];
    out[f] = Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  return out;
}

function framePositions(t: TrajectoryData, f: number, into?: MinimizeRmsd.Positions) {
  const n = t.nRes;
  into ??= MinimizeRmsd.Positions.empty(n);
  const x = into.x as Float64Array,
    y = into.y as Float64Array,
    z = into.z as Float64Array;
  const o = f * n * 3;
  for (let i = 0; i < n; i++) {
    x[i] = t.xyz[o + i * 3];
    y[i] = t.xyz[o + i * 3 + 1];
    z[i] = t.xyz[o + i * 3 + 2];
  }
  return into;
}

/** RMSD (nm) of every frame to frame `refFrame` after optimal superposition of the representative atoms. */
export function rmsdSeries(t: TrajectoryData, refFrame = 0) {
  const out = new Float32Array(t.nFrames);
  const a = framePositions(t, refFrame);
  const b = MinimizeRmsd.Positions.empty(t.nRes);
  for (let f = 0; f < t.nFrames; f++) {
    framePositions(t, f, b);
    out[f] = MinimizeRmsd.compute({ a, b }).rmsd;
  }
  return out;
}

/** Radius of gyration (nm) of the representative atoms in every frame. */
export function gyrationSeries(t: TrajectoryData) {
  const out = new Float32Array(t.nFrames);
  const n = t.nRes;
  for (let f = 0; f < t.nFrames; f++) {
    const o = f * n * 3;
    let cx = 0,
      cy = 0,
      cz = 0;
    for (let i = 0; i < n; i++) {
      cx += t.xyz[o + i * 3];
      cy += t.xyz[o + i * 3 + 1];
      cz += t.xyz[o + i * 3 + 2];
    }
    cx /= n;
    cy /= n;
    cz /= n;
    let s = 0;
    for (let i = 0; i < n; i++) {
      const dx = t.xyz[o + i * 3] - cx,
        dy = t.xyz[o + i * 3 + 1] - cy,
        dz = t.xyz[o + i * 3 + 2] - cz;
      s += dx * dx + dy * dy + dz * dz;
    }
    out[f] = Math.sqrt(s / n);
  }
  return out;
}

/** Per-residue RMSF (nm) over the frames from `first` on: frames fitted to the first one, then to their mean. */
export function rmsf(t: TrajectoryData, first = 0) {
  const n = t.nRes;
  const out = new Float32Array(n);
  if (t.nFrames - first < 2) return out;
  let ref = framePositions(t, first);
  const b = MinimizeRmsd.Positions.empty(n);
  const sum = new Float64Array(n * 3),
    sumSq = new Float64Array(n);
  for (let pass = 0; pass < 2; pass++) {
    sum.fill(0);
    sumSq.fill(0);
    let count = 0;
    for (let f = first; f < t.nFrames; f++) {
      framePositions(t, f, b);
      const m = MinimizeRmsd.compute({ a: ref, b }).bTransform;
      for (let i = 0; i < n; i++) {
        const x = b.x[i],
          y = b.y[i],
          z = b.z[i];
        const tx = m[0] * x + m[4] * y + m[8] * z + m[12],
          ty = m[1] * x + m[5] * y + m[9] * z + m[13],
          tz = m[2] * x + m[6] * y + m[10] * z + m[14];
        sum[i * 3] += tx;
        sum[i * 3 + 1] += ty;
        sum[i * 3 + 2] += tz;
        sumSq[i] += tx * tx + ty * ty + tz * tz;
      }
      count++;
    }
    const mean = MinimizeRmsd.Positions.empty(n);
    for (let i = 0; i < n; i++) {
      const mx = sum[i * 3] / count,
        my = sum[i * 3 + 1] / count,
        mz = sum[i * 3 + 2] / count;
      (mean.x as Float64Array)[i] = mx;
      (mean.y as Float64Array)[i] = my;
      (mean.z as Float64Array)[i] = mz;
      out[i] = Math.sqrt(Math.max(0, sumSq[i] / count - (mx * mx + my * my + mz * mz)));
    }
    ref = mean;
  }
  return out;
}

/** Pairs in contact in the reference frame (closer than `cutoff`, at least `minSep` apart in sequence). */
export function nativeContacts(t: TrajectoryData, refFrame: number, cutoff: number, minSep = 4) {
  const d = frameDistances(t, refFrame);
  const n = t.nRes;
  const pairs: number[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + minSep; j < n; j++) if (d[i * n + j] < cutoff) pairs.push(i, j);
  }
  return Int32Array.from(pairs);
}

/** Fraction of the native contacts kept in every frame (a pair counts while closer than `cutoff * tolerance`). */
export function nativeContactSeries(t: TrajectoryData, pairs: Int32Array, cutoff: number, tolerance = 1.2) {
  const out = new Float32Array(t.nFrames);
  const n = t.nRes,
    lim2 = (cutoff * tolerance) ** 2,
    np = pairs.length / 2;
  for (let f = 0; f < t.nFrames; f++) {
    const o = f * n * 3;
    let kept = 0;
    for (let p = 0; p < np; p++) {
      const a = o + pairs[p * 2] * 3,
        b = o + pairs[p * 2 + 1] * 3;
      const dx = t.xyz[a] - t.xyz[b],
        dy = t.xyz[a + 1] - t.xyz[b + 1],
        dz = t.xyz[a + 2] - t.xyz[b + 2];
      if (dx * dx + dy * dy + dz * dz < lim2) kept++;
    }
    out[f] = np ? kept / np : 1;
  }
  return out;
}

function meanOf(a: ArrayLike<number>, from: number, to: number) {
  let s = 0;
  for (let i = from; i < to; i++) s += a[i];
  return to > from ? s / (to - from) : 0;
}

/** Least-squares slope of `y` against `x` over `[from, to)`. */
function slope(x: ArrayLike<number>, y: ArrayLike<number>, from: number, to: number) {
  const n = to - from;
  if (n < 2) return 0;
  const mx = meanOf(x, from, to),
    my = meanOf(y, from, to);
  let sxy = 0,
    sxx = 0;
  for (let i = from; i < to; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
  }
  return sxx > 0 ? sxy / sxx : 0;
}

/**
 * Two questions, kept apart: how far from the start the run ended up, and whether it is still moving.
 * `holds`: near the start and flat. `shifted`: flat now, but settled away from the start (it moved early on).
 * `drifting`: still moving over the second half. `departed`: the starting shape is lost. `short`: too few frames.
 */
export type StabilityVerdict = 'holds' | 'shifted' | 'drifting' | 'departed' | 'short';

/** fewer analysed frames than this cannot be judged */
export const MinFramesForVerdict = 20;

export interface StabilityThresholds {
  /** nm: dRMSD drift over the second half of the run below this counts as a plateau */
  drift: number;
  /** nm: late dRMSD to the start above this counts as having settled somewhere else than the start */
  shifted: number;
  /** nm: late dRMSD to the start above this counts as having left the starting shape */
  departed: number;
  /** native-contact fraction a held shape keeps */
  contactsHeld: number;
  /** native-contact fraction below which the shape is lost */
  contactsLost: number;
}

export const DefaultStabilityThresholds: StabilityThresholds = {
  drift: 0.03,
  shifted: 0.2,
  departed: 0.4,
  contactsHeld: 0.85,
  contactsLost: 0.7,
};

export interface StabilityResult {
  /** superposed RMSD to the first analysed frame, nm */
  rmsd: Float32Array;
  /** superposition-free dRMSD to the first analysed frame, nm */
  drmsd: Float32Array;
  rg: Float32Array;
  /** fraction of native contacts */
  q: Float32Array;
  rmsf: Float32Array;
  /** first analysed frame (frames before it are equilibration) */
  first: number;
  nativePairs: number;
  /** means over the last fifth of the run */
  rmsdLate: number;
  drmsdLate: number;
  qLate: number;
  rgLate: number;
  /** dRMSD change over the second half of the run (fitted slope times its duration), nm */
  drift: number;
  /** |mean of the last quarter - mean of the third quarter| of dRMSD, nm */
  blockShift: number;
  verdict: StabilityVerdict;
  reasons: string[];
}

/**
 * Stability of one replica against its own start: does the shape it began in hold? The verdict is a heuristic on the
 * superposition-free numbers (dRMSD and native contacts), so a rigid-body tumble or a floppy tail weighs little.
 */
export function evaluateStability(
  t: TrajectoryData,
  options?: { skipPs?: number; cutoff?: number; thresholds?: StabilityThresholds },
): StabilityResult {
  const th = options?.thresholds ?? DefaultStabilityThresholds;
  const cutoff = options?.cutoff ?? 0.8;
  let first = 0;
  const skip = options?.skipPs ?? 0;
  if (skip > 0 && t.nFrames > 1) {
    const t0 = t.times[0];
    while (first < t.nFrames - 2 && t.times[first] - t0 < skip) first++;
  }
  const rmsd = rmsdSeries(t, first);
  const ref = frameDistances(t, first);
  const drmsd = drmsdSeries(t, ref, t.nRes);
  const rg = gyrationSeries(t);
  const pairs = nativeContacts(t, first, cutoff);
  const q = nativeContactSeries(t, pairs, cutoff);
  const fluct = rmsf(t, first);

  const n = t.nFrames,
    span = n - first;
  const late = Math.max(first, n - Math.max(1, Math.round(span / 5)));
  const half = first + Math.floor(span / 2);
  const q3 = first + Math.floor(span / 2),
    q4 = first + Math.floor((span * 3) / 4);
  const drift = slope(t.times, drmsd, half, n) * (t.times[n - 1] - t.times[half]);
  const blockShift = Math.abs(meanOf(drmsd, q4, n) - meanOf(drmsd, q3, q4));
  const drmsdLate = meanOf(drmsd, late, n),
    qLate = meanOf(q, late, n);

  const reasons: string[] = [];
  let verdict: StabilityVerdict = 'holds';
  if (drmsdLate > th.departed || qLate < th.contactsLost) {
    verdict = 'departed';
    if (drmsdLate > th.departed) reasons.push(`late dRMSD to the start ${drmsdLate.toFixed(2)} nm > ${th.departed} nm`);
    if (qLate < th.contactsLost) reasons.push(`only ${(qLate * 100).toFixed(0)}% of the starting contacts kept`);
  } else if (Math.abs(drift) > th.drift || blockShift > th.drift) {
    verdict = 'drifting';
    if (Math.abs(drift) > th.drift)
      reasons.push(`dRMSD still moving: ${drift >= 0 ? '+' : ''}${drift.toFixed(3)} nm over the second half`);
    if (blockShift > th.drift) reasons.push(`last two quarters differ by ${blockShift.toFixed(3)} nm`);
    if (drmsdLate > th.shifted)
      reasons.push(`already ${drmsdLate.toFixed(2)} nm from the start (more than ${th.shifted} nm)`);
    if (qLate < th.contactsHeld) reasons.push(`${(qLate * 100).toFixed(0)}% of the starting contacts kept`);
  } else if (drmsdLate > th.shifted || qLate < th.contactsHeld) {
    // flat now, but not where it started: it moved early and settled; a flat second half alone is not "holds"
    verdict = 'shifted';
    if (drmsdLate > th.shifted)
      reasons.push(
        `settled ${drmsdLate.toFixed(2)} nm dRMSD from the start (more than ${th.shifted} nm), flat over the second half`,
      );
    if (qLate < th.contactsHeld) reasons.push(`${(qLate * 100).toFixed(0)}% of the starting contacts kept`);
  } else {
    reasons.push(
      `dRMSD plateau at ${drmsdLate.toFixed(2)} nm, ${(qLate * 100).toFixed(0)}% of the starting contacts kept`,
    );
  }
  if (span < MinFramesForVerdict) {
    verdict = 'short';
    reasons.length = 0;
    reasons.push(`only ${span} frame${span === 1 ? '' : 's'} to analyse: too short to judge`);
  }

  return {
    rmsd,
    drmsd,
    rg,
    q,
    rmsf: fluct,
    first,
    nativePairs: pairs.length / 2,
    rmsdLate: meanOf(rmsd, late, n),
    drmsdLate,
    qLate,
    rgLate: meanOf(rg, late, n),
    drift,
    blockShift,
    verdict,
    reasons,
  };
}

export interface ReplicaAgreement {
  n: number;
  /** std of the replica means of every distance, nm: where the replicas sat in different shapes */
  betweenStd: Float32Array;
  /** dRMSD between the mean maps of every replica pair, nm (`r x r`) */
  pairDrmsd: Float32Array;
  maxPairDrmsd: number;
  /** share of residue pairs whose replica means differ by more than the spread inside the replicas */
  disagreeingShare: number;
  agree: boolean;
}

/** Do the replicas of one simulation agree? Compares their mean distance maps. */
export function replicaAgreement(stats: MapStats[], tolerance = 0.1): ReplicaAgreement {
  const n = stats[0].n,
    r = stats.length;
  const betweenStd = new Float32Array(n * n);
  let disagreeing = 0,
    pairs = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const k = i * n + j;
      let m = 0;
      for (const s of stats) m += s.mean[k];
      m /= r;
      let v = 0,
        within = 0;
      for (const s of stats) {
        v += (s.mean[k] - m) ** 2;
        within += s.std[k] * s.std[k];
      }
      const b = Math.sqrt(v / Math.max(1, r - 1));
      betweenStd[k] = betweenStd[j * n + i] = b;
      if (j - i >= 3) {
        pairs++;
        if (b > Math.max(tolerance, 2 * Math.sqrt(within / r))) disagreeing++;
      }
    }
  }
  const pairDrmsd = new Float32Array(r * r);
  let maxPairDrmsd = 0;
  for (let a = 0; a < r; a++) {
    for (let b = a + 1; b < r; b++) {
      const d = mapDrmsd(stats[a].mean, stats[b].mean, n);
      pairDrmsd[a * r + b] = pairDrmsd[b * r + a] = d;
      if (d > maxPairDrmsd) maxPairDrmsd = d;
    }
  }
  const disagreeingShare = pairs ? disagreeing / pairs : 0;
  return {
    n,
    betweenStd,
    pairDrmsd,
    maxPairDrmsd,
    disagreeingShare,
    agree: maxPairDrmsd < tolerance && disagreeingShare < 0.01,
  };
}

export interface SwitchPair {
  i: number;
  j: number;
  meanA: number;
  meanB: number;
  stdA: number;
  stdB: number;
  /** |meanA - meanB|, nm */
  delta: number;
  /** 0..1: the share of this distance's variance that is between the two groups */
  betweenShare: number;
  /** delta * betweenShare */
  score: number;
  /** in contact in one group and not in the other */
  contactSwitch: boolean;
}

/** The residue pairs whose distance tells two groups of frames (conformers, states) apart, best first. */
export function switchPairs(a: MapStats, b: MapStats, options?: { top?: number; minSep?: number; cutoff?: number }) {
  const n = a.n,
    minSep = options?.minSep ?? 3,
    cutoff = options?.cutoff ?? a.cutoff;
  const out: SwitchPair[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = i + minSep; j < n; j++) {
      const k = i * n + j;
      const delta = Math.abs(a.mean[k] - b.mean[k]);
      if (delta < 0.05) continue;
      const between = (delta * delta) / 4;
      const within = (a.std[k] * a.std[k] + b.std[k] * b.std[k]) / 2;
      const betweenShare = between / (between + within || 1);
      out.push({
        i,
        j,
        meanA: a.mean[k],
        meanB: b.mean[k],
        stdA: a.std[k],
        stdB: b.std[k],
        delta,
        betweenShare,
        score: delta * betweenShare,
        contactSwitch: a.mean[k] < cutoff !== b.mean[k] < cutoff,
      });
    }
  }
  out.sort((x, y) => y.score - x.score);
  return out.slice(0, options?.top ?? 50);
}

/** Per-pair separation of two groups as a map: upper values `betweenShare`, for display. */
export function separationMap(a: MapStats, b: MapStats) {
  const n = a.n;
  const out = new Float32Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const k = i * n + j;
      const delta = a.mean[k] - b.mean[k];
      const between = (delta * delta) / 4;
      const within = (a.std[k] * a.std[k] + b.std[k] * b.std[k]) / 2;
      out[k] = out[j * n + i] = between / (between + within || 1);
    }
  }
  return out;
}

// conformers

/** A frame of one of the trajectories under analysis. */
export interface FrameRef {
  traj: number;
  frame: number;
}

export interface ConformerOptions {
  /** frames sampled per trajectory for the clustering itself (the rest are assigned afterwards) */
  framesPerTrajectory: number;
  /** residues used for the fingerprint (more is slower; the rest are skipped evenly) */
  maxResidues: number;
  /** 0: automatic (every split of the tree that has a gap between its two groups) */
  nClusters: number;
  /** ps skipped at the start of every trajectory */
  skipPs: number;
  /** contact cutoff of the per-conformer maps, nm */
  cutoff: number;
  /**
   * Automatic mode: a split of the tree separates two conformers when few frames lie between the two groups. This is
   * the most frames allowed half-way between them, relative to the frames at the (less populated) group itself.
   */
  maxDip: number;
  /**
   * Automatic mode, the scale of a conformer: two conformers are at least this many times further apart than the
   * structure moves from one sampled stretch of a trajectory to the next (the median of that, `ConformerResult.step`).
   * A floppy protein thus needs a bigger change than a stiff one.
   */
  stepFactor: number;
  /** automatic mode: and never closer than this; nm of dRMSD (3Di: 0.6 of it, as share of residues) */
  minSeparation: number;
  /** automatic mode: a group smaller than this share of one trajectory's sampled frames is not a conformer of its own */
  minShare: number;
  /** automatic mode: a group without a gap to a bigger one is in transit to it when it is below this share of both */
  transitShare: number;
  /** automatic mode: never more conformers than this */
  maxClusters: number;
  /** what frames are compared by: their Cα distance maps (dRMSD), or their 3Di structural letters (share that differ) */
  fingerprint: 'distance' | '3di';
}

export const DefaultConformerOptions: ConformerOptions = {
  framesPerTrajectory: 150,
  maxResidues: 110,
  nClusters: 0,
  skipPs: 0,
  cutoff: 0.8,
  maxDip: 0.5,
  stepFactor: 2.25,
  minSeparation: 0.05,
  minShare: 0.2,
  transitShare: 0.2,
  maxClusters: 8,
  fingerprint: 'distance',
};

export interface Conformer {
  id: number;
  /** frames assigned, of all trajectories */
  size: number;
  /** frames assigned per trajectory */
  perTrajectory: number[];
  /** the most typical sampled frame */
  medoid: FrameRef;
  /** mean dRMSD of the sampled members to the medoid, nm */
  spread: number;
  stats: MapStats;
}

export interface ConformerResult {
  options: ConformerOptions;
  /** how many conformers the automatic search finds (also reported when a number was asked for) */
  suggested: number;
  /** the splits of the tree the automatic search looked at, from the top down */
  splits: SplitTest[];
  /** median distance between one sampled stretch of a trajectory and the next: how much the structure moves anyway */
  step: number;
  /** the smallest distance between two groups that counts as two conformers (automatic mode) */
  scale: number;
  /** share of the sampled frames on the way between conformers (no gap to where they lead): not a conformer */
  transit: number;
  /** how many conformers the automatic search finds at a finer and at a coarser scale: how firm `suggested` is */
  nearby: { scale: number; count: number }[];
  /** highest dendrogram merges, nm, largest first */
  merges: number[];
  conformers: Conformer[];
  /** conformer of every frame of every trajectory (-1: skipped as equilibration) */
  assignment: Int16Array[];
  /** the sampled frames, their conformer and a 2D map of them (classical MDS of the dRMSD matrix) */
  sample: { ref: FrameRef; conformer: number; x: number; y: number }[];
  /** dRMSD between conformer means, nm (`k x k`) */
  betweenDrmsd: Float32Array;
  /** the top of the clustering tree, for an icicle plot: `tree[0]` is the root */
  tree: TreeNode[];
  /**
   * the merge height at which the tree was cut into the conformers (between the last merge kept apart and the next);
   * undefined when the automatic search kept splits at different heights
   */
  cutHeight?: number;
  /** what the merge heights measure: nm of dRMSD, or the share of residues whose 3Di letter differs */
  unit: 'nm' | 'letters';
}

/** Two groups of sampled frames, as the automatic search judged them. */
export interface SplitTest {
  /** sampled frames in the two groups */
  sizes: [number, number];
  /** distance between the two group centres: nm of dRMSD, or the share of residues with another 3Di letter */
  separation: number;
  /**
   * frames at the emptiest place between the groups relative to frames at the densest place of the less populated
   * group: 0 a clean gap, 1 a continuum
   */
  dip: number;
  /**
   * - `apart`: two conformers (far enough apart, with a gap)
   * - `joined`: far enough apart, but frames fill the space between: one motion, counted as one conformer
   * - `transit`: as `joined`, and the first group is so thin that its frames are on the way to the second
   * - `close`: closer than the scale of a conformer, not looked at further
   */
  outcome: 'apart' | 'joined' | 'transit' | 'close';
  /** `outcome === 'apart'` */
  accepted: boolean;
}

/** A cluster of the tree: it exists between its own merge height (where its two children join) and its parent's. */
export interface TreeNode {
  id: number;
  parent: number;
  /** ids of the two clusters it was merged from; none for a cluster not opened further */
  children?: [number, number];
  /** merge height of its children (0 for a single frame) */
  height: number;
  /** merge height of its parent (the root: its own) */
  top: number;
  /** horizontal extent, 0..1, by the number of sampled frames */
  x0: number;
  x1: number;
  /** indices into `sample` */
  members: Int32Array;
  /** index into `sample` of its most typical frame */
  medoid: number;
  /** sampled frames per trajectory */
  perTrajectory: number[];
  /** the conformer all its frames belong to, or -1 when it spans several */
  conformer: number;
}

interface FullNode {
  height: number;
  size: number;
  first: number;
  left?: FullNode;
  right?: FullNode;
}

/** The top `maxLeaves` clusters of the tree (opening the highest merge first), laid out left to right. */
function buildTree(
  merges: { height: number; members: number[][] }[],
  n: number,
  dist: Float32Array,
  label: ArrayLike<number>,
  sample: FrameRef[],
  nTraj: number,
  maxLeaves = 28,
): TreeNode[] {
  const byFirst = new Map<number, FullNode>();
  const leaf = (i: number): FullNode => ({ height: 0, size: 1, first: i });
  let root: FullNode = leaf(0);
  for (const m of merges) {
    const [a, b] = m.members;
    const left = byFirst.get(a[0]) ?? leaf(a[0]),
      right = byFirst.get(b[0]) ?? leaf(b[0]);
    root = { height: m.height, size: left.size + right.size, first: a[0], left, right };
    byFirst.set(a[0], root);
    byFirst.delete(b[0]);
  }
  // open the highest merges
  const open = new Set<FullNode>();
  const frontier = [root];
  while (frontier.length < maxLeaves) {
    let best = -1;
    for (let i = 0; i < frontier.length; i++) {
      if (frontier[i].left && (best < 0 || frontier[i].height > frontier[best].height)) best = i;
    }
    if (best < 0) break;
    const node = frontier[best];
    open.add(node);
    frontier.splice(best, 1, node.left!, node.right!);
  }
  const collect = (node: FullNode, into: number[]) => {
    const stack = [node];
    while (stack.length) {
      const c = stack.pop()!;
      if (c.left) stack.push(c.right!, c.left);
      else into.push(c.first);
    }
    return into;
  };
  const out: TreeNode[] = [];
  const visit = (node: FullNode, parent: number, top: number, x0: number): number => {
    const id = out.length;
    const members = Int32Array.from(collect(node, []));
    let medoid = members[0],
      bestSum = Infinity;
    // the member closest to all others (a subsample of the others keeps this linear for big clusters)
    const step = Math.max(1, Math.floor(members.length / 200));
    for (const a of members) {
      let sum = 0;
      for (let k = 0; k < members.length; k += step) sum += dist[a * n + members[k]];
      if (sum < bestSum) {
        bestSum = sum;
        medoid = a;
      }
    }
    const perTrajectory = new Array<number>(nTraj).fill(0);
    let conformer = label[members[0]];
    for (const m of members) {
      perTrajectory[sample[m].traj]++;
      if (label[m] !== conformer) conformer = -1;
    }
    const x1 = x0 + node.size / n;
    const t: TreeNode = { id, parent, height: node.height, top, x0, x1, members, medoid, perTrajectory, conformer };
    out.push(t);
    if (open.has(node)) {
      const l = visit(node.left!, id, node.height, x0);
      const r = visit(node.right!, id, node.height, x0 + node.left!.size / n);
      t.children = [l, r];
    }
    return id;
  };
  visit(root, -1, root.height, 0);
  return out;
}

/** Mean, std and contact frequency of the distance map over some frames of several trajectories (a tree cluster). */
export function statsOfFrames(
  trajectories: TrajectoryData[],
  frames: FrameRef[],
  cutoff: number,
  maps?: (ArrayLike<number> | undefined)[],
) {
  const n = maps?.[0] ? maps[0].length : trajectories[0].nRes;
  const acc = new MapAccumulator(n, cutoff);
  for (const f of frames) acc.add(trajectories[f.traj], f.frame, maps?.[f.traj]);
  return acc.finish();
}

function fingerprintResidues(n: number, max: number) {
  if (n <= max) return Int32Array.from({ length: n }, (_, i) => i);
  const out = new Int32Array(max);
  for (let i = 0; i < max; i++) out[i] = Math.round((i * (n - 1)) / (max - 1));
  return out;
}

function fingerprint(
  t: TrajectoryData,
  frame: number,
  res: ArrayLike<number>,
  map: ArrayLike<number> | undefined,
  out: Float32Array,
) {
  const m = res.length,
    o = frame * t.nRes * 3;
  let p = 0;
  for (let a = 0; a < m; a++) {
    const ia = o + (map ? map[res[a]] : res[a]) * 3;
    const x = t.xyz[ia],
      y = t.xyz[ia + 1],
      z = t.xyz[ia + 2];
    for (let b = a + 1; b < m; b++) {
      const ib = o + (map ? map[res[b]] : res[b]) * 3;
      const dx = t.xyz[ib] - x,
        dy = t.xyz[ib + 1] - y,
        dz = t.xyz[ib + 2] - z;
      out[p++] = Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
}

/** Average-linkage agglomerative clustering of a full distance matrix; returns the merges, lowest first. */
function averageLinkage(dist: Float32Array, n: number) {
  const d = Float32Array.from(dist);
  const size = new Int32Array(n).fill(1);
  const active = new Uint8Array(n).fill(1);
  const members: number[][] = Array.from({ length: n }, (_, i) => [i]);
  const merges: { height: number; members: number[][] }[] = [];
  // nearest active neighbour of every cluster
  const nn = new Int32Array(n),
    nnD = new Float32Array(n);
  const findNN = (i: number) => {
    let best = -1,
      bd = Infinity;
    for (let j = 0; j < n; j++) {
      if (j === i || !active[j]) continue;
      const v = d[i * n + j];
      if (v < bd) {
        bd = v;
        best = j;
      }
    }
    nn[i] = best;
    nnD[i] = bd;
  };
  for (let i = 0; i < n; i++) findNN(i);
  for (let step = 0; step < n - 1; step++) {
    let a = -1,
      bd = Infinity;
    for (let i = 0; i < n; i++) {
      if (active[i] && nnD[i] < bd) {
        bd = nnD[i];
        a = i;
      }
    }
    const b = nn[a];
    merges.push({ height: bd, members: [members[a], members[b]] });
    // merge b into a
    const sa = size[a],
      sb = size[b];
    for (let k = 0; k < n; k++) {
      if (!active[k] || k === a || k === b) continue;
      const v = (d[a * n + k] * sa + d[b * n + k] * sb) / (sa + sb);
      d[a * n + k] = d[k * n + a] = v;
    }
    active[b] = 0;
    size[a] = sa + sb;
    members[a] = members[a].concat(members[b]);
    members[b] = [];
    findNN(a);
    for (let k = 0; k < n; k++) {
      if (!active[k] || k === a) continue;
      if (nn[k] === a || nn[k] === b) findNN(k);
      else if (d[k * n + a] < nnD[k]) {
        nn[k] = a;
        nnD[k] = d[k * n + a];
      }
    }
  }
  return merges;
}

interface LinkNode {
  members: number[];
  height: number;
  left?: LinkNode;
  right?: LinkNode;
}

/**
 * Is there a gap between two groups of sampled frames? Every frame is projected on the line between the two group
 * centres (0: at the first centre, 1: at the second), from the pair distances alone. Separate conformers leave part of
 * the middle of that line nearly empty; one continuous motion cut in two fills it as densely as the groups themselves.
 *
 * A frame is left out of its own group's centre: in many dimensions any two groups of noise look separated along the
 * line between their own means, because every frame pulls its centre towards itself.
 *
 * The densities are counts in a window sliding along the line: of each group where it is densest (a group made of
 * several sub-states is not dense at its centre), and of all frames where the middle is emptiest. One frame is added
 * to both, so a handful of frames is not taken for a gap.
 */
function splitTest(a: number[], b: number[], sq: (i: number, j: number) => number): SplitTest {
  const na = a.length,
    nb = b.length;
  const rowA = new Float64Array(na + nb),
    rowB = new Float64Array(na + nb);
  const all = a.concat(b);
  let sumA = 0,
    sumB = 0,
    cross = 0;
  for (let p = 0; p < all.length; p++) {
    const i = all[p];
    let ra = 0,
      rb = 0;
    for (let q = 0; q < na; q++) ra += sq(i, a[q]);
    for (let q = 0; q < nb; q++) rb += sq(i, b[q]);
    rowA[p] = ra;
    rowB[p] = rb;
    if (p < na) {
      sumA += ra;
      cross += rb;
    } else sumB += rb;
  }
  // squared radius of a group: half its mean squared pair distance
  const vA = sumA / (2 * na * na),
    vB = sumB / (2 * nb * nb);
  const d2 = cross / (na * nb) - vA - vB;
  const result: SplitTest = {
    sizes: [na, nb],
    separation: Math.sqrt(Math.max(0, d2)),
    dip: 1,
    outcome: 'joined',
    accepted: false,
  };
  if (!(d2 > 1e-12) || na < 2 || nb < 2) return result;
  const t = new Float64Array(na + nb);
  for (let p = 0; p < all.length; p++) {
    let da: number, db: number, dd: number;
    if (p < na) {
      // centre of the first group without this frame
      const v = (sumA - 2 * rowA[p]) / (2 * (na - 1) * (na - 1));
      da = rowA[p] / (na - 1) - v;
      db = rowB[p] / nb - vB;
      dd = (cross - rowB[p]) / ((na - 1) * nb) - v - vB;
    } else {
      const v = (sumB - 2 * rowB[p]) / (2 * (nb - 1) * (nb - 1));
      da = rowA[p] / na - vA;
      db = rowB[p] / (nb - 1) - v;
      dd = (cross - rowA[p]) / (na * (nb - 1)) - vA - v;
    }
    t[p] = dd > 1e-12 ? (da - db + dd) / (2 * dd) : 0.5;
  }
  const W = 0.15;
  // most frames of a (sorted) group within any window
  const peak = (ts: Float64Array) => {
    ts.sort();
    let best = 0;
    for (let lo = 0, hi = 0; hi < ts.length; hi++) {
      while (ts[hi] - ts[lo] > 2 * W) lo++;
      best = Math.max(best, hi - lo + 1);
    }
    return best;
  };
  let valley = Infinity;
  for (let c = 0.3; c < 0.701; c += 0.05) {
    let count = 0;
    for (let p = 0; p < t.length; p++) if (Math.abs(t[p] - c) <= W) count++;
    valley = Math.min(valley, count);
  }
  result.dip = (valley + 1) / (Math.min(peak(t.slice(0, na)), peak(t.slice(na))) + 1);
  return result;
}

/**
 * The automatic conformers: groups of sampled frames that are, pair by pair, far enough apart (`minSeparation`) with a
 * gap between them (`splitTest`), and populated enough (`minSize`).
 *
 * Proteins move on every scale: below any two shapes there are smaller sub-states, also separated by gaps, down to
 * single side chains. `minSeparation` is the scale asked about. So, in two steps:
 *
 * 1. Pieces: the clustering tree is opened from the top down to groups that cannot be halved into two of `minSize`.
 *    (Judging the splits of the tree itself is not enough: one side of a high split is a mixture of shapes, whose
 *    centre is closer to the other side than any of the shapes are, and one unclear split would hide all below it.)
 * 2. The two closest pieces that are not two conformers are joined, again and again, until every pair left is far
 *    enough apart and has a gap. Pieces closer than the scale are one conformer. Pieces further apart but without a
 *    gap are one motion; when one of them is much smaller, its frames are in `transit` instead (all runs start in one
 *    shape and relax into others) and join whatever conformer each is nearest to.
 */
function autoClusters(
  merges: { height: number; members: number[][] }[],
  n: number,
  sq: (i: number, j: number) => number,
  o: { maxDip: number; minSeparation: number; minSize: number; maxClusters: number; transitShare: number },
) {
  const byFirst = new Map<number, LinkNode>();
  let root: LinkNode = { members: [0], height: 0 };
  for (const m of merges) {
    const [a, b] = m.members;
    const left = byFirst.get(a[0]) ?? { members: a, height: 0 },
      right = byFirst.get(b[0]) ?? { members: b, height: 0 };
    root = { members: a.concat(b), height: m.height, left, right };
    byFirst.set(a[0], root);
    byFirst.delete(b[0]);
  }
  interface Group {
    id: number;
    members: number[];
  }
  const strays: number[] = [],
    splits: SplitTest[] = [];
  const groups: Group[] = [];
  let nextId = 0;
  const stack = [root];
  while (stack.length) {
    let node = stack.pop()!;
    for (;;) {
      const { left, right } = node;
      const [small, big] =
        left && right ? (left.members.length < right.members.length ? [left, right] : [right, left]) : [];
      if (!small || !big || big.members.length < o.minSize) {
        groups.push({ id: nextId++, members: node.members });
        break;
      }
      if (small.members.length >= o.minSize) {
        stack.push(right!, left!);
        break;
      }
      // strays: peel them off and look at the rest
      for (const i of small.members) strays.push(i);
      node = big;
    }
  }

  // sums of squared distances between the pieces: they add up when pieces are joined, and give the centre distances
  let cross = groups.map((x) =>
    groups.map((y) => {
      let sum = 0;
      for (const i of x.members) for (const j of y.members) sum += sq(i, j);
      return sum;
    }),
  );
  const separation2 = (a: number, b: number) => {
    const na = groups[a].members.length,
      nb = groups[b].members.length;
    return cross[a][b] / (na * nb) - cross[a][a] / (2 * na * na) - cross[b][b] / (2 * nb * nb);
  };
  const tested = new Map<string, SplitTest>();
  const testOf = (x: Group, y: Group) => {
    const key = `${x.id} ${y.id}`;
    let test = tested.get(key);
    if (!test) {
      test = splitTest(x.members, y.members, sq);
      tested.set(key, test);
    }
    return test;
  };
  const remove = (a: number) => {
    groups.splice(a, 1);
    cross.splice(a, 1);
    for (const row of cross) row.splice(a, 1);
  };
  const min2 = o.minSeparation * o.minSeparation;
  let transit = 0;
  for (let changed = true; changed && groups.length > 1; ) {
    changed = false;
    const pairs: [number, number, number][] = [];
    for (let a = 0; a < groups.length; a++)
      for (let b = a + 1; b < groups.length; b++) pairs.push([separation2(a, b), a, b]);
    pairs.sort((x, y) => x[0] - y[0]);
    for (const [d2, a, b] of pairs) {
      // the smaller of the two first
      const [x, y] = groups[a].members.length <= groups[b].members.length ? [a, b] : [b, a];
      const far = d2 >= min2;
      const test = testOf(groups[x], groups[y]);
      if (far && test.dip <= o.maxDip) continue;
      const nx = groups[x].members.length,
        ny = groups[y].members.length;
      if (far && nx < o.transitShare * (nx + ny)) {
        test.outcome = 'transit';
        transit += nx;
        for (const m of groups[x].members) strays.push(m);
      } else {
        test.outcome = far ? 'joined' : 'close';
        groups[y] = { id: nextId++, members: groups[y].members.concat(groups[x].members) };
        const cyy = cross[y][y] + cross[x][x] + 2 * cross[x][y];
        for (let k = 0; k < groups.length; k++) cross[y][k] = cross[k][y] = cross[y][k] + cross[x][k];
        cross[y][y] = cyy;
      }
      remove(x);
      splits.push(test);
      changed = true;
      break;
    }
  }
  // the cap keeps the most populated
  const order = groups.map((_, i) => i).sort((x, y) => groups[y].members.length - groups[x].members.length);
  const kept = order.slice(0, o.maxClusters).map((i) => groups[i]);
  for (const i of order.slice(o.maxClusters)) for (const m of groups[i].members) strays.push(m);
  cross = [];
  // every pair left is apart
  const apart: SplitTest[] = [];
  for (let i = 0; i < kept.length; i++) {
    for (let j = i + 1; j < kept.length; j++) {
      const test = tested.get(`${kept[j].id} ${kept[i].id}`) ?? testOf(kept[i], kept[j]);
      test.outcome = 'apart';
      test.accepted = true;
      apart.push(test);
    }
  }
  apart.sort((x, y) => x.separation - y.separation);
  const label = new Int32Array(n).fill(-1);
  kept.forEach((g, c) => {
    for (const i of g.members) label[i] = c;
  });
  // strays and frames in transit join the conformer whose centre they are closest to
  if (strays.length && kept.length > 1) {
    const radius = kept.map((g) => {
      let sum = 0;
      for (const i of g.members) for (const j of g.members) sum += sq(i, j);
      return sum / (2 * g.members.length * g.members.length);
    });
    for (const i of strays) {
      let best = 0,
        bd = Infinity;
      kept.forEach((g, c) => {
        let sum = 0;
        for (const j of g.members) sum += sq(i, j);
        const d = sum / g.members.length - radius[c];
        if (d < bd) {
          bd = d;
          best = c;
        }
      });
      label[i] = best;
    }
  } else for (const i of strays) label[i] = 0;
  return { label, count: kept.length, splits: apart.concat(splits), transit: kept.length > 1 ? transit / n : 0 };
}

function cutTree(merges: { height: number; members: number[][] }[], n: number, k: number) {
  const label = new Int32Array(n);
  for (let i = 0; i < n; i++) label[i] = i;
  // apply all but the last k-1 merges
  const upTo = Math.max(0, merges.length - (k - 1));
  for (let m = 0; m < upTo; m++) {
    const [a, b] = merges[m].members;
    const target = label[a[0]];
    for (const i of b) label[i] = target;
    for (const i of a) label[i] = target;
  }
  // relabel 0..k-1 by size, largest first
  const counts = new Map<number, number>();
  for (let i = 0; i < n; i++) counts.set(label[i], (counts.get(label[i]) ?? 0) + 1);
  const order = [...counts.entries()].sort((x, y) => y[1] - x[1]).map((e) => e[0]);
  const rank = new Map(order.map((l, i) => [l, i]));
  for (let i = 0; i < n; i++) label[i] = rank.get(label[i])!;
  return label;
}

/** Top two coordinates of classical MDS, by power iteration on the double-centred squared distances. */
function mds2(dist: Float32Array, n: number) {
  const b = new Float64Array(n * n);
  const rowMean = new Float64Array(n);
  let total = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const v = dist[i * n + j] ** 2;
      b[i * n + j] = v;
      rowMean[i] += v / n;
      total += v / (n * n);
    }
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) b[i * n + j] = -0.5 * (b[i * n + j] - rowMean[i] - rowMean[j] + total);
  }
  const coords: Float64Array[] = [];
  for (let c = 0; c < 2; c++) {
    let v = new Float64Array(n);
    for (let i = 0; i < n; i++) v[i] = Math.sin(i * (c + 1) * 12.9898) + 0.1;
    let lambda = 0;
    for (let it = 0; it < 80; it++) {
      const w = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        let s = 0;
        for (let j = 0; j < n; j++) s += b[i * n + j] * v[j];
        w[i] = s;
      }
      for (const prev of coords) {
        let dot = 0,
          pp = 0;
        for (let i = 0; i < n; i++) {
          dot += w[i] * prev[i];
          pp += prev[i] * prev[i];
        }
        if (pp > 0) for (let i = 0; i < n; i++) w[i] -= (dot / pp) * prev[i];
      }
      let norm = 0;
      for (let i = 0; i < n; i++) norm += w[i] * w[i];
      norm = Math.sqrt(norm);
      if (norm < 1e-12) break;
      lambda = norm;
      for (let i = 0; i < n; i++) w[i] /= norm;
      v = w;
    }
    const scale = Math.sqrt(Math.max(0, lambda));
    coords.push(v.map((x) => x * scale));
  }
  return coords;
}

/**
 * Finds conformers in one or more trajectories of the same residues: frames are fingerprinted by their distance map,
 * a sample of them is clustered by dRMSD (average linkage), and every frame is then assigned to the nearest conformer.
 *
 * `maps[t]`: for trajectory `t`, the residue index of each common residue (when the trajectories differ in residues).
 */
export async function findConformers(
  trajectories: TrajectoryData[],
  options: Partial<ConformerOptions> = {},
  hooks?: {
    maps?: (ArrayLike<number> | undefined)[];
    /** 3Di letters of every trajectory, `letters[frame * nRes + residue]`, for the '3di' fingerprint */
    letters?: (Uint8Array | undefined)[];
    progress?: (message: string, fraction: number) => Promise<void> | void;
  },
): Promise<ConformerResult> {
  const o = { ...DefaultConformerOptions, ...options };
  const maps = hooks?.maps ?? trajectories.map(() => undefined);
  const nCommon = maps[0] ? maps[0].length : trajectories[0].nRes;
  const res = fingerprintResidues(nCommon, o.maxResidues);
  const letters = hooks?.letters;
  const use3Di = o.fingerprint === '3di';
  if (use3Di && !(letters && trajectories.every((_, i) => letters[i])))
    throw new Error('The 3Di fingerprint needs 3Di letters for every simulation compared.');
  // a frame's fingerprint: its distances between the chosen residues, or its 3Di letter at every residue
  const nPairs = use3Di ? nCommon : (res.length * (res.length - 1)) / 2;
  const fill = (ti: number, frame: number, out: Float32Array) => {
    const t = trajectories[ti];
    if (!use3Di) return fingerprint(t, frame, res, maps[ti], out);
    const l = letters![ti]!,
      map = maps[ti],
      o3 = frame * t.nRes;
    for (let r = 0; r < nCommon; r++) out[r] = l[o3 + (map ? map[r] : r)];
  };
  const tick = async (m: string, f: number) => {
    if (hooks?.progress) await hooks.progress(m, f);
  };

  // first analysed frame of each trajectory, and the sample
  const firsts = trajectories.map((t) => {
    let f = 0;
    if (o.skipPs > 0) while (f < t.nFrames - 1 && t.times[f] - t.times[0] < o.skipPs) f++;
    return f;
  });
  const sample: FrameRef[] = [];
  trajectories.forEach((t, ti) => {
    const avail = t.nFrames - firsts[ti];
    const count = Math.min(avail, o.framesPerTrajectory);
    for (let s = 0; s < count; s++) {
      const frame = firsts[ti] + (count > 1 ? Math.round((s * (avail - 1)) / (count - 1)) : 0);
      sample.push({ traj: ti, frame });
    }
  });
  const ns = sample.length;
  if (ns < 2) throw new Error('Not enough frames to look for conformers.');

  // A sampled frame stands for the stretch of the trajectory around it: its fingerprint is the average over that
  // stretch (3Di: the most common letter). Conformers live much longer than one frame, thermal noise does not, so this
  // keeps what tells conformers apart and removes most of what blurs them.
  await tick('Fingerprinting frames', 0);
  const fp = new Float32Array(ns * nPairs);
  {
    const MAX_PER_BLOCK = 16;
    const one = new Float32Array(nPairs);
    const votes = use3Di ? new Uint16Array(nPairs * 128) : undefined;
    let s = 0,
      lastTick = Date.now();
    for (let ti = 0; ti < trajectories.length; ti++) {
      const t = trajectories[ti];
      const avail = t.nFrames - firsts[ti];
      const count = Math.min(avail, o.framesPerTrajectory);
      for (let b = 0; b < count; b++, s++) {
        const centre = sample[s].frame;
        const half = count > 1 ? (avail - 1) / (count - 1) / 2 : 0;
        const from = Math.max(firsts[ti], Math.ceil(centre - half)),
          to = Math.min(t.nFrames - 1, Math.floor(centre + half));
        const len = Math.max(1, to - from + 1);
        const take = Math.min(len, MAX_PER_BLOCK);
        const out = fp.subarray(s * nPairs, (s + 1) * nPairs);
        if (take === 1) {
          fill(ti, centre, out);
          continue;
        }
        votes?.fill(0);
        for (let q = 0; q < take; q++) {
          const frame = from + Math.round((q * (len - 1)) / (take - 1));
          fill(ti, frame, one);
          if (votes) for (let p = 0; p < nPairs; p++) votes[p * 128 + (one[p] & 127)]++;
          else for (let p = 0; p < nPairs; p++) out[p] += one[p];
        }
        if (votes) {
          for (let p = 0; p < nPairs; p++) {
            let best = 0;
            for (let l = 1; l < 128; l++) if (votes[p * 128 + l] > votes[p * 128 + best]) best = l;
            out[p] = best;
          }
        } else for (let p = 0; p < nPairs; p++) out[p] /= take;
        if (Date.now() - lastTick > 60) {
          await tick('Fingerprinting frames', s / ns);
          lastTick = Date.now();
        }
      }
    }
  }

  const dist = new Float32Array(ns * ns);
  let last = Date.now();
  for (let a = 0; a < ns; a++) {
    const oa = a * nPairs;
    for (let b = a + 1; b < ns; b++) {
      const ob = b * nPairs;
      let s = 0;
      if (use3Di) {
        // the share of residues whose letter differs
        for (let p = 0; p < nPairs; p++) if (fp[oa + p] !== fp[ob + p]) s++;
        dist[a * ns + b] = dist[b * ns + a] = s / nPairs;
        continue;
      }
      for (let p = 0; p < nPairs; p++) {
        const d = fp[oa + p] - fp[ob + p];
        s += d * d;
      }
      dist[a * ns + b] = dist[b * ns + a] = Math.sqrt(s / nPairs);
    }
    if (Date.now() - last > 60) {
      await tick('Comparing frames', (a / ns) * (2 - a / ns));
      last = Date.now();
    }
  }

  await tick('Clustering', 1);
  const merges = averageLinkage(dist, ns);
  const heights = merges.map((m) => m.height).reverse();
  // squared distances in the space the fingerprints live in (the share of differing letters already is one)
  const sq = use3Di ? (i: number, j: number) => dist[i * ns + j] : (i: number, j: number) => dist[i * ns + j] ** 2;
  // how far the structure moves from one sampled stretch to the next, within a trajectory
  const steps: number[] = [];
  for (let a = 1; a < ns; a++) if (sample[a].traj === sample[a - 1].traj) steps.push(dist[a * ns + a - 1]);
  steps.sort((x, y) => x - y);
  const step = steps.length ? steps[steps.length >> 1] : 0;
  const scale = Math.max(use3Di ? o.minSeparation * 0.6 : o.minSeparation, o.stepFactor * step);
  const search = (at: number) =>
    autoClusters(merges, ns, sq, {
      maxDip: o.maxDip,
      // the share of differing letters is a squared distance: `separation` is its root
      minSeparation: use3Di ? Math.sqrt(at) : at,
      minSize: Math.max(8, Math.ceil((ns / trajectories.length) * o.minShare)),
      maxClusters: o.maxClusters,
      transitShare: o.transitShare,
    });
  const auto = search(scale);
  const nearby = [0.8, 1.25].map((f) => ({ scale: scale * f, count: search(scale * f).count }));
  if (use3Di) for (const t of auto.splits) t.separation = t.separation ** 2;
  const suggested = auto.count;
  const forced = o.nClusters > 0 && o.nClusters !== suggested;
  const k = forced ? Math.max(1, Math.min(o.nClusters, ns)) : suggested;
  const label = forced ? cutTree(merges, ns, k) : auto.label;
  // one height cuts the tree into the conformers only when they are what a flat cut gives
  let flat = forced;
  if (!forced && k > 1) {
    const cut = cutTree(merges, ns, k);
    const pairing = new Map<number, number>();
    flat = true;
    for (let s = 0; s < ns && flat; s++) {
      if (!pairing.has(label[s])) pairing.set(label[s], cut[s]);
      flat = pairing.get(label[s]) === cut[s];
    }
    flat &&= new Set(pairing.values()).size === k;
  }

  // conformer mean fingerprints, medoids
  const centre = new Float32Array(k * nPairs);
  const counts = new Int32Array(k);
  for (let s = 0; s < ns; s++) {
    const c = label[s];
    counts[c]++;
    for (let p = 0; p < nPairs; p++) centre[c * nPairs + p] += fp[s * nPairs + p];
  }
  for (let c = 0; c < k; c++) for (let p = 0; p < nPairs; p++) centre[c * nPairs + p] /= Math.max(1, counts[c]);
  if (use3Di) {
    // letters have no mean: a conformer's centre is its most common letter at every residue
    const votes = new Int32Array(128);
    for (let c = 0; c < k; c++) {
      for (let p = 0; p < nPairs; p++) {
        votes.fill(0);
        let best = 0;
        for (let q = 0; q < ns; q++) {
          if (label[q] !== c) continue;
          const l = fp[q * nPairs + p] & 127;
          if (++votes[l] > votes[best]) best = l;
        }
        centre[c * nPairs + p] = best;
      }
    }
  }
  const medoids: number[] = [],
    spreads: number[] = [];
  for (let c = 0; c < k; c++) {
    let best = -1,
      bs = Infinity;
    for (let a = 0; a < ns; a++) {
      if (label[a] !== c) continue;
      let s = 0;
      for (let b = 0; b < ns; b++) if (label[b] === c) s += dist[a * ns + b];
      if (s < bs) {
        bs = s;
        best = a;
      }
    }
    medoids.push(best);
    spreads.push(counts[c] > 1 ? bs / (counts[c] - 1) : 0);
  }

  // assign every frame to the nearest conformer centre, and accumulate full-resolution maps per conformer
  await tick('Assigning all frames', 1);
  const accs = Array.from({ length: k }, () => new MapAccumulator(nCommon, o.cutoff));
  const assignment: Int16Array[] = [];
  const perTraj = Array.from({ length: k }, () => new Array<number>(trajectories.length).fill(0));
  const buf = new Float32Array(nPairs);
  for (let ti = 0; ti < trajectories.length; ti++) {
    const t = trajectories[ti];
    const a = new Int16Array(t.nFrames).fill(-1);
    for (let f = firsts[ti]; f < t.nFrames; f++) {
      let c = 0;
      if (k > 1) {
        fill(ti, f, buf);
        let bd = Infinity;
        for (let cc = 0; cc < k; cc++) {
          let s = 0;
          const oc = cc * nPairs;
          if (use3Di) {
            for (let p = 0; p < nPairs; p++) if (buf[p] !== centre[oc + p]) s++;
          } else {
            for (let p = 0; p < nPairs; p++) {
              const d = buf[p] - centre[oc + p];
              s += d * d;
            }
          }
          if (s < bd) {
            bd = s;
            c = cc;
          }
        }
      }
      a[f] = c;
      perTraj[c][ti]++;
      accs[c].add(t, f, maps[ti]);
      if (Date.now() - last > 60) {
        await tick('Assigning all frames', (ti + f / t.nFrames) / trajectories.length);
        last = Date.now();
      }
    }
    assignment.push(a);
  }

  const conformers: Conformer[] = accs.map((acc, c) => ({
    id: c,
    size: acc.count,
    perTrajectory: perTraj[c],
    medoid: sample[medoids[c]],
    spread: spreads[c],
    stats: acc.finish(),
  }));
  const betweenDrmsd = new Float32Array(k * k);
  for (let a = 0; a < k; a++) {
    for (let b = a + 1; b < k; b++) {
      betweenDrmsd[a * k + b] = betweenDrmsd[b * k + a] = mapDrmsd(
        conformers[a].stats.mean,
        conformers[b].stats.mean,
        nCommon,
      );
    }
  }
  const [mx, my] = mds2(dist, ns);
  return {
    options: o,
    suggested,
    merges: heights.slice(0, 8),
    splits: auto.splits,
    step,
    scale,
    transit: forced ? 0 : auto.transit,
    nearby,
    conformers,
    assignment,
    sample: sample.map((ref, s) => ({ ref, conformer: label[s], x: mx[s], y: my[s] })),
    betweenDrmsd,
    tree: buildTree(merges, ns, dist, label, sample, trajectories.length),
    // heights[k - 2] is the last merge kept apart, heights[k - 1] the first one made
    cutHeight: k === 1 ? heights[0] * 1.02 : flat ? (heights[k - 2] + (heights[k - 1] ?? 0)) / 2 : undefined,
    unit: use3Di ? 'letters' : 'nm',
  };
}

/** Residues present in every trajectory (by chain and number), with their index in each. */
export function commonResidues(trajectories: TrajectoryData[]) {
  const first = trajectories[0];
  const indexMaps = trajectories.map((t) => new Map(t.residues.map((r, i) => [residueKey(r), i])));
  const keep: number[] = [];
  first.residues.forEach((r, i) => {
    if (indexMaps.every((m) => m.has(residueKey(r)))) keep.push(i);
  });
  const identical = trajectories.every((t) => t.nRes === keep.length);
  const maps = trajectories.map((_, ti) =>
    Int32Array.from(keep, (i) => indexMaps[ti].get(residueKey(first.residues[i]))!),
  );
  return { residues: keep.map((i) => first.residues[i]), maps, identical };
}
