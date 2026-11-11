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

/** `short`: too few frames to judge */
export type StabilityVerdict = 'holds' | 'drifting' | 'departed' | 'short';

/** fewer analysed frames than this cannot be judged */
export const MinFramesForVerdict = 20;

export interface StabilityThresholds {
  /** nm: dRMSD drift over the second half of the run below this counts as a plateau */
  drift: number;
  /** nm: late dRMSD to the start above this counts as having left the starting shape */
  departed: number;
  /** native-contact fraction a held shape keeps */
  contactsHeld: number;
  /** native-contact fraction below which the shape is lost */
  contactsLost: number;
}

export const DefaultStabilityThresholds: StabilityThresholds = {
  drift: 0.03,
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
  } else if (Math.abs(drift) > th.drift || blockShift > th.drift || qLate < th.contactsHeld) {
    verdict = 'drifting';
    if (Math.abs(drift) > th.drift)
      reasons.push(`dRMSD still moving: ${drift >= 0 ? '+' : ''}${drift.toFixed(3)} nm over the second half`);
    if (blockShift > th.drift) reasons.push(`last two quarters differ by ${blockShift.toFixed(3)} nm`);
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
  /** 0: automatic from the merge heights */
  nClusters: number;
  /** ps skipped at the start of every trajectory */
  skipPs: number;
  /** contact cutoff of the per-conformer maps, nm */
  cutoff: number;
  /** a merge this many times higher than the next one separates conformers (automatic mode) */
  gapRatio: number;
  /** what frames are compared by: their Cα distance maps (dRMSD), or their 3Di structural letters (share that differ) */
  fingerprint: 'distance' | '3di';
}

export const DefaultConformerOptions: ConformerOptions = {
  framesPerTrajectory: 150,
  maxResidues: 110,
  nClusters: 0,
  skipPs: 0,
  cutoff: 0.8,
  gapRatio: 1.6,
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
  /** how many conformers the merge heights suggest */
  suggested: number;
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
  /** the merge height at which the tree was cut into the conformers (between the last merge kept apart and the next) */
  cutHeight: number;
  /** what the merge heights measure: nm of dRMSD, or the share of residues whose 3Di letter differs */
  unit: 'nm' | 'letters';
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

  await tick('Fingerprinting frames', 0);
  const fp = new Float32Array(ns * nPairs);
  for (let s = 0; s < ns; s++) {
    const { traj, frame } = sample[s];
    fill(traj, frame, fp.subarray(s * nPairs, (s + 1) * nPairs));
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
  // automatic count: the largest relative gap among the top merges
  let suggested = 1,
    bestRatio = o.gapRatio;
  for (let i = 0; i < Math.min(7, heights.length - 1); i++) {
    const ratio = heights[i] / Math.max(1e-6, heights[i + 1]);
    if (ratio > bestRatio) {
      bestRatio = ratio;
      suggested = i + 2;
    }
  }
  const k = Math.max(1, Math.min(o.nClusters > 0 ? o.nClusters : suggested, ns));
  const label = cutTree(merges, ns, k);

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
    conformers,
    assignment,
    sample: sample.map((ref, s) => ({ ref, conformer: label[s], x: mx[s], y: my[s] })),
    betweenDrmsd,
    tree: buildTree(merges, ns, dist, label, sample, trajectories.length),
    // heights[k - 2] is the last merge kept apart, heights[k - 1] the first one made
    cutHeight: k > 1 ? (heights[k - 2] + (heights[k - 1] ?? 0)) / 2 : heights[0] * 1.02,
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
