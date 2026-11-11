/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

/**
 * A trajectory in the 3Di structural alphabet (Foldseek): one letter per residue and frame for the residue's local
 * 3D environment. The frames are then rows of an alignment, and alignment measures apply: how conserved a column is
 * over time, which columns differ between two sets of frames. The letters come from Foldseek itself (the companion
 * server runs it); nothing here encodes structures.
 */

export interface ThreeDi {
  nFrames: number;
  nRes: number;
  /** character codes, `letters[frame * nRes + residue]` */
  letters: Uint8Array;
}

export interface ThreeDiSummary {
  /** Shannon entropy of every column over the analysed frames, bits (0: one letter throughout) */
  entropy: Float32Array;
  /** the most common letter of every column (character code), and the share of frames that have it */
  consensus: Uint8Array;
  consensusShare: Float32Array;
  /** share of residues whose letter differs from the first analysed frame's, per frame */
  changed: Float32Array;
  first: number;
}

/**
 * Reads the server's 3Di file: one line of letters per frame. `pick` selects the residues analysed out of the
 * residues Foldseek saw (all of them when the counts agree).
 */
export function parseThreeDi(
  text: string,
  nFrames: number,
  nRes: number,
  totalResidues: number,
  allIndex: ArrayLike<number>,
): ThreeDi {
  const lines = text.split(/\r?\n/).filter((l) => l.length > 0);
  if (lines.length !== nFrames) throw new Error(`The 3Di file has ${lines.length} frames, the trajectory ${nFrames}.`);
  const len = lines[0].length;
  let pick: ArrayLike<number>;
  if (len === nRes) pick = Int32Array.from({ length: nRes }, (_, i) => i);
  else if (len === totalResidues) pick = allIndex;
  else
    throw new Error(
      `The 3Di strings have ${len} letters, the structure ${nRes} residues with a Cα (${totalResidues} in all).`,
    );
  const letters = new Uint8Array(nFrames * nRes);
  for (let f = 0; f < nFrames; f++) {
    const l = lines[f];
    if (l.length !== len) throw new Error(`Frame ${f + 1} of the 3Di file has ${l.length} letters, not ${len}.`);
    for (let r = 0; r < nRes; r++) letters[f * nRes + r] = l.charCodeAt(pick[r]) & 127;
  }
  return { nFrames, nRes, letters };
}

export function summarizeThreeDi(t: ThreeDi, first = 0): ThreeDiSummary {
  const { nFrames, nRes, letters } = t;
  const entropy = new Float32Array(nRes),
    consensus = new Uint8Array(nRes),
    consensusShare = new Float32Array(nRes),
    changed = new Float32Array(nFrames);
  const counts = new Int32Array(128);
  const n = Math.max(1, nFrames - first);
  for (let r = 0; r < nRes; r++) {
    counts.fill(0);
    for (let f = first; f < nFrames; f++) counts[letters[f * nRes + r]]++;
    let h = 0,
      best = 0;
    for (let c = 0; c < 128; c++) {
      if (!counts[c]) continue;
      const p = counts[c] / n;
      h -= p * Math.log2(p);
      if (counts[c] > counts[best]) best = c;
    }
    entropy[r] = h;
    consensus[r] = best;
    consensusShare[r] = counts[best] / n;
  }
  for (let f = 0; f < nFrames; f++) {
    let d = 0;
    for (let r = 0; r < nRes; r++) if (letters[f * nRes + r] !== letters[first * nRes + r]) d++;
    changed[f] = d / nRes;
  }
  return { entropy, consensus, consensusShare, changed, first };
}

export interface LocalSwitch {
  /** residue index (in the common residues of the frames compared) */
  residue: number;
  /** most common letter in each set and the share of its frames that have it */
  letterA: string;
  shareA: number;
  letterB: string;
  shareB: number;
  /** Jensen-Shannon divergence of the two letter distributions, 0 (same) to 1 (no letter in common) */
  divergence: number;
}

/**
 * Residues whose local structure differs between two sets of frames: the columns whose letter distributions are
 * furthest apart. `frames`: for each set, the letters of each frame as (3Di, frame index, residue map) triples.
 */
export function localSwitches(
  a: { t: ThreeDi; frame: number; map?: ArrayLike<number> }[],
  b: { t: ThreeDi; frame: number; map?: ArrayLike<number> }[],
  nRes: number,
  top = 15,
): LocalSwitch[] {
  const ca = new Float64Array(128),
    cb = new Float64Array(128);
  const out: LocalSwitch[] = [];
  if (!a.length || !b.length) return out;
  for (let r = 0; r < nRes; r++) {
    ca.fill(0);
    cb.fill(0);
    for (const f of a) ca[f.t.letters[f.frame * f.t.nRes + (f.map ? f.map[r] : r)]] += 1 / a.length;
    for (const f of b) cb[f.t.letters[f.frame * f.t.nRes + (f.map ? f.map[r] : r)]] += 1 / b.length;
    let js = 0,
      ba = 0,
      bb = 0;
    for (let c = 0; c < 128; c++) {
      const p = ca[c],
        q = cb[c];
      if (p > ca[ba]) ba = c;
      if (q > cb[bb]) bb = c;
      const m = (p + q) / 2;
      if (p > 0) js += 0.5 * p * Math.log2(p / m);
      if (q > 0) js += 0.5 * q * Math.log2(q / m);
    }
    if (js < 0.05) continue;
    out.push({
      residue: r,
      letterA: String.fromCharCode(ba),
      shareA: ca[ba],
      letterB: String.fromCharCode(bb),
      shareB: cb[bb],
      divergence: js,
    });
  }
  out.sort((x, y) => y.divergence - x.divergence);
  return out.slice(0, top);
}
