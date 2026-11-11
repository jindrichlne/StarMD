/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useMemo, useState } from 'react';
import { oneLetter } from '../align.js';
import { type MapStats, residueLabel } from '../analysis.js';
import { cssGradient, type ScaleName, Scales } from '../colors.js';
import type { MdContactsSession, Simulation } from '../session.js';

/** contacts closer in sequence than this are the chain itself, not counted */
const MIN_SEP = 4;

function contactCounts(map: ArrayLike<number>, n: number, cutoff: number) {
  const out = new Uint16Array(n);
  for (let i = 0; i < n; i++) {
    for (let j = i + MIN_SEP; j < n; j++) {
      if (map[i * n + j] < cutoff) {
        out[i]++;
        out[j]++;
      }
    }
  }
  return out;
}

function shade(scale: ScaleName, t: number) {
  const lut = Scales[scale];
  const k = Math.round(Math.max(0, Math.min(1, t)) * 255) * 3;
  const [r, g, b] = [lut[k], lut[k + 1], lut[k + 2]];
  // dark fills take white letters
  const dark = 0.299 * r + 0.587 * g + 0.114 * b < 140;
  return { background: `rgb(${r},${g},${b})`, color: dark ? '#fff' : '#0b0b0b' };
}

type Mode = 'now' | 'change' | '3di' | '3di-var';

/**
 * The sequence with the contacts drawn on it: every residue shaded by how many contacts it has in the frame shown (or
 * how many it gained or lost since the start); pin a residue to see its partners along the sequence.
 */
export function SequenceView({
  session,
  sim,
  frameMap,
  startMap,
  stats,
}: {
  session: MdContactsSession;
  sim: Simulation;
  frameMap: ArrayLike<number>;
  startMap: ArrayLike<number>;
  stats: MapStats;
}) {
  const n = sim.data.nRes;
  const cutoff = stats.cutoff;
  const [mode, setMode] = useState<Mode>('now');
  const [pinnedKey, setPinned] = useState<number | undefined>();
  const [hover, setHover] = useState<number | undefined>();
  const pinned = pinnedKey !== undefined && pinnedKey < n ? pinnedKey : undefined;
  const pair = session.pair.value;

  const now = useMemo(() => contactCounts(frameMap, n, cutoff), [frameMap, n, cutoff]);
  const start = useMemo(() => contactCounts(startMap, n, cutoff), [startMap, n, cutoff]);
  const maxCount = useMemo(() => Math.max(4, ...start), [start]);
  const residues = sim.data.residues;
  const td = sim.threeDi,
    tds = sim.threeDiSummary;
  const frame = Math.min(session.frame.value, sim.data.nFrames - 1);
  const mode3 = td && tds && (mode === '3di' || mode === '3di-var') ? mode : undefined;

  const cell = (i: number) => {
    if (pinned !== undefined) {
      if (i === pinned) return { background: '#0b0b0b', color: '#fff' };
      if (Math.abs(i - pinned) < MIN_SEP) return { background: '#ecebe6', color: '#85847f' };
      return shade('share', stats.contact[pinned * n + i]);
    }
    if (mode3 === '3di') {
      // the letter now; dark where it is not the residue's usual one
      const same = td!.letters[frame * n + i] === tds!.consensus[i];
      return same ? { background: '#eef4fc', color: '#0b0b0b' } : { background: '#963b15', color: '#fff' };
    }
    if (mode3 === '3di-var') return shade('spread', tds!.entropy[i] / 2);
    if (mode === 'change') return shade('diverging', 0.5 + (now[i] - start[i]) / 8);
    return shade('share', now[i] / maxCount);
  };

  const blocks: number[][] = [];
  for (let i = 0; i < n; i += 10) blocks.push(Array.from({ length: Math.min(10, n - i) }, (_, k) => i + k));

  const h = hover;
  let lost = 0,
    gained = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + MIN_SEP; j < n; j++) {
      const a = startMap[i * n + j] < cutoff,
        b = frameMap[i * n + j] < cutoff;
      if (a && !b) lost++;
      else if (!a && b) gained++;
    }
  }

  return (
    <div className="mdc-card">
      <div className="mdc-card-head">
        <b>Sequence and contacts</b>
        <span className="mdc-row">
          {pinned === undefined ? (
            <select aria-label="Shade the sequence by" value={mode} onChange={(e) => setMode(e.target.value as Mode)}>
              <option value="now">Contacts in this frame</option>
              <option value="change">Gained or lost since the start</option>
              {td && <option value="3di">3Di letters of this frame</option>}
              {td && <option value="3di-var">3Di variability over time</option>}
            </select>
          ) : (
            <>
              <span>
                partners of <b>{residueLabel(residues[pinned])}</b>
              </span>
              <button type="button" className="mdc-btn" onClick={() => setPinned(undefined)}>
                Unpin
              </button>
            </>
          )}
        </span>
      </div>
      <div
        className="mdc-seq"
        onMouseLeave={() => (setHover(undefined), session.highlight(pinned !== undefined ? [pinned] : undefined))}
      >
        {blocks.map((block) => (
          <span className="mdc-seq-block" key={block[0]}>
            <span className="mdc-seq-num">{residues[block[0]].seq}</span>
            <span className="mdc-seq-row">
              {block.map((i) => {
                const touching =
                  pinned !== undefined &&
                  i !== pinned &&
                  Math.abs(i - pinned) >= MIN_SEP &&
                  frameMap[pinned * n + i] < cutoff;
                const selected = !!pair && (pair.i === i || pair.j === i);
                return (
                  <span
                    key={i}
                    className={`mdc-seq-res${touching ? ' mdc-seq-touch' : ''}${selected ? ' mdc-seq-selected' : ''}`}
                    style={cell(i)}
                    onMouseEnter={() => {
                      setHover(i);
                      session.highlight(pinned !== undefined && pinned !== i ? [pinned, i] : [i]);
                    }}
                    onClick={() => {
                      if (pinned === undefined) setPinned(i);
                      else if (pinned === i) setPinned(undefined);
                      else void session.selectPair({ i: pinned, j: i });
                    }}
                  >
                    {mode3 === '3di' && pinned === undefined
                      ? String.fromCharCode(td!.letters[frame * n + i]).toLowerCase()
                      : oneLetter(residues[i].name)}
                  </span>
                );
              })}
            </span>
          </span>
        ))}
      </div>
      <div className="mdc-legend-bar">
        {pinned !== undefined ? (
          <>
            <span>never in contact</span>
            <div style={{ background: cssGradient('share') }} />
            <span>always · outlined: touching now</span>
          </>
        ) : mode3 === '3di' ? (
          <span>3Di structural letters (lower case). Dark: not the letter this residue has most of the time.</span>
        ) : mode3 === '3di-var' ? (
          <>
            <span>one letter throughout</span>
            <div style={{ background: cssGradient('spread') }} />
            <span>2 bits of entropy and more</span>
          </>
        ) : mode === 'change' ? (
          <>
            <span>4 lost</span>
            <div style={{ background: cssGradient('diverging') }} />
            <span>4 gained</span>
          </>
        ) : (
          <>
            <span>0</span>
            <div style={{ background: cssGradient('share') }} />
            <span>{maxCount} contacts</span>
          </>
        )}
      </div>
      <div className="mdc-readout">
        {h !== undefined ? (
          <>
            <span>
              <b>{residueLabel(residues[h])}</b>
            </span>
            <span>
              {now[h]} contacts now, {start[h]} at the start
            </span>
            {td && tds && (
              <span>
                3Di {String.fromCharCode(td.letters[frame * n + h])} now, usually{' '}
                {String.fromCharCode(tds.consensus[h])} ({(tds.consensusShare[h] * 100).toFixed(0)}% of the time)
              </span>
            )}
            {pinned !== undefined && h !== pinned && (
              <span>
                with {residueLabel(residues[pinned])}: {frameMap[pinned * n + h].toFixed(2)} nm now, in contact{' '}
                {(stats.contact[pinned * n + h] * 100).toFixed(0)}% of the time
              </span>
            )}
          </>
        ) : (
          <span className="mdc-muted">
            Since the start this frame lost {lost} and gained {gained} contacts (Cα closer than {cutoff} nm, at least{' '}
            {MIN_SEP} residues apart).{' '}
            {pinned === undefined
              ? 'Click a residue to see its partners.'
              : 'Click a partner to follow that pair; click the black residue to unpin.'}
          </span>
        )}
      </div>
    </div>
  );
}
