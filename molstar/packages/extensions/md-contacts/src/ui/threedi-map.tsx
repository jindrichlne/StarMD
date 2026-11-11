/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useEffect, useRef, useState } from 'react';
import { residueLabel } from '../analysis.js';
import type { MdContactsSession, Simulation } from '../session.js';
import { ns } from './common.js';

/**
 * The trajectory as an alignment of its frames in 3Di letters, drawn small: time runs left to right, residues top to
 * bottom, and a cell is dark where the residue's letter is not its usual one. Vertical bands are moments when much of
 * the local structure changed; horizontal bands are residues that switch.
 */
export function ThreeDiMap({ session, sim, frame }: { session: MdContactsSession; sim: Simulation; frame: number }) {
  const canvas = useRef<HTMLCanvasElement>(null);
  const [hover, setHover] = useState<{ f: number; r: number } | undefined>();
  const td = sim.threeDi,
    sum = sim.threeDiSummary;
  useEffect(() => {
    const c = canvas.current;
    if (!c || !td || !sum) return;
    const { nFrames, nRes, letters } = td;
    c.width = nFrames;
    c.height = nRes;
    const ctx = c.getContext('2d')!;
    const img = ctx.createImageData(nFrames, nRes);
    for (let r = 0; r < nRes; r++) {
      const usual = sum.consensus[r];
      for (let f = 0; f < nFrames; f++) {
        const o = (r * nFrames + f) * 4;
        const differs = letters[f * nRes + r] !== usual;
        // skipped equilibration frames are paler
        const pale = f < sum.first;
        img.data[o] = differs ? (pale ? 0xd9 : 0x96) : 0xf4;
        img.data[o + 1] = differs ? (pale ? 0xb5 : 0x3b) : 0xf6;
        img.data[o + 2] = differs ? (pale ? 0xa4 : 0x15) : 0xfb;
        img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }, [td, sum]);
  if (!td || !sum) return null;
  const locate = (ev: React.MouseEvent) => {
    const b = ev.currentTarget.getBoundingClientRect();
    return {
      f: Math.max(0, Math.min(td.nFrames - 1, Math.floor(((ev.clientX - b.left) / b.width) * td.nFrames))),
      r: Math.max(0, Math.min(td.nRes - 1, Math.floor(((ev.clientY - b.top) / b.height) * td.nRes))),
    };
  };
  return (
    <div className="mdc-card">
      <div className="mdc-card-head">
        <b>Local structure over time (3Di) · {sim.label}</b>
      </div>
      <div
        className="mdc-3di-map"
        role="img"
        aria-label="Residues by time: dark where a residue's 3Di letter is not its usual one"
        onMouseMove={(ev) => {
          const p = locate(ev);
          setHover(p);
          session.highlight([p.r]);
        }}
        onMouseLeave={() => {
          setHover(undefined);
          session.highlight(undefined);
        }}
        onClick={(ev) => session.setFrame(locate(ev).f)}
      >
        <canvas ref={canvas} />
        <i style={{ left: `${((frame + 0.5) / td.nFrames) * 100}%` }} />
      </div>
      <div className="mdc-readout">
        {hover ? (
          <>
            <span>
              <b>{residueLabel(sim.data.residues[hover.r])}</b> at {ns(sim.data.times[hover.f], 1)} ns
            </span>
            <span>
              letter {String.fromCharCode(td.letters[hover.f * td.nRes + hover.r])}, usually{' '}
              {String.fromCharCode(sum.consensus[hover.r])} ({(sum.consensusShare[hover.r] * 100).toFixed(0)}% of the
              time)
            </span>
          </>
        ) : (
          <span className="mdc-muted">
            Time left to right, residues top ({sim.data.residues[0].seq}) to bottom (
            {sim.data.residues[td.nRes - 1].seq}). Dark: the residue's 3Di letter is not its usual one. Click to show
            that frame.
          </span>
        )}
      </div>
    </div>
  );
}
