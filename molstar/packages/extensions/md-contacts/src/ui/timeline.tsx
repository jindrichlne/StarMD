/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import type { MdContactsSession } from '../session.js';
import { ns, shortGroup, simColor, useSubject } from './common.js';

/** Which simulation is shown in 3D, and which of its frames: always visible above the tabs. */
export function Timeline({ session }: { session: MdContactsSession }) {
  const sims = useSubject(session.simulations);
  const active = useSubject(session.active);
  const frame = useSubject(session.frame);
  const playing = useSubject(session.playing);
  const smooth = useSubject(session.smooth);
  if (!active) return <div className="mdc-timeline mdc-muted">Load a simulation in the Data tab to begin.</div>;
  const t = active.data.times;
  const groups = session.groups;
  return (
    <div className="mdc-timeline">
      <div className="mdc-row">
        <span className="mdc-dot" style={{ background: simColor(active) }} />
        <select
          aria-label="Simulation shown in 3D"
          value={active.id}
          onChange={(e) => {
            const sim = sims.find((s) => s.id === e.target.value);
            if (sim) void session.setActive(sim);
          }}
        >
          {groups.map((g) => (
            <optgroup key={g} label={shortGroup(g)}>
              {session.simsOf(g).map((s) => (
                <option key={s.id} value={s.id}>
                  {groups.length > 1 ? `${shortGroup(g)} / ${s.label}` : s.label}
                </option>
              ))}
            </optgroup>
          ))}
        </select>
        <button
          type="button"
          className="mdc-btn mdc-btn-icon"
          aria-label={playing ? 'Pause' : 'Play'}
          onClick={() => (playing ? session.pause() : session.play())}
        >
          {playing ? '❚❚' : '▶'}
        </button>
        <select
          aria-label="Playback step"
          defaultValue={String(session.playStep)}
          onChange={(e) => {
            session.playStep = Number(e.target.value);
          }}
          title="Frames advanced per step while playing"
        >
          {[1, 2, 5, 10, 25].map((s) => (
            <option key={s} value={s}>
              ×{s}
            </option>
          ))}
        </select>
        <select
          aria-label="Smoothing"
          value={smooth}
          onChange={(e) => void session.setSmooth(Number(e.target.value))}
          title="Draw interpolated frames between the saved ones: slower, smoother motion. The analysis always uses the saved frames."
        >
          <option value={1}>Smooth: off</option>
          {[2, 4, 8, 16].map((s) => (
            <option key={s} value={s}>
              Smooth: {s - 1} in between
            </option>
          ))}
        </select>
        <span className="mdc-time">
          <b>{ns(t[frame] ?? 0, 2)} ns</b> · frame {frame + 1}/{active.data.nFrames}
        </span>
      </div>
      <input
        type="range"
        aria-label="Frame"
        min={0}
        max={active.data.nFrames - 1}
        value={frame}
        onChange={(e) => session.setFrame(Number(e.target.value))}
      />
    </div>
  );
}
