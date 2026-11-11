/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useMemo, useState } from 'react';
import { oneLetter, superpose } from '../align.js';
import { residueLabel } from '../analysis.js';
import { cssGradient, Scales } from '../colors.js';
import type { MdContactsSession } from '../session.js';
import { ns, shortGroup, simColor, useSubject } from './common.js';
import { LinePlot } from './plot.js';

function shade(t: number) {
  const lut = Scales.spread;
  const k = Math.round(Math.max(0, Math.min(1, t)) * 255) * 3;
  const [r, g, b] = [lut[k], lut[k + 1], lut[k + 2]];
  return { background: `rgb(${r},${g},${b})`, color: 0.299 * r + 0.587 * g + 0.114 * b < 140 ? '#fff' : '#0b0b0b' };
}

/**
 * The frame shown against a frame of another simulation: their sequences aligned, the two structures superposed on
 * the aligned residues, and every column shaded by how far the two residues are apart after the fit.
 */
export function AlignTab({ session }: { session: MdContactsSession }) {
  const sims = useSubject(session.simulations);
  const active = useSubject(session.active);
  const frame = useSubject(session.frame);
  const reference = useSubject(session.reference);
  const run = useSubject(session.conformers);
  const [hover, setHover] = useState<number | undefined>();

  const others = sims.filter((s) => s !== active);
  const ref = reference && sims.includes(reference.sim) && reference.sim !== active ? reference : undefined;

  const alignment = useMemo(
    () => (active && ref ? session.alignment(active, ref.sim) : undefined),
    [active, ref?.sim, session],
  );
  const fit = useMemo(
    () =>
      active && ref && alignment
        ? superpose(active.data, Math.min(frame, active.data.nFrames - 1), ref.sim.data, ref.frame, alignment.pairs)
        : undefined,
    [active, ref?.sim, ref?.frame, alignment, frame],
  );
  // displacement by column, and the scale's top: the 95th percentile, at least 0.2 nm
  const view = useMemo(() => {
    if (!alignment || !fit) return undefined;
    const byColumn = new Float32Array(alignment.columns.length).fill(Number.NaN);
    let k = 0;
    alignment.columns.forEach((c, i) => {
      if (c.a >= 0 && c.b >= 0) byColumn[i] = fit.displacement[k++];
    });
    const sorted = Float32Array.from(fit.displacement).sort();
    const p95 = sorted.length ? sorted[Math.floor((sorted.length - 1) * 0.95)] : 0;
    return { byColumn, max: Math.max(0.2, Math.ceil(p95 * 10) / 10) };
  }, [alignment, fit]);

  if (!active) return <div className="mdc-empty">Load a simulation in the Data tab first.</div>;
  if (!others.length)
    return (
      <div className="mdc-empty">Load a second simulation (another replica or state) to align it with this one.</div>
    );

  const pick = (id: string) => {
    const sim = sims.find((s) => s.id === id);
    void session.setReference(sim, sim ? sim.stability.first : 0, ref?.show ?? true);
  };
  const medoids = run ? run.result.conformers.filter((c) => run.sims[c.medoid.traj] === ref?.sim) : [];

  return (
    <div className="mdc-tab">
      <div className="mdc-card">
        <div className="mdc-controls">
          <label>
            <span>Compare the frame shown with</span>
            <select value={ref?.sim.id ?? ''} onChange={(e) => pick(e.target.value)}>
              <option value="">nothing</option>
              {session.groups.map((g) => (
                <optgroup key={g} label={shortGroup(g)}>
                  {session
                    .simsOf(g)
                    .filter((s) => s !== active)
                    .map((s) => (
                      <option key={s.id} value={s.id}>
                        {shortGroup(g)} / {s.label}
                      </option>
                    ))}
                </optgroup>
              ))}
            </select>
          </label>
          {ref && (
            <>
              <label style={{ flex: 1, minWidth: 160 }}>
                <span>
                  its frame at <b>{ns(ref.sim.data.times[ref.frame], 1)} ns</b>
                </span>
                <input
                  type="range"
                  min={0}
                  max={ref.sim.data.nFrames - 1}
                  value={ref.frame}
                  style={{ width: '100%' }}
                  onChange={(e) => void session.setReference(ref.sim, Number(e.target.value), ref.show)}
                />
              </label>
              <button
                type="button"
                className="mdc-btn"
                onClick={() => void session.setReference(ref.sim, ref.sim.stability.first, ref.show)}
              >
                Start
              </button>
              <button
                type="button"
                className="mdc-btn"
                onClick={() => void session.setReference(ref.sim, ref.sim.data.nFrames - 1, ref.show)}
              >
                End
              </button>
              {medoids.map((c) => (
                <button
                  key={c.id}
                  type="button"
                  className="mdc-btn"
                  title="The most typical frame of this conformer"
                  onClick={() => void session.setReference(ref.sim, c.medoid.frame, ref.show)}
                >
                  Medoid C{c.id + 1}
                </button>
              ))}
              <label className="mdc-check">
                <input
                  type="checkbox"
                  checked={ref.show}
                  onChange={(e) => void session.setReference(ref.sim, ref.frame, e.target.checked)}
                />
                <span>Show superposed in 3D (grey)</span>
              </label>
            </>
          )}
        </div>
        {!ref && (
          <div className="mdc-note">
            Pick another replica or state. Its residues are aligned with this simulation's, its frame is superposed on
            the frame shown, and you see residue by residue where the two differ.
          </div>
        )}
      </div>

      {ref && alignment && fit && view && (
        <>
          <div className="mdc-card">
            <div className="mdc-stats">
              <span>
                <span className="mdc-dot" style={{ background: simColor(active) }} />
                <b>{active.label}</b> at {ns(active.data.times[frame] ?? 0, 1)} ns
              </span>
              <span>
                <span className="mdc-dot" style={{ background: '#8f8e88' }} />
                <b>
                  {shortGroup(ref.sim.group)} / {ref.sim.label}
                </b>{' '}
                at {ns(ref.sim.data.times[ref.frame], 1)} ns
              </span>
              <span>
                Cα RMSD <b>{fit.rmsd.toFixed(2)} nm</b> over {alignment.pairs.length / 2} aligned residues
              </span>
              <span>
                sequence identity <b>{(alignment.identity * 100).toFixed(0)}%</b> ({alignment.method})
              </span>
            </div>
          </div>

          <div className="mdc-card">
            <div className="mdc-card-head">
              <b>Aligned sequences, shaded by distance after superposition</b>
            </div>
            <div className="mdc-seq" onMouseLeave={() => (setHover(undefined), session.highlight(undefined))}>
              {Array.from({ length: Math.ceil(alignment.columns.length / 10) }, (_, bi) => {
                const cols = alignment.columns.slice(bi * 10, bi * 10 + 10);
                const first = cols.find((c) => c.a >= 0);
                return (
                  <span className="mdc-seq-block" key={bi}>
                    <span className="mdc-seq-num">{first ? active.data.residues[first.a].seq : ' '}</span>
                    {(['a', 'b'] as const).map((row) => (
                      <span className="mdc-seq-row" key={row}>
                        {cols.map((c, ci) => {
                          const col = bi * 10 + ci;
                          const idx = c[row];
                          const res = idx >= 0 ? (row === 'a' ? active : ref.sim).data.residues[idx] : undefined;
                          const d = view.byColumn[col];
                          const differs =
                            c.a >= 0 &&
                            c.b >= 0 &&
                            oneLetter(active.data.residues[c.a].name) !== oneLetter(ref.sim.data.residues[c.b].name);
                          return (
                            <span
                              key={ci}
                              className={`mdc-seq-res${differs ? ' mdc-seq-differs' : ''}${hover === col ? ' mdc-seq-hover' : ''}`}
                              style={
                                !res
                                  ? { background: 'transparent', color: '#85847f' }
                                  : Number.isNaN(d)
                                    ? { background: '#ecebe6', color: '#52514e' }
                                    : shade(d / view.max)
                              }
                              onMouseEnter={() => {
                                setHover(col);
                                session.highlight(c.a >= 0 ? [c.a] : undefined);
                              }}
                              onClick={() => c.a >= 0 && session.focus([c.a])}
                            >
                              {res ? oneLetter(res.name) : '-'}
                            </span>
                          );
                        })}
                      </span>
                    ))}
                  </span>
                );
              })}
            </div>
            <div className="mdc-legend-bar">
              <span>0</span>
              <div style={{ background: cssGradient('spread') }} />
              <span>{view.max.toFixed(1)} nm apart and more</span>
            </div>
            <div className="mdc-readout">
              {hover !== undefined ? (
                (() => {
                  const c = alignment.columns[hover];
                  return (
                    <>
                      <span>
                        <b>{c.a >= 0 ? residueLabel(active.data.residues[c.a]) : 'gap'}</b> ↔{' '}
                        <b>{c.b >= 0 ? residueLabel(ref.sim.data.residues[c.b]) : 'gap'}</b>
                      </span>
                      {!Number.isNaN(view.byColumn[hover]) && (
                        <span>{view.byColumn[hover].toFixed(2)} nm apart after the fit</span>
                      )}
                    </>
                  );
                })()
              ) : (
                <span className="mdc-muted">
                  Top row: the simulation shown; bottom row: the reference. Underlined: different residue types. Hover
                  to highlight in 3D, click to zoom.
                </span>
              )}
            </div>
          </div>

          <LinePlot
            title="Distance between aligned residues after superposition"
            series={[
              {
                key: 'd',
                label: 'displacement',
                color: '#eb6834',
                x: Array.from(
                  { length: alignment.pairs.length / 2 },
                  (_, k) => active.data.residues[alignment.pairs[2 * k]].seq,
                ),
                y: fit.displacement,
              },
            ]}
            xLabel="residue"
            yUnit="nm"
            yMin={0}
            xFormat={(x) => String(Math.round(x))}
            onHoverX={(x) => {
              if (x === undefined) return session.highlight(undefined);
              let best = 0;
              const res = active.data.residues;
              for (let i = 1; i < res.length; i++) if (Math.abs(res[i].seq - x) < Math.abs(res[best].seq - x)) best = i;
              session.highlight([best]);
            }}
          />
          <div className="mdc-note">
            Furthest apart:{' '}
            {Array.from(fit.displacement, (v, k) => ({ v, k }))
              .sort((a, b) => b.v - a.v)
              .slice(0, 6)
              .map(({ v, k }) => `${residueLabel(active.data.residues[alignment.pairs[2 * k]])} (${v.toFixed(2)} nm)`)
              .join(', ')}
            . The numbers follow the frame slider and playback; the grey structure in 3D is refitted when the frame
            rests. The fit uses all aligned Cα atoms, so a moving domain shows as displacement on both sides of its
            hinge.
          </div>
        </>
      )}
    </div>
  );
}
