/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { residueKey, residueLabel, separationMap, statsOfFrames, switchPairs, type TreeNode } from '../analysis.js';
import { localSwitches } from '../threedi.js';
import { seriesColor } from '../colors.js';
import type { ConformerRun, MdContactsSession, Simulation } from '../session.js';
import { Help, ns, Sections, shortGroup, simColor, useSection, useSubject } from './common.js';
import { Heatmap } from './heatmap.js';
import { Icicle } from './icicle.js';

const confColor = (c: number) => (c < 0 ? '#d9d8d4' : seriesColor(c));
const confName = (c: number) => `C${c + 1}`;
/** a simulation's name in a run: with its group when the run spans several */
const simName = (run: ConformerRun, sim: Simulation) =>
  new Set(run.sims.map((s) => s.group)).size > 1 ? `${shortGroup(sim.group)} / ${sim.label}` : sim.label;

function Strip({
  run,
  index,
  sim,
  cursor,
  longest,
  onPick,
}: {
  run: ConformerRun;
  index: number;
  sim: Simulation;
  cursor?: number;
  /** ps: the longest simulation of the run, which spans the full width */
  longest: number;
  onPick: (frame: number) => void;
}) {
  const canvas = useRef<HTMLCanvasElement>(null);
  const a = run.result.assignment[index];
  useEffect(() => {
    const c = canvas.current;
    if (!c) return;
    c.width = a.length;
    c.height = 1;
    const ctx = c.getContext('2d')!;
    for (let f = 0; f < a.length; f++) {
      ctx.fillStyle = confColor(a[f]);
      ctx.fillRect(f, 0, 1, 1);
    }
  }, [a]);
  return (
    <div className="mdc-strip-row">
      <span className="mdc-strip-label" title={sim.id}>
        <span className="mdc-dot" style={{ background: simColor(sim) }} />
        {simName(run, sim)}
      </span>
      <div
        className="mdc-strip"
        style={{
          width: `${Math.max(1, ((sim.data.times[a.length - 1] - sim.data.times[0]) / (longest || 1)) * 100)}%`,
        }}
        onClick={(ev) => {
          const r = ev.currentTarget.getBoundingClientRect();
          onPick(Math.min(a.length - 1, Math.max(0, Math.floor(((ev.clientX - r.left) / r.width) * a.length))));
        }}
        title="Click to show that frame"
      >
        <canvas ref={canvas} />
        {cursor !== undefined && <i style={{ left: `${((cursor + 0.5) / a.length) * 100}%` }} />}
      </div>
    </div>
  );
}

function Scatter({ run, onPick }: { run: ConformerRun; onPick: (traj: number, frame: number) => void }) {
  const pts = run.result.sample;
  const dom = useMemo(() => {
    let x0 = Infinity,
      x1 = -Infinity,
      y0 = Infinity,
      y1 = -Infinity;
    for (const p of pts) {
      x0 = Math.min(x0, p.x);
      x1 = Math.max(x1, p.x);
      y0 = Math.min(y0, p.y);
      y1 = Math.max(y1, p.y);
    }
    // one scale for both axes: distances on the map mean the same in every direction
    const span = Math.max(x1 - x0, y1 - y0, 1e-6);
    return { cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, span };
  }, [pts]);
  const [hover, setHover] = useState<number | undefined>();
  const S = 300,
    P = 14;
  const sx = (x: number) => S / 2 + ((x - dom.cx) / dom.span) * (S - 2 * P),
    sy = (y: number) => S / 2 - ((y - dom.cy) / dom.span) * (S - 2 * P);
  const medoids = run.result.conformers.map((c) => c.medoid);
  const h = hover !== undefined ? pts[hover] : undefined;
  return (
    <div className="mdc-scatter">
      <svg
        viewBox={`0 0 ${S} ${S}`}
        width="100%"
        role="img"
        aria-label="Map of sampled frames by distance-map similarity"
      >
        <rect x={0.5} y={0.5} width={S - 1} height={S - 1} className="mdc-frame" />
        {pts.map((p, i) => (
          <circle
            key={i}
            cx={sx(p.x)}
            cy={sy(p.y)}
            r={hover === i ? 5 : 3}
            fill={confColor(p.conformer)}
            fillOpacity={0.75}
            onMouseEnter={() => setHover(i)}
            onMouseLeave={() => setHover(undefined)}
            onClick={() => onPick(p.ref.traj, p.ref.frame)}
            style={{ cursor: 'pointer' }}
          />
        ))}
        {medoids.map((m, c) => {
          const p = pts.find((q) => q.ref.traj === m.traj && q.ref.frame === m.frame);
          if (!p) return null;
          return (
            <g key={c} pointerEvents="none">
              <circle cx={sx(p.x)} cy={sy(p.y)} r={7} fill="none" stroke="#0b0b0b" strokeWidth={1.5} />
              <text x={sx(p.x) + 10} y={sy(p.y) + 4} className="mdc-scatter-label">
                {confName(c)}
              </text>
            </g>
          );
        })}
      </svg>
      <div className="mdc-readout">
        {h ? (
          <span>
            {simName(run, run.sims[h.ref.traj])}, {ns(run.sims[h.ref.traj].data.times[h.ref.frame], 1)} ns ·{' '}
            {confName(h.conformer)}
          </span>
        ) : (
          <span className="mdc-muted">
            Each dot is a sampled frame; close dots have similar distance maps. Rings mark the medoids. Click a dot to
            show it.
          </span>
        )}
      </div>
    </div>
  );
}

/** Conformers: frames grouped by distance-map similarity, when each replica visits them, and what tells them apart. */
export function ConformersTab({ session }: { session: MdContactsSession }) {
  const sims = useSubject(session.simulations);
  const active = useSubject(session.active);
  const frame = useSubject(session.frame);
  const run = useSubject(session.conformers);
  const pair = useSubject(session.pair);
  const busy = useSubject(session.busy);
  const overlay = useSubject(session.overlay);
  const [scope, setScope] = useState<string>('');
  const [k, setK] = useState(0);
  const [stepFactor, setStepFactor] = useState(2.25);
  const [perTraj, setPerTraj] = useState(150);
  const [error, setError] = useState('');
  const [a, setA] = useState(0);
  const [b, setB] = useState(1);
  const [fingerprint, setFingerprint] = useState<'distance' | '3di'>('distance');
  /** tree blocks picked in the icicle plot for comparison (ids, at most two) */
  const [blocks, setBlocks] = useState<number[]>([]);

  const groups = session.groups;
  const scopeValue = scope && (scope === '*' || groups.includes(scope)) ? scope : (active?.group ?? groups[0] ?? '');

  const scopeSims = scopeValue === '*' ? sims : session.simsOf(scopeValue);
  const has3Di = scopeSims.length > 0 && scopeSims.every((s) => s.threeDi);
  const fp = fingerprint === '3di' && has3Di ? '3di' : 'distance';
  const start = async () => {
    setError('');
    try {
      await session.runConformers(scopeValue, {
        nClusters: k,
        stepFactor,
        framesPerTrajectory: perTraj,
        fingerprint: fp,
      });
      setA(0);
      setB(1);
      setBlocks([]);
    } catch (e) {
      setError(String((e as Error).message ?? e));
    }
  };

  const nConf = run?.result.conformers.length ?? 0;
  const ia = Math.min(a, Math.max(0, nConf - 1)),
    ib = Math.min(b, Math.max(0, nConf - 1));
  const tree = run?.result.tree ?? [];
  const picked =
    blocks.length === 2 && blocks.every((id) => id < tree.length) ? [tree[blocks[0]], tree[blocks[1]]] : undefined;
  // the two sides compared: two blocks of the tree when picked, else two conformers
  const sides = useMemo(() => {
    if (!run) return undefined;
    const frames = (members: ArrayLike<number>) => Array.from(members, (m) => run.result.sample[m].ref);
    if (picked) {
      const data = run.sims.map((s) => s.data);
      const cutoff = run.result.options.cutoff;
      return picked.map((n, i) => ({
        name: `block ${'AB'[i]}`,
        frames: frames(n.members),
        stats: statsOfFrames(data, frames(n.members), cutoff, run.maps),
      }));
    }
    if (nConf < 2 || ia === ib) return undefined;
    return [ia, ib].map((c) => ({
      name: confName(c),
      frames: run.result.sample.filter((p) => p.conformer === c).map((p) => p.ref),
      stats: run.result.conformers[c].stats,
    }));
  }, [run, picked?.[0], picked?.[1], ia, ib, nConf]);
  const compare = useMemo(() => {
    if (!sides) return undefined;
    const A = sides[0].stats,
      B = sides[1].stats;
    const n = A.n;
    const diff = new Float32Array(n * n);
    for (let q = 0; q < n * n; q++) diff[q] = A.mean[q] - B.mean[q];
    let max = 0.2;
    for (let q = 0; q < n * n; q++) if (Math.abs(diff[q]) > max) max = Math.abs(diff[q]);
    return { n, diff, sep: separationMap(A, B), pairs: switchPairs(A, B, { top: 15 }), max: Math.min(max, 3) };
  }, [sides]);
  // residues whose 3Di letter differs between the two sides (local structure, where the pairs above are global)
  const local = useMemo(() => {
    if (!run || !sides || !run.sims.every((s) => s.threeDi)) return undefined;
    const of = (side: (typeof sides)[0]) =>
      side.frames.map((f) => ({ t: run.sims[f.traj].threeDi!, frame: f.frame, map: run.maps[f.traj] }));
    return localSwitches(of(sides[0]), of(sides[1]), run.residues.length, 12);
  }, [run, sides]);
  const nameA = sides?.[0].name ?? confName(ia),
    nameB = sides?.[1].name ?? confName(ib);
  const unit = run?.result.unit ?? 'nm';
  const fmtHeight = (h: number) => (unit === 'nm' ? h.toFixed(h < 0.1 ? 3 : 2) : `${(h * 100).toFixed(0)}%`);
  const unitLabel = unit === 'nm' ? 'nm dRMSD' : 'of residues with another 3Di letter';

  /** index in the active simulation of a residue of the run */
  const toActive = useCallback(
    (i: number) => {
      if (!run || !active) return -1;
      const idx = run.sims.indexOf(active);
      if (idx >= 0) return run.maps[idx][i];
      const key = residueKey(run.residues[i]);
      return active.data.residues.findIndex((r) => residueKey(r) === key);
    },
    [run, active],
  );
  const fromActive = useCallback(
    (i: number) => {
      if (!run || !active) return -1;
      const key = residueKey(active.data.residues[i]);
      return run.residues.findIndex((r) => residueKey(r) === key);
    },
    [run, active],
  );
  const label = useCallback((i: number) => (run ? residueLabel(run.residues[i]) : ''), [run]);
  const tick = useCallback((i: number) => (run ? run.residues[i].seq : i), [run]);
  const selected = useMemo(() => {
    if (!pair || !run) return undefined;
    const i = fromActive(pair.i),
      j = fromActive(pair.j);
    return i >= 0 && j >= 0 ? { i, j } : undefined;
  }, [pair, run, fromActive]);

  const sectionOptions = [
    ['overview', 'Overview'],
    ['tree', 'Tree'],
    ['frames', 'Frame map'],
    ...(nConf > 1 || picked
      ? ([
          ['diff', 'Difference map'],
          ['switch', 'Switch residues'],
        ] as const)
      : []),
  ] as const;
  const [section, setSection] = useSection(
    'conformers',
    sectionOptions.map((o) => o[0]),
  );
  if (!sims.length) return <div className="mdc-empty">Load a simulation in the Data tab to look for conformers.</div>;

  const total = run ? run.result.conformers.reduce((s, c) => s + c.size, 0) : 0;
  const activeIdx = run && active ? run.sims.indexOf(active) : -1;
  const current = run && activeIdx >= 0 ? run.result.assignment[activeIdx][frame] : undefined;

  return (
    <div className="mdc-tab">
      <div className="mdc-controls">
        <label>
          <span>Look in</span>
          <select value={scopeValue} onChange={(e) => setScope(e.target.value)}>
            {groups.map((g) => (
              <option key={g} value={g}>
                {shortGroup(g)} ({session.simsOf(g).length} replica{session.simsOf(g).length > 1 ? 's' : ''})
              </option>
            ))}
            {groups.length > 1 && <option value="*">all loaded simulations together</option>}
          </select>
        </label>
        <label>
          <span>Conformers</span>
          <select value={k} onChange={(e) => setK(Number(e.target.value))}>
            <option value={0}>automatic</option>
            {[2, 3, 4, 5, 6, 7, 8].map((v) => (
              <option key={v} value={v}>
                {v}
              </option>
            ))}
          </select>
        </label>
        {k === 0 && (
          <label title="How different two groups of frames must be to count as two conformers, relative to how much the structure moves anyway between neighbouring stretches of a run. Proteins have sub-states at every scale: finer finds more of them.">
            <span>Scale of a conformer</span>
            <select value={stepFactor} onChange={(e) => setStepFactor(Number(e.target.value))}>
              <option value={3}>coarse: big changes only</option>
              <option value={2.25}>default</option>
              <option value={1.75}>fine: also sub-states</option>
            </select>
          </label>
        )}
        <label>
          <span>Frames clustered per replica</span>
          <select value={perTraj} onChange={(e) => setPerTraj(Number(e.target.value))}>
            {[75, 150, 300, 500].map((v) => (
              <option key={v} value={v}>
                {v}
              </option>
            ))}
          </select>
        </label>
        <label title="What frames are compared by. 3Di letters describe each residue's local environment (Foldseek); they need the server to have Foldseek.">
          <span>Compare frames by</span>
          <select value={fp} onChange={(e) => setFingerprint(e.target.value as 'distance' | '3di')}>
            <option value="distance">Cα distance maps</option>
            <option value="3di" disabled={!has3Di}>
              3Di letters{has3Di ? '' : ' (not loaded)'}
            </option>
          </select>
        </label>
        <button type="button" className="mdc-btn mdc-btn-primary" disabled={!!busy} onClick={() => void start()}>
          {busy || (run ? 'Run again' : 'Find conformers')}
        </button>
      </div>
      {error && <div className="mdc-error">{error}</div>}
      {!run && !error && (
        <div className="mdc-note">
          Frames are fingerprinted by their Cα distance map (averaged over short stretches, which removes thermal
          noise), clustered by dRMSD, and every frame is then assigned to the nearest conformer. Automatic mode keeps a
          split of the tree when the two groups are far enough apart, each hold enough frames and have a gap between
          them. No superposition is involved, so rigid-body motion does not count.
        </div>
      )}

      {run && (
        <>
          <Sections label="Conformer sections" options={sectionOptions} value={section} onChange={setSection}>
            {current !== undefined && current >= 0 && (
              <span>
                frame in 3D: <span className="mdc-chip" style={{ background: confColor(current) }} />{' '}
                {confName(current)}
              </span>
            )}
          </Sections>
          {section === 'overview' && (
            <>
              <div className="mdc-card">
                <div className="mdc-card-head">
                  <b>
                    {nConf} conformer{nConf > 1 ? 's' : ''} in{' '}
                    {run.scope === '*' ? 'all simulations' : shortGroup(run.scope)}
                  </b>
                </div>
                <div className="mdc-reason">
                  {(() => {
                    const r = run.result;
                    const fmtSep = (v: number) => (unit === 'nm' ? `${v.toFixed(2)} nm` : `${(v * 100).toFixed(0)}%`);
                    const pct = (v: number) => `${Math.round(v * 100)}%`;
                    const apart = r.splits.filter((t) => t.outcome === 'apart');
                    const joined = r.splits.filter((t) => t.outcome === 'joined' || t.outcome === 'transit');
                    // sub-states: a gap, but closer than the scale
                    const finer = r.splits.filter((t) => t.outcome === 'close' && t.dip <= r.options.maxDip);
                    const scale = `${fmtSep(r.scale)} (${r.options.stepFactor}× the ${fmtSep(r.step)} the structure moves between neighbouring stretches of a run)`;
                    const rule = `To count, two groups must be at least ${scale} apart, with few frames between them, and each must hold at least ${pct(r.options.minShare)} of one run's frames.`;
                    const closest = apart[0];
                    const weakest = apart.slice().sort((x, y) => y.dip - x.dip)[0];
                    let found =
                      r.suggested === 1
                        ? `No separate conformers found at this scale. ${rule}`
                        : `${r.suggested} separate conformers found. The two closest are ${fmtSep(closest.separation)} apart${unit === 'nm' ? ' (dRMSD)' : ' (residues with another 3Di letter)'}; the least clear gap has ${pct(weakest.dip)} as many frames in the middle as in the groups (${pct(r.options.maxDip)} allowed). ${rule}`;
                    if (joined.length) {
                      const widest = joined.slice().sort((x, y) => y.separation - x.separation)[0];
                      found += ` ${joined.length} group${joined.length > 1 ? 's were' : ' was'} far enough from a neighbour (up to ${fmtSep(widest.separation)}) but frames fill the space between them: one continuous motion, not counted.`;
                    }
                    if (finer.length) {
                      const widest = finer.slice().sort((x, y) => y.separation - x.separation)[0];
                      found += ` Below this scale there are sub-states with a gap (up to ${fmtSep(widest.separation)} apart).`;
                    }
                    // how firm the number is
                    const [fine, coarse] = r.nearby;
                    found +=
                      fine.count === r.suggested && coarse.count === r.suggested
                        ? ` The number is the same at a finer (${fmtSep(fine.scale)}) and a coarser (${fmtSep(coarse.scale)}) scale.`
                        : ` The number depends on the scale: ${fine.count} at ${fmtSep(fine.scale)}, ${coarse.count} at ${fmtSep(coarse.scale)}.`;
                    const asked =
                      r.options.nClusters > 0 && r.options.nClusters !== r.suggested
                        ? r.options.nClusters > r.suggested
                          ? ` You asked for ${r.options.nClusters}: the extra groups are closer than that scale, too small or have no gap, so read them as regions of one motion.`
                          : ` You asked for ${r.options.nClusters}, so separate conformers are shown merged.`
                        : '';
                    const transit =
                      !asked && r.transit > 0.005
                        ? ` ${(r.transit * 100).toFixed(0)}% of the frames lie beside a conformer without a gap to it (typically on the way there from the start) and are shown with the nearest one.`
                        : '';
                    return found + transit + asked;
                  })()}
                </div>
                <table className="mdc-table">
                  <thead>
                    <tr>
                      <th>Conformer</th>
                      <th>Frames</th>
                      {run.sims.map((s) => (
                        <th key={s.id} title={s.id}>
                          {simName(run, s)}
                        </th>
                      ))}
                      <th title="Mean dRMSD of the members to the medoid">Spread</th>
                      <th>Most typical frame</th>
                    </tr>
                  </thead>
                  <tbody>
                    {run.result.conformers.map((c) => {
                      const ms = run.sims[c.medoid.traj];
                      return (
                        <tr key={c.id}>
                          <td>
                            <span className="mdc-chip" style={{ background: confColor(c.id) }} /> {confName(c.id)}
                          </td>
                          <td>
                            {c.size} ({((c.size / Math.max(1, total)) * 100).toFixed(0)}%)
                          </td>
                          {c.perTrajectory.map((v, i) => (
                            <td key={i}>
                              {(
                                (v / Math.max(1, run.sims[i].data.nFrames - run.sims[i].stability.first)) *
                                100
                              ).toFixed(0)}
                              %
                            </td>
                          ))}
                          <td>{c.spread.toFixed(3)} nm</td>
                          <td>
                            <button
                              type="button"
                              className="mdc-btn"
                              onClick={() => void session.show(ms, c.medoid.frame)}
                            >
                              {simName(run, ms)}, {ns(ms.data.times[c.medoid.frame], 1)} ns
                            </button>
                          </td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
                <div className="mdc-row">
                  <button
                    type="button"
                    className={overlay ? 'mdc-btn' : 'mdc-btn mdc-btn-primary'}
                    disabled={!!busy}
                    onClick={() => void session.showConformers(!overlay)}
                  >
                    {overlay
                      ? 'Back to the trajectory'
                      : `Show ${nConf > 1 ? `all ${nConf} conformers` : 'the conformer'} superposed in 3D`}
                  </button>
                  {overlay && (
                    <span className="mdc-plot-legend">
                      {overlay.map((o) => (
                        <span
                          key={o.conformer}
                          title={`${simName(run, o.sim)}, ${ns(o.sim.data.times[o.frame], 1)} ns`}
                        >
                          <i style={{ background: confColor(o.conformer), height: 10 }} />
                          {confName(o.conformer)}
                          {o.conformer > 0 ? ` · ${o.rmsd.toFixed(2)} nm from C1` : ''}
                        </span>
                      ))}
                    </span>
                  )}
                </div>
                <Help>
                  Shares are of simulated frames, not populations: they depend on where the runs started and how long
                  they are. The superposed view shows the most typical frame of every conformer, each in its colour,
                  fitted on C1 over all aligned Cα atoms; the simulations are hidden meanwhile, and moving the frame
                  slider, playing or picking a frame goes back to the trajectory.
                </Help>
              </div>

              <div className="mdc-card">
                <div className="mdc-card-head">
                  <b title="One row per replica, time left to right on a shared scale. Click to show that frame.">
                    When each replica is in which conformer
                  </b>
                  <span className="mdc-plot-legend">
                    {run.result.conformers.map((c) => (
                      <span key={c.id}>
                        <i style={{ background: confColor(c.id) }} />
                        {confName(c.id)}
                      </span>
                    ))}
                  </span>
                </div>
                {run.sims.map((s, i) => (
                  <Strip
                    key={s.id}
                    run={run}
                    index={i}
                    sim={s}
                    longest={Math.max(...run.sims.map((x) => x.data.times[x.data.nFrames - 1] - x.data.times[0]))}
                    cursor={s === active ? frame : undefined}
                    onPick={(f) => void session.show(s, f)}
                  />
                ))}
              </div>
            </>
          )}

          {section === 'tree' && (
            <div className="mdc-card">
              <div className="mdc-card-head">
                <b>The clustering tree</b>
                {blocks.length > 0 && (
                  <button type="button" className="mdc-btn" onClick={() => setBlocks([])}>
                    Clear picked blocks
                  </button>
                )}
              </div>
              <Icicle
                tree={tree}
                cutHeight={run.result.cutHeight}
                format={fmtHeight}
                unitLabel={unitLabel}
                color={confColor}
                selected={blocks}
                total={run.result.sample.length}
                describe={(n: TreeNode) =>
                  n.perTrajectory
                    .map((v, i) =>
                      v ? `${simName(run, run.sims[i])} ${((v / n.members.length) * 100).toFixed(0)}%` : '',
                    )
                    .filter(Boolean)
                    .join(', ')
                }
                onPick={(n: TreeNode) => {
                  const m = run.result.sample[n.medoid].ref;
                  void session.show(run.sims[m.traj], m.frame);
                  setBlocks((old) => (old.includes(n.id) ? old.filter((x) => x !== n.id) : [...old, n.id].slice(-2)));
                }}
              />
            </div>
          )}
          {section === 'frames' && (
            <div className="mdc-card">
              <div className="mdc-card-head">
                <b>Map of the sampled frames</b>
              </div>
              <Scatter run={run} onPick={(t, f) => void session.show(run.sims[t], f)} />
            </div>
          )}

          {(section === 'diff' || section === 'switch') && (nConf > 1 || picked) && (
            <div className="mdc-card">
              <div className="mdc-card-head">
                <b>{picked ? 'What tells the two picked blocks apart' : 'What tells two conformers apart'}</b>
                <span className="mdc-row" style={{ display: picked ? 'none' : undefined }}>
                  <select aria-label="First conformer" value={ia} onChange={(e) => setA(Number(e.target.value))}>
                    {run.result.conformers.map((c) => (
                      <option key={c.id} value={c.id}>
                        {confName(c.id)}
                      </option>
                    ))}
                  </select>
                  vs
                  <select aria-label="Second conformer" value={ib} onChange={(e) => setB(Number(e.target.value))}>
                    {run.result.conformers.map((c) => (
                      <option key={c.id} value={c.id}>
                        {confName(c.id)}
                      </option>
                    ))}
                  </select>
                  <span className="mdc-muted">dRMSD {run.result.betweenDrmsd[ia * nConf + ib].toFixed(3)} nm</span>
                </span>
              </div>
              {compare ? (
                <>
                  {section === 'diff' && (
                    <div className="mdc-col">
                      <Heatmap
                        n={compare.n}
                        upper={{
                          label: `Mean ${nameA} − mean ${nameB}`,
                          unit: 'nm',
                          values: compare.diff,
                          scale: 'diverging',
                          min: -compare.max,
                          max: compare.max,
                        }}
                        lower={{
                          label: 'Share of the variance between the two',
                          unit: '',
                          values: compare.sep,
                          scale: 'share',
                          min: 0,
                          max: 1,
                        }}
                        label={label}
                        tick={tick}
                        selected={selected}
                        marks={compare.pairs.slice(0, 8)}
                        onHover={(p) => {
                          if (!p) return session.highlight(undefined);
                          const r = [toActive(p.i), toActive(p.j)].filter((x) => x >= 0);
                          session.highlight(r);
                        }}
                        onSelect={(p) => {
                          const i = toActive(p.i),
                            j = toActive(p.j);
                          if (i >= 0 && j >= 0) void session.selectPair({ i, j });
                        }}
                      />
                      <div className="mdc-note">
                        Blue: closer in {nameA}; red: closer in {nameB}. Rings mark the top switch pairs. Click a cell
                        to select the pair; it is followed in 3D and in the Contacts tab.
                      </div>
                    </div>
                  )}
                  {section === 'switch' && (
                    <div className="mdc-pair">
                      <div className="mdc-col">
                        <table className="mdc-table">
                          <thead>
                            <tr>
                              <th>Switch pair</th>
                              <th>{nameA}</th>
                              <th>{nameB}</th>
                              <th>Change</th>
                              <th title="Share of this distance's variance that is between the two conformers">
                                Between share
                              </th>
                              <th title="In contact in one conformer and not in the other">Contact switch</th>
                            </tr>
                          </thead>
                          <tbody>
                            {compare.pairs.map((p) => (
                              <tr
                                key={`${p.i}-${p.j}`}
                                className={
                                  selected && selected.i === p.i && selected.j === p.j ? 'mdc-active' : undefined
                                }
                                onClick={() => {
                                  const i = toActive(p.i),
                                    j = toActive(p.j);
                                  if (i >= 0 && j >= 0) {
                                    void session.selectPair({ i, j });
                                    session.focus([i, j]);
                                  }
                                }}
                                onMouseEnter={() =>
                                  session.highlight([toActive(p.i), toActive(p.j)].filter((x) => x >= 0))
                                }
                                onMouseLeave={() => session.highlight(undefined)}
                              >
                                <td>
                                  {label(p.i)} – {label(p.j)}
                                </td>
                                <td>
                                  {p.meanA.toFixed(2)} ± {p.stdA.toFixed(2)}
                                </td>
                                <td>
                                  {p.meanB.toFixed(2)} ± {p.stdB.toFixed(2)}
                                </td>
                                <td>{p.delta.toFixed(2)} nm</td>
                                <td>{p.betweenShare.toFixed(2)}</td>
                                <td>{p.contactSwitch ? 'yes' : '–'}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                      <div className="mdc-col">
                        {local ? (
                          <>
                            <div
                              className="mdc-note"
                              title="Residues whose structural letter differs between the two: a change of the local fold, where the switch pairs also catch rigid parts moving against each other. Divergence 1 means no letter in common."
                            >
                              <b>Local switch residues (3Di)</b>
                            </div>
                            {local.length ? (
                              <table className="mdc-table">
                                <thead>
                                  <tr>
                                    <th>Residue</th>
                                    <th>{nameA}</th>
                                    <th>{nameB}</th>
                                    <th>Divergence</th>
                                  </tr>
                                </thead>
                                <tbody>
                                  {local.map((l) => (
                                    <tr
                                      key={l.residue}
                                      onClick={() => {
                                        const i = toActive(l.residue);
                                        if (i >= 0) session.focus([i]);
                                      }}
                                      onMouseEnter={() =>
                                        session.highlight([toActive(l.residue)].filter((x) => x >= 0))
                                      }
                                      onMouseLeave={() => session.highlight(undefined)}
                                    >
                                      <td>{label(l.residue)}</td>
                                      <td>
                                        {l.letterA} ({(l.shareA * 100).toFixed(0)}%)
                                      </td>
                                      <td>
                                        {l.letterB} ({(l.shareB * 100).toFixed(0)}%)
                                      </td>
                                      <td>{l.divergence.toFixed(2)}</td>
                                    </tr>
                                  ))}
                                </tbody>
                              </table>
                            ) : (
                              <div className="mdc-muted">
                                No residue changes its 3Di letter between the two: the difference is rigid-body.
                              </div>
                            )}
                          </>
                        ) : (
                          <div className="mdc-note">
                            Local switch residues need 3Di letters for every simulation compared (the server makes them
                            when it has Foldseek).
                          </div>
                        )}
                      </div>
                    </div>
                  )}
                </>
              ) : (
                <div className="mdc-muted">Pick two different conformers.</div>
              )}
            </div>
          )}
        </>
      )}
    </div>
  );
}
