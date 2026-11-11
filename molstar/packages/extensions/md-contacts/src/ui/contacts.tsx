/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useCallback, useMemo, useState } from 'react';
import { frameDistances, pairSeries, residueKey, residueLabel } from '../analysis.js';
import type { MdContactsSession, Simulation } from '../session.js';
import { ns, shortGroup, simColor, timesNs, useSubject } from './common.js';
import { Heatmap, type MapLayer } from './heatmap.js';
import { LinePlot, type Series } from './plot.js';
import { SequenceView } from './sequence.js';

type LayerKind = 'frame' | 'mean' | 'std' | 'contact' | 'dmean' | 'dstart' | 'gmx' | 'repstd' | `dgroup:${string}`;

function sameResidues(a: Simulation, b: Simulation) {
  if (a.data.nRes !== b.data.nRes) return false;
  for (let i = 0; i < a.data.nRes; i++)
    if (residueKey(a.data.residues[i]) !== residueKey(b.data.residues[i])) return false;
  return true;
}

function subtract(a: ArrayLike<number>, b: ArrayLike<number>) {
  const out = new Float32Array(a.length);
  for (let k = 0; k < out.length; k++) out[k] = a[k] - b[k];
  return out;
}

/** The contact map of the frame shown in 3D beside its statistics; hover highlights the pair, click follows it. */
export function ContactsTab({ session }: { session: MdContactsSession }) {
  const sim = useSubject(session.active);
  const frame = useSubject(session.frame);
  const pair = useSubject(session.pair);
  const settings = useSubject(session.settings);
  const revision = useSubject(session.revision);
  const [upper, setUpper] = useState<LayerKind>('frame');
  const [lower, setLower] = useState<LayerKind>('mean');
  const [scope, setScope] = useState<'replica' | 'group'>('replica');
  const [dmax, setDmax] = useState(3);
  const [drange, setDrange] = useState(0.5);

  const summary = useMemo(() => (sim ? session.groupSummary(sim.group) : undefined), [sim, revision, session]);
  const pooled = scope === 'group' && summary && sim && summary.sims.includes(sim) && summary.sims.length > 1;
  const stats = pooled ? summary!.stats : sim?.stats;
  const otherGroups = useMemo(() => {
    if (!sim) return [];
    return session.groups.filter((g) => {
      if (g === sim.group) return false;
      const s = session.groupSummary(g);
      return !!s && sameResidues(s.sims[0], sim);
    });
  }, [sim, revision, session]);

  const frameMap = useMemo(
    () => (sim ? frameDistances(sim.data, Math.min(frame, sim.data.nFrames - 1)) : undefined),
    [sim, frame],
  );
  const startMap = useMemo(() => (sim ? frameDistances(sim.data, sim.stability.first) : undefined), [sim, revision]);

  const build = useCallback(
    (kind: LayerKind): MapLayer => {
      const nm = 'nm';
      if (!sim || !stats || !frameMap || !startMap)
        return { label: '', unit: '', values: undefined, scale: 'distance', min: 0, max: 1 };
      const where = pooled ? `${summary!.sims.length} replicas` : sim.label;
      switch (kind) {
        case 'frame':
          return { label: 'This frame', unit: nm, values: frameMap, scale: 'distance', min: 0, max: dmax };
        case 'mean':
          return { label: `Mean (${where})`, unit: nm, values: stats.mean, scale: 'distance', min: 0, max: dmax };
        case 'std':
          return {
            label: `Std over time (${where})`,
            unit: nm,
            values: stats.std,
            scale: 'spread',
            min: 0,
            max: drange,
          };
        case 'contact':
          return {
            label: `Time in contact, < ${stats.cutoff} nm (${where})`,
            unit: '',
            values: stats.contact,
            scale: 'share',
            min: 0,
            max: 1,
          };
        case 'dmean':
          return {
            label: 'This frame − mean',
            unit: nm,
            values: subtract(frameMap, stats.mean),
            scale: 'diverging',
            min: -drange,
            max: drange,
          };
        case 'dstart':
          return {
            label: 'This frame − start',
            unit: nm,
            values: subtract(frameMap, startMap),
            scale: 'diverging',
            min: -drange,
            max: drange,
          };
        case 'gmx':
          return {
            label: 'gmx mdmat mean smallest distance',
            unit: nm,
            values: sim.gmxMean,
            scale: 'distance',
            min: 0,
            max: Math.min(dmax, sim.gmxInfo?.max ?? dmax),
          };
        case 'repstd':
          return {
            label: 'Disagreement between replicas (std of their means)',
            unit: nm,
            values: summary?.agreement?.betweenStd,
            scale: 'spread',
            min: 0,
            max: drange,
          };
        default: {
          const g = kind.slice('dgroup:'.length);
          const other = session.groupSummary(g);
          return {
            label: `Mean − mean of ${shortGroup(g)}`,
            unit: nm,
            values: other && other.stats.n === stats.n ? subtract(stats.mean, other.stats.mean) : undefined,
            scale: 'diverging',
            min: -drange * 2,
            max: drange * 2,
          };
        }
      }
    },
    [sim, stats, frameMap, startMap, dmax, drange, pooled, summary, session],
  );

  const label = useCallback((i: number) => (sim ? residueLabel(sim.data.residues[i]) : ''), [sim]);
  const tick = useCallback((i: number) => (sim ? sim.data.residues[i].seq : i), [sim]);

  const pairPlot = useMemo(() => {
    if (!sim || !pair || pair.j >= sim.data.nRes) return undefined;
    const sims = session.simsOf(sim.group).filter((s) => sameResidues(s, sim));
    const series: Series[] = sims.map((s) => ({
      key: s.id,
      label: s.label,
      color: simColor(s),
      x: timesNs(s.data.times),
      y: pairSeries(s.data, pair.i, pair.j),
      emphasis: s === sim,
    }));
    return series;
  }, [sim, pair, revision, session]);

  if (!sim || !stats) return <div className="mdc-empty">Load a simulation in the Data tab to see its contact map.</div>;

  const options: [LayerKind, string][] = [
    ['frame', 'This frame'],
    ['mean', 'Mean over time'],
    ['std', 'Std over time'],
    ['contact', 'Time in contact'],
    ['dmean', 'This frame − mean'],
    ['dstart', 'This frame − start'],
  ];
  if (sim.gmxMean) options.push(['gmx', 'gmx mdmat mean']);
  if (summary?.agreement && summary.sims.includes(sim)) options.push(['repstd', 'Replica disagreement']);
  for (const g of otherGroups) options.push([`dgroup:${g}`, `Mean − ${shortGroup(g)}`]);
  const kinds = options.map((o) => o[0]);
  const up = kinds.includes(upper) ? upper : 'frame',
    lo = kinds.includes(lower) ? lower : 'mean';

  const upperLayer = build(up);
  const lowerLayer = build(lo);
  const n = sim.data.nRes;
  const k = pair && pair.j < n ? pair.i * n + pair.j : -1;

  return (
    <div className="mdc-tab">
      <div className="mdc-controls">
        <label>
          <span>◥ Upper triangle</span>
          <select value={up} onChange={(e) => setUpper(e.target.value as LayerKind)}>
            {options.map(([v, l]) => (
              <option key={v} value={v}>
                {l}
              </option>
            ))}
          </select>
        </label>
        <label>
          <span>◣ Lower triangle</span>
          <select value={lo} onChange={(e) => setLower(e.target.value as LayerKind)}>
            {options.map(([v, l]) => (
              <option key={v} value={v}>
                {l}
              </option>
            ))}
          </select>
        </label>
        {summary && summary.sims.length > 1 && summary.sims.includes(sim) && (
          <label>
            <span>Statistics over</span>
            <select value={scope} onChange={(e) => setScope(e.target.value as 'replica' | 'group')}>
              <option value="replica">this replica</option>
              <option value="group">all {summary.sims.length} replicas</option>
            </select>
          </label>
        )}
        <label>
          <span>Distance scale, to {dmax.toFixed(1)} nm</span>
          <input
            type="range"
            min={1}
            max={8}
            step={0.5}
            value={dmax}
            onChange={(e) => setDmax(Number(e.target.value))}
          />
        </label>
        <label>
          <span>Std / difference scale, {drange.toFixed(2)} nm</span>
          <input
            type="range"
            min={0.1}
            max={2}
            step={0.05}
            value={drange}
            onChange={(e) => setDrange(Number(e.target.value))}
          />
        </label>
      </div>

      <Heatmap
        n={n}
        upper={upperLayer}
        lower={lowerLayer}
        label={label}
        tick={tick}
        selected={pair && pair.j < n ? pair : undefined}
        onHover={(p) => session.highlight(p ? [p.i, p.j] : undefined)}
        onSelect={(p) => void session.selectPair(p)}
      />
      <div className="mdc-note">
        {n} residues, Cα–Cα distances, {stats.count} frames in the statistics
        {settings.skipPs > 0 ? ` (first ${ns(settings.skipPs)} ns left out)` : ''}. Hover the map to highlight a pair in
        3D, click to follow its distance.
      </div>

      <SequenceView session={session} sim={sim} frameMap={frameMap!} startMap={startMap!} stats={stats} />

      {pair && pairPlot && k >= 0 ? (
        <div className="mdc-card">
          <div className="mdc-card-head">
            <span>
              <b>{label(pair.i)}</b> – <b>{label(pair.j)}</b>
            </span>
            <span className="mdc-row">
              <button type="button" className="mdc-btn" onClick={() => session.focus([pair.i, pair.j])}>
                Focus in 3D
              </button>
              <button type="button" className="mdc-btn" onClick={() => void session.selectPair(undefined)}>
                Clear
              </button>
            </span>
          </div>
          <div className="mdc-stats">
            <span>
              now <b>{frameMap![k].toFixed(2)} nm</b>
            </span>
            <span>
              mean <b>{stats.mean[k].toFixed(2)}</b> ± {stats.std[k].toFixed(2)} nm
            </span>
            <span>
              in contact <b>{(stats.contact[k] * 100).toFixed(0)}%</b> of the time
            </span>
            {sim.gmxMean && (
              <span>
                gmx mdmat <b>{sim.gmxMean[k].toFixed(2)} nm</b>
              </span>
            )}
          </div>
          <LinePlot
            title="Distance over time"
            series={pairPlot}
            xLabel="ns"
            yUnit="nm"
            cursor={sim.data.times[frame] / 1000}
            shadeTo={settings.skipPs > 0 ? (sim.data.times[0] + settings.skipPs) / 1000 : undefined}
            guides={[{ y: settings.cutoff, label: `contact < ${settings.cutoff} nm` }]}
            xFormat={(x) => x.toFixed(x < 10 ? 1 : 0)}
            onPick={(_, nearest) => nearest && session.setFrame(nearest.index)}
          />
        </div>
      ) : (
        <div className="mdc-card mdc-muted">
          No pair selected. Click a cell of the map to plot that distance over time.
        </div>
      )}
    </div>
  );
}
