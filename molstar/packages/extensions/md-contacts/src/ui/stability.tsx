/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useMemo } from 'react';
import { DefaultStabilityThresholds, residueLabel, type StabilityVerdict } from '../analysis.js';
import { cssGradient, seriesColor, STATUS } from '../colors.js';
import type { MdContactsSession, Simulation } from '../session.js';
import { Help, ns, Sections, shortGroup, simColor, timesNs, useSection, useSubject } from './common.js';
import { LinePlot, type Series } from './plot.js';
import { ThreeDiMap } from './threedi-map.js';

const VERDICT: Record<StabilityVerdict, { label: string; icon: string; color: string; text: string }> = {
  holds: { label: 'Holds', icon: '✓', color: STATUS.good, text: 'The starting shape is kept.' },
  shifted: {
    label: 'Moved, settled',
    icon: '→',
    color: STATUS.serious,
    text: 'Left the starting shape early and settled somewhere else.',
  },
  drifting: {
    label: 'Drifting',
    icon: '≈',
    color: STATUS.warning,
    text: 'Still changing: extend the run before concluding.',
  },
  departed: { label: 'Departed', icon: '✕', color: STATUS.critical, text: 'The starting shape was left.' },
  short: { label: 'Too short', icon: '?', color: '#85847f', text: 'Not judged.' },
};

export function VerdictBadge({ verdict }: { verdict: StabilityVerdict }) {
  const v = VERDICT[verdict];
  return (
    <span className="mdc-badge" style={{ borderColor: v.color }}>
      <i style={{ background: v.color }}>{v.icon}</i>
      {v.label}
    </span>
  );
}

function groupVerdict(all: Simulation[]) {
  const sims = all.filter((s) => s.stability.verdict !== 'short');
  const skipped = all.length - sims.length;
  const note = skipped ? ` ${skipped} too short to judge.` : '';
  if (!sims.length) return 'Too short to judge.';
  return groupVerdictOf(sims) + note;
}

function groupVerdictOf(sims: Simulation[]) {
  const count = (v: StabilityVerdict) => sims.filter((s) => s.stability.verdict === v).length;
  const held = count('holds'),
    shifted = count('shifted'),
    drifting = count('drifting'),
    departed = count('departed');
  const of = `${sims.length} replica${sims.length > 1 ? 's' : ''}`;
  if (held === sims.length) return `Held in ${held} of ${of}: a plausible stable state.`;
  if (departed === sims.length) return `Left in all ${of}: not stable with this force field.`;
  if (held === 0 && drifting === 0) return `Moved away from the start in all ${of}: the starting shape is not kept.`;
  const parts = [
    held && `held in ${held}`,
    shifted && `moved and settled in ${shifted}`,
    drifting && `still drifting in ${drifting}`,
    departed && `left in ${departed}`,
  ].filter(Boolean);
  const text = `${parts.join(', ')} of ${of}`;
  return `${text[0].toUpperCase()}${text.slice(1)}: ${held === 0 ? 'the starting shape is not kept.' : 'not settled yet.'}`;
}

/** Does each replica keep the shape it started in, and do the replicas agree? */
export function StabilityTab({ session }: { session: MdContactsSession }) {
  const sims = useSubject(session.simulations);
  const active = useSubject(session.active);
  const frame = useSubject(session.frame);
  const settings = useSubject(session.settings);
  const revision = useSubject(session.revision);
  const coloring = useSubject(session.coloring);
  const rmsfMax = useSubject(session.rmsfMax);

  const group = active?.group;
  const groupSims = useMemo(() => (group ? session.simsOf(group) : []), [group, sims, session]);

  const metric = (pick: (s: Simulation) => ArrayLike<number>): Series[] =>
    groupSims.map((s) => ({
      key: s.id,
      label: s.label,
      color: simColor(s),
      x: timesNs(s.data.times),
      y: pick(s),
      emphasis: s === active,
    }));

  const rmsfSeries: Series[] = useMemo(
    () =>
      groupSims.map((s) => ({
        key: s.id,
        label: s.label,
        color: simColor(s),
        x: s.data.residues.map((r) => r.seq),
        y: s.stability.rmsf,
        emphasis: s === active,
      })),
    [groupSims, active, revision],
  );

  const stateSeries: Series[] = useMemo(() => {
    if (!active || session.groups.length < 2) return [];
    return session.stateDistances(active).map((d) => {
      const first = session.simsOf(d.group)[0];
      return {
        key: d.group,
        label: `to ${shortGroup(d.group)}${d.group === active.group ? ' (own)' : ''}`,
        color: first ? seriesColor(first.slot) : '#888',
        x: timesNs(active.data.times),
        y: d.series,
        emphasis: d.group === active.group,
        dashed: d.group !== active.group,
      };
    });
  }, [active, revision, session]);

  const has3Di = groupSims.some((x) => x.threeDiSummary);
  const sectionOptions = [
    ['verdict', 'Verdict'],
    ['time', 'Over time'],
    ['flex', 'Flexibility'],
    ...(has3Di ? ([['3di', 'Local structure (3Di)']] as const) : []),
  ] as const;
  const [section, setSection] = useSection(
    'stability',
    sectionOptions.map((o) => o[0]),
  );
  if (!active) return <div className="mdc-empty">Load a simulation in the Data tab to evaluate its stability.</div>;

  const cursor = active.data.times[frame] / 1000;
  const shadeTo = settings.skipPs > 0 ? (active.data.times[0] + settings.skipPs) / 1000 : undefined;
  const pick = (_: number, nearest: { index: number } | undefined) => nearest && session.setFrame(nearest.index);
  const xf = (x: number) => x.toFixed(x < 10 ? 1 : 0);
  const th = DefaultStabilityThresholds;
  const lastStates = stateSeries.map((s) => {
    const y = s.y;
    const from = Math.floor(y.length * 0.8);
    let m = 0;
    for (let i = from; i < y.length; i++) m += y[i];
    return { key: s.key, mean: m / Math.max(1, y.length - from) };
  });
  const closest = lastStates.length ? lastStates.reduce((a, b) => (b.mean < a.mean ? b : a)) : undefined;

  return (
    <div className="mdc-tab">
      <Sections label="Stability sections" options={sectionOptions} value={section} onChange={setSection}>
        {section !== 'verdict' && (
          <span className="mdc-muted">
            replicas of <b>{shortGroup(active.group)}</b>, heavier line: in 3D
          </span>
        )}
      </Sections>
      {section === 'verdict' &&
        session.groups.map((g) => {
          const gs = session.simsOf(g);
          const sum = session.groupSummary(g);
          return (
            <div className="mdc-card" key={g}>
              <div className="mdc-card-head">
                <b>{shortGroup(g)}</b>
                <span className="mdc-muted">{groupVerdict(gs)}</span>
              </div>
              <table className="mdc-table">
                <thead>
                  <tr>
                    <th>Replica</th>
                    <th>Verdict</th>
                    <th title="RMSD after superposition to the first analysed frame, mean of the last fifth">RMSD</th>
                    <th title="Superposition-free distance-map difference to the first analysed frame, mean of the last fifth">
                      dRMSD
                    </th>
                    <th title="Share of the starting Cα contacts still present, mean of the last fifth">
                      Contacts kept
                    </th>
                    <th title="Change of dRMSD over the second half of the run">Drift</th>
                    <th>Length</th>
                  </tr>
                </thead>
                <tbody>
                  {gs.map((s) => (
                    <tr
                      key={s.id}
                      className={s === active ? 'mdc-active' : undefined}
                      onClick={() => void session.setActive(s)}
                      title={s.stability.reasons.join('; ')}
                    >
                      <td>
                        <span className="mdc-dot" style={{ background: simColor(s) }} />
                        {s.label}
                      </td>
                      <td>
                        <VerdictBadge verdict={s.stability.verdict} />
                      </td>
                      <td>{s.stability.rmsdLate.toFixed(2)} nm</td>
                      <td>{s.stability.drmsdLate.toFixed(2)} nm</td>
                      <td>{(s.stability.qLate * 100).toFixed(0)}%</td>
                      <td>
                        {s.stability.drift >= 0 ? '+' : ''}
                        {s.stability.drift.toFixed(3)} nm
                      </td>
                      <td>{ns(s.data.times[s.data.nFrames - 1] - s.data.times[0], 0)} ns</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {sum?.agreement ? (
                <div className="mdc-reason">
                  <VerdictInline ok={sum.agreement.agree} />
                  Replicas {sum.agreement.agree ? 'agree' : 'disagree'}: their mean distance maps differ by up to{' '}
                  <b>{sum.agreement.maxPairDrmsd.toFixed(3)} nm</b> (dRMSD), and{' '}
                  <b>{(sum.agreement.disagreeingShare * 100).toFixed(1)}%</b> of the residue pairs sit at different
                  distances in different replicas.{' '}
                  {sum.agreement.agree
                    ? ''
                    : 'See "Replica disagreement" in the Contacts tab for where; extend the runs before concluding.'}
                </div>
              ) : (
                <div className="mdc-reason mdc-muted">
                  Fewer than two replicas long enough to compare: "stable" cannot be told from "stuck". Load at least
                  three replicas of this system.
                </div>
              )}
              <Help title="Why these verdicts">
                {gs.map((s) => (
                  <div className="mdc-reason" key={s.id}>
                    <b>{s.label}:</b> {VERDICT[s.stability.verdict].text} {s.stability.reasons.join('; ')}.
                  </div>
                ))}
              </Help>
            </div>
          );
        })}

      {section === 'verdict' && (
        <Help title="How the verdict is made">
          Click a replica to show it in 3D. Two things are judged apart. Is it still moving: the dRMSD changes by more
          than {th.drift} nm over the second half (drifting). How far from the start did it end up: the dRMSD of the
          last fifth is above {th.shifted} nm or fewer than {th.contactsHeld * 100}% of the starting contacts are kept
          (moved, settled), above {th.departed} nm or below {th.contactsLost * 100}% contacts (departed). Holds means
          neither. It is a heuristic: read the plots under "Over time".
        </Help>
      )}

      {section === 'time' && (
        <div className="mdc-plots">
          {stateSeries.length > 1 && (
            <div className="mdc-card">
              <LinePlot
                title={`Which state is ${active.label} in? dRMSD to each state's mean map`}
                series={stateSeries}
                xLabel="ns"
                yUnit="nm"
                cursor={cursor}
                shadeTo={shadeTo}
                yMin={0}
                xFormat={xf}
                onPick={(x) => {
                  const t = active.data.times;
                  let best = 0;
                  for (let i = 1; i < t.length; i++)
                    if (Math.abs(t[i] / 1000 - x) < Math.abs(t[best] / 1000 - x)) best = i;
                  session.setFrame(best);
                }}
              />
              {closest && (
                <div className="mdc-reason">
                  Over its last fifth, {active.label} is closest to <b>{shortGroup(closest.key)}</b>
                  {closest.key === active.group
                    ? ' (its own state): it did not cross over.'
                    : ': it moved towards another state.'}
                </div>
              )}
            </div>
          )}

          <LinePlot
            title="dRMSD to the start (superposition-free)"
            series={metric((s) => s.stability.drmsd)}
            xLabel="ns"
            yUnit="nm"
            cursor={cursor}
            shadeTo={shadeTo}
            yMin={0}
            xFormat={xf}
            onPick={pick}
          />
          <LinePlot
            title="Starting contacts kept"
            series={metric((s) => s.stability.q)}
            xLabel="ns"
            yUnit=""
            cursor={cursor}
            shadeTo={shadeTo}
            yMax={1}
            xFormat={xf}
            onPick={pick}
            guides={[{ y: th.contactsHeld, label: 'held' }]}
          />
          <LinePlot
            title="Cα RMSD to the start (after superposition)"
            series={metric((s) => s.stability.rmsd)}
            xLabel="ns"
            yUnit="nm"
            cursor={cursor}
            shadeTo={shadeTo}
            yMin={0}
            xFormat={xf}
            onPick={pick}
          />
          <LinePlot
            title="Radius of gyration (Cα)"
            series={metric((s) => s.stability.rg)}
            xLabel="ns"
            yUnit="nm"
            cursor={cursor}
            shadeTo={shadeTo}
            xFormat={xf}
            onPick={pick}
          />
        </div>
      )}
      {section === 'flex' && (
        <>
          <div className="mdc-card">
            <div className="mdc-card-head">
              <label className="mdc-check">
                <input
                  type="checkbox"
                  checked={coloring === 'rmsf'}
                  onChange={(e) => void session.setColoring(e.target.checked ? 'rmsf' : 'default')}
                />
                <b>Colour the 3D structure by flexibility (RMSF)</b>
              </label>
              {coloring === 'rmsf' && (
                <div className="mdc-legend-bar" style={{ minWidth: 220 }}>
                  <span>rigid 0</span>
                  <div style={{ background: cssGradient('spread') }} />
                  <span>{rmsfMax.toFixed(2)} nm and above</span>
                </div>
              )}
            </div>
          </div>
          <LinePlot
            title="Flexibility per residue (Cα RMSF)"
            series={rmsfSeries}
            xLabel="residue"
            yUnit="nm"
            yMin={0}
            xFormat={(x) => String(Math.round(x))}
            onHoverX={(x) => {
              if (x === undefined) return session.highlight(undefined);
              const res = active.data.residues;
              let best = 0;
              for (let i = 1; i < res.length; i++) if (Math.abs(res[i].seq - x) < Math.abs(res[best].seq - x)) best = i;
              session.highlight([best]);
            }}
            onPick={(x) => {
              const res = active.data.residues;
              let best = 0;
              for (let i = 1; i < res.length; i++) if (Math.abs(res[i].seq - x) < Math.abs(res[best].seq - x)) best = i;
              session.focus([best]);
            }}
          />
          <div className="mdc-note">
            Most flexible in {active.label}:{' '}
            {[...active.stability.rmsf]
              .map((v, i) => ({ v, i }))
              .sort((a, b) => b.v - a.v)
              .slice(0, 5)
              .map(({ v, i }) => `${residueLabel(active.data.residues[i])} (${v.toFixed(2)} nm)`)
              .join(', ')}
            . Hover the plot to highlight a residue in 3D, click to zoom to it. In 3D each replica is coloured by its
            own RMSF on one shared scale; residues without a Cα stay grey.
          </div>
        </>
      )}
      {section === '3di' && (
        <>
          <div className="mdc-plots">
            <LinePlot
              title="Local structure changed since the start (3Di letters)"
              series={groupSims
                .filter((x) => x.threeDiSummary)
                .map((x) => ({
                  key: x.id,
                  label: x.label,
                  color: simColor(x),
                  x: timesNs(x.data.times),
                  y: x.threeDiSummary!.changed,
                  emphasis: x === active,
                }))}
              xLabel="ns"
              yUnit="share of residues"
              cursor={cursor}
              shadeTo={shadeTo}
              yMin={0}
              xFormat={xf}
              onPick={pick}
            />
            <LinePlot
              title="Local variability per residue (3Di column entropy)"
              series={groupSims
                .filter((x) => x.threeDiSummary)
                .map((x) => ({
                  key: x.id,
                  label: x.label,
                  color: simColor(x),
                  x: x.data.residues.map((r) => r.seq),
                  y: x.threeDiSummary!.entropy,
                  emphasis: x === active,
                }))}
              xLabel="residue"
              yUnit="bits"
              yMin={0}
              xFormat={(x) => String(Math.round(x))}
            />
          </div>
          <ThreeDiMap session={session} sim={active} frame={frame} />
          <Help title="What 3Di letters are">
            3Di (Foldseek) gives every residue a letter for its local 3D environment, so these measures need no
            superposition and see local refolding. A rigid domain moving as a whole changes few letters: read them
            beside the distance-based plots, not instead of them.
          </Help>
        </>
      )}
      {!has3Di && section === 'flex' && (
        <div className="mdc-note">
          3Di structural letters are not loaded: the server makes them when it has Foldseek (see the README).
        </div>
      )}
    </div>
  );
}

function VerdictInline({ ok }: { ok: boolean }) {
  return (
    <span className="mdc-inline-icon" style={{ background: ok ? STATUS.good : STATUS.warning }}>
      {ok ? '✓' : '!'}
    </span>
  );
}
