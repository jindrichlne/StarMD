/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { FoundSimulation, Job, MdServer, ProcessParams, ScanResult, ServerInfo } from '../server-api.js';
import type { MdContactsSession, Simulation, SimulationSource } from '../session.js';
import { parseXpm } from '../xpm.js';
import { bytes, Help, ns, Sections, shortGroup, simColor, useSection, useSubject } from './common.js';
import { VerdictBadge } from './stability.js';

interface RowState {
  job?: Job;
  status: 'idle' | 'working' | 'loading' | 'loaded' | 'error';
  message?: string;
}

const STRUCT = /\.(pdb|gro)$/i;
const TRAJ = /\.(xtc|trr|dcd)$/i;
const ext = (name: string) => name.slice(name.lastIndexOf('.') + 1).toLowerCase();
const stem = (name: string) => name.slice(0, name.lastIndexOf('.'));

/** Finding and loading simulations: through the companion server (which runs GROMACS), or from prepared files. */
export function DataTab({ session, server }: { session: MdContactsSession; server?: MdServer }) {
  const sims = useSubject(session.simulations);
  const active = useSubject(session.active);
  const settings = useSubject(session.settings);
  useSubject(session.revision);

  const [section, setSection] = useSection('data', ['find', 'loaded'] as const);
  const [info, setInfo] = useState<ServerInfo | undefined>();
  const [serverError, setServerError] = useState('');
  const [dir, setDir] = useState('');
  const [scan, setScan] = useState<ScanResult | undefined>();
  const [scanning, setScanning] = useState(false);
  const [checked, setChecked] = useState<Set<string>>(new Set());
  const [showEq, setShowEq] = useState(false);
  const [rows, setRows] = useState<Record<string, RowState>>({});
  const [params, setParams] = useState<Partial<ProcessParams>>({});
  const [advanced, setAdvanced] = useState(false);
  const [dropMessage, setDropMessage] = useState('');
  const [dragging, setDragging] = useState(false);
  const [skipNs, setSkipNs] = useState(String(settings.skipPs / 1000));
  const [cutoff, setCutoff] = useState(String(settings.cutoff));
  const addQueue = useRef<Promise<unknown>>(Promise.resolve());
  const alive = useRef(true);
  useEffect(() => () => void (alive.current = false), []);

  const setRow = useCallback((id: string, patch: Partial<RowState>) => {
    setRows((r) => ({ ...r, [id]: { ...(r[id] ?? { status: 'idle' }), ...patch } }));
  }, []);

  const doScan = useCallback(
    async (d: string) => {
      if (!server) return;
      setScanning(true);
      setServerError('');
      try {
        const res = await server.scan(d);
        setScan(res);
        // production runs with a topology, pre-selected
        setChecked(new Set(res.simulations.filter((s) => !s.equilibration && (s.tpr || s.structure)).map((s) => s.id)));
      } catch (e) {
        setServerError(String((e as Error).message ?? e));
      } finally {
        setScanning(false);
      }
    },
    [server],
  );

  useEffect(() => {
    if (!server) return;
    server
      .info()
      .then((i) => {
        setInfo(i);
        setParams(i.defaults);
        void doScan('');
      })
      .catch(() => setServerError('The companion server is not reachable: only dropped files can be loaded.'));
  }, [server, doScan]);

  const enqueueAdd = useCallback(
    (source: SimulationSource) => {
      const p = addQueue.current.then(() => session.add(source));
      addQueue.current = p.catch(() => undefined);
      return p;
    },
    [session],
  );

  const labelOf = useCallback((s: FoundSimulation, selected: FoundSimulation[]) => {
    const sameDir = selected.filter((o) => o.group === s.group && o.replica === s.replica);
    return sameDir.length > 1 && s.replica !== s.stage ? `${s.replica}/${s.stage}` : s.replica;
  }, []);

  const attachMap = useCallback(
    async (sim: Simulation, job: Job) => {
      if (!server || !job.files.mdmat) return;
      const text = await server.text(server.url(job.files.mdmat));
      const xpm = parseXpm(text)[0];
      if (!xpm) throw new Error('The gmx mdmat output could not be read.');
      session.attachGmxMap(sim, xpm, {
        xpmUrl: server.url(job.files.mdmat),
        epsUrl: job.files.mdmatEps ? server.url(job.files.mdmatEps) : undefined,
      });
    },
    [server, session],
  );

  const loadOne = useCallback(
    async (s: FoundSimulation, label: string) => {
      if (!server) return;
      setRow(s.id, { status: 'working', message: 'queued', job: undefined });
      try {
        let job = await server.process(s, params);
        let sim: Simulation | undefined;
        let adding: Promise<void> | undefined;
        for (;;) {
          if (!alive.current) return;
          setRow(s.id, {
            job,
            status: sim ? 'loaded' : 'working',
            message:
              job.status === 'done'
                ? undefined
                : `${job.stage}${job.progress > 0 && job.progress < 1 ? ` ${(job.progress * 100).toFixed(0)}%` : ''}`,
          });
          if (job.status === 'error') throw new Error(job.error ?? 'GROMACS failed.');
          if (job.trajectoryReady && !adding) {
            const files = job.files;
            adding = enqueueAdd({
              id: s.id,
              group: s.group,
              label,
              structure: { url: server.url(files.structure!), format: 'pdb' },
              trajectory: { url: server.url(files.trajectory!), format: 'xtc' },
            }).then((added) => {
              sim = added;
            });
            // surface a load failure without waiting for mdmat
            adding.catch(() => undefined);
          }
          if (job.status === 'done') break;
          await new Promise((r) => setTimeout(r, 700));
          job = await server.job(job.id);
        }
        await adding;
        if (sim && job.files.threeDi) session.attachThreeDi(sim, await server.text(server.url(job.files.threeDi)));
        if (sim && job.files.mdmat) await attachMap(sim, job);
        setRow(s.id, { job, status: 'loaded', message: undefined });
      } catch (e) {
        setRow(s.id, { status: 'error', message: String((e as Error).message ?? e) });
      }
    },
    [server, params, enqueueAdd, attachMap, setRow, session],
  );

  const loadAsIs = useCallback(
    async (s: FoundSimulation, label: string) => {
      if (!server || !s.structure) return;
      setRow(s.id, { status: 'loading', message: 'reading files' });
      try {
        await enqueueAdd({
          id: s.id,
          group: s.group,
          label,
          structure: { url: server.fileUrl(s.structure), format: ext(s.structure) as 'pdb' | 'gro' },
          trajectory: { url: server.fileUrl(s.trajectory), format: ext(s.trajectory) as 'xtc' | 'trr' },
        });
        setRow(s.id, { status: 'loaded', message: undefined });
      } catch (e) {
        setRow(s.id, { status: 'error', message: String((e as Error).message ?? e) });
      }
    },
    [server, enqueueAdd, setRow],
  );

  const visible = useMemo(() => (scan?.simulations ?? []).filter((s) => showEq || !s.equilibration), [scan, showEq]);
  const byGroup = useMemo(() => {
    const m = new Map<string, FoundSimulation[]>();
    for (const s of visible) m.set(s.group, [...(m.get(s.group) ?? []), s]);
    return [...m.entries()];
  }, [visible]);
  const selected = visible.filter((s) => checked.has(s.id));
  const loadedIds = new Set(sims.map((s) => s.id));
  const pending = selected.filter(
    (s) => !loadedIds.has(s.id) && !['working', 'loading'].includes(rows[s.id]?.status ?? 'idle'),
  );

  const handleFiles = useCallback(
    async (files: File[]) => {
      setDropMessage('');
      const structures = files.filter((f) => STRUCT.test(f.name)),
        trajs = files.filter((f) => TRAJ.test(f.name)),
        xpms = files.filter((f) => /\.xpm$/i.test(f.name)),
        threeDis = files.filter((f) => /3di.*\.(txt|fasta)$|\.3di$/i.test(f.name));
      const notes: string[] = [];
      for (const t of trajs) {
        const s =
          structures.find((f) => stem(f.name) === stem(t.name)) ??
          (structures.length === 1 ? structures[0] : undefined);
        if (!s) {
          notes.push(`${t.name}: drop its structure (.pdb or .gro with the same atoms) together with it.`);
          continue;
        }
        const rel = (t as File & { webkitRelativePath?: string }).webkitRelativePath || '';
        const parts = rel.split('/').filter(Boolean);
        try {
          await enqueueAdd({
            id: `dropped:${rel || t.name}:${t.size}`,
            group: parts.length > 2 ? parts.slice(0, -2).join('/') : 'dropped files',
            label: parts.length > 1 ? parts[parts.length - 2] : stem(t.name),
            structure: { data: await s.text(), format: ext(s.name) as 'pdb' | 'gro' },
            trajectory: { data: new Uint8Array(await t.arrayBuffer()), format: ext(t.name) as 'xtc' | 'trr' | 'dcd' },
          });
        } catch (e) {
          notes.push(`${t.name}: ${String((e as Error).message ?? e)}`);
        }
      }
      for (const x of xpms) {
        const target = session.active.value;
        if (!target) {
          notes.push(`${x.name}: load a trajectory first, the map is attached to the simulation shown in 3D.`);
          continue;
        }
        try {
          const m = parseXpm(await x.text())[0];
          if (!m) throw new Error('not a GROMACS XPM matrix');
          session.attachGmxMap(target, m);
          notes.push(`${x.name}: attached to ${target.label} (see "gmx mdmat mean" in the Contacts tab).`);
        } catch (e) {
          notes.push(`${x.name}: ${String((e as Error).message ?? e)}`);
        }
      }
      for (const x of threeDis) {
        const target = session.active.value;
        if (!target) {
          notes.push(`${x.name}: load a trajectory first, the letters are attached to the simulation shown in 3D.`);
          continue;
        }
        try {
          // one line of letters per frame; FASTA headers are skipped
          const text = (await x.text())
            .split(/\r?\n/)
            .filter((l) => l && !l.startsWith('>'))
            .join('\n');
          session.attachThreeDi(target, text);
          notes.push(`${x.name}: 3Di letters attached to ${target.label}.`);
        } catch (e) {
          notes.push(`${x.name}: ${String((e as Error).message ?? e)}`);
        }
      }
      if (!trajs.length && !xpms.length && !threeDis.length)
        notes.push(
          'Nothing to load: drop a trajectory (.xtc, .trr, .dcd) with its structure (.pdb, .gro), or an .xpm map.',
        );
      setDropMessage(notes.join('\n'));
    },
    [enqueueAdd, session],
  );

  const num = (k: keyof ProcessParams) => (e: React.ChangeEvent<HTMLInputElement>) =>
    setParams((p) => ({ ...p, [k]: e.target.value === '' ? (k === 'end' ? -1 : 0) : Number(e.target.value) }));

  return (
    <div className="mdc-tab">
      <Sections
        label="Data sections"
        options={[
          ['find', 'Add simulations'],
          ['loaded', `Loaded (${sims.length}) and settings`],
        ]}
        value={section}
        onChange={setSection}
      />
      {section === 'find' && server && (
        <div className="mdc-card">
          <div className="mdc-card-head">
            <b>Simulations on this computer</b>
            {info && (
              <span className="mdc-muted">
                {info.gmx ? `GROMACS ${info.gmx.version}` : 'GROMACS not found: prepared files only'}
                {info.foldseek ? ' · Foldseek (3Di)' : ' · no Foldseek: no 3Di letters'}
              </span>
            )}
          </div>
          {serverError && <div className="mdc-error">{serverError}</div>}
          {info && (
            <>
              <div className="mdc-row">
                <span className="mdc-path" title={info.root}>
                  {info.root}/
                </span>
                <input
                  type="text"
                  aria-label="Subfolder to scan"
                  placeholder="subfolder (empty: everything)"
                  value={dir}
                  onChange={(e) => setDir(e.target.value)}
                  onKeyDown={(e) => e.key === 'Enter' && void doScan(dir)}
                />
                <button type="button" className="mdc-btn" disabled={scanning} onClick={() => void doScan(dir)}>
                  {scanning ? 'Scanning' : 'Scan'}
                </button>
              </div>

              <div className="mdc-controls">
                <label title="Index group written out and analysed (gmx default groups of the .tpr, or of your index file)">
                  <span>Group</span>
                  <input
                    type="text"
                    value={params.group ?? ''}
                    onChange={(e) => setParams((p) => ({ ...p, group: e.target.value }))}
                  />
                </label>
                <label title="Time between the frames kept. 100 ps gives 1000 frames per 100 ns.">
                  <span>Frame every (ps)</span>
                  <input type="number" min={0} step={10} value={params.dt ?? 0} onChange={num('dt')} />
                </label>
                <label title="gmx mdmat cuts every distance above this value (-t). Its default of 1.5 nm hides large motions.">
                  <span>mdmat cut-off (nm)</span>
                  <input type="number" min={0.5} step={0.5} value={params.truncate ?? 4} onChange={num('truncate')} />
                </label>
                <button type="button" className="mdc-btn" onClick={() => setAdvanced(!advanced)}>
                  {advanced ? 'Fewer options' : 'More options'}
                </button>
              </div>
              {advanced && (
                <div className="mdc-controls">
                  <label>
                    <span>From (ps)</span>
                    <input type="number" min={0} value={params.begin ?? 0} onChange={num('begin')} />
                  </label>
                  <label>
                    <span>To (ps, empty: end)</span>
                    <input
                      type="number"
                      min={0}
                      value={params.end !== undefined && params.end >= 0 ? params.end : ''}
                      onChange={num('end')}
                    />
                  </label>
                  <label title="How trjconv makes the molecules whole (-pbc)">
                    <span>Periodic images</span>
                    <select
                      value={params.pbc ?? 'mol'}
                      onChange={(e) => setParams((p) => ({ ...p, pbc: e.target.value as ProcessParams['pbc'] }))}
                    >
                      {['mol', 'whole', 'cluster', 'nojump', 'none'].map((v) => (
                        <option key={v}>{v}</option>
                      ))}
                    </select>
                  </label>
                  <label title="Index file relative to the folder above, when the group is not a default one">
                    <span>Index file</span>
                    <input
                      type="text"
                      placeholder="e.g. state1/index.ndx"
                      value={params.index ?? ''}
                      onChange={(e) => setParams((p) => ({ ...p, index: e.target.value }))}
                    />
                  </label>
                  <label className="mdc-check">
                    <input
                      type="checkbox"
                      checked={params.fit ?? true}
                      onChange={(e) => setParams((p) => ({ ...p, fit: e.target.checked }))}
                    />
                    <span>Fit to the first frame (backbone)</span>
                  </label>
                  <label className="mdc-check">
                    <input
                      type="checkbox"
                      checked={params.mdmat ?? true}
                      onChange={(e) => setParams((p) => ({ ...p, mdmat: e.target.checked }))}
                    />
                    <span>Run gmx mdmat and xpm2ps</span>
                  </label>
                </div>
              )}

              {scan && !scan.simulations.length && (
                <div className="mdc-muted">No trajectories (.xtc, .trr) found here.</div>
              )}
              {byGroup.map(([group, list]) => (
                <div className="mdc-group" key={group}>
                  <div className="mdc-group-head">
                    <label className="mdc-check">
                      <input
                        type="checkbox"
                        aria-label={`Select all of ${group}`}
                        checked={list.every((s) => checked.has(s.id))}
                        onChange={(e) => {
                          const next = new Set(checked);
                          for (const s of list) e.target.checked ? next.add(s.id) : next.delete(s.id);
                          setChecked(next);
                        }}
                      />
                      <b title={group}>{group}</b>
                    </label>
                    <span className="mdc-muted">
                      {list.length} trajector{list.length > 1 ? 'ies' : 'y'}
                    </span>
                  </div>
                  {list.map((s) => {
                    const row = rows[s.id];
                    const loaded = loadedIds.has(s.id);
                    const label = labelOf(s, selected.includes(s) ? selected : [...selected, s]);
                    const canProcess = !!info.gmx && !!(s.tpr || s.structure);
                    return (
                      <div className="mdc-sim-row" key={s.id}>
                        <label className="mdc-check" title={s.trajectory}>
                          <input
                            type="checkbox"
                            checked={checked.has(s.id)}
                            onChange={(e) => {
                              const next = new Set(checked);
                              e.target.checked ? next.add(s.id) : next.delete(s.id);
                              setChecked(next);
                            }}
                          />
                          <span>
                            {s.replica === s.stage ? s.stage : `${s.replica} / ${s.stage}`}
                            <span className="mdc-muted"> · {bytes(s.size)}</span>
                          </span>
                        </label>
                        <span className="mdc-sim-status">
                          {row?.status === 'error' ? (
                            <span className="mdc-error-inline" title={row.message}>
                              ✕ {row.message}
                            </span>
                          ) : loaded ? (
                            <span>✓ loaded{row?.message ? ` · ${row.message}` : ''}</span>
                          ) : row?.status === 'working' || row?.status === 'loading' ? (
                            <span>
                              <span className="mdc-spinner" /> {row.message}
                            </span>
                          ) : !s.tpr && !s.structure ? (
                            <span className="mdc-muted">no .tpr or structure beside it</span>
                          ) : (
                            <>
                              {!s.tpr && <span className="mdc-muted">no .tpr</span>}
                              {s.structure && s.size < 4e8 && (
                                <button
                                  type="button"
                                  className="mdc-btn"
                                  title={`Load ${s.trajectory} with ${s.structure} without running GROMACS (already prepared: protein only, whole, fitted)`}
                                  onClick={() => void loadAsIs(s, label)}
                                >
                                  Load as is
                                </button>
                              )}
                              {canProcess && (
                                <button type="button" className="mdc-btn" onClick={() => void loadOne(s, label)}>
                                  Prepare & load
                                </button>
                              )}
                            </>
                          )}
                        </span>
                      </div>
                    );
                  })}
                </div>
              ))}
              {scan && (
                <div className="mdc-row mdc-row-between">
                  <label className="mdc-check">
                    <input type="checkbox" checked={showEq} onChange={(e) => setShowEq(e.target.checked)} />
                    <span>Show equilibration stages ({scan.simulations.filter((s) => s.equilibration).length})</span>
                  </label>
                  <button
                    type="button"
                    className="mdc-btn mdc-btn-primary"
                    disabled={!pending.length || !info.gmx}
                    onClick={() => {
                      for (const s of pending) void loadOne(s, labelOf(s, selected));
                    }}
                  >
                    Prepare & load {pending.length || ''} selected
                  </button>
                </div>
              )}
              <Help title="What Prepare does">
                "Prepare" runs <code>gmx trjconv</code> (group only, molecules whole, centred, fitted, thinned), then{' '}
                <code>gmx mdmat</code> and <code>gmx xpm2ps</code>. Results are cached in <code>{info.work}</code>; the
                originals are never changed.
              </Help>
            </>
          )}
        </div>
      )}

      {section === 'find' && (
        <div
          className={`mdc-card mdc-drop${dragging ? ' mdc-drop-over' : ''}`}
          onDragOver={(e) => {
            e.preventDefault();
            setDragging(true);
          }}
          onDragLeave={() => setDragging(false)}
          onDrop={(e) => {
            e.preventDefault();
            setDragging(false);
            void handleFiles([...e.dataTransfer.files]);
          }}
        >
          <div className="mdc-card-head">
            <b>Already prepared files</b>
            <label className="mdc-btn">
              Choose files
              <input
                type="file"
                multiple
                accept=".pdb,.gro,.xtc,.trr,.dcd,.xpm,.txt,.fasta,.3di"
                style={{ display: 'none' }}
                onChange={(e) => {
                  void handleFiles([...(e.target.files ?? [])]);
                  e.target.value = '';
                }}
              />
            </label>
          </div>
          <div className="mdc-muted">
            Drop trajectories (.xtc, .trr, .dcd) with their structure (.pdb or .gro of the same atoms; same file name,
            or one structure for all), and <code>gmx mdmat</code> maps (.xpm) for the simulation shown in 3D. Nothing
            leaves your browser.
          </div>
          {dropMessage && <pre className="mdc-drop-message">{dropMessage}</pre>}
        </div>
      )}

      {section === 'loaded' && (
        <div className="mdc-card">
          <div className="mdc-card-head">
            <b>Analysis settings</b>
          </div>
          <div className="mdc-controls">
            <label title="Left out of the mean maps, the stability measures and the conformers">
              <span>Skip at the start (ns)</span>
              <input type="number" min={0} step={1} value={skipNs} onChange={(e) => setSkipNs(e.target.value)} />
            </label>
            <label title="Two residues are in contact when their Cα atoms are closer than this">
              <span>Contact cut-off (nm)</span>
              <input type="number" min={0.3} step={0.05} value={cutoff} onChange={(e) => setCutoff(e.target.value)} />
            </label>
            <button
              type="button"
              className="mdc-btn"
              disabled={Number(skipNs) * 1000 === settings.skipPs && Number(cutoff) === settings.cutoff}
              onClick={() =>
                session.updateSettings({
                  skipPs: Math.max(0, Number(skipNs) || 0) * 1000,
                  cutoff: Math.max(0.1, Number(cutoff) || 0.8),
                })
              }
            >
              Apply
            </button>
          </div>
        </div>
      )}

      {section === 'loaded' && sims.length > 0 && (
        <div className="mdc-card">
          <div className="mdc-card-head">
            <b>Loaded</b>
            <span className="mdc-muted">
              {sims.length} simulation{sims.length > 1 ? 's' : ''}
            </span>
          </div>
          <table className="mdc-table">
            <tbody>
              {sims.map((s) => (
                <tr
                  key={s.id}
                  className={s === active ? 'mdc-active' : undefined}
                  onClick={() => void session.setActive(s)}
                >
                  <td title={s.id}>
                    <span className="mdc-dot" style={{ background: simColor(s) }} />
                    {shortGroup(s.group)} / {s.label}
                  </td>
                  <td>
                    {s.data.nFrames} frames, {ns(s.data.times[s.data.nFrames - 1] - s.data.times[0], 0)} ns,{' '}
                    {s.data.nRes} residues
                  </td>
                  <td>
                    <VerdictBadge verdict={s.stability.verdict} />
                  </td>
                  <td onClick={(e) => e.stopPropagation()}>
                    {s.threeDi && <span title="3Di letters of every frame are loaded">3Di</span>}{' '}
                    {s.gmxInfo?.epsUrl && (
                      <a
                        className="mdc-link"
                        href={s.gmxInfo.epsUrl}
                        download
                        title="gmx xpm2ps output of the mdmat mean map"
                      >
                        EPS
                      </a>
                    )}{' '}
                    {s.gmxInfo?.xpmUrl && (
                      <a className="mdc-link" href={s.gmxInfo.xpmUrl} download title="gmx mdmat mean map">
                        XPM
                      </a>
                    )}{' '}
                    <button
                      type="button"
                      className="mdc-btn"
                      aria-label={`Remove ${s.label}`}
                      onClick={() => void session.remove(s)}
                    >
                      Remove
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
