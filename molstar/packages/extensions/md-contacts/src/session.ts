/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { BehaviorSubject, type Subscription } from 'rxjs';
import { Interval, OrderedSet, SortedArray } from '@molstar/core/data/int';
import type { StateObjectSelector } from '@molstar/core/state';
import { EmptyLoci } from '@molstar/model/model/loci';
import type { Coordinates } from '@molstar/model/model/structure/coordinates';
import type { Model } from '@molstar/model/model/structure';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import type { PluginContext } from '@molstar/plugin/context';
import {
  type ConformerOptions,
  type ConformerResult,
  commonResidues,
  drmsdSeries,
  evaluateStability,
  findConformers,
  type MapStats,
  MinFramesForVerdict,
  mapStats,
  poolStats,
  type ReplicaAgreement,
  type ResidueInfo,
  replicaAgreement,
  type StabilityResult,
  type TrajectoryData,
} from './analysis.js';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Color } from '@molstar/core/util/color';
import { alignResidues, type SequenceAlignment, superpose } from './align.js';
import { seriesColor } from './colors.js';
import { SmoothTrajectory } from './interpolate.js';
import { MdContactsRmsfThemeName, MdContactsRmsfThemeProvider, ResidueValues } from './theme.js';
import { parseThreeDi, summarizeThreeDi, type ThreeDi, type ThreeDiSummary } from './threedi.js';
import type { XpmMatrix } from './xpm.js';

/** Mol*'s computed secondary structure property (`SecondaryStructureProvider.descriptor.name`) */
const SecondaryStructureName = 'molstar_computed_secondary_structure';

/** Where the bytes of a simulation come from: URLs (the companion server) or files the user dropped. */
export interface SimulationSource {
  id: string;
  /** simulations of one group are replicas of the same system (a state) */
  group: string;
  label: string;
  structure: { url?: string; data?: string; format: 'pdb' | 'gro' };
  trajectory: { url?: string; data?: Uint8Array<ArrayBuffer>; format: 'xtc' | 'trr' | 'dcd' };
}

export interface Simulation {
  id: string;
  group: string;
  label: string;
  /** categorical colour slot */
  slot: number;
  data: TrajectoryData;
  /** number of residues of the structure, including the ones without a CA (what `gmx mdmat` counts) */
  totalResidues: number;
  /** index among all residues of each analysed residue */
  allIndex: Int32Array;
  /** model residue index of each analysed residue */
  modelResidue: Int32Array;
  stats: MapStats;
  stability: StabilityResult;
  /** `gmx mdmat` mean smallest-distance map, subset to the analysed residues (nm) */
  gmxMean?: Float32Array;
  gmxInfo?: { title: string; max: number; epsUrl?: string; xpmUrl?: string };
  /** the trajectory in the 3Di structural alphabet (Foldseek), when the server made it or a file was dropped */
  threeDi?: ThreeDi;
  threeDiSummary?: ThreeDiSummary;
  refs: { model: string; coords: string; trajectory: string; frameModel: string; structure: string };
}

export interface Settings {
  /** contact cutoff, nm (representative atoms) */
  cutoff: number;
  /** ps left out at the start of every trajectory as equilibration */
  skipPs: number;
}

export interface GroupSummary {
  group: string;
  sims: Simulation[];
  stats: MapStats;
  agreement?: ReplicaAgreement;
}

export interface ResiduePair {
  i: number;
  j: number;
}

export interface ConformerRun {
  scope: string;
  sims: Simulation[];
  /** residues shared by the simulations, and each simulation's index of them */
  residues: ResidueInfo[];
  maps: Int32Array[];
  result: ConformerResult;
}

function extractTrajectory(model: Model, coords: Coordinates) {
  const h = model.atomicHierarchy;
  const nRes = h.residues._rowCount;
  const offsets = h.residueAtomSegments.offsets;
  const chainIndex = h.chainAtomSegments.index;
  const atomName = h.atoms.label_atom_id;
  const compId = h.atoms.label_comp_id;
  const all: (ResidueInfo & { model: number; all: number })[] = [];
  for (let r = 0; r < nRes; r++) {
    let ca = -1;
    for (let a = offsets[r]; a < offsets[r + 1]; a++) {
      if (atomName.value(a) === 'CA') {
        ca = a;
        break;
      }
    }
    const first = offsets[r];
    all.push({
      seq: h.residues.auth_seq_id.value(r),
      name: compId.value(first),
      chain: h.chains.auth_asym_id.value(chainIndex[first]),
      atom: h.atomSourceIndex.value(ca >= 0 ? ca : first),
      hasCA: ca >= 0,
      model: r,
      all: r,
    });
  }
  // residues in input order (what GROMACS counts), CA-bearing ones analysed
  all.sort((a, b) => a.atom - b.atom);
  all.forEach((r, i) => (r.all = i));
  let kept = all.filter((r) => r.hasCA);
  if (kept.length < 3) kept = all;
  const residues: ResidueInfo[] = kept.map(({ seq, name, chain, atom, hasCA }) => ({ seq, name, chain, atom, hasCA }));

  const frames = coords.frames;
  const nFrames = frames.length;
  const n = residues.length;
  const xyz = new Float32Array(nFrames * n * 3);
  const times = new Float32Array(nFrames);
  let maxAtom = 0;
  for (const r of residues) if (r.atom > maxAtom) maxAtom = r.atom;
  for (let f = 0; f < nFrames; f++) {
    const fr = frames[f];
    if (fr.elementCount <= maxAtom)
      throw new Error(
        `The trajectory has ${fr.elementCount} atoms per frame, the structure needs at least ${maxAtom + 1}: they do not belong together.`,
      );
    const o = f * n * 3;
    for (let i = 0; i < n; i++) {
      const a = residues[i].atom;
      // Mol* keeps Å, GROMACS and this analysis nm
      xyz[o + i * 3] = fr.x[a] * 0.1;
      xyz[o + i * 3 + 1] = fr.y[a] * 0.1;
      xyz[o + i * 3 + 2] = fr.z[a] * 0.1;
    }
    // GROMACS readers give ps here (the XTC adapter labels the unit 'step', the value is the file's time)
    times[f] = fr.time.value;
  }
  // a file without times (or one frame per step): fall back to the frame index
  if (nFrames > 1 && !(times[nFrames - 1] > times[0])) for (let f = 0; f < nFrames; f++) times[f] = f;
  const data: TrajectoryData = { nFrames, nRes: n, times, xyz, residues };
  return {
    data,
    totalResidues: all.length,
    allIndex: Int32Array.from(kept, (r) => r.all),
    modelResidue: Int32Array.from(kept, (r) => r.model),
  };
}

/** The state of one analysis session: the loaded simulations, the frame shown, the pair selected, and the results. */
export class MdContactsSession {
  readonly simulations = new BehaviorSubject<Simulation[]>([]);
  readonly active = new BehaviorSubject<Simulation | undefined>(undefined);
  readonly frame = new BehaviorSubject<number>(0);
  readonly playing = new BehaviorSubject<boolean>(false);
  readonly pair = new BehaviorSubject<ResiduePair | undefined>(undefined);
  readonly settings = new BehaviorSubject<Settings>({ cutoff: 0.8, skipPs: 0 });
  readonly conformers = new BehaviorSubject<ConformerRun | undefined>(undefined);
  /** bumps whenever derived results change */
  readonly revision = new BehaviorSubject<number>(0);
  readonly busy = new BehaviorSubject<string>('');
  /** how the 3D structures are coloured: Mol*'s default, or by each simulation's per-residue RMSF */
  readonly coloring = new BehaviorSubject<'default' | 'rmsf'>('default');
  /** nm: the RMSF drawn darkest, shared by all simulations so their colours compare */
  readonly rmsfMax = new BehaviorSubject<number>(0.4);

  private subs: Subscription[] = [];
  private groupCache = new Map<string, GroupSummary>();
  private stateDistanceCache = new Map<string, { group: string; series: Float32Array }[]>();
  private nextSlot = 0;
  /** position to show, in steps of the interpolated trajectory (`frame * smooth` and the in-between ones) */
  private pendingFrame: number | undefined;
  private position = 0;
  /** frames shown per frame of the file: 1 is off, more draws linearly interpolated frames in between */
  readonly smooth = new BehaviorSubject<number>(1);
  /** a frame of another simulation compared with the one shown, and drawn superposed on it in 3D when `show` */
  readonly reference = new BehaviorSubject<{ sim: Simulation; frame: number; show: boolean } | undefined>(undefined);
  /** the conformers shown together in 3D (their medoid frames, superposed on the first), when that view is on */
  readonly overlay = new BehaviorSubject<
    { conformer: number; sim: Simulation; frame: number; rmsd: number }[] | undefined
  >(undefined);
  private alignments = new Map<string, SequenceAlignment>();
  private refitTimer: ReturnType<typeof setTimeout> | undefined;
  private committing = false;
  private timer: ReturnType<typeof setTimeout> | undefined;
  private measurement: StateObjectSelector | undefined;
  /** steps advanced per tick while playing (frames of the file when not smoothing) */
  playStep = 1;
  playIntervalMs = 40;

  constructor(readonly plugin: PluginContext) {
    // the superposed reference follows the frame shown, once it rests (refitting rebuilds the reference's shapes)
    const scheduleRefit = () => {
      if (this.refitTimer) clearTimeout(this.refitTimer);
      if (!this.reference.value?.show || this.playing.value) return;
      this.refitTimer = setTimeout(() => void this.refitReference(), 200);
    };
    this.subs.push(this.frame.subscribe(scheduleRefit), this.playing.subscribe(scheduleRefit));
    // the view of all conformers belongs to one run: a new or dropped run ends it
    this.subs.push(
      this.conformers.subscribe(() => {
        if (this.overlay.value || this.overlayRefs.length) void this.showConformers(false);
      }),
    );
    const themes = plugin.representation.structure.themes.colorThemeRegistry;
    if (!themes.has(MdContactsRmsfThemeProvider)) themes.add(MdContactsRmsfThemeProvider);
    // follow frame changes made by Mol*'s own trajectory controls
    this.subs.push(
      plugin.state.data.events.changed.subscribe(() => {
        const sim = this.active.value;
        if (!sim || this.committing) return;
        const cell = plugin.state.data.cells.get(sim.refs.frameModel);
        const idx = (cell?.transform.params as { modelIndex?: number } | undefined)?.modelIndex;
        if (typeof idx !== 'number' || this.pendingFrame !== undefined || idx === this.position) return;
        this.position = idx;
        const f = Math.min(sim.data.nFrames - 1, Math.round(idx / this.smooth.value));
        if (f !== this.frame.value) this.frame.next(idx);
      }),
    );
  }

  dispose() {
    if (this.refitTimer) clearTimeout(this.refitTimer);
    this.pause();
    for (const s of this.subs) s.unsubscribe();
  }

  private bump() {
    this.groupCache.clear();
    this.stateDistanceCache.clear();
    this.revision.next(this.revision.value + 1);
  }

  get groups() {
    const out: string[] = [];
    for (const s of this.simulations.value) if (!out.includes(s.group)) out.push(s.group);
    return out;
  }

  simsOf(group: string) {
    return this.simulations.value.filter((s) => s.group === group);
  }

  private analyse(data: TrajectoryData) {
    const { cutoff, skipPs } = this.settings.value;
    let first = 0;
    if (skipPs > 0) while (first < data.nFrames - 2 && data.times[first] - data.times[0] < skipPs) first++;
    return {
      stats: mapStats(data, cutoff, { first }),
      stability: evaluateStability(data, { skipPs, cutoff }),
    };
  }

  async add(source: SimulationSource) {
    const plugin = this.plugin;
    await plugin.initialized;
    if (this.simulations.value.some((s) => s.id === source.id)) throw new Error(`${source.label} is already loaded.`);
    this.busy.next(`Loading ${source.label}`);
    const created: string[] = [];
    try {
      const modelData =
        source.structure.data !== undefined
          ? await plugin.builders.data.rawData({ data: source.structure.data, label: `${source.label} structure` })
          : await plugin.builders.data.download({ url: source.structure.url!, label: `${source.label} structure` });
      created.push(modelData.ref);
      const topoTraj = await plugin.builders.structure.parseTrajectory(modelData, source.structure.format);
      const topoModel = await plugin.builders.structure.createModel(topoTraj);

      const coordData =
        source.trajectory.data !== undefined
          ? await plugin.builders.data.rawData({ data: source.trajectory.data, label: `${source.label} trajectory` })
          : await plugin.builders.data.download({
              url: source.trajectory.url!,
              isBinary: true,
              label: `${source.label} trajectory`,
            });
      created.push(coordData.ref);
      const provider = plugin.dataFormats.get(source.trajectory.format);
      if (!provider) throw new Error(`No reader for .${source.trajectory.format} trajectories.`);
      const coords = await provider.parse(plugin, coordData);

      const model = topoModel.cell?.obj?.data as Model | undefined;
      const coordinates = plugin.state.data.cells.get(coords.ref)?.obj?.data as Coordinates | undefined;
      if (!model || !coordinates) throw new Error('The structure or the trajectory could not be read.');
      if (!coordinates.frames.length) throw new Error('The trajectory has no frames.');
      const extracted = extractTrajectory(model, coordinates);

      const trajectory = await plugin
        .build()
        .toRoot()
        .apply(
          SmoothTrajectory,
          { modelRef: topoModel.ref, coordinatesRef: coords.ref, steps: this.smooth.value },
          { dependsOn: [topoModel.ref, coords.ref] },
        )
        .commit();
      created.push(trajectory.ref);
      const preset = await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
        showUnitcell: false,
      });
      if (!preset?.model || !preset.structure) throw new Error('The structure could not be shown.');
      await this.useModelSecondaryStructure(preset.structureProperties?.ref);

      await new Promise((r) => setTimeout(r, 0));
      const { stats, stability } = this.analyse(extracted.data);
      const sim: Simulation = {
        id: source.id,
        group: source.group,
        label: source.label,
        slot: this.nextSlot++,
        ...extracted,
        stats,
        stability,
        refs: {
          model: modelData.ref,
          coords: coordData.ref,
          trajectory: trajectory.ref,
          frameModel: preset.model.ref,
          structure: preset.structure.ref,
        },
      };
      this.simulations.next(
        [...this.simulations.value, sim].sort((x, y) => x.id.localeCompare(y.id, undefined, { numeric: true })),
      );
      this.conformers.next(undefined);
      this.bump();
      if (!this.active.value) await this.setActive(sim);
      else setSubtreeVisibility(plugin.state.data, sim.refs.trajectory, true);
      if (this.coloring.value === 'rmsf') await this.setColoring('rmsf');
      return sim;
    } catch (e) {
      const b = plugin.build();
      for (const ref of created) if (plugin.state.data.cells.has(ref)) b.delete(ref);
      await b.commit();
      throw e;
    } finally {
      this.busy.next('');
    }
  }

  async remove(sim: Simulation) {
    const b = this.plugin.build();
    for (const ref of [sim.refs.trajectory, sim.refs.model, sim.refs.coords])
      if (this.plugin.state.data.cells.has(ref)) b.delete(ref);
    await b.commit();
    const rest = this.simulations.value.filter((s) => s !== sim);
    ResidueValues.delete(sim.id);
    if (this.reference.value?.sim === sim) this.reference.next(undefined);
    this.simulations.next(rest);
    this.conformers.next(undefined);
    this.bump();
    if (this.active.value === sim) {
      this.active.next(undefined);
      this.pair.next(undefined);
      if (rest.length) await this.setActive(rest[0]);
    }
  }

  /** Attaches a `gmx mdmat` map (all residues of the group, in order) to a simulation. */
  attachGmxMap(sim: Simulation, xpm: XpmMatrix, info?: { epsUrl?: string; xpmUrl?: string }) {
    const n = sim.data.nRes;
    let pick: ArrayLike<number> | undefined;
    if (xpm.width === sim.totalResidues && xpm.height === sim.totalResidues) pick = sim.allIndex;
    else if (xpm.width === n && xpm.height === n) pick = Int32Array.from({ length: n }, (_, i) => i);
    if (!pick)
      throw new Error(
        `The XPM map is ${xpm.width} x ${xpm.height}, the structure has ${sim.totalResidues} residues (${n} with a CA): it belongs to another selection.`,
      );
    const out = new Float32Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) out[i * n + j] = xpm.values[pick[i] * xpm.width + pick[j]];
    }
    sim.gmxMean = out;
    sim.gmxInfo = { title: xpm.title || 'gmx mdmat', max: xpm.max, ...info };
    this.revision.next(this.revision.value + 1);
  }

  /** Attaches the 3Di letters of a simulation's frames (one line per frame, as the server writes them). */
  attachThreeDi(sim: Simulation, text: string) {
    sim.threeDi = parseThreeDi(text, sim.data.nFrames, sim.data.nRes, sim.totalResidues, sim.allIndex);
    sim.threeDiSummary = summarizeThreeDi(sim.threeDi, sim.stability.first);
    this.revision.next(this.revision.value + 1);
  }

  updateSettings(patch: Partial<Settings>) {
    const next = { ...this.settings.value, ...patch };
    this.settings.next(next);
    for (const sim of this.simulations.value) {
      const { stats, stability } = this.analyse(sim.data);
      sim.stats = stats;
      sim.stability = stability;
      if (sim.threeDi) sim.threeDiSummary = summarizeThreeDi(sim.threeDi, stability.first);
    }
    this.conformers.next(undefined);
    this.bump();
    if (this.coloring.value === 'rmsf') void this.setColoring('rmsf');
  }

  groupSummary(group: string): GroupSummary | undefined {
    const cached = this.groupCache.get(group);
    if (cached) return cached;
    const sims = this.simsOf(group);
    if (!sims.length) return undefined;
    const same = sims.every((s) => s.data.nRes === sims[0].data.nRes);
    const usable = same ? sims : [sims[0]];
    // replicas too short to judge do not count for the agreement
    const judged = usable.filter((s) => s.stats.count >= MinFramesForVerdict);
    const summary: GroupSummary = {
      group,
      sims: usable,
      stats: usable.length > 1 ? poolStats(usable.map((s) => s.stats)) : usable[0].stats,
      agreement: judged.length > 1 ? replicaAgreement(judged.map((s) => s.stats)) : undefined,
    };
    this.groupCache.set(group, summary);
    return summary;
  }

  /** dRMSD of every frame of `sim` to the mean map of every group (on the residues they share). */
  stateDistances(sim: Simulation) {
    const cached = this.stateDistanceCache.get(sim.id);
    if (cached) return cached;
    const out: { group: string; series: Float32Array }[] = [];
    for (const group of this.groups) {
      const summary = this.groupSummary(group);
      if (!summary) continue;
      const other = summary.sims[0];
      const { residues, maps } = commonResidues([sim.data, other.data]);
      if (residues.length < 10) continue;
      const m = residues.length,
        n = other.data.nRes;
      const ref = new Float32Array(m * m);
      for (let i = 0; i < m; i++)
        for (let j = 0; j < m; j++) ref[i * m + j] = summary.stats.mean[maps[1][i] * n + maps[1][j]];
      out.push({ group, series: drmsdSeries(sim.data, ref, m, maps[0]) });
    }
    this.stateDistanceCache.set(sim.id, out);
    return out;
  }

  // the 3D view

  async setActive(sim: Simulation) {
    const state = this.plugin.state.data;
    if (this.reference.value?.sim === sim) await this.setReference(undefined);
    if (this.overlay.value) await this.showConformers(false);
    const prev = this.active.value;
    this.active.next(sim);
    this.syncVisibility();
    if (prev && prev.data.nRes !== sim.data.nRes) this.pair.next(undefined);
    const cell = state.cells.get(sim.refs.frameModel);
    const idx = (cell?.transform.params as { modelIndex?: number } | undefined)?.modelIndex ?? 0;
    this.position = idx;
    this.frame.next(Math.min(Math.round(idx / this.smooth.value), sim.data.nFrames - 1));
    await this.syncMeasurement();
    if (this.reference.value) await this.refitReference();
    this.plugin.managers.camera.reset();
  }

  setFrame(frame: number) {
    const sim = this.active.value;
    if (!sim) return;
    // moving through the trajectory ends the view of all conformers
    if (this.overlay.value) void this.showConformers(false);
    const f = Math.max(0, Math.min(sim.data.nFrames - 1, Math.round(frame)));
    this.setPosition(f * this.smooth.value);
  }

  /** `position`: in steps of the interpolated trajectory; the frame of the analysis is the nearest one of the file. */
  private setPosition(position: number) {
    const sim = this.active.value;
    if (!sim) return;
    const k = this.smooth.value;
    const p = Math.max(0, Math.min((sim.data.nFrames - 1) * k, Math.round(position)));
    if (p === this.position && this.pendingFrame === undefined) return;
    this.position = p;
    const f = Math.round(p / k);
    if (f !== this.frame.value) this.frame.next(f);
    this.pendingFrame = p;
    void this.flushFrame();
  }

  /** Draws `steps - 1` interpolated frames between the frames of the file (1: off), for every loaded simulation. */
  async setSmooth(steps: number) {
    const k = Math.max(1, Math.round(steps));
    const old = this.smooth.value;
    if (k === old) return;
    const wasPlaying = this.playing.value;
    this.pause();
    while (this.committing) await new Promise((r) => setTimeout(r, 10));
    this.committing = true;
    try {
      const state = this.plugin.state.data;
      const b = this.plugin.build();
      for (const sim of this.simulations.value) {
        const idx =
          (state.cells.get(sim.refs.frameModel)?.transform.params as { modelIndex?: number } | undefined)?.modelIndex ??
          0;
        const ref = this.reference.value;
        const frame =
          sim === this.active.value ? this.frame.value : ref?.sim === sim ? ref.frame : Math.round(idx / old);
        b.to(sim.refs.trajectory).update((p: any) => ({ ...p, steps: k }));
        b.to(sim.refs.frameModel).update({ modelIndex: Math.min(frame, sim.data.nFrames - 1) * k });
      }
      await b.commit();
      this.smooth.next(k);
      this.position = this.frame.value * k;
      // the structures were made again: hidden ones stay hidden, the measurement and colours follow
      this.syncVisibility();
    } finally {
      this.committing = false;
    }
    await this.syncMeasurement();
    if (this.coloring.value === 'rmsf') await this.setColoring('rmsf');
    if (wasPlaying) this.play();
  }

  private async flushFrame() {
    if (this.committing) return;
    this.committing = true;
    try {
      while (this.pendingFrame !== undefined) {
        const sim = this.active.value;
        const f = this.pendingFrame;
        this.pendingFrame = undefined;
        if (!sim) break;
        await this.plugin.build().to(sim.refs.frameModel).update({ modelIndex: f }).commit();
      }
    } finally {
      this.committing = false;
    }
  }

  /** Shows a frame of any loaded simulation (activating it when needed). */
  async show(sim: Simulation, frame: number) {
    if (this.overlay.value) await this.showConformers(false);
    if (this.active.value !== sim) await this.setActive(sim);
    this.setFrame(frame);
  }

  play() {
    if (this.playing.value) return;
    if (this.overlay.value) void this.showConformers(false);
    this.playing.next(true);
    const tick = async () => {
      if (!this.playing.value) return;
      const sim = this.active.value;
      if (sim) {
        let next = this.position + this.playStep;
        if (next > (sim.data.nFrames - 1) * this.smooth.value) next = 0;
        this.setPosition(next);
        await this.flushFrame();
      }
      this.timer = setTimeout(tick, this.playIntervalMs);
    };
    void tick();
  }

  pause() {
    this.playing.next(false);
    if (this.timer) clearTimeout(this.timer);
    this.timer = undefined;
  }

  private structureOf(sim: Simulation) {
    return this.plugin.state.data.cells.get(sim.refs.structure)?.obj?.data as Structure | undefined;
  }

  /** Loci of residues (whole residues) or of their representative atoms in the current frame's structure. */
  loci(sim: Simulation, residues: number[], atomOnly = false): StructureElement.Loci | undefined {
    const structure = this.structureOf(sim);
    if (!structure) return undefined;
    const elements: {
      unit: StructureElement.Loci['elements'][0]['unit'];
      indices: OrderedSet<StructureElement.UnitIndex>;
    }[] = [];
    for (const unit of structure.units) {
      const h = unit.model.atomicHierarchy;
      const offsets = h.residueAtomSegments.offsets;
      const idx: number[] = [];
      for (const r of residues) {
        if (r < 0 || r >= sim.data.nRes) continue;
        const mr = sim.modelResidue[r];
        if (atomOnly) {
          // the representative atom: by its source index
          for (let a = offsets[mr]; a < offsets[mr + 1]; a++) {
            if (h.atomSourceIndex.value(a) === sim.data.residues[r].atom) {
              const k = SortedArray.indexOf(unit.elements, a);
              if (k >= 0) idx.push(k);
              break;
            }
          }
        } else {
          const range = SortedArray.findRange(unit.elements, offsets[mr], offsets[mr + 1] - 1);
          for (let k = Interval.start(range); k < Interval.end(range); k++) idx.push(k);
        }
      }
      if (idx.length) {
        idx.sort((a, b) => a - b);
        elements.push({
          unit,
          indices: OrderedSet.ofSortedArray(idx as unknown as ArrayLike<StructureElement.UnitIndex>),
        });
      }
    }
    return elements.length ? StructureElement.Loci(structure, elements) : undefined;
  }

  highlight(residues: number[] | undefined) {
    const sim = this.active.value;
    const loci = sim && residues?.length ? this.loci(sim, residues) : undefined;
    if (!loci) this.plugin.managers.interactivity.lociHighlights.clearHighlights();
    else this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci }, false);
  }

  focus(residues: number[]) {
    const sim = this.active.value;
    const loci = sim ? this.loci(sim, residues) : undefined;
    if (loci) this.plugin.managers.camera.focusLoci(loci);
  }

  async selectPair(pair: ResiduePair | undefined) {
    if (pair && pair.i > pair.j) pair = { i: pair.j, j: pair.i };
    if (pair && pair.i === pair.j) pair = undefined;
    this.pair.next(pair);
    await this.syncMeasurement();
  }

  private async syncMeasurement() {
    const plugin = this.plugin;
    if (this.measurement) {
      const ref = this.measurement.ref;
      this.measurement = undefined;
      if (plugin.state.data.cells.has(ref)) await plugin.build().delete(ref).commit();
    }
    const sim = this.active.value,
      pair = this.pair.value;
    if (!sim || !pair) return;
    const a = this.loci(sim, [pair.i], true),
      b = this.loci(sim, [pair.j], true);
    if (!a || !b) return;
    const added = await plugin.managers.structure.measurement.addDistance(a, b);
    this.measurement = added?.selection;
  }

  clearHighlight() {
    this.plugin.managers.interactivity.lociHighlights.highlightOnly({ loci: EmptyLoci }, false);
  }

  // colouring

  private themeRevision = 0;

  /** Colours every loaded structure by its simulation's RMSF, or back by Mol*'s default. */
  async setColoring(mode: 'default' | 'rmsf') {
    this.coloring.next(mode);
    const sims = this.simulations.value;
    if (mode === 'rmsf') {
      // one scale for all: the 98th percentile of every RMSF value, rounded up to 0.05 nm
      const all: number[] = [];
      for (const sim of sims) {
        if (sim.stability.verdict === 'short') continue;
        for (const v of sim.stability.rmsf) all.push(v);
      }
      all.sort((a, b) => a - b);
      const p98 = all.length ? all[Math.floor((all.length - 1) * 0.98)] : 0.4;
      this.rmsfMax.next(Math.max(0.1, Math.ceil(p98 / 0.05) * 0.05));
      for (const sim of sims) {
        const model = this.structureOf(sim)?.models[0];
        if (!model) continue;
        const values = new Float32Array(model.atomicHierarchy.residues._rowCount).fill(Number.NaN);
        for (let i = 0; i < sim.data.nRes; i++) values[sim.modelResidue[i]] = sim.stability.rmsf[i];
        ResidueValues.set(sim.id, { values });
      }
    }
    this.themeRevision++;
    for (const sim of sims) await this.applyTheme(sim);
  }

  /** The colours of one simulation: grey as the superposed reference, else by RMSF or Mol*'s default. */
  private async applyTheme(sim: Simulation) {
    const ref = this.plugin.managers.structure.hierarchy.current.structures.find(
      (s) => s.cell.transform.ref === sim.refs.structure,
    );
    if (!ref) return;
    const params =
      this.reference.value?.sim === sim
        ? { color: 'uniform', colorParams: { value: Color(0x8f8e88) } }
        : this.coloring.value === 'rmsf'
          ? {
              color: MdContactsRmsfThemeName,
              colorParams: { sim: sim.id, max: this.rmsfMax.value, revision: this.themeRevision },
            }
          : { color: 'default' };
    await this.plugin.managers.structure.component.updateRepresentationsTheme(ref.components, params as any);
  }

  // comparing with a reference

  private syncVisibility() {
    const active = this.active.value,
      ref = this.reference.value,
      overlay = !!this.overlay.value;
    // while the conformers are shown together, the simulations themselves are out of the way
    for (const s of this.simulations.value)
      setSubtreeVisibility(
        this.plugin.state.data,
        s.refs.trajectory,
        overlay || !(s === active || (ref?.sim === s && ref.show)),
      );
  }

  /** The frames share one secondary structure (set on the model by the trajectory): take it, not DSSP per frame. */
  private async useModelSecondaryStructure(propertiesRef: string | undefined) {
    if (!propertiesRef) return;
    await this.plugin
      .build()
      .to(propertiesRef)
      .update((old: any) => ({
        ...old,
        autoAttach: old.autoAttach.includes(SecondaryStructureName)
          ? old.autoAttach
          : [...old.autoAttach, SecondaryStructureName],
        properties: { ...old.properties, [SecondaryStructureName]: { type: { name: 'model', params: {} } } },
      }))
      .commit();
  }

  // all conformers at once

  private overlayRefs: string[] = [];

  /**
   * Shows the most typical frame of every conformer of the current run together, each in its conformer's colour and
   * superposed on the first conformer's over the aligned Cα atoms; the simulations are hidden meanwhile. `false`
   * goes back to the trajectory.
   */
  async showConformers(on: boolean) {
    const plugin = this.plugin;
    const cells = plugin.state.data.cells;
    if (this.overlayRefs.length) {
      const b = plugin.build();
      for (const ref of this.overlayRefs) if (cells.has(ref)) b.delete(ref);
      this.overlayRefs = [];
      await b.commit();
    }
    const run = this.conformers.value;
    if (!on || !run) {
      if (this.overlay.value) this.overlay.next(undefined);
      this.syncVisibility();
      return;
    }
    this.pause();
    this.busy.next('Superposing the conformers');
    try {
      const medoids = run.result.conformers.map((c) => ({ sim: run.sims[c.medoid.traj], frame: c.medoid.frame }));
      const target = medoids[0];
      const shown: { conformer: number; sim: Simulation; frame: number; rmsd: number }[] = [];
      for (let c = 0; c < medoids.length; c++) {
        const { sim, frame } = medoids[c];
        const source = cells.get(sim.refs.trajectory)?.transform.params as
          | { modelRef: string; coordinatesRef: string }
          | undefined;
        if (!source) continue;
        const { pairs } = this.alignment(target.sim, sim);
        const fit = superpose(target.sim.data, target.frame, sim.data, frame, pairs);
        const m = Mat4.clone(fit.transform);
        m[12] *= 10;
        m[13] *= 10;
        m[14] *= 10;
        // its own copy of the simulation's trajectory (frames are made on demand, so this holds no coordinates)
        const trajectory = await plugin
          .build()
          .toRoot()
          .apply(
            SmoothTrajectory,
            { modelRef: source.modelRef, coordinatesRef: source.coordinatesRef, steps: 1, transform: m },
            { dependsOn: [source.modelRef, source.coordinatesRef] },
          )
          .commit();
        this.overlayRefs.push(trajectory.ref);
        const preset = await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
          showUnitcell: false,
        });
        if (!preset?.model || !preset.structure) continue;
        await plugin.build().to(preset.model.ref).update({ modelIndex: frame }).commit();
        await this.useModelSecondaryStructure(preset.structureProperties?.ref);
        const ref = plugin.managers.structure.hierarchy.current.structures.find(
          (x) => x.cell.transform.ref === preset.structure!.ref,
        );
        if (ref)
          await plugin.managers.structure.component.updateRepresentationsTheme(ref.components, {
            color: 'uniform',
            colorParams: { value: Color(Number.parseInt(seriesColor(c).slice(1), 16)) },
          } as any);
        shown.push({ conformer: c, sim, frame, rmsd: fit.rmsd });
      }
      this.overlay.next(shown);
      this.syncVisibility();
      plugin.managers.camera.reset();
    } finally {
      this.busy.next('');
    }
  }

  /** The residue alignment of two simulations (kept). */
  alignment(a: Simulation, b: Simulation) {
    const key = `${a.id}\n${b.id}`;
    let out = this.alignments.get(key);
    if (!out) {
      out = alignResidues(a.data.residues, b.data.residues);
      this.alignments.set(key, out);
    }
    return out;
  }

  /**
   * Compares the simulation shown with a frame of another one; with `show`, that frame is drawn in grey, superposed
   * on the frame shown over the aligned residues. `undefined` ends the comparison.
   */
  async setReference(sim: Simulation | undefined, frame = 0, show = true) {
    const prev = this.reference.value;
    if (sim && sim === this.active.value) sim = undefined;
    const next = sim ? { sim, frame: Math.max(0, Math.min(sim.data.nFrames - 1, Math.round(frame))), show } : undefined;
    this.reference.next(next);
    if (prev && prev.sim !== next?.sim && this.simulations.value.includes(prev.sim)) await this.applyTheme(prev.sim);
    if (next) {
      if (prev?.sim !== next.sim) await this.applyTheme(next.sim);
      await this.refitReference();
    }
    this.syncVisibility();
  }

  /** Puts the reference at its frame, superposed on the frame shown. */
  private async refitReference() {
    const ref = this.reference.value,
      active = this.active.value;
    if (!ref || !active) return;
    const cells = this.plugin.state.data.cells;
    if (!cells.has(ref.sim.refs.trajectory)) return;
    const { pairs } = this.alignment(active, ref.sim);
    const fit = superpose(active.data, this.frame.value, ref.sim.data, ref.frame, pairs);
    // nm to Å: the rotation stays, the translation scales
    const m = Mat4.clone(fit.transform);
    m[12] *= 10;
    m[13] *= 10;
    m[14] *= 10;
    const b = this.plugin.build();
    const old = cells.get(ref.sim.refs.trajectory)?.transform.params as { transform?: Mat4 } | undefined;
    if (ref.show && !(old?.transform && Mat4.areEqual(old.transform, m, 1e-6)))
      b.to(ref.sim.refs.trajectory).update((p: any) => ({ ...p, transform: m }));
    b.to(ref.sim.refs.frameModel).update({ modelIndex: ref.frame * this.smooth.value });
    await b.commit();
    this.syncVisibility();
  }

  // conformers

  async runConformers(scope: string, options: Partial<ConformerOptions>) {
    const sims = scope === '*' ? this.simulations.value : this.simsOf(scope);
    if (!sims.length) throw new Error('Nothing loaded to analyse.');
    const { residues, maps, identical } = commonResidues(sims.map((s) => s.data));
    if (residues.length < 10) throw new Error('The simulations share too few residues to compare.');
    this.busy.next('Looking for conformers');
    try {
      const result = await findConformers(
        sims.map((s) => s.data),
        { cutoff: this.settings.value.cutoff, skipPs: this.settings.value.skipPs, ...options },
        {
          maps: identical ? undefined : maps,
          letters: sims.map((s) => s.threeDi?.letters),
          progress: async (message, fraction) => {
            this.busy.next(`${message} ${(fraction * 100).toFixed(0)}%`);
            await new Promise((r) => setTimeout(r, 0));
          },
        },
      );
      const run: ConformerRun = { scope, sims, residues, maps, result };
      this.conformers.next(run);
      return run;
    } finally {
      this.busy.next('');
    }
  }
}
