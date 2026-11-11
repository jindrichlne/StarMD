/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { StateTransformer, type StateTransform } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Model } from '@molstar/model/model/structure';
import { Coordinates, type Frame, Time } from '@molstar/model/model/structure/coordinates';
import type { Trajectory } from '@molstar/model/model/structure/trajectory';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { setTrajectorySecondaryStructure } from './secondary-structure.js';

/**
 * A trajectory with `steps - 1` linearly interpolated frames between every two frames of the file, made when asked
 * for (and a few kept): playback moves smoothly without holding `steps` times the coordinates. Frame `i` of the file
 * is index `i * steps` here. Meant for trajectories made whole and fitted: atoms that jump across the periodic box
 * between two frames would be drawn flying through it.
 */
export class InterpolatedTrajectory implements Trajectory {
  readonly frameCount: number;
  readonly duration: number;
  readonly representative: Model;
  private cache = new Map<number, Model>();

  constructor(
    private model: Model,
    private coordinates: Coordinates,
    readonly steps: number,
    /** a rigid move of every frame (Å): a reference superposed on another simulation */
    private transform?: Mat4,
  ) {
    if (transform && Mat4.isIdentity(transform)) this.transform = undefined;
    const n = coordinates.frames.length;
    this.frameCount = n > 0 ? (n - 1) * steps + 1 : 0;
    this.duration = this.frameCount;
    this.representative = this.getFrameAtIndex(0);
  }

  private frame(index: number): Frame {
    const { frames } = this.coordinates;
    const a = Math.min(frames.length - 1, Math.floor(index / this.steps));
    const t = (index - a * this.steps) / this.steps;
    const fa = frames[a];
    const m = this.transform;
    const between = t !== 0 && a + 1 < frames.length;
    if (!between && !m) return fa;
    const fb = between ? frames[a + 1] : fa;
    const n = fa.elementCount;
    const x = new Float32Array(n),
      y = new Float32Array(n),
      z = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const px = fa.x[i] + (fb.x[i] - fa.x[i]) * t,
        py = fa.y[i] + (fb.y[i] - fa.y[i]) * t,
        pz = fa.z[i] + (fb.z[i] - fa.z[i]) * t;
      if (m) {
        x[i] = m[0] * px + m[4] * py + m[8] * pz + m[12];
        y[i] = m[1] * px + m[5] * py + m[9] * pz + m[13];
        z[i] = m[2] * px + m[6] * py + m[10] * pz + m[14];
      } else {
        x[i] = px;
        y[i] = py;
        z[i] = pz;
      }
    }
    return {
      elementCount: n,
      time: Time(fa.time.value + (fb.time.value - fa.time.value) * t, fa.time.unit),
      x,
      y,
      z,
      // a moved frame is no longer in its box
      cell: m ? undefined : fa.cell,
      box: m ? undefined : fa.box,
      xyzOrdering: fa.xyzOrdering,
    };
  }

  getFrameAtIndex(i: number): Model {
    const index = Math.max(0, Math.min(this.frameCount - 1, Math.round(i)));
    const cached = this.cache.get(index);
    if (cached) return cached;
    // the model's own way of making a frame (cell, periodic box, bonds), for this one frame
    const single = Coordinates.create([this.frame(index)], this.coordinates.deltaTime, this.coordinates.timeOffset);
    const made = Model.trajectoryFromModelAndCoordinates(this.model, single).representative;
    if (this.cache.size >= 48) this.cache.delete(this.cache.keys().next().value!);
    this.cache.set(index, made);
    return made;
  }
}

const Factory = StateTransformer.builderFactory('md-contacts-extension');

/** `TrajectoryFromModelAndCoordinates` with in-between frames: `steps` 1 is the file's frames as they are. */
export const SmoothTrajectory = Factory({
  name: 'smooth-trajectory',
  display: {
    name: 'Trajectory with in-between frames',
    description: 'Topology and coordinates, linearly interpolated between frames.',
  },
  from: SO.Root,
  to: SO.Molecule.Trajectory,
  params: {
    modelRef: PD.Text('', { isHidden: true }),
    coordinatesRef: PD.Text('', { isHidden: true }),
    steps: PD.Numeric(1, { min: 1, max: 16, step: 1 }, { description: 'Frames shown per frame of the file' }),
    transform: PD.Mat4(Mat4.identity(), { isHidden: true }),
  },
})({
  getDependencies: ({ modelRef, coordinatesRef }: { modelRef: string; coordinatesRef: string }) => {
    const deps: StateTransform.Ref[] = [];
    if (modelRef) deps.push(modelRef as StateTransform.Ref);
    if (coordinatesRef) deps.push(coordinatesRef as StateTransform.Ref);
    return deps;
  },
  apply({ params, dependencies }) {
    return Task.create('Trajectory with in-between frames', async () => {
      const model = dependencies![params.modelRef].data as Model;
      const coordinates = dependencies![params.coordinatesRef].data as Coordinates;
      // one secondary structure for every frame, or the cartoon flickers (see secondary-structure.ts)
      await setTrajectorySecondaryStructure(model, coordinates);
      const trajectory = new InterpolatedTrajectory(
        model,
        coordinates,
        Math.max(1, Math.round(params.steps)),
        params.transform,
      );
      return new SO.Molecule.Trajectory(trajectory, {
        label: 'Trajectory',
        description: `${coordinates.frames.length} frames${params.steps > 1 ? `, ${params.steps - 1} in between` : ''}`,
      });
    });
  },
});
