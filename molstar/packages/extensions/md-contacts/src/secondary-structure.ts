/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { ModelSecondaryStructure } from '@molstar/model/formats/structure/property/secondary-structure';
import { Model, type ResidueIndex, Structure, Unit } from '@molstar/model/model/structure';
import { Coordinates } from '@molstar/model/model/structure/coordinates';
import { SecondaryStructure } from '@molstar/model/model/structure/model/properties/secondary-structure';
import type { SecondaryStructureType } from '@molstar/model/model/structure/model/types';
import { computeUnitDSSP, DefaultDSSPComputationProps } from '@molstar/model/props/computed/secondary-structure/dssp';

/**
 * One secondary structure for a whole trajectory. Mol* assigns it (DSSP) to every frame on its own when the structure
 * file has none, so a helix turn whose hydrogen bond breathes switches between helix and coil from frame to frame and
 * the cartoon flickers. Here DSSP runs on a few frames spread over the trajectory, every residue takes the assignment
 * it has most often, and the result is put on the model that all frames share.
 */

const done = new WeakSet<Model>();

// `SecondaryStructureType.Flag` (a const enum: its values, for the element kinds)
const HELIX = 0x2,
  BETA = 0x4,
  BEND = 0x8,
  TURN = 0x10;

function element(type: number): SecondaryStructure.Element {
  const flags = type as SecondaryStructureType;
  if (type & HELIX) return { kind: 'helix', flags } as SecondaryStructure.Helix;
  if (type & BETA) return { kind: 'sheet', flags } as SecondaryStructure.Sheet;
  if (type & (BEND | TURN)) return { kind: 'turn', flags };
  return { kind: 'none' };
}

/** DSSP of one frame, by model residue (0 where there is none). */
async function frameTypes(model: Model, into: Uint32Array) {
  const structure = Structure.ofModel(model);
  for (const group of structure.unitSymmetryGroups) {
    const unit = group.units[0];
    if (!Unit.isAtomic(unit) || Model.isCoarseGrained(model)) continue;
    const ss = await computeUnitDSSP(unit, DefaultDSSPComputationProps);
    let last = -1;
    for (let i = 0; i < unit.elements.length; i++) {
      const rI = unit.residueIndex[unit.elements[i]];
      if (rI === last) continue;
      last = rI;
      const k = ss.getIndex(rI as ResidueIndex);
      if (k >= 0) into[rI] = ss.type[k] as unknown as number;
    }
  }
}

/** Sets the consensus secondary structure of `samples` frames on `model` (once; not when the file brought its own). */
export async function setTrajectorySecondaryStructure(model: Model, coordinates: Coordinates, samples = 9) {
  if (done.has(model) || Model.hasSecondaryStructure(model)) return;
  done.add(model);
  const n = coordinates.frames.length;
  const nRes = model.atomicHierarchy.residues._rowCount;
  if (!n || !nRes) return;
  const count = Math.min(samples, n);
  const types: Uint32Array[] = [];
  for (let s = 0; s < count; s++) {
    const f = count > 1 ? Math.round((s * (n - 1)) / (count - 1)) : 0;
    const single = Coordinates.create([coordinates.frames[f]], coordinates.deltaTime, coordinates.timeOffset);
    const frame = Model.trajectoryFromModelAndCoordinates(model, single).representative;
    const t = new Uint32Array(nRes);
    await frameTypes(frame, t);
    types.push(t);
  }

  // the assignment each residue has most often (ties: the earlier frame's)
  const type = new Uint32Array(nRes);
  const votes = new Map<number, number>();
  for (let r = 0; r < nRes; r++) {
    votes.clear();
    let best = 0,
      bestCount = 0;
    for (const t of types) {
      const c = (votes.get(t[r]) ?? 0) + 1;
      votes.set(t[r], c);
      if (c > bestCount) {
        bestCount = c;
        best = t[r];
      }
    }
    type[r] = best;
  }

  // elements: runs of one assignment inside a chain
  const chainIndex = model.atomicHierarchy.chainAtomSegments.index;
  const offsets = model.atomicHierarchy.residueAtomSegments.offsets;
  const key = new Int32Array(nRes);
  const elements: SecondaryStructure.Element[] = [];
  for (let r = 0; r < nRes; r++) {
    const newRun = r === 0 || type[r] !== type[r - 1] || chainIndex[offsets[r]] !== chainIndex[offsets[r - 1]];
    if (newRun) elements.push(element(type[r]));
    key[r] = elements.length - 1;
  }
  ModelSecondaryStructure.Provider.set(
    model,
    SecondaryStructure(type as unknown as SecondaryStructureType[], key, elements, (rI) => rI),
  );
}
