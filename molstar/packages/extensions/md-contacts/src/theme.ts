/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import { ColorThemeCategory } from '@molstar/graphics/theme/color/categories';
import type { ThemeDataContext } from '@molstar/graphics/theme/theme';
import type { Location } from '@molstar/model/model/location';
import { Bond, type ElementIndex, StructureElement, Unit } from '@molstar/model/model/structure';
import { Scales } from './colors.js';

/**
 * Colours residues by a per-residue value of a loaded simulation (its RMSF). The values are not model data: the
 * session puts them here by simulation id, and the theme's `revision` parameter changes when they do.
 */
export const ResidueValues = new Map<
  string,
  { /** by model residue index; NaN where there is none */ values: Float32Array }
>();

export const MdContactsRmsfThemeName = 'md-contacts-rmsf';
const Missing = Color(0xcfcec9);

export const MdContactsRmsfThemeParams = {
  sim: PD.Text('', { isHidden: true }),
  /** nm: the value drawn darkest */
  max: PD.Numeric(0.4, { min: 0.05, max: 5, step: 0.05 }),
  revision: PD.Numeric(0, {}, { isHidden: true }),
};
export type MdContactsRmsfThemeParams = typeof MdContactsRmsfThemeParams;

export function MdContactsRmsfTheme(
  ctx: ThemeDataContext,
  props: PD.Values<MdContactsRmsfThemeParams>,
): ColorTheme<MdContactsRmsfThemeParams> {
  const values = ResidueValues.get(props.sim)?.values;
  const lut = Scales.spread;
  const of = (unit: Unit, element: ElementIndex) => {
    if (!values || !Unit.isAtomic(unit)) return Missing;
    const v = values[unit.model.atomicHierarchy.residueAtomSegments.index[element]];
    if (!(v >= 0)) return Missing;
    const k = Math.round(Math.min(1, v / props.max) * 255) * 3;
    return Color.fromRgb(lut[k], lut[k + 1], lut[k + 2]);
  };
  function color(location: Location): Color {
    if (StructureElement.Location.is(location)) return of(location.unit, location.element);
    if (Bond.isLocation(location)) return of(location.aUnit, location.aUnit.elements[location.aIndex]);
    return Missing;
  }
  return {
    factory: MdContactsRmsfTheme,
    granularity: 'group',
    preferSmoothing: true,
    color,
    props,
    description: 'Per-residue Cα RMSF of the simulation (nm): pale is rigid, dark orange is flexible.',
  };
}

export const MdContactsRmsfThemeProvider: ColorTheme.Provider<
  MdContactsRmsfThemeParams,
  typeof MdContactsRmsfThemeName
> = {
  name: MdContactsRmsfThemeName,
  label: 'MD flexibility (RMSF)',
  category: ColorThemeCategory.Residue,
  factory: MdContactsRmsfTheme,
  getParams: () => MdContactsRmsfThemeParams,
  defaultValues: PD.getDefaultValues(MdContactsRmsfThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.structure,
};
