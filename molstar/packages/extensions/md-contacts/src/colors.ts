/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

/** Colour scales of the maps and plots: one-hue sequential ramps, a diverging pair with a neutral middle, fixed series. */

type RGB = [number, number, number];

function hex(c: string): RGB {
  const v = Number.parseInt(c.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

function ramp(stops: string[]) {
  const rgb = stops.map(hex);
  // 256-entry lookup, t in [0, 1]
  const lut = new Uint8ClampedArray(256 * 3);
  for (let k = 0; k < 256; k++) {
    const t = (k / 255) * (rgb.length - 1);
    const a = Math.min(rgb.length - 2, Math.floor(t)),
      f = t - a;
    for (let c = 0; c < 3; c++) lut[k * 3 + c] = rgb[a][c] + (rgb[a + 1][c] - rgb[a][c]) * f;
  }
  return lut;
}

const BLUE = ['#f1f6fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
const ORANGE = ['#fdf3ee', '#fbd9c9', '#f6b497', '#f08e63', '#eb6834', '#c44f1f', '#963b15', '#6b290d'];

export const Scales = {
  /** distance: close is dark */
  distance: ramp([...BLUE].reverse()),
  /** magnitude from zero: zero recedes */
  spread: ramp(ORANGE),
  /** frequency / share, 0..1 */
  share: ramp(BLUE),
  /** signed difference: blue below zero, red above, grey at zero */
  diverging: ramp(['#104281', '#2a78d6', '#86b6ef', '#f0efec', '#f0a3a2', '#e34948', '#8e1f1f']),
};
export type ScaleName = keyof typeof Scales;

export function cssGradient(scale: ScaleName) {
  const lut = Scales[scale];
  const stops: string[] = [];
  for (let k = 0; k <= 8; k++) {
    const i = Math.round((k / 8) * 255) * 3;
    stops.push(`rgb(${lut[i]},${lut[i + 1]},${lut[i + 2]}) ${(k / 8) * 100}%`);
  }
  return `linear-gradient(to right, ${stops.join(', ')})`;
}

/** Fixed categorical order (never cycled by rank): a simulation keeps its slot for the session. */
export const SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];

export function seriesColor(slot: number) {
  return SERIES[slot % SERIES.length];
}

export const STATUS = {
  good: '#0ca30c',
  warning: '#fab219',
  serious: '#ec835a',
  critical: '#d03b3b',
};
