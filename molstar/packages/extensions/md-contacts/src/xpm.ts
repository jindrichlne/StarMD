/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

/** One matrix of a GROMACS XPM file (`gmx mdmat`, `gmx xpm2ps`), rows flipped so row 0 is the first y value. */
export interface XpmMatrix {
  title: string;
  legend: string;
  xLabel: string;
  yLabel: string;
  /** 'Continuous' or 'Discrete' */
  type: string;
  width: number;
  height: number;
  /** row-major, `values[y * width + x]`; NaN where the colour has no numeric level */
  values: Float32Array;
  xAxis: number[];
  yAxis: number[];
  min: number;
  max: number;
}

function unquote(s: string) {
  const a = s.indexOf('"'),
    b = s.lastIndexOf('"');
  return a >= 0 && b > a ? s.slice(a + 1, b) : s.trim();
}

/** Parses every matrix of a GROMACS XPM file (a `-frames` file of `gmx mdmat` holds one per frame). */
export function parseXpm(text: string): XpmMatrix[] {
  const lines = text.split(/\r?\n/);
  const out: XpmMatrix[] = [];
  let i = 0;
  while (i < lines.length) {
    const meta: Record<string, string> = {};
    const xAxis: number[] = [],
      yAxis: number[] = [];
    // header comments up to the "w h ncolors cpp" line
    let dims: number[] | undefined;
    for (; i < lines.length; i++) {
      const l = lines[i].trim();
      const m = /^\/\*\s*([\w-]+):\s*(.*?)\s*\*\/$/.exec(l);
      if (m) meta[m[1].toLowerCase()] = unquote(m[2]);
      else if (l.startsWith('"')) {
        const d = unquote(l).trim().split(/\s+/).map(Number);
        if (d.length >= 4 && d.every((x) => Number.isFinite(x))) {
          dims = d;
          i++;
          break;
        }
      }
    }
    if (!dims) break;
    const [width, height, nColors, cpp] = dims;
    const levels = new Map<string, number>();
    for (let c = 0; c < nColors && i < lines.length; i++) {
      const l = lines[i];
      if (!l.trim().startsWith('"')) continue;
      const q = l.indexOf('"');
      const code = l.slice(q + 1, q + 1 + cpp);
      const cm = /\/\*\s*"(.*?)"\s*\*\//.exec(l);
      const v = cm ? Number(cm[1]) : Number.NaN;
      levels.set(code, Number.isFinite(v) ? v : c);
      c++;
    }
    const values = new Float32Array(width * height);
    let row = 0,
      min = Infinity,
      max = -Infinity;
    for (; i < lines.length && row < height; i++) {
      const l = lines[i].trim();
      if (l.startsWith('/*')) {
        const am = /^\/\*\s*([xy])-axis:\s*(.*?)\s*\*\/$/.exec(l);
        if (am) {
          const nums = am[2].split(/\s+/).map(Number);
          (am[1] === 'x' ? xAxis : yAxis).push(...nums);
        }
        continue;
      }
      if (!l.startsWith('"')) continue;
      const q = lines[i].indexOf('"');
      const data = lines[i].slice(q + 1);
      // XPM rows go top to bottom, the y axis of GROMACS bottom to top
      const y = height - 1 - row;
      for (let x = 0; x < width; x++) {
        const v = levels.get(data.substr(x * cpp, cpp)) ?? Number.NaN;
        values[y * width + x] = v;
        if (v < min) min = v;
        if (v > max) max = v;
      }
      row++;
    }
    if (row < height) break;
    out.push({
      title: meta.title ?? '',
      legend: meta.legend ?? '',
      xLabel: meta['x-label'] ?? '',
      yLabel: meta['y-label'] ?? '',
      type: meta.type ?? '',
      width,
      height,
      values,
      xAxis,
      yAxis,
      min,
      max,
    });
  }
  return out;
}

export interface XvgData {
  title: string;
  xLabel: string;
  yLabel: string;
  legends: string[];
  /** columns[0] is x */
  columns: Float64Array[];
}

/** Parses a GROMACS XVG file (`gmx rms`, `gmx gyrate`, `gmx distance`, ...). */
export function parseXvg(text: string): XvgData {
  const legends: string[] = [];
  let title = '',
    xLabel = '',
    yLabel = '';
  const rows: number[][] = [];
  for (const raw of text.split(/\r?\n/)) {
    const l = raw.trim();
    if (!l || l.startsWith('#')) continue;
    if (l.startsWith('@')) {
      let m: RegExpExecArray | null;
      if ((m = /^@\s*title\s+"(.*)"/.exec(l))) title = m[1];
      else if ((m = /^@\s*xaxis\s+label\s+"(.*)"/.exec(l))) xLabel = m[1];
      else if ((m = /^@\s*yaxis\s+label\s+"(.*)"/.exec(l))) yLabel = m[1];
      else if ((m = /^@\s*s(\d+)\s+legend\s+"(.*)"/.exec(l))) legends[Number(m[1])] = m[2];
      continue;
    }
    if (l.startsWith('&')) break;
    const nums = l.split(/\s+/).map(Number);
    if (nums.every((x) => Number.isFinite(x))) rows.push(nums);
  }
  const nCols = rows.length ? Math.min(...rows.map((r) => r.length)) : 0;
  const columns: Float64Array[] = [];
  for (let c = 0; c < nCols; c++) {
    const col = new Float64Array(rows.length);
    for (let r = 0; r < rows.length; r++) col[r] = rows[r][c];
    columns.push(col);
  }
  return { title, xLabel, yLabel, legends, columns };
}
