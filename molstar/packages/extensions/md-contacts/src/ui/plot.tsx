/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useMemo, useRef, useState } from 'react';

export interface Series {
  key: string;
  label: string;
  color: string;
  x: ArrayLike<number>;
  y: ArrayLike<number>;
  /** drawn heavier (the simulation shown in 3D) */
  emphasis?: boolean;
  dashed?: boolean;
}

export interface LinePlotProps {
  title: string;
  series: Series[];
  xLabel: string;
  yUnit: string;
  height?: number;
  /** x position of the frame shown in 3D */
  cursor?: number;
  /** shaded x range (equilibration left out of the statistics) */
  shadeTo?: number;
  yMin?: number;
  yMax?: number;
  /** horizontal reference lines */
  guides?: { y: number; label: string }[];
  xFormat?: (x: number) => string;
  onPick?: (x: number, nearest: { series: Series; index: number } | undefined) => void;
  onHoverX?: (x: number | undefined) => void;
}

const W = 520;
const PAD = { l: 44, r: 10, t: 8, b: 24 };

function niceTicks(min: number, max: number, count: number) {
  const span = max - min || 1;
  const raw = span / count;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= count) ?? raw;
  const out: number[] = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-6; v += step) out.push(Number(v.toPrecision(10)));
  return out;
}

function path(s: Series, sx: (x: number) => number, sy: (y: number) => number, maxPoints: number) {
  const n = Math.min(s.x.length, s.y.length);
  if (!n) return '';
  // min/max decimation per pixel column keeps spikes
  const bucket = Math.max(1, Math.floor(n / maxPoints));
  let d = '';
  for (let i = 0; i < n; i += bucket) {
    let lo = Infinity,
      hi = -Infinity;
    const end = Math.min(n, i + bucket);
    for (let k = i; k < end; k++) {
      const v = s.y[k];
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    const x = sx(s.x[i]);
    d += `${d ? 'L' : 'M'}${x.toFixed(1)},${sy(lo).toFixed(1)}`;
    if (hi !== lo) d += `L${x.toFixed(1)},${sy(hi).toFixed(1)}`;
  }
  return d;
}

function nearestIndex(x: ArrayLike<number>, v: number) {
  let lo = 0,
    hi = x.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (x[m] < v) lo = m;
    else hi = m;
  }
  return Math.abs(x[lo] - v) <= Math.abs(x[hi] - v) ? lo : hi;
}

/** Lines over a shared x axis with a crosshair read-out; a click picks the x position (a frame). */
export function LinePlot(props: LinePlotProps) {
  const { title, series, xLabel, yUnit, cursor, shadeTo, guides, onPick, onHoverX } = props;
  const H = props.height ?? 150;
  const svg = useRef<SVGSVGElement>(null);
  const [hoverX, setHoverX] = useState<number | undefined>();

  const dom = useMemo(() => {
    let x0 = Infinity,
      x1 = -Infinity,
      y0 = Infinity,
      y1 = -Infinity;
    for (const s of series) {
      const n = Math.min(s.x.length, s.y.length);
      for (let i = 0; i < n; i++) {
        const x = s.x[i],
          y = s.y[i];
        if (x < x0) x0 = x;
        if (x > x1) x1 = x;
        if (y < y0) y0 = y;
        if (y > y1) y1 = y;
      }
    }
    for (const g of guides ?? []) {
      if (g.y < y0) y0 = g.y;
      if (g.y > y1) y1 = g.y;
    }
    if (!Number.isFinite(x0)) return undefined;
    if (props.yMin !== undefined) y0 = props.yMin;
    if (props.yMax !== undefined) y1 = props.yMax;
    if (y1 - y0 < 1e-9) y1 = y0 + 1;
    const padY = props.yMax === undefined ? (y1 - y0) * 0.06 : 0;
    return { x0, x1: x1 > x0 ? x1 : x0 + 1, y0: props.yMin === undefined ? y0 - padY : y0, y1: y1 + padY };
  }, [series, guides, props.yMin, props.yMax]);

  const paths = useMemo(() => {
    if (!dom) return [];
    const sx = (x: number) => PAD.l + ((x - dom.x0) / (dom.x1 - dom.x0)) * (W - PAD.l - PAD.r);
    const sy = (y: number) => H - PAD.b - ((y - dom.y0) / (dom.y1 - dom.y0)) * (H - PAD.t - PAD.b);
    return series.map((s) => ({ s, d: path(s, sx, sy, W - PAD.l - PAD.r) }));
  }, [series, dom, H]);

  if (!dom) {
    return (
      <div className="mdc-plot">
        <div className="mdc-plot-title">{title}</div>
        <div className="mdc-muted">no data</div>
      </div>
    );
  }
  const sx = (x: number) => PAD.l + ((x - dom.x0) / (dom.x1 - dom.x0)) * (W - PAD.l - PAD.r);
  const sy = (y: number) => H - PAD.b - ((y - dom.y0) / (dom.y1 - dom.y0)) * (H - PAD.t - PAD.b);
  const xf = props.xFormat ?? ((x: number) => String(x));
  const toX = (ev: React.MouseEvent) => {
    const r = svg.current!.getBoundingClientRect();
    const px = ((ev.clientX - r.left) / r.width) * W;
    const t = (px - PAD.l) / (W - PAD.l - PAD.r);
    return dom.x0 + Math.max(0, Math.min(1, t)) * (dom.x1 - dom.x0);
  };
  const readout =
    hoverX === undefined
      ? []
      : series
          .filter((s) => s.x.length)
          .map((s) => {
            const i = nearestIndex(s.x, hoverX);
            return { s, i, y: s.y[i] };
          });
  const yTicks = niceTicks(dom.y0, dom.y1, 4);
  const xTicks = niceTicks(dom.x0, dom.x1, 6);

  return (
    <div className="mdc-plot">
      <div className="mdc-plot-head">
        <span className="mdc-plot-title">
          {title}
          {yUnit && <span className="mdc-plot-unit"> ({yUnit})</span>}
        </span>
        {series.length > 1 && (
          <span className="mdc-plot-legend">
            {series.map((s) => (
              <span key={s.key} className={s.emphasis ? 'mdc-emph' : undefined}>
                <i style={{ background: s.color }} />
                {s.label}
              </span>
            ))}
          </span>
        )}
      </div>
      <svg
        ref={svg}
        viewBox={`0 0 ${W} ${H}`}
        width="100%"
        role="img"
        aria-label={title}
        onMouseMove={(ev) => {
          const x = toX(ev);
          setHoverX(x);
          onHoverX?.(x);
        }}
        onMouseLeave={() => {
          setHoverX(undefined);
          onHoverX?.(undefined);
        }}
        onClick={(ev) => {
          if (!onPick) return;
          const x = toX(ev);
          const target = series.find((s) => s.emphasis) ?? series[0];
          onPick(x, target ? { series: target, index: nearestIndex(target.x, x) } : undefined);
        }}
        style={{ cursor: onPick ? 'pointer' : 'default' }}
      >
        {shadeTo !== undefined && shadeTo > dom.x0 && (
          <rect
            x={PAD.l}
            y={PAD.t}
            width={Math.max(0, sx(Math.min(shadeTo, dom.x1)) - PAD.l)}
            height={H - PAD.t - PAD.b}
            className="mdc-shade"
          />
        )}
        {yTicks.map((t) => (
          <g key={`y${t}`}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="mdc-grid" />
            <text x={PAD.l - 5} y={sy(t) + 3.5} textAnchor="end" className="mdc-tick">
              {t}
            </text>
          </g>
        ))}
        {xTicks.map((t) => (
          <text key={`x${t}`} x={sx(t)} y={H - PAD.b + 13} textAnchor="middle" className="mdc-tick">
            {xf(t)}
          </text>
        ))}
        <text x={W - PAD.r} y={H - 1} textAnchor="end" className="mdc-axis-label">
          {xLabel}
        </text>
        {guides?.map((g) => (
          <g key={g.label}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(g.y)} y2={sy(g.y)} className="mdc-guide" />
            <text x={W - PAD.r - 2} y={sy(g.y) - 3} textAnchor="end" className="mdc-axis-label">
              {g.label}
            </text>
          </g>
        ))}
        {paths.map(({ s, d }) => (
          <path
            key={s.key}
            d={d}
            fill="none"
            stroke={s.color}
            strokeWidth={s.emphasis ? 2 : 1.2}
            strokeOpacity={s.emphasis || series.length === 1 ? 1 : 0.7}
            strokeDasharray={s.dashed ? '5 3' : undefined}
            strokeLinejoin="round"
          />
        ))}
        {cursor !== undefined && cursor >= dom.x0 && cursor <= dom.x1 && (
          <line x1={sx(cursor)} x2={sx(cursor)} y1={PAD.t} y2={H - PAD.b} className="mdc-cursor" />
        )}
        {hoverX !== undefined && (
          <>
            <line x1={sx(hoverX)} x2={sx(hoverX)} y1={PAD.t} y2={H - PAD.b} className="mdc-crosshair" />
            {readout.map(({ s, i, y }) => (
              <circle key={s.key} cx={sx(s.x[i])} cy={sy(y)} r={3.5} fill={s.color} stroke="#fff" strokeWidth={1.5} />
            ))}
          </>
        )}
      </svg>
      <div className="mdc-readout">
        {hoverX !== undefined ? (
          <>
            <span>
              {xf(hoverX)} {xLabel}
            </span>
            {readout.map(({ s, y }) => (
              <span key={s.key}>
                <i style={{ background: s.color }} />
                {series.length > 1 ? `${s.label}: ` : ''}
                {Number(y).toFixed(3)} {yUnit}
              </span>
            ))}
          </>
        ) : (
          <span className="mdc-muted">
            {onPick ? 'Hover to read values, click to show that frame' : 'Hover to read values'}
          </span>
        )}
      </div>
    </div>
  );
}
