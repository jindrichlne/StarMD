/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useEffect, useMemo, useRef, useState } from 'react';
import { cssGradient, type ScaleName, Scales } from '../colors.js';

/** One triangle of a residue-residue map: `values[i * n + j]`, coloured from `min` to `max`. */
export interface MapLayer {
  label: string;
  unit: string;
  values: ArrayLike<number> | undefined;
  scale: ScaleName;
  min: number;
  max: number;
}

export interface HeatmapProps {
  n: number;
  /** above the diagonal (column > row) */
  upper: MapLayer;
  /** below the diagonal */
  lower: MapLayer;
  /** residue label of an index, for the tooltip and the axes */
  label: (i: number) => string;
  /** residue number of an index, for the axis ticks */
  tick: (i: number) => number;
  selected?: { i: number; j: number };
  /** extra markers (switch pairs), drawn as rings */
  marks?: { i: number; j: number }[];
  onHover?: (pair: { i: number; j: number } | undefined) => void;
  onSelect?: (pair: { i: number; j: number }) => void;
}

function paint(ctx: CanvasRenderingContext2D, n: number, upper: MapLayer, lower: MapLayer) {
  const img = ctx.createImageData(n, n);
  const px = img.data;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const o = (i * n + j) * 4;
      const layer = j > i ? upper : lower;
      const v = i === j || !layer.values ? Number.NaN : layer.values[i * n + j];
      if (Number.isNaN(v)) {
        px[o] = px[o + 1] = px[o + 2] = i === j ? 255 : 244;
        px[o + 3] = 255;
        continue;
      }
      const lut = Scales[layer.scale];
      let t = (v - layer.min) / (layer.max - layer.min || 1);
      t = t < 0 ? 0 : t > 1 ? 1 : t;
      const k = Math.round(t * 255) * 3;
      px[o] = lut[k];
      px[o + 1] = lut[k + 1];
      px[o + 2] = lut[k + 2];
      px[o + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
}

function fmt(v: number) {
  return Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);
}

function Legend({ layer, corner }: { layer: MapLayer; corner: string }) {
  return (
    <div className="mdc-legend">
      <div className="mdc-legend-title">
        <span className="mdc-corner">{corner}</span> {layer.label}
      </div>
      {layer.values ? (
        <div className="mdc-legend-bar">
          <span>{fmt(layer.min)}</span>
          <div style={{ background: cssGradient(layer.scale) }} />
          <span>
            {fmt(layer.max)} {layer.unit}
          </span>
        </div>
      ) : (
        <div className="mdc-muted">not available</div>
      )}
    </div>
  );
}

/** A square residue-residue map with a different quantity in each triangle, hover read-out and pair selection. */
export function Heatmap(props: HeatmapProps) {
  const { n, upper, lower, label, tick, selected, marks, onHover, onSelect } = props;
  const canvas = useRef<HTMLCanvasElement>(null);
  const [hover, setHover] = useState<{ i: number; j: number; x: number; y: number } | undefined>();

  useEffect(() => {
    const c = canvas.current;
    if (!c || !n) return;
    if (c.width !== n) c.width = c.height = n;
    const ctx = c.getContext('2d');
    if (ctx) paint(ctx, n, upper, lower);
  }, [n, upper, lower]);

  const ticks = useMemo(() => {
    if (!n) return [];
    const first = tick(0),
      last = tick(n - 1);
    const span = Math.max(1, last - first);
    const step = [10, 20, 25, 50, 100, 200, 250, 500, 1000].find((s) => span / s <= 8) ?? 1000;
    const out: { at: number; text: number }[] = [];
    for (let i = 0; i < n; i++) {
      const t = tick(i);
      if (t % step === 0 && (!out.length || out[out.length - 1].text !== t)) out.push({ at: (i + 0.5) / n, text: t });
    }
    return out;
  }, [n, tick]);

  const locate = (ev: React.MouseEvent) => {
    const r = ev.currentTarget.getBoundingClientRect();
    const j = Math.floor(((ev.clientX - r.left) / r.width) * n),
      i = Math.floor(((ev.clientY - r.top) / r.height) * n);
    if (i < 0 || j < 0 || i >= n || j >= n) return undefined;
    return { i, j, x: ev.clientX - r.left, y: ev.clientY - r.top };
  };

  const pct = (i: number) => `${((i + 0.5) / n) * 100}%`;
  const hv = hover ? (hover.j > hover.i ? upper : lower) : undefined;
  const hValue = hover && hv?.values && hover.i !== hover.j ? hv.values[hover.i * n + hover.j] : undefined;

  return (
    <div className="mdc-heatmap">
      <div className="mdc-heatmap-grid">
        <div className="mdc-axis-y">
          {ticks.map((t) => (
            <span key={t.text} style={{ top: `${t.at * 100}%` }}>
              {t.text}
            </span>
          ))}
        </div>
        <div
          className="mdc-heatmap-plot"
          role="img"
          aria-label={`Residue-residue map: ${upper.label} above the diagonal, ${lower.label} below`}
          onMouseMove={(ev) => {
            const p = locate(ev);
            setHover(p);
            onHover?.(p ? { i: p.i, j: p.j } : undefined);
          }}
          onMouseLeave={() => {
            setHover(undefined);
            onHover?.(undefined);
          }}
          onClick={(ev) => {
            const p = locate(ev);
            if (p && p.i !== p.j) onSelect?.({ i: p.i, j: p.j });
          }}
        >
          <canvas ref={canvas} />
          <span className="mdc-corner-tag mdc-corner-tr">{upper.label}</span>
          <span className="mdc-corner-tag mdc-corner-bl">{lower.label}</span>
          {hover && (
            <>
              <div className="mdc-cross-h" style={{ top: pct(hover.i) }} />
              <div className="mdc-cross-v" style={{ left: pct(hover.j) }} />
            </>
          )}
          {marks?.map((m, k) => (
            <div key={k} className="mdc-mark" style={{ left: pct(Math.max(m.i, m.j)), top: pct(Math.min(m.i, m.j)) }} />
          ))}
          {selected && (
            <>
              <div className="mdc-selected" style={{ left: pct(selected.j), top: pct(selected.i) }} />
              <div className="mdc-selected" style={{ left: pct(selected.i), top: pct(selected.j) }} />
            </>
          )}
          {hover && (
            <div
              className="mdc-tooltip"
              style={{
                left: hover.x > 190 ? undefined : hover.x + 14,
                right: hover.x > 190 ? `calc(100% - ${hover.x - 14}px)` : undefined,
                top: hover.y + 14,
              }}
            >
              <b>{label(hover.i)}</b> – <b>{label(hover.j)}</b>
              {hValue !== undefined && !Number.isNaN(hValue) && (
                <div>
                  {hv!.label}: {fmt(hValue)} {hv!.unit}
                </div>
              )}
            </div>
          )}
        </div>
        <div />
        <div className="mdc-axis-x">
          {ticks.map((t) => (
            <span key={t.text} style={{ left: `${t.at * 100}%` }}>
              {t.text}
            </span>
          ))}
        </div>
      </div>
      <div className="mdc-legends">
        <Legend layer={upper} corner="◥" />
        <Legend layer={lower} corner="◣" />
      </div>
    </div>
  );
}
