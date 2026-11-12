/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useState } from 'react';
import type { TreeNode } from '../analysis.js';

export interface IcicleProps {
  tree: TreeNode[];
  /** merge height at which the tree was cut into conformers */
  cutHeight?: number;
  /** formats a merge height */
  format: (height: number) => string;
  unitLabel: string;
  color: (conformer: number) => string;
  /** up to two blocks picked for comparison */
  selected: number[];
  /** total sampled frames */
  total: number;
  describe: (node: TreeNode) => string;
  onPick: (node: TreeNode) => void;
}

const W = 520,
  H = 230,
  PAD = { l: 44, r: 6, t: 8, b: 6 };

/**
 * The clustering tree as an icicle plot: every block is a cluster, as wide as its share of the frames, hanging from
 * the height at which it joins its sibling down to the height at which it splits itself. One tall split at the top
 * means separate states; a staircase of slivers means one continuous motion.
 */
export function Icicle({ tree, cutHeight, format, unitLabel, color, selected, total, describe, onPick }: IcicleProps) {
  const [hover, setHover] = useState<number | undefined>();
  if (!tree.length) return null;
  const top = Math.max(tree[0].height * 1.08, 1e-6);
  const sx = (x: number) => PAD.l + x * (W - PAD.l - PAD.r),
    sy = (h: number) => PAD.t + (1 - h / top) * (H - PAD.t - PAD.b);
  const ticks: number[] = [];
  const step = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2].find((s) => top / s <= 6) ?? 5;
  for (let v = 0; v <= top; v += step) ticks.push(Number(v.toFixed(4)));
  const h = hover !== undefined ? tree[hover] : undefined;
  return (
    <div className="mdc-plot" style={{ border: 0, padding: 0 }}>
      <svg
        viewBox={`0 0 ${W} ${H}`}
        width="100%"
        role="img"
        aria-label="Clustering tree of the sampled frames as an icicle plot"
      >
        {ticks.map((t) => (
          <g key={t}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="mdc-grid" />
            <text x={PAD.l - 5} y={sy(t) + 3.5} textAnchor="end" className="mdc-tick">
              {format(t)}
            </text>
          </g>
        ))}
        {tree.map((n) => {
          const y0 = sy(n.parent < 0 ? top : n.top),
            y1 = sy(n.height);
          const picked = selected.indexOf(n.id);
          const w = sx(n.x1) - sx(n.x0);
          return (
            <g key={n.id}>
              <rect
                x={sx(n.x0)}
                y={y0}
                width={Math.max(0.5, w)}
                height={Math.max(0.5, y1 - y0)}
                fill={n.conformer >= 0 ? color(n.conformer) : '#c9c8c2'}
                fillOpacity={hover === n.id ? 1 : 0.85}
                stroke={picked >= 0 ? '#0b0b0b' : '#fcfcfb'}
                strokeWidth={picked >= 0 ? 2 : 1}
                onMouseEnter={() => setHover(n.id)}
                onMouseLeave={() => setHover(undefined)}
                onClick={() => onPick(n)}
                style={{ cursor: 'pointer' }}
              />
              {picked >= 0 && w > 12 && (
                <text
                  x={sx(n.x0) + w / 2}
                  y={y0 + 13}
                  textAnchor="middle"
                  className="mdc-scatter-label"
                  pointerEvents="none"
                >
                  {picked === 0 ? 'A' : 'B'}
                </text>
              )}
            </g>
          );
        })}
        {cutHeight !== undefined && (
          <>
            <line
              x1={PAD.l}
              x2={W - PAD.r}
              y1={sy(cutHeight)}
              y2={sy(cutHeight)}
              className="mdc-guide"
              pointerEvents="none"
            />
            <text
              x={W - PAD.r - 3}
              y={sy(cutHeight) - 3}
              textAnchor="end"
              className="mdc-axis-label"
              pointerEvents="none"
            >
              cut into conformers
            </text>
          </>
        )}
      </svg>
      <div className="mdc-readout">
        {h ? (
          <>
            <span>
              <b>
                {h.members.length} frames ({((h.members.length / Math.max(1, total)) * 100).toFixed(0)}%)
              </b>
            </span>
            <span>
              joins its sibling at {format(h.parent < 0 ? h.height : h.top)} {unitLabel}
              {h.children ? `, splits at ${format(h.height)}` : ''}
            </span>
            <span>{describe(h)}</span>
          </>
        ) : (
          <span className="mdc-muted">
            Width: share of the sampled frames. Height: {unitLabel} at which clusters join. Grey blocks span several
            conformers. Click a block to show its most typical frame; click two to compare them below.
          </span>
        )}
      </div>
    </div>
  );
}
