/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { type ReactNode, useEffect, useState } from 'react';
import type { BehaviorSubject } from 'rxjs';
import { seriesColor } from '../colors.js';
import type { Simulation } from '../session.js';

export function useSubject<T>(subject: BehaviorSubject<T>): T {
  const [value, setValue] = useState(subject.value);
  useEffect(() => {
    setValue(subject.value);
    const sub = subject.subscribe(setValue);
    return () => sub.unsubscribe();
  }, [subject]);
  return value;
}

export function ns(ps: number, digits = 1) {
  return (ps / 1000).toFixed(digits);
}

export function simColor(sim: Simulation) {
  return seriesColor(sim.slot);
}

export function shortGroup(group: string) {
  const parts = group.split('/');
  return parts[parts.length - 1] || group;
}

export function bytes(n: number) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(2)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(0)} MB`;
  return `${Math.max(1, Math.round(n / 1e3))} kB`;
}

/** x axis of a simulation in ns */
const nsCache = new WeakMap<Float32Array, Float32Array>();
export function timesNs(times: Float32Array) {
  let out = nsCache.get(times);
  if (!out) {
    out = times.map((t) => t / 1000);
    nsCache.set(times, out);
  }
  return out;
}

const sectionMemory = new Map<string, string>();
/** Which section of a tab is open; remembered while the page lives, so changing tabs does not reset it. */
export function useSection<T extends string>(key: string, available: readonly T[]): [T, (v: T) => void] {
  const [value, set] = useState<string>(sectionMemory.get(key) ?? available[0]);
  const current = (available as readonly string[]).includes(value) ? (value as T) : available[0];
  return [
    current,
    (v: T) => {
      sectionMemory.set(key, v);
      set(v);
    },
  ];
}

/** Second-level navigation inside a tab: one section on screen at a time instead of one long page. */
export function Sections<T extends string>({
  label,
  options,
  value,
  onChange,
  children,
}: {
  label: string;
  options: readonly (readonly [T, string])[];
  value: T;
  onChange: (v: T) => void;
  children?: ReactNode;
}) {
  return (
    <div className="mdc-subnav">
      <div className="mdc-seg" role="tablist" aria-label={label}>
        {options.map(([id, name]) => (
          <button
            key={id}
            type="button"
            role="tab"
            aria-selected={id === value}
            className={id === value ? 'mdc-seg-on' : undefined}
            onClick={() => onChange(id)}
          >
            {name}
          </button>
        ))}
      </div>
      {children}
    </div>
  );
}

/** Explanations that are worth having but not worth the room: closed until asked for. */
export function Help({ title = 'How to read this', children }: { title?: string; children: ReactNode }) {
  return (
    <details className="mdc-help">
      <summary>{title}</summary>
      <div>{children}</div>
    </details>
  );
}
