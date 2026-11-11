/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useEffect, useState } from 'react';
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
