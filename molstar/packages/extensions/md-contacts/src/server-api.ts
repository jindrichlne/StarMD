/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

/** Client of the companion server (`apps/md-contacts/server/server.mjs`) that finds simulations and runs GROMACS. */

export interface ServerInfo {
  root: string;
  work: string;
  gmx: { path: string; version: string } | null;
  /** encodes frames in the 3Di structural alphabet; optional */
  foldseek?: { path: string; version: string } | null;
  defaults: ProcessParams;
  maxJobs: number;
}

export interface ProcessParams {
  /** index group written out and analysed */
  group: string;
  fitGroup: string;
  /** ps between the frames kept (0: every frame) */
  dt: number;
  /** ps; -1 for the end of the file */
  begin: number;
  end: number;
  pbc: 'mol' | 'whole' | 'cluster' | 'nojump' | 'none';
  center: boolean;
  fit: boolean;
  mdmat: boolean;
  /** encode the frames in 3Di when the server has Foldseek */
  threeDi: boolean;
  /** `gmx mdmat -t`, nm: distances above it are cut */
  truncate: number;
  nlevels: number;
  /** optional index file, relative to the root */
  index: string;
}

export interface FoundSimulation {
  id: string;
  group: string;
  replica: string;
  stage: string;
  trajectory: string;
  size: number;
  mtime: number;
  tpr?: string;
  structure?: string;
  equilibration: boolean;
}

export interface ScanResult {
  root: string;
  dir: string;
  simulations: FoundSimulation[];
  analyses: { path: string; size: number }[];
}

export interface Job {
  id: string;
  trajectory: string;
  params: ProcessParams;
  status: 'queued' | 'running' | 'done' | 'error';
  stage: string;
  progress: number;
  trajectoryReady: boolean;
  files: { structure?: string; trajectory?: string; mdmat?: string; mdmatEps?: string; threeDi?: string };
  log: string[];
  error?: string;
}

export class MdServer {
  constructor(readonly base = '') {}

  private async json<T>(path: string, init?: RequestInit): Promise<T> {
    const res = await fetch(this.base + path, init);
    const data = await res.json().catch(() => ({ error: `Server answered ${res.status}` }));
    if (!res.ok) throw new Error(data.error ?? `Server answered ${res.status}`);
    return data as T;
  }

  info() {
    return this.json<ServerInfo>('/api/info');
  }

  scan(dir = '') {
    return this.json<ScanResult>(`/api/scan?dir=${encodeURIComponent(dir)}`);
  }

  process(sim: FoundSimulation, params: Partial<ProcessParams>) {
    return this.json<Job>('/api/process', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ trajectory: sim.trajectory, tpr: sim.tpr, structure: sim.structure, params }),
    });
  }

  job(id: string) {
    return this.json<Job>(`/api/jobs/${id}`);
  }

  url(path: string) {
    return this.base + path;
  }

  fileUrl(relPath: string) {
    return `${this.base}/api/file?path=${encodeURIComponent(relPath)}`;
  }

  async text(url: string) {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Could not read ${url} (${res.status})`);
    return res.text();
  }
}
