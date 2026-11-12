#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 *
 * Local companion server of the Mol* MD contacts app. It
 *   - finds simulations (.xtc/.trr with a .tpr or structure) under a root folder,
 *   - runs GROMACS on them (trjconv: protein only, whole, centred, fitted, thinned; mdmat; xpm2ps),
 *   - serves the results and the built app.
 * No dependencies beyond Node. It listens on localhost only and never writes outside its work folder.
 *
 *   node server.mjs --root <folder with simulations> [--work <cache folder>] [--port 8765] [--gmx <gmx binary>]
 *
 * The cache is deleted when the last page of the app is closed and when the server stops (Ctrl+C, closed terminal);
 * --keep-cache keeps it for the next run.
 */

import { spawn, spawnSync } from 'node:child_process';
import crypto from 'node:crypto';
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const argv = process.argv.slice(2);
const arg = (flag, fallback) => {
  const i = argv.indexOf(flag);
  return i >= 0 && argv[i + 1] ? argv[i + 1] : fallback;
};
if (argv.includes('--help') || argv.includes('-h')) {
  console.log(
    'Usage: node server.mjs --root <simulations folder> [--work <cache folder>] [--port 8765] [--gmx <binary>] [--app <built app folder>] [--jobs <n>] [--keep-cache]',
  );
  process.exit(0);
}

const root = path.resolve(arg('--root', process.cwd()));
const work = path.resolve(arg('--work', path.join(process.cwd(), 'md-contacts-work')));
const port = Number(arg('--port', '8765'));
const appDir = path.resolve(arg('--app', path.join(here, '../build')));
const maxJobs = Number(arg('--jobs', String(Math.max(1, Math.min(4, Math.floor(os.cpus().length / 2))))));
const keepCache = argv.includes('--keep-cache');
fs.mkdirSync(work, { recursive: true });

// GROMACS

function findGmx() {
  const candidates = [arg('--gmx'), process.env.GMX, 'gmx', 'gmx_mpi', 'gmx_d'].filter(Boolean);
  // a GROMACS built by scripts/build-gromacs.sh beside the work folder
  for (const base of [
    path.join(work, '..', '.tools', 'gromacs', 'bin', 'gmx'),
    path.join(here, 'gromacs', 'bin', 'gmx'),
  ])
    candidates.push(base);
  for (const c of candidates) {
    try {
      const r = spawnSync(c, ['--version'], { encoding: 'utf8', timeout: 15000 });
      if (r.status === 0) {
        const m = /GROMACS version:\s*(.+)/.exec(r.stdout);
        return { path: c, version: m ? m[1].trim() : 'unknown' };
      }
    } catch {
      // try the next one
    }
  }
  return null;
}
const gmx = findGmx();

// Foldseek (optional): encodes the frames in the 3Di structural alphabet
function findFoldseek() {
  const candidates = [
    arg('--foldseek'),
    process.env.FOLDSEEK,
    'foldseek',
    path.join(work, '..', '.tools', 'foldseek', 'bin', 'foldseek'),
  ].filter(Boolean);
  for (const c of candidates) {
    try {
      const r = spawnSync(c, ['version'], { encoding: 'utf8', timeout: 15000 });
      if (r.status === 0) return { path: c, version: r.stdout.trim().split(/\s+/)[0] || 'unknown' };
    } catch {
      // try the next one
    }
  }
  return null;
}
const foldseek = findFoldseek();

function runTool(job, bin, args, cwd) {
  return new Promise((resolve, reject) => {
    job.log.push(`$ ${path.basename(bin)} ${args.join(' ')}`);
    const child = spawn(bin, args, { cwd });
    job.child = child;
    let tail = '';
    const keep = (d) => {
      tail = (tail + d.toString()).slice(-4000);
    };
    child.stderr.on('data', keep);
    child.stdout.on('data', keep);
    child.on('error', reject);
    child.on('close', (code) => {
      job.child = undefined;
      if (code === 0) resolve(tail);
      else {
        job.log.push(
          ...tail
            .split(/\r?\n/)
            .filter((l) => l.trim())
            .slice(-8),
        );
        reject(new Error(`${path.basename(bin)} ${args[0]} failed (exit ${code})`));
      }
    });
  });
}

// biome-ignore format: a table
const STANDARD_RESIDUE = {
  HID: 'HIS', HIE: 'HIS', HIP: 'HIS', HISD: 'HIS', HISE: 'HIS', HISH: 'HIS', HSD: 'HIS', HSE: 'HIS', HSP: 'HIS',
  CYX: 'CYS', CYM: 'CYS', CYS2: 'CYS', ASH: 'ASP', ASPH: 'ASP', GLH: 'GLU', GLUH: 'GLU', LYN: 'LYS', LYSH: 'LYS',
};

/** The trajectory in 3Di letters: one line per frame, written to `3di.txt` in the job's folder. */
async function encode3Di(job, dir) {
  const frames = path.join(dir, 'frames3di'),
    db = path.join(dir, 'db3di');
  fs.rmSync(frames, { recursive: true, force: true });
  fs.rmSync(db, { recursive: true, force: true });
  fs.mkdirSync(frames);
  fs.mkdirSync(db);
  try {
    // backbone and CB of every frame as its own PDB file (all Foldseek reads)
    await runGmx(
      job,
      ['trjconv', '-s', 'ref.pdb', '-f', 'traj.xtc', '-o', path.join('frames3di', 'f.pdb'), '-sep'],
      'MainChain+Cb\n',
      dir,
    );
    // force-field residue names (HISE, CYX, ...) are unknown to Foldseek, which would drop them and split the chain
    for (const f of fs.readdirSync(frames)) {
      const file = path.join(frames, f);
      const text = fs.readFileSync(file, 'latin1');
      const fixed = text.replace(/^(ATOM  .{11})(.{4})/gm, (m, head, name) => {
        const std = STANDARD_RESIDUE[name.trim()];
        return std ? `${head}${std} ` : m;
      });
      if (fixed !== text) fs.writeFileSync(file, fixed, 'latin1');
    }
    await runTool(
      job,
      foldseek.path,
      ['createdb', 'frames3di', path.join('db3di', 'frames'), '--threads', '2', '-v', '1'],
      dir,
    );
    await runTool(
      job,
      foldseek.path,
      ['lndb', path.join('db3di', 'frames_h'), path.join('db3di', 'frames_ss_h'), '-v', '1'],
      dir,
    );
    await runTool(
      job,
      foldseek.path,
      ['convert2fasta', path.join('db3di', 'frames_ss'), path.join('db3di', '3di.fasta'), '-v', '1'],
      dir,
    );
    // entries are named after the files (f<frame>, with a chain suffix when there are several): chains joined per frame
    const byFrame = new Map();
    let current;
    for (const line of fs.readFileSync(path.join(db, '3di.fasta'), 'utf8').split(/\r?\n/)) {
      if (line.startsWith('>')) {
        const m = /^>f(\d+)/.exec(line);
        current = m ? Number(m[1]) : undefined;
        if (current !== undefined && !byFrame.has(current)) byFrame.set(current, '');
      } else if (current !== undefined && line) byFrame.set(current, byFrame.get(current) + line.trim());
    }
    const order = [...byFrame.keys()].sort((a, b) => a - b);
    if (!order.length) throw new Error('Foldseek returned no 3Di sequences.');
    fs.writeFileSync(path.join(dir, '3di.txt'), order.map((f) => byFrame.get(f)).join('\n') + '\n');
  } finally {
    fs.rmSync(frames, { recursive: true, force: true });
    fs.rmSync(db, { recursive: true, force: true });
  }
}

function runGmx(job, args, stdin, cwd, onStderr) {
  return new Promise((resolve, reject) => {
    job.log.push(`$ gmx ${args.join(' ')}${stdin ? `   (stdin: ${stdin.trim().split('\n').join(', ')})` : ''}`);
    const child = spawn(gmx.path, args, { cwd, env: { ...process.env, GMX_MAXBACKUP: '-1' } });
    job.child = child;
    let tail = '';
    child.stderr.on('data', (d) => {
      const s = d.toString();
      tail = (tail + s).slice(-6000);
      onStderr?.(s);
    });
    child.stdout.on('data', () => {});
    child.on('error', reject);
    child.on('close', (code) => {
      job.child = undefined;
      if (code === 0) resolve(tail);
      else {
        const lines = tail.split(/\r?\n/).filter((l) => l.trim());
        const err = lines.filter((l) => /error|fatal|not found|cannot/i.test(l)).slice(-4);
        job.log.push(...lines.slice(-12));
        reject(new Error(`gmx ${args[0]} failed (exit ${code})${err.length ? `: ${err.join(' ')}` : ''}`));
      }
    });
    if (stdin) child.stdin.write(stdin);
    child.stdin.end();
  });
}

// scanning

const TRAJ_EXT = new Set(['.xtc', '.trr']);
const STRUCT_EXT = new Set(['.gro', '.pdb']);
const SKIP_DIRS = new Set(['node_modules', '.git', path.basename(work)]);

function inside(base, p) {
  const rel = path.relative(base, p);
  return !rel.startsWith('..') && !path.isAbsolute(rel);
}

function resolveInRoot(rel) {
  const p = path.resolve(root, rel || '.');
  if (!inside(root, p)) throw Object.assign(new Error('Path is outside the root folder.'), { status: 400 });
  return p;
}

function scan(startRel, maxDepth = 7) {
  const start = resolveInRoot(startRel);
  const sims = [];
  const analyses = [];
  const walk = (dir, depth) => {
    let entries;
    try {
      entries = fs.readdirSync(dir, { withFileTypes: true });
    } catch {
      return;
    }
    const files = entries.filter((e) => e.isFile() && !e.name.startsWith('#') && !e.name.startsWith('.'));
    const byExt = (set) => files.filter((f) => set.has(path.extname(f.name).toLowerCase())).map((f) => f.name);
    const trajs = byExt(TRAJ_EXT),
      tprs = byExt(new Set(['.tpr'])),
      structs = byExt(STRUCT_EXT);
    for (const f of byExt(new Set(['.xpm'])))
      analyses.push({ path: path.relative(root, path.join(dir, f)), size: fs.statSync(path.join(dir, f)).size });
    const stem = (f) => f.slice(0, -path.extname(f).length);
    const newest = (list) =>
      list.map((f) => [f, fs.statSync(path.join(dir, f)).mtimeMs]).sort((a, b) => b[1] - a[1])[0]?.[0];
    for (const t of trajs) {
      const s = stem(t);
      const tpr = tprs.find((f) => stem(f) === s) ?? (tprs.length === 1 ? tprs[0] : undefined);
      const structure = structs.find((f) => stem(f) === s) ?? (structs.length === 1 ? structs[0] : newest(structs));
      const dirName = path.basename(dir);
      const isReplicaDir = /^(rep|replica|run|r|seed|traj|copy)[_-]?\d+$/i.test(dirName);
      const groupDir = isReplicaDir ? path.dirname(dir) : dir;
      const equilibration = /^(nvt|npt|em|min|eq|equil|heat|relax|step[0-6]|anneal)/i.test(s);
      const full = path.join(dir, t);
      const st = fs.statSync(full);
      sims.push({
        id: path.relative(root, full),
        group: path.relative(root, groupDir) || path.basename(root),
        replica: isReplicaDir ? dirName : s,
        stage: s,
        trajectory: path.relative(root, full),
        size: st.size,
        mtime: st.mtimeMs,
        tpr: tpr ? path.relative(root, path.join(dir, tpr)) : undefined,
        structure: structure ? path.relative(root, path.join(dir, structure)) : undefined,
        equilibration,
      });
    }
    if (depth >= maxDepth) return;
    for (const e of entries) {
      if (!e.isDirectory() || e.name.startsWith('.') || SKIP_DIRS.has(e.name)) continue;
      if (path.resolve(dir, e.name) === work) continue;
      walk(path.join(dir, e.name), depth + 1);
    }
  };
  walk(start, 0);
  sims.sort((a, b) => a.id.localeCompare(b.id, undefined, { numeric: true }));
  // the stage most groups treat as production: prefer it in the UI
  return { root, dir: path.relative(root, start), simulations: sims, analyses };
}

// processing

const jobs = new Map();
const queue = [];
let running = 0;

const DEFAULTS = {
  group: 'Protein',
  fitGroup: 'Backbone',
  dt: 100,
  begin: 0,
  end: -1,
  pbc: 'mol',
  center: true,
  fit: true,
  mdmat: true,
  threeDi: true,
  truncate: 4,
  nlevels: 80,
  index: '',
};

function jobKey(sim, p) {
  const st = fs.statSync(sim.trajectoryAbs);
  const h = crypto.createHash('sha1');
  h.update(
    JSON.stringify([
      sim.trajectoryAbs,
      st.size,
      st.mtimeMs,
      sim.topologyAbs,
      p.group,
      p.fitGroup,
      p.dt,
      p.begin,
      p.end,
      p.pbc,
      p.center,
      p.fit,
      p.index,
    ]),
  );
  return h.digest('hex').slice(0, 16);
}

function publicJob(job) {
  const { child, sim, ...rest } = job;
  return { ...rest, log: job.log.slice(-40) };
}

function startProcess(body) {
  if (!gmx)
    throw Object.assign(new Error('GROMACS (gmx) was not found. Start the server with --gmx <binary>.'), {
      status: 409,
    });
  const p = { ...DEFAULTS, ...(body.params ?? {}) };
  for (const k of ['dt', 'begin', 'end', 'truncate', 'nlevels']) p[k] = Number(p[k]);
  for (const k of ['group', 'fitGroup']) {
    if (!/^[\w .+-]{1,64}$/.test(String(p[k])))
      throw Object.assign(new Error(`Bad group name: ${p[k]}`), { status: 400 });
  }
  if (!['mol', 'whole', 'cluster', 'nojump', 'none'].includes(p.pbc))
    throw Object.assign(new Error('Bad pbc mode.'), { status: 400 });
  const trajectoryAbs = resolveInRoot(body.trajectory);
  const topologyRel = body.tpr || body.structure;
  if (!topologyRel)
    throw Object.assign(new Error('The trajectory needs a .tpr (or a structure file) beside it.'), { status: 400 });
  const topologyAbs = resolveInRoot(topologyRel);
  for (const f of [trajectoryAbs, topologyAbs]) {
    if (!fs.existsSync(f)) throw Object.assign(new Error(`Not found: ${path.relative(root, f)}`), { status: 404 });
  }
  const indexAbs = p.index ? resolveInRoot(p.index) : undefined;
  const sim = { trajectoryAbs, topologyAbs, indexAbs };
  const key = jobKey(sim, p);
  const mdmatName = `mdmat_mean_t${p.truncate}_n${p.nlevels}`;
  const existing = jobs.get(key);
  if (existing && existing.status !== 'error') {
    // the trajectory is cached; a new mdmat truncation only reruns mdmat
    if (existing.status === 'done' && p.mdmat && !existing.files.mdmat?.includes(mdmatName)) {
      existing.params = p;
      existing.status = 'queued';
      existing.stage = 'queued';
      enqueue(existing);
    }
    return existing;
  }
  const dir = path.join(work, key);
  fs.mkdirSync(dir, { recursive: true });
  const job = {
    id: key,
    trajectory: body.trajectory,
    params: p,
    status: 'queued',
    stage: 'queued',
    progress: 0,
    trajectoryReady: false,
    files: {},
    log: [],
    error: undefined,
    sim,
    dir,
  };
  jobs.set(key, job);
  enqueue(job);
  return job;
}

function enqueue(job) {
  queue.push(job);
  pump();
}

function pump() {
  while (running < maxJobs && queue.length) {
    const job = queue.shift();
    running++;
    runJob(job)
      .catch((e) => {
        job.status = 'error';
        job.error = String(e.message ?? e);
        job.log.push(job.error);
      })
      .finally(() => {
        running--;
        pump();
      });
  }
}

const url = (job, name) => `/work/${job.id}/${name}`;

async function runJob(job) {
  const { sim, params: p, dir } = job;
  const has = (f) => fs.existsSync(path.join(dir, f)) && fs.statSync(path.join(dir, f)).size > 0;
  const ndx = sim.indexAbs ? ['-n', sim.indexAbs] : [];
  job.status = 'running';
  const pbcArgs = p.pbc === 'none' ? [] : ['-pbc', p.pbc];
  const centerArgs = p.center ? ['-center'] : [];
  const select = (p.center ? `${p.group}\n` : '') + `${p.group}\n`;
  const range = [...(p.begin > 0 ? ['-b', String(p.begin)] : []), ...(p.end >= 0 ? ['-e', String(p.end)] : [])];

  if (!has('traj.xtc') || !has('ref.pdb')) {
    job.stage = 'reference structure';
    await runGmx(
      job,
      [
        'trjconv',
        '-s',
        sim.topologyAbs,
        '-f',
        sim.trajectoryAbs,
        '-o',
        'ref.pdb',
        ...pbcArgs,
        ...centerArgs,
        '-dump',
        String(Math.max(0, p.begin)),
        ...ndx,
      ],
      select,
      dir,
    );
    job.stage = 'extracting trajectory';
    job.progress = 0;
    const size = fs.statSync(sim.trajectoryAbs).size;
    await runGmx(
      job,
      [
        'trjconv',
        '-s',
        sim.topologyAbs,
        '-f',
        sim.trajectoryAbs,
        '-o',
        'whole.xtc',
        ...pbcArgs,
        ...centerArgs,
        ...(p.dt > 0 ? ['-dt', String(p.dt)] : []),
        ...range,
        ...ndx,
      ],
      select,
      dir,
      () => {
        // trjconv prints frame counters only; estimate from the output file's growth is not possible, so pulse
        job.progress = Math.min(0.95, job.progress + 0.002 * (2e9 / Math.max(size, 1e7)));
      },
    );
    if (p.fit) {
      job.stage = 'fitting';
      try {
        await runGmx(
          job,
          ['trjconv', '-s', 'ref.pdb', '-f', 'whole.xtc', '-o', 'traj.xtc', '-fit', 'rot+trans'],
          `${p.fitGroup}\nSystem\n`,
          dir,
        );
        fs.rmSync(path.join(dir, 'whole.xtc'), { force: true });
      } catch (e) {
        job.log.push(`Fitting skipped: ${e.message}`);
        fs.renameSync(path.join(dir, 'whole.xtc'), path.join(dir, 'traj.xtc'));
      }
    } else fs.renameSync(path.join(dir, 'whole.xtc'), path.join(dir, 'traj.xtc'));
  }
  job.files.structure = url(job, 'ref.pdb');
  job.files.trajectory = url(job, 'traj.xtc');
  job.trajectoryReady = true;

  if (foldseek && p.threeDi) {
    if (!has('3di.txt')) {
      job.stage = '3Di (Foldseek)';
      try {
        await encode3Di(job, dir);
      } catch (e) {
        job.log.push(`3Di skipped: ${e.message}`);
      }
    }
    if (has('3di.txt')) job.files.threeDi = url(job, '3di.txt');
  }
  job.progress = 0;

  if (p.mdmat) {
    const name = `mdmat_mean_t${p.truncate}_n${p.nlevels}`;
    if (!has(`${name}.xpm`)) {
      job.stage = 'gmx mdmat';
      await runGmx(
        job,
        [
          'mdmat',
          '-s',
          'ref.pdb',
          '-f',
          'traj.xtc',
          '-mean',
          `${name}.xpm`,
          '-t',
          String(p.truncate),
          '-nlevels',
          String(p.nlevels),
        ],
        'System\n',
        dir,
        (s) => {
          const m = [...s.matchAll(/(\d+)%/g)].pop();
          if (m) job.progress = Number(m[1]) / 100;
        },
      );
    }
    job.files.mdmat = url(job, `${name}.xpm`);
    if (!has(`${name}.eps`)) {
      job.stage = 'gmx xpm2ps';
      try {
        await runGmx(job, ['xpm2ps', '-f', `${name}.xpm`, '-o', `${name}.eps`, '-rainbow', 'blue'], '', dir);
      } catch (e) {
        job.log.push(`xpm2ps failed: ${e.message}`);
      }
    }
    if (has(`${name}.eps`)) job.files.mdmatEps = url(job, `${name}.eps`);
  }
  job.stage = 'done';
  job.progress = 1;
  job.status = 'done';
}

// http

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
  '.pdb': 'text/plain; charset=utf-8',
  '.gro': 'text/plain; charset=utf-8',
  '.xpm': 'text/plain; charset=utf-8',
  '.xvg': 'text/plain; charset=utf-8',
  '.txt': 'text/plain; charset=utf-8',
  '.eps': 'application/postscript',
  '.map': 'application/json',
};

function sendJson(res, status, data) {
  const body = JSON.stringify(data);
  res.writeHead(status, { 'Content-Type': 'application/json', 'Cache-Control': 'no-store' });
  res.end(body);
}

function sendFile(res, file, download) {
  fs.stat(file, (err, st) => {
    if (err || !st.isFile()) {
      res.writeHead(404, { 'Content-Type': 'text/plain' });
      res.end('Not found');
      return;
    }
    res.writeHead(200, {
      'Content-Type': MIME[path.extname(file).toLowerCase()] ?? 'application/octet-stream',
      'Content-Length': st.size,
      'Cache-Control': 'no-store',
      ...(download ? { 'Content-Disposition': `attachment; filename="${path.basename(file)}"` } : {}),
    });
    fs.createReadStream(file).pipe(res);
  });
}

function readBody(req) {
  return new Promise((resolve, reject) => {
    let data = '';
    req.on('data', (c) => {
      data += c;
      if (data.length > 1e6) reject(Object.assign(new Error('Body too large'), { status: 413 }));
    });
    req.on('end', () => {
      try {
        resolve(data ? JSON.parse(data) : {});
      } catch {
        reject(Object.assign(new Error('Bad JSON'), { status: 400 }));
      }
    });
    req.on('error', reject);
  });
}

const DIRECT_EXT = new Set(['.xtc', '.trr', '.gro', '.pdb', '.xpm', '.xvg']);

const server = http.createServer(async (req, res) => {
  try {
    const u = new URL(req.url, 'http://localhost');
    const p = decodeURIComponent(u.pathname);
    if (p === '/api/alive') return pageOpened(req, res);
    if (p === '/api/info') {
      return sendJson(res, 200, { root, work, gmx, foldseek, defaults: DEFAULTS, maxJobs });
    }
    if (p === '/api/scan') {
      return sendJson(res, 200, scan(u.searchParams.get('dir') ?? ''));
    }
    if (p === '/api/process' && req.method === 'POST') {
      return sendJson(res, 200, publicJob(startProcess(await readBody(req))));
    }
    if (p.startsWith('/api/jobs/')) {
      const job = jobs.get(p.slice('/api/jobs/'.length));
      if (!job) return sendJson(res, 404, { error: 'Unknown job' });
      if (req.method === 'DELETE') {
        job.child?.kill();
        const qi = queue.indexOf(job);
        if (qi >= 0) queue.splice(qi, 1);
        job.status = 'error';
        job.error = 'Cancelled';
        jobs.delete(job.id);
      }
      return sendJson(res, 200, publicJob(job));
    }
    if (p === '/api/file') {
      // a file of the root folder as it is (already prepared trajectories, structures, XPM/XVG analyses)
      const file = resolveInRoot(u.searchParams.get('path') ?? '');
      if (!DIRECT_EXT.has(path.extname(file).toLowerCase()))
        return sendJson(res, 400, { error: 'File type not served.' });
      return sendFile(res, file);
    }
    if (p.startsWith('/work/')) {
      const file = path.resolve(work, p.slice('/work/'.length));
      if (!inside(work, file)) return sendJson(res, 400, { error: 'Bad path' });
      return sendFile(res, file, file.endsWith('.eps'));
    }
    // the built app
    const rel = p === '/' ? 'index.html' : p.slice(1);
    const file = path.resolve(appDir, rel);
    if (!inside(appDir, file)) return sendJson(res, 400, { error: 'Bad path' });
    return sendFile(res, file);
  } catch (e) {
    sendJson(res, e.status ?? 500, { error: String(e.message ?? e) });
  }
});

// Closing: stop running tools and delete the cache (only the job folders this server writes, never anything else).

const JOB_DIR = /^[0-9a-f]{12,40}$/;
function clearCache() {
  let removed = 0;
  let entries = [];
  try {
    entries = fs.readdirSync(work, { withFileTypes: true });
  } catch {
    return removed;
  }
  for (const e of entries) {
    if (!e.isDirectory() || !JOB_DIR.test(e.name)) continue;
    try {
      fs.rmSync(path.join(work, e.name), { recursive: true, force: true });
      removed++;
    } catch {}
  }
  try {
    fs.rmdirSync(work); // only succeeds when nothing else is left in it
  } catch {}
  return removed;
}

function stopJobs() {
  queue.length = 0;
  for (const job of jobs.values()) {
    try {
      job.child?.kill('SIGKILL');
    } catch {}
  }
  jobs.clear();
}

// Open pages hold /api/alive; when the last one is gone (tab closed, not just reloaded) the cache is deleted too.
const TAB_GRACE_MS = 15000;
let pages = 0;
let tabTimer;
function pageOpened(req, res) {
  res.writeHead(200, {
    'Content-Type': 'text/event-stream',
    'Cache-Control': 'no-store',
    Connection: 'keep-alive',
  });
  res.write(': connected\n\n');
  const beat = setInterval(() => res.write(': alive\n\n'), 20000);
  pages++;
  clearTimeout(tabTimer);
  req.on('close', () => {
    clearInterval(beat);
    pages--;
    if (pages > 0 || keepCache || closing) return;
    clearTimeout(tabTimer);
    tabTimer = setTimeout(() => {
      if (pages > 0) return;
      stopJobs();
      const removed = clearCache();
      console.log(`Page closed; deleted ${removed} cached simulation${removed === 1 ? '' : 's'} from ${work}`);
    }, TAB_GRACE_MS);
  });
}

let closing = false;
function shutdown(signal) {
  if (closing) return;
  closing = true;
  stopJobs();
  if (keepCache) {
    console.log(`\nStopped (${signal}); cache kept in ${work}`);
  } else {
    const removed = clearCache();
    console.log(`\nStopped (${signal}); deleted ${removed} cached simulation${removed === 1 ? '' : 's'} from ${work}`);
  }
  process.exit(0);
}
for (const signal of ['SIGINT', 'SIGTERM', 'SIGHUP']) process.on(signal, () => shutdown(signal));

server.listen(port, '127.0.0.1', () => {
  console.log(`Mol* MD contacts`);
  console.log(`  simulations: ${root}`);
  console.log(`  work folder: ${work}${keepCache ? ' (kept on exit)' : ' (deleted on exit; --keep-cache keeps it)'}`);
  console.log(
    `  GROMACS:     ${gmx ? `${gmx.path} (${gmx.version})` : 'not found: only prepared files can be loaded'}`,
    `  Foldseek:    ${foldseek ? `${foldseek.path} (${foldseek.version})` : 'not found: no 3Di letters (--foldseek <binary>)'}`,
  );
  console.log(`  open         http://localhost:${port}/`);
  if (!fs.existsSync(path.join(appDir, 'index.html'))) console.log(`  (app not built yet: ${appDir})`);
});
