#!/usr/bin/env node
// Records GIFs of the StarMD app, mostly of the Conformers tab: drives the app in headless Chrome (Playwright from
// molstar/node_modules, the Chrome installed on this computer), draws a cursor and captions over it, captures the
// screen with the DevTools screencast and turns every scene into a GIF with ffmpeg.
//
//   source scripts/env.sh && node scripts/record-gifs.mjs                 # every scene, simulations under ./MD
//   node scripts/record-gifs.mjs --only find,tree --out docs/gifs        # some scenes
//   node scripts/record-gifs.mjs --url http://localhost:8765/            # a server that is already running
//
// Without --url it starts the app server on --port (8790) with --keep-cache, so md-contacts-work/ is reused by the
// next recording (delete it by hand, or let a normal ./run-md-contacts.sh session clean it up). Frames are written
// to .tools/gif-frames/ and deleted after encoding unless --keep-frames. Needs Google Chrome and ffmpeg.

import { spawn, spawnSync } from 'node:child_process';
import fs from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const repo = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const { chromium } = createRequire(path.join(repo, 'molstar/package.json'))('playwright');

// ---------------------------------------------------------------------------------------------------------------
// options

const argv = process.argv.slice(2);
const opt = (name, fallback) => {
  const i = argv.indexOf(`--${name}`);
  return i >= 0 && argv[i + 1] !== undefined ? argv[i + 1] : fallback;
};
const flag = (name) => argv.includes(`--${name}`);
const SCENES = ['load', 'find', 'timeline', 'superposed', 'tree', 'frame-map', 'switch', 'scale'];
const O = {
  url: opt('url'),
  root: path.resolve(opt('root', path.join(repo, 'MD'))),
  port: Number(opt('port', '8790')),
  out: path.resolve(opt('out', path.join(repo, 'docs/gifs'))),
  only: opt('only', SCENES.join(',')).split(','),
  /** groups unticked in the Data tab before loading (regular expression on the group name) */
  exclude: opt('exclude', 'negative_control'),
  width: Number(opt('width', '1440')),
  height: Number(opt('height', '900')),
  panel: opt('panel', '50%'),
  gifWidth: Number(opt('gif-width', '1200')),
  fps: Number(opt('fps', '12')),
  combined: !flag('no-combined'),
  keepFrames: flag('keep-frames'),
  headed: flag('headed'),
};
for (const s of O.only) if (!SCENES.includes(s)) throw new Error(`Unknown scene ${s}; scenes: ${SCENES.join(', ')}`);
if (spawnSync('ffmpeg', ['-version']).status !== 0) throw new Error('ffmpeg is needed on the PATH.');

const framesRoot = path.join(repo, '.tools/gif-frames');
fs.rmSync(framesRoot, { recursive: true, force: true });
fs.mkdirSync(framesRoot, { recursive: true });
fs.mkdirSync(O.out, { recursive: true });

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const log = (...a) => console.log(`[${new Date().toISOString().slice(11, 19)}]`, ...a);

// ---------------------------------------------------------------------------------------------------------------
// server

let server;
async function startServer() {
  if (O.url) return O.url;
  const script = path.join(repo, 'molstar/apps/md-contacts/server/server.mjs');
  const work = path.join(repo, 'md-contacts-work');
  server = spawn(process.execPath, [script, '--root', O.root, '--work', work, '--port', String(O.port), '--keep-cache'], {
    stdio: ['ignore', 'pipe', 'inherit'],
  });
  server.stdout.on('data', () => {});
  const url = `http://localhost:${O.port}/`;
  for (let i = 0; i < 100; i++) {
    try {
      if ((await fetch(`${url}api/info`)).ok) return url;
    } catch {}
    await sleep(200);
  }
  throw new Error(`The app server did not start on port ${O.port}.`);
}

// ---------------------------------------------------------------------------------------------------------------
// what is drawn over the page: cursor, click ripple, captions, highlight boxes, a stand-in for native select menus

function overlay() {
  const install = () => {
    if (document.getElementById('rec-cursor')) return;
    const style = document.createElement('style');
    style.textContent = `
      #rec-cursor{position:fixed;left:0;top:0;width:28px;height:28px;z-index:2147483647;pointer-events:none;
        margin:-3px 0 0 -4px;filter:drop-shadow(0 1px 2px rgba(0,0,0,.4));will-change:transform}
      .rec-ripple{position:fixed;width:14px;height:14px;margin:-7px 0 0 -7px;border-radius:50%;
        border:3px solid #f39c12;z-index:2147483646;pointer-events:none;animation:rec-ripple .55s ease-out forwards}
      @keyframes rec-ripple{from{transform:scale(.5);opacity:1}to{transform:scale(3.6);opacity:0}}
      #rec-caption{position:fixed;left:16px;top:16px;z-index:2147483640;max-width:var(--rec-cap-w,620px);
        background:rgba(18,21,27,.94);color:#fff;border-radius:12px;padding:12px 16px 14px;pointer-events:none;
        font:15px/1.42 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
        box-shadow:0 8px 28px rgba(0,0,0,.28);opacity:0;transform:translateY(-10px);transition:opacity .35s,transform .35s}
      #rec-caption.on{opacity:1;transform:none}
      #rec-caption .k{font-size:11.5px;font-weight:650;letter-spacing:.09em;text-transform:uppercase;color:#f5b041}
      #rec-caption .h{display:block;font-size:19px;font-weight:650;margin:1px 0 3px}
      #rec-caption .t{color:#d9dce2}
      .rec-spot{position:fixed;border:3px solid #f39c12;border-radius:9px;pointer-events:none;z-index:2147483639;
        box-shadow:0 0 0 5px rgba(243,156,18,.22);opacity:0;transition:opacity .3s}
      .rec-spot.on{opacity:1}
      .rec-menu{position:fixed;z-index:2147483641;background:#fff;border:1px solid #c9c8c2;border-radius:6px;
        box-shadow:0 8px 24px rgba(0,0,0,.18);padding:4px 0;pointer-events:none;
        font:13px/1.2 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;color:#1b1b1b}
      .rec-menu div{padding:5px 12px;white-space:nowrap}
      .rec-menu div.on{background:#2a78d6;color:#fff}`;
    document.head.appendChild(style);
    const cursor = document.createElement('div');
    cursor.id = 'rec-cursor';
    cursor.innerHTML =
      '<svg viewBox="0 0 28 28" width="28" height="28"><path d="M4 3 L4 22 L9 17.2 L12.6 25 L16 23.5 L12.5 16 L19.5 16 Z" fill="#111" stroke="#fff" stroke-width="1.6" stroke-linejoin="round"/></svg>';
    cursor.style.transform = 'translate(-100px,-100px)';
    document.body.appendChild(cursor);
    const cap = document.createElement('div');
    cap.id = 'rec-caption';
    document.body.appendChild(cap);
    const move = (e) => {
      cursor.style.transform = `translate(${e.clientX}px,${e.clientY}px)`;
    };
    window.addEventListener('mousemove', move, true);
    window.addEventListener('pointermove', move, true);
    window.addEventListener('mousedown', (e) => window.__rec.ripple(e.clientX, e.clientY), true);
  };
  window.__rec = {
    ripple(x, y) {
      const d = document.createElement('div');
      d.className = 'rec-ripple';
      d.style.left = `${x}px`;
      d.style.top = `${y}px`;
      document.body.appendChild(d);
      setTimeout(() => d.remove(), 700);
    },
    caption(kicker, head, text) {
      const c = document.getElementById('rec-caption');
      if (!head) return void c.classList.remove('on');
      const esc = (s) => String(s).replace(/[&<>]/g, (m) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' })[m]);
      c.innerHTML = `${kicker ? `<span class="k">${esc(kicker)}</span>` : ''}<span class="h">${esc(head)}</span>${text ? `<span class="t">${esc(text)}</span>` : ''}`;
      c.classList.add('on');
    },
    spot(el, ms, pad = 6) {
      const d = document.createElement('div');
      d.className = 'rec-spot';
      document.body.appendChild(d);
      let live = true;
      const place = () => {
        if (!live) return;
        const r = el.getBoundingClientRect();
        Object.assign(d.style, {
          left: `${r.x - pad}px`,
          top: `${r.y - pad}px`,
          width: `${r.width + 2 * pad}px`,
          height: `${r.height + 2 * pad}px`,
        });
        requestAnimationFrame(place);
      };
      place();
      requestAnimationFrame(() => d.classList.add('on'));
      setTimeout(() => {
        d.classList.remove('on');
        setTimeout(() => {
          live = false;
          d.remove();
        }, 350);
      }, ms);
    },
    menu(r, items, pick, ms) {
      const m = document.createElement('div');
      m.className = 'rec-menu';
      Object.assign(m.style, { left: `${r.x}px`, top: `${r.y + r.height + 2}px`, minWidth: `${r.width}px` });
      m.innerHTML = items.map((t) => `<div>${t}</div>`).join('');
      document.body.appendChild(m);
      setTimeout(() => m.children[pick]?.classList.add('on'), ms * 0.45);
      setTimeout(() => m.remove(), ms);
    },
    clearSpots() {
      for (const d of document.querySelectorAll('.rec-spot,.rec-menu')) d.remove();
    },
  };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', install);
  else install();
}

// ---------------------------------------------------------------------------------------------------------------
// recorder: the page, the screencast, and the moves a scene is made of

class Recorder {
  constructor(page, cdp) {
    this.page = page;
    this.cdp = cdp;
    this.pos = { x: O.width * 0.75, y: O.height * 0.5 };
    this.factor = 1;
    this.scene = undefined;
    this.done = [];
    cdp.on('Page.screencastFrame', ({ data, metadata, sessionId }) => {
      cdp.send('Page.screencastFrameAck', { sessionId }).catch(() => {});
      const s = this.scene;
      if (!s) return;
      const file = path.join(s.dir, `${String(s.frames.length).padStart(5, '0')}.png`);
      fs.writeFileSync(file, Buffer.from(data, 'base64'));
      s.frames.push({ file, t: metadata.timestamp ?? Date.now() / 1000, factor: this.factor });
    });
  }

  async start(name) {
    const dir = path.join(framesRoot, name);
    fs.mkdirSync(dir, { recursive: true });
    this.factor = 1;
    this.scene = { name, dir, frames: [] };
    await this.cdp.send('Page.startScreencast', { format: 'png', everyNthFrame: 1 });
    // a first frame even when nothing moves
    await this.page.mouse.move(this.pos.x + 1, this.pos.y);
    await this.page.mouse.move(this.pos.x, this.pos.y);
    await sleep(300);
  }

  async stop() {
    const end = Date.now() / 1000;
    await this.cdp.send('Page.stopScreencast');
    const s = this.scene;
    this.scene = undefined;
    const frames = s.frames;
    for (let i = 0; i < frames.length; i++) {
      const next = i + 1 < frames.length ? frames[i + 1].t : Math.max(end, frames[i].t + 0.1);
      frames[i].dur = Math.max(0.001, (next - frames[i].t) / frames[i].factor);
    }
    this.done.push(s);
    await this.caption();
    await this.page.evaluate(() => window.__rec.clearSpots());
    return s;
  }

  /** playback speed of what follows (2: twice as fast in the GIF) */
  speed(f) {
    this.factor = f;
  }

  wait(ms) {
    return sleep(ms);
  }

  caption(kicker, head, text) {
    return this.page.evaluate(([k, h, t]) => window.__rec.caption(k, h, t), [kicker, head, text]);
  }

  /** scrolls the element into view (smoothly, so the GIF shows it) and returns its box */
  async box(target) {
    if (typeof target.x === 'number' && typeof target.y === 'number' && !target.boundingBox) return target;
    const scrolled = await target.evaluate((el) => {
      const r = el.getBoundingClientRect();
      if (r.top >= 0 && r.bottom <= innerHeight) return false;
      el.scrollIntoView({ behavior: 'smooth', block: r.height > innerHeight * 0.7 ? 'start' : 'nearest' });
      return true;
    });
    if (scrolled) await sleep(700);
    const b = await target.boundingBox();
    if (!b || !b.width || !b.height) throw new Error(`Not on screen: ${target}`);
    return b;
  }

  async point(target, fx = 0.5, fy = 0.5) {
    const b = await this.box(target);
    if (b.width === undefined) return b;
    return { x: b.x + b.width * fx, y: b.y + b.height * fy };
  }

  async move(target, { ms, fx, fy } = {}) {
    const p = await this.point(target, fx, fy);
    const x0 = this.pos.x,
      y0 = this.pos.y;
    const dist = Math.hypot(p.x - x0, p.y - y0);
    const dur = ms ?? Math.min(900, 250 + dist * 0.9);
    const n = Math.max(6, Math.round(dur / 16));
    for (let i = 1; i <= n; i++) {
      const t = i / n;
      const e = t < 0.5 ? 4 * t * t * t : 1 - (-2 * t + 2) ** 3 / 2;
      await this.page.mouse.move(x0 + (p.x - x0) * e, y0 + (p.y - y0) * e);
      await sleep(dur / n);
    }
    this.pos = p;
  }

  async click(target, { fx, fy, pause = 180, after = 350 } = {}) {
    await this.move(target, { fx, fy });
    await sleep(pause);
    await this.page.mouse.down();
    await sleep(70);
    await this.page.mouse.up();
    await sleep(after);
  }

  /** a native select menu does not show in the screencast: draws a stand-in, then picks the value */
  async choose(select, value, { ms = 1100 } = {}) {
    await this.move(select);
    await sleep(150);
    const b = await select.boundingBox();
    const { items, pick } = await select.evaluate(
      (el, v) => ({
        items: Array.from(el.options, (o) => o.textContent),
        pick: Array.from(el.options).findIndex((o) => o.value === String(v)),
      }),
      value,
    );
    await this.page.evaluate(({ x, y }) => window.__rec.ripple(x, y), this.pos);
    await this.page.evaluate(([r, it, p, t]) => window.__rec.menu(r, it, p, t), [b, items, pick, ms]);
    await sleep(ms * 0.75);
    // menu: 2 px below the select, 1 px border, 4 px padding, rows of 13 px × 1.2 + 10 px
    await this.page.evaluate(({ x, y }) => window.__rec.ripple(x, y), {
      x: b.x + Math.min(b.width - 20, 60),
      y: b.y + b.height + 7 + (pick + 0.5) * 25.6,
    });
    await sleep(ms * 0.25);
    await select.selectOption(String(value));
    await sleep(300);
  }

  /** a box drawn around an element for a while; skipped when the element is not there */
  async spot(target, ms = 2000, pad = 6) {
    if (!(await target.count())) return log('nothing to highlight:', String(target));
    await this.box(target);
    await target.first().evaluate((el, [t, p]) => window.__rec.spot(el, t, p), [ms, pad]);
  }

  until(fn, arg, timeout = 600000) {
    return this.page.waitForFunction(fn, arg, { timeout, polling: 100 });
  }
}

// ---------------------------------------------------------------------------------------------------------------
// encoding

function encode(segments, file, { width = O.gifWidth, fps = O.fps } = {}) {
  const frames = segments.flat();
  const list = path.join(framesRoot, `${path.basename(file, '.gif')}.ffconcat`);
  let txt = 'ffconcat version 1.0\n';
  for (const f of frames) txt += `file '${f.file}'\nduration ${f.dur.toFixed(4)}\n`;
  txt += `file '${frames.at(-1).file}'\n`;
  fs.writeFileSync(list, txt);
  const vf = [
    `fps=${fps}`,
    `scale=${width}:-1:flags=lanczos`,
    'split[a][b]',
    '[a]palettegen=max_colors=256:stats_mode=full[p]',
    '[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle',
  ].join(',');
  const r = spawnSync(
    'ffmpeg',
    ['-v', 'error', '-y', '-f', 'concat', '-safe', '0', '-i', list, '-filter_complex', vf, '-loop', '0', file],
    { stdio: 'inherit' },
  );
  if (r.status !== 0) throw new Error(`ffmpeg failed for ${file}`);
  const seconds = frames.reduce((s, f) => s + f.dur, 0);
  log(`${path.relative(repo, file)}: ${(fs.statSync(file).size / 1e6).toFixed(1)} MB, ${seconds.toFixed(1)} s`);
}

// ---------------------------------------------------------------------------------------------------------------
// the app: helpers that read its state

const KICK = 'StarMD · Conformers';

async function main() {
  const url = await startServer();
  const browser = await chromium.launch({ channel: 'chrome', headless: !O.headed });
  const context = await browser.newContext({ viewport: { width: O.width, height: O.height },
    deviceScaleFactor: 1,
    locale: 'en-US',
  });
  await context.addInitScript(overlay);
  const page = await context.newPage();
  page.on('pageerror', (e) => log('page error:', e.message));
  const cdp = await context.newCDPSession(page);
  const r = new Recorder(page, cdp);

  const tab = (name) => page.locator('.mdc-tabs').getByRole('tab', { name, exact: true });
  // the tab shown beside the Data tab (which stays mounted, hidden)
  const view = page.locator('.mdc-body > .mdc-tab');
  const section = (name) => view.locator('.mdc-subnav').getByRole('tab', { name, exact: true });
  const control = (label) => view.locator('.mdc-controls label', { hasText: label }).locator('select');
  const js = (fn, arg) => page.evaluate(fn, arg);
  const busyDone = () => r.until(() => !window.mdContacts.session.busy.value);
  const sideView = (ms = 0) =>
    js((ms) => {
      // the membrane normal (z of the simulation box) up, looking along y
      const { plugin } = window.mdContacts.viewer;
      const { center, radius } = plugin.canvas3d.boundingSphereVisible;
      const cam = plugin.canvas3d.camera;
      cam.setState(cam.getInvariantFocus([...center], radius * 0.9, [0, 0, 1], [0, 1, 0]), ms);
    }, ms);
  /** rocks the 3D view about the vertical axis (a full spin makes a GIF several times larger); 0 stops it */
  const rock = (speed, angle = 35) =>
    js(
      ([speed, angle]) => {
        const { plugin } = window.mdContacts.viewer;
        plugin.canvas3d.setProps({
          trackball: {
            animate: speed ? { name: 'rock', params: { speed, angle, axis: [0, -1, 0] } } : { name: 'off', params: {} },
          },
        });
      },
      [speed, angle],
    );
  /** what the last conformer run found, in words: which runs hold each conformer */
  const describeRun = () =>
    js(() => {
      const run = window.mdContacts.session.conformers.value;
      const short = (g) => g.split('/').pop().replace(/^state\d+_/, '');
      const conf = run.result.conformers.map((c, ci) => {
        const holders = run.sims
          .map((s, i) => ({ s, share: c.perTrajectory[i] / Math.max(1, s.data.nFrames - s.stability.first) }))
          .filter((x) => x.share >= 0.2);
        const byGroup = new Map();
        for (const { s } of holders) byGroup.set(short(s.group), [...(byGroup.get(short(s.group)) ?? []), s.label]);
        const who = [...byGroup].map(([g, reps]) => `${g} ${reps.join(', ')}`).join(' + ');
        return `C${ci + 1}: ${who || 'transit frames'}`;
      });
      const r = run.result;
      // runs that stay in one conformer (95% of their frames), and the others
      const stay = run.sims.filter((s, i) => {
        const a = r.assignment[i];
        const counts = new Map();
        for (const c of a) counts.set(c, (counts.get(c) ?? 0) + 1);
        return Math.max(...counts.values()) >= 0.95 * a.length;
      }).length;
      return {
        n: r.conformers.length,
        conf,
        stay,
        fine: r.nearby[0],
        coarse: r.nearby[1],
        scale: r.scale,
        nSims: run.sims.length,
        nFrames: r.conformers.reduce((s, c) => s + c.size, 0),
      };
    });

  log('opening', url);
  await page.goto(url);
  await page.waitForSelector('.mdc-panel');
  await js((panel) => {
    document.getElementById('layout').style.setProperty('--panel', panel);
    window.mdContacts.viewer.handleResize();
    // captions sit over the 3D view
    const w = document.getElementById('viewer').getBoundingClientRect().width;
    document.documentElement.style.setProperty('--rec-cap-w', `${Math.round(w - 100)}px`);
  }, O.panel);
  await page.getByRole('button', { name: /Prepare & load \d+ selected/ }).waitFor({ timeout: 60000 });

  const want = new Set(O.only);
  let loaded = false;
  const names = {
    load: 'starmd-load',
    find: 'conformers-1-find',
    timeline: 'conformers-2-timeline',
    superposed: 'conformers-3-superposed',
    tree: 'conformers-4-tree',
    'frame-map': 'conformers-5-frame-map',
    switch: 'conformers-6-switch-pairs',
    scale: 'conformers-7-scale',
  };
  const failed = [];
  /** runs a scene: recorded and encoded when asked for, else only for the state it leaves (a failure skips the GIF) */
  const record = async (name, body) => {
    if (!want.has(name)) return body(false);
    log('recording', name);
    await r.start(name);
    try {
      await body(true);
    } catch (e) {
      if (name === 'load') throw e;
      log(`scene ${name} failed, no GIF:`, String(e.message ?? e).split('\n')[0]);
      failed.push(name);
      await r.stop();
      r.done.pop();
      return;
    }
    const s = await r.stop();
    encode([s.frames], path.join(O.out, `${names[s.name]}.gif`), { fps: name === 'superposed' ? 10 : O.fps });
  };

  // ---- 1. load -------------------------------------------------------------------------------------------------
  await record('load', async () => {
    await r.caption('StarMD · Data', 'Load the simulations', 'The server lists every GROMACS run under the folder, grouped by system; production runs are pre-selected.');
    await r.wait(2200);
    const excluded = page.locator('label.mdc-check input[aria-label^="Select all of "]');
    for (const box of await excluded.all()) {
      const label = await box.getAttribute('aria-label');
      if (new RegExp(O.exclude).test(label) && (await box.isChecked())) await r.click(box);
    }
    const n = Number((await page.getByRole('button', { name: /Prepare & load \d+ selected/ }).textContent()).match(/\d+/)[0]);
    await r.caption('StarMD · Data', `Prepare & load ${n} runs`, 'gmx trjconv, gmx mdmat and Foldseek run once per trajectory; the results are cached, the trajectories stream into the 3D view.');
    await r.click(page.getByRole('button', { name: /Prepare & load \d+ selected/ }), { after: 600 });
    r.speed(2);
    await r.until((n) => {
      const s = window.mdContacts.session.simulations.value;
      return s.length === n && s.every((x) => x.gmxMean && x.threeDi);
    }, n);
    r.speed(1);
    await sideView(700);
    await busyDone();
    await r.wait(1200);
    await r.click(page.locator('.mdc-seg').getByRole('tab', { name: /Loaded/ }));
    await r.caption('StarMD · Data', `${n} replicas loaded`, 'Two states of the SLC35A2 transporter, three 100 ns replicas each. Next: which conformations do they visit?');
    await r.wait(2600);
    await r.click(tab('Conformers'), { after: 900 });
    loaded = true;
  });
  if (!loaded) throw new Error('Nothing loaded.');
  await js(() => window.mdContacts.session.setFrame(0));

  // ---- 2. find conformers ------------------------------------------------------------------------------------------
  await record('find', async (rec) => {
    if (!rec) {
      if (!(await tab('Conformers').getAttribute('aria-selected')).includes('true')) await tab('Conformers').click();
      await control('Look in').selectOption('*');
      await page.getByRole('button', { name: /Find conformers|Run again/ }).click();
      await r.until(() => window.mdContacts.session.conformers.value);
      await busyDone();
      return;
    }
    await sideView();
    await r.caption(KICK, 'Find conformers', 'Each frame is fingerprinted by its Cα distance map and the frames are clustered by dRMSD: no superposition, so rigid-body motion does not count.');
    await r.spot(view.locator(':scope > .mdc-note').first(), 2600);
    await r.wait(2800);
    await r.choose(control('Look in'), '*');
    await r.caption(KICK, 'All six runs together', 'Clustering both states at once shows whether their replicas end up in shared or separate conformations.');
    await r.wait(1400);
    await r.click(page.getByRole('button', { name: 'Find conformers' }));
    r.speed(1.5);
    await r.until(() => window.mdContacts.session.conformers.value);
    await busyDone();
    r.speed(1);
    await r.wait(500);
    const d = await describeRun();
    await r.caption(KICK, `${d.n} conformers in ${d.nSims} runs`, d.conf.join(' · '));
    await r.spot(view.locator('.mdc-card').first().locator('.mdc-card-head'), 2400);
    await r.move(view.locator('.mdc-card').first().locator('.mdc-reason'), { fx: 0.4 });
    await r.wait(2600);
    await r.spot(view.locator('.mdc-card').nth(1), 3200);
    await r.move(view.locator('.mdc-strip').nth(4), { fx: 0.55 });
    const moving = d.nSims - d.stay;
    await r.caption(KICK, 'When each replica is in which conformer', `One row per replica on a shared time axis: ${d.stay} of ${d.nSims} runs stay in one conformer${moving ? `, ${moving === 1 ? 'one switches' : `${moving} switch`} part-way` : ''}.`);
    await r.wait(3800);
  });

  // ---- 3. timeline: strips and playback ------------------------------------------------------------------------------
  await record('timeline', async (rec) => {
    // the replica that changes conformer most often
    const pick = await js(() => {
      const run = window.mdContacts.session.conformers.value;
      let best = 0,
        most = -1;
      run.result.assignment.forEach((a, i) => {
        let k = 0;
        for (let f = 1; f < a.length; f++) if (a[f] !== a[f - 1]) k++;
        if (k > most) (most = k), (best = i);
      });
      const a = run.result.assignment[best];
      // first frame of the conformer it ends in, and one well before
      const last = a[a.length - 1];
      let from = a.length - 1;
      while (from > 0 && a[from - 1] === last) from--;
      return { index: best, n: a.length, early: Math.round(from * 0.4), late: Math.round((from + a.length) / 2), first: a[0], last };
    });
    if (!rec) return;
    await section('Overview').click();
    const strip = view.locator('.mdc-strip').nth(pick.index);
    await r.caption(KICK, 'Click a stripe to see that frame', 'The 3D view jumps to the frame; the chip at the top right names its conformer.');
    await r.click(strip, { fx: (pick.early + 0.5) / pick.n });
    await r.spot(view.locator('.mdc-subnav > span'), 1800, 4);
    await r.wait(1900);
    await r.click(strip, { fx: (pick.late + 0.5) / pick.n });
    await r.spot(view.locator('.mdc-subnav > span'), 1800, 4);
    await r.wait(2000);
    await r.click(strip, { fx: 0.004 });
    await r.wait(500);
    await r.caption(KICK, 'Play through the switch', `The cursor runs along the stripe while the trajectory plays; the chip turns from C${pick.first + 1} to C${pick.last + 1}.`);
    await r.choose(page.getByLabel('Playback step'), 5, { ms: 900 });
    await r.click(page.getByRole('button', { name: 'Play' }), { after: 0 });
    await r.move(strip, { fx: 0.5, fy: 3 });
    r.speed(1.6);
    await r.until(() => window.mdContacts.session.frame.value > 0.9 * window.mdContacts.session.active.value.data.nFrames, null, 120000);
    await r.wait(500);
    r.speed(1);
    await page.getByRole('button', { name: 'Pause' }).click();
    await r.wait(1500);
    await page.getByLabel('Playback step').selectOption('1');
  });

  // ---- 4. superposed medoids -------------------------------------------------------------------------------------------
  await record('superposed', async (rec) => {
    if (!rec) return;
    await section('Overview').click();
    await sideView();
    await r.caption(KICK, 'The conformers in 3D', 'The most typical frame (medoid) of every conformer, each in its colour, superposed on C1 over the aligned Cα atoms.');
    await r.click(page.getByRole('button', { name: /superposed in 3D/ }));
    await r.until(() => window.mdContacts.session.overlay.value);
    await busyDone();
    await sideView(600);
    await r.wait(700);
    await r.spot(view.locator('.mdc-card').first().locator('.mdc-plot-legend'), 2600, 5);
    await r.wait(1200);
    await r.move({ x: O.width * 0.36, y: O.height * 0.94 });
    // one swing each way and back
    await rock(0.2);
    await r.wait(5000);
    await rock(0);
    await r.wait(600);
    await r.click(page.getByRole('button', { name: 'Back to the trajectory' }));
    await r.wait(1000);
  });
  await js(() => window.mdContacts.session.overlay.value && window.mdContacts.session.showConformers(false));

  // ---- 5. the clustering tree --------------------------------------------------------------------------------------------
  await record('tree', async (rec) => {
    if (!rec) return;
    await sideView();
    await r.click(section('Tree'));
    await r.caption(KICK, 'The clustering tree', 'Every block is a cluster: as wide as its share of the frames, hanging from the dRMSD at which it joins its sibling.');
    const rects = page.locator('svg[aria-label^="Clustering tree"] rect');
    const nodes = await js(() => {
      const r = window.mdContacts.session.conformers.value.result;
      // the largest block of each conformer
      return r.conformers.map((c) => {
        let best = -1;
        r.tree.forEach((n, i) => {
          if (n.conformer === c.id && (best < 0 || n.members.length > r.tree[best].members.length)) best = i;
        });
        return best;
      });
    });
    await r.move(rects.nth(nodes[1] ?? nodes[0]), { fy: 0.3 });
    await r.wait(1500);
    await r.move(rects.nth(0), { fx: 0.75, fy: 0.5 });
    await r.wait(1300);
    await r.move(rects.nth(nodes[0]), { fy: 0.3 });
    await r.wait(1200);
    await r.spot(view.locator('svg[aria-label^="Clustering tree"] .mdc-axis-label'), 2200, 4);
    await r.caption(KICK, 'Reading the tree', 'One tall split means separate states; a staircase of slivers is one continuous motion. Grey blocks span several conformers.');
    await r.wait(2600);
    const last = nodes[nodes.length - 1];
    await r.caption(KICK, 'Click two blocks to compare them', 'Each click shows the block’s most typical frame in 3D; two blocks get a difference map and switch pairs of their own.');
    await r.click(rects.nth(nodes[0]), { fy: 0.3 });
    await r.wait(1300);
    await r.click(rects.nth(last), { fy: 0.3 });
    await r.wait(1500);
    await r.click(section('Difference map'));
    await r.caption(KICK, `Block A vs block B`, 'Upper triangle: mean distance A − B (blue: closer in A). Lower: the share of each distance’s variance that lies between the two.');
    await r.move(view.locator('.mdc-mark').first(), { ms: 700 });
    await r.wait(3500);
  });
  // the next scenes compare conformers, not picked blocks
  if (await section('Tree').count()) {
    await section('Tree').click();
    const clear = view.getByRole('button', { name: 'Clear picked blocks' });
    if (await clear.count()) await clear.click();
  }

  // ---- 6. the frame map ------------------------------------------------------------------------------------------------------
  await record('frame-map', async (rec) => {
    if (!rec) return;
    await sideView();
    await r.click(section('Frame map'));
    await r.caption(KICK, 'Map of the sampled frames', 'Each dot is a frame, placed so that close dots have similar distance maps; rings mark the medoids.');
    const dots = page.locator('svg[aria-label^="Map of sampled frames"] circle');
    // a path through every conformer: a few dots near each medoid, the medoid last
    const plan = await js(() => {
      const r = window.mdContacts.session.conformers.value.result;
      return r.conformers.map((c) => {
        const m = r.sample.findIndex((p) => p.ref.traj === c.medoid.traj && p.ref.frame === c.medoid.frame);
        // dots whose frame is also assigned to that conformer (the chip then names the same one)
        const mine = r.sample
          .map((p, i) => ({ p, i }))
          .filter(({ p }) => p.conformer === c.id && r.assignment[p.ref.traj][p.ref.frame] === c.id);
        const mp = r.sample[m];
        mine.sort((a, b) => Math.hypot(b.p.x - mp.x, b.p.y - mp.y) - Math.hypot(a.p.x - mp.x, a.p.y - mp.y));
        return { medoid: m, far: [mine[Math.floor(mine.length * 0.05)].i, mine[Math.floor(mine.length * 0.3)].i] };
      });
    });
    for (const c of plan) {
      for (const i of c.far) {
        await r.move(dots.nth(i), { ms: 520 });
        await r.wait(650);
      }
    }
    const lastC = plan[plan.length - 1];
    await r.caption(KICK, 'Click a dot to show the frame', 'The 3D view and the frame slider jump to it; the chip names its conformer.');
    await r.click(dots.nth(lastC.far[1]));
    await r.spot(view.locator('.mdc-subnav > span'), 1800, 4);
    await r.wait(2000);
    await r.click(dots.nth(plan[0].far[1]));
    await r.spot(view.locator('.mdc-subnav > span'), 1800, 4);
    await r.wait(2200);
  });

  // ---- 7. what tells two conformers apart ---------------------------------------------------------------------------------
  await record('switch', async (rec) => {
    if (!rec) return;
    await sideView();
    await r.click(section('Difference map'));
    const names = await js(() => {
      const sel = document.querySelectorAll('.mdc-card-head select');
      return Array.from(sel, (s) => s.selectedOptions[0]?.textContent);
    });
    await r.caption(KICK, `What tells ${names[0]} and ${names[1]} apart`, `Upper triangle: mean Cα distance in ${names[0]} minus ${names[1]}; lower: the share of each distance’s variance between the two. Rings mark the top switch pairs.`);
    await r.wait(1800);
    const marks = view.locator('.mdc-mark');
    const nMarks = Math.min(3, await marks.count());
    for (let i = 0; i < nMarks; i++) {
      await r.move(marks.nth(i), { ms: 600 });
      await r.wait(1200);
    }
    await r.click(marks.nth(0));
    await r.caption(KICK, 'Click a pair to follow it in 3D', 'The two residues are highlighted and their distance measured in the 3D view, on every frame.');
    await r.wait(2600);
    await r.click(section('Switch residues'));
    await r.caption(KICK, 'Switch pairs, ranked', 'Distances that change most between the two, with the share of their variance that is between them; below, residues whose 3Di letter changes.');
    const rows = view.locator('.mdc-pair table').first().locator('tbody tr');
    for (let i = 0; i < 4; i++) {
      await r.move(rows.nth(i), { fx: 0.2, ms: 450 });
      await r.wait(600);
    }
    await r.click(rows.nth(1), { fx: 0.2 });
    await r.wait(2200);
    await r.spot(view.locator('.mdc-pair .mdc-col').nth(1), 2600);
    await r.move(view.locator('.mdc-pair .mdc-col').nth(1).locator('tbody tr').first(), { fx: 0.2 });
    await r.wait(2800);
  });
  await js(() => window.mdContacts.session.selectPair(undefined));

  // ---- 8. scale of a conformer ---------------------------------------------------------------------------------------------
  await record('scale', async (rec) => {
    if (!rec) return;
    await sideView();
    await r.click(section('Overview'));
    const base = await describeRun();
    await r.caption(KICK, 'How firm is the number?', `Proteins have sub-states at every scale. The default counts groups at least ${base.scale.toFixed(2)} nm apart: ${base.n} conformers.`);
    await r.spot(view.locator('.mdc-card').first().locator('.mdc-card-head'), 2200);
    await r.wait(2600);
    for (const [value, word] of [
      [1.75, 'fine'],
      [3, 'coarse'],
    ]) {
      await r.choose(control('Scale of a conformer'), value);
      await r.click(page.getByRole('button', { name: 'Run again' }));
      r.speed(1.5);
      await r.until(() => !window.mdContacts.session.busy.value);
      r.speed(1);
      await r.wait(300);
      const d = await describeRun();
      const text =
        d.n > base.n
          ? `Groups just ${d.scale.toFixed(2)} nm apart count now, so sub-states inside the runs come out as conformers of their own.`
          : `Only groups at least ${d.scale.toFixed(2)} nm apart count: ${d.conf.join(' · ')}`;
      await r.caption(KICK, `${word[0].toUpperCase()}${word.slice(1)}: ${d.n} conformer${d.n > 1 ? 's' : ''}`, text);
      await r.spot(view.locator('.mdc-card').first().locator('.mdc-card-head'), 1800);
      await r.wait(1200);
      await r.spot(view.locator('.mdc-card').nth(1), 2200);
      await r.wait(2600);
    }
    await r.choose(control('Scale of a conformer'), 2.25);
    await r.click(page.getByRole('button', { name: 'Run again' }));
    await r.until(() => !window.mdContacts.session.busy.value);
    await r.caption(KICK, `Default: ${base.n} conformers`, 'Each run says what a finer and a coarser scale would give: when these differ, the number is a matter of scale, not a property of the simulation.');
    await r.move(view.locator('.mdc-card').first().locator('.mdc-reason'), { fx: 0.5, fy: 0.85 });
    await r.wait(3800);
  });

  await browser.close();

  // one GIF with the whole workflow, smaller and a little faster
  const story = ['load', 'find', 'superposed', 'tree', 'switch'].map((n) => r.done.find((s) => s.name === n));
  if (O.combined && story.every(Boolean)) {
    const fast = story.map((s) => s.frames.map((f) => ({ ...f, dur: f.dur / 1.4 })));
    encode(fast, path.join(O.out, 'conformer-workflow.gif'), { width: 960, fps: 10 });
  }
  if (!O.keepFrames) fs.rmSync(framesRoot, { recursive: true, force: true });
  if (failed.length) {
    log(`failed: ${failed.join(', ')}`);
    process.exitCode = 1;
  }
}

try {
  await main();
} finally {
  server?.kill();
}
