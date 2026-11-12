/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { createRoot } from 'react-dom/client';
import { MdServer } from '@molstar/md-contacts-extension/server-api';
import { MdContactsSession } from '@molstar/md-contacts-extension/session';
import { MdContactsPanel } from '@molstar/md-contacts-extension/ui/panel';
import { Viewer } from '@molstar/viewer/app';
import './index.html';
import '@molstar/plugin-ui/skin/light.scss';
import '@molstar/md-contacts-extension/ui/style.css';

async function init() {
  const viewer = await Viewer.create('viewer', {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLeftPanel: false,
    layoutShowLog: false,
    layoutShowSequence: true,
    layoutShowRemoteState: false,
    viewportShowExpand: false,
    viewportShowTrajectoryControls: false,
    viewportShowAnimation: false,
    volumeStreamingDisabled: true,
  });
  const session = new MdContactsSession(viewer.plugin);
  // served by the companion server: same origin; opened as a file: prepared files only
  const server = location.protocol.startsWith('http') ? new MdServer('') : undefined;
  // tells the server this page is open: it deletes its cache once the last page is closed
  const alive = server ? new EventSource('/api/alive') : undefined;
  window.addEventListener('pagehide', () => alive?.close());
  createRoot(document.getElementById('panel')!).render(<MdContactsPanel session={session} server={server} />);

  // draggable split between the 3D view and the panel
  const layout = document.getElementById('layout')!,
    splitter = document.getElementById('splitter')!;
  splitter.addEventListener('pointerdown', (down) => {
    splitter.setPointerCapture(down.pointerId);
    const move = (ev: PointerEvent) => {
      const pct = Math.max(25, Math.min(70, ((window.innerWidth - ev.clientX) / window.innerWidth) * 100));
      layout.style.setProperty('--panel', `${pct}%`);
      viewer.handleResize();
    };
    const up = () => {
      splitter.removeEventListener('pointermove', move);
      splitter.removeEventListener('pointerup', up);
    };
    splitter.addEventListener('pointermove', move);
    splitter.addEventListener('pointerup', up);
  });

  Object.assign(window, { mdContacts: { viewer, session, server } });
}

void init();
