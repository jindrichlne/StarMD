/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jindrich Lnenicka
 */

import { useState } from 'react';
import type { MdServer } from '../server-api.js';
import type { MdContactsSession } from '../session.js';
import { AlignTab } from './align.js';
import { useSubject } from './common.js';
import { ConformersTab } from './conformers.js';
import { ContactsTab } from './contacts.js';
import { DataTab } from './data.js';
import { StabilityTab } from './stability.js';
import { Timeline } from './timeline.js';

type Tab = 'data' | 'contacts' | 'stability' | 'conformers' | 'align';
const TABS: [Tab, string][] = [
  ['data', 'Data'],
  ['contacts', 'Contacts'],
  ['stability', 'Stability'],
  ['conformers', 'Conformers'],
  ['align', 'Align'],
];

/** The MD contacts panel: load simulations, then read contact maps, stability and conformers beside the 3D view. */
export function MdContactsPanel({ session, server }: { session: MdContactsSession; server?: MdServer }) {
  const [tab, setTab] = useState<Tab>('data');
  const sims = useSubject(session.simulations);
  const busy = useSubject(session.busy);
  return (
    <div className="mdc-panel">
      <div className="mdc-header">
        <div className="mdc-tabs" role="tablist">
          {TABS.map(([id, name]) => (
            <button
              key={id}
              type="button"
              role="tab"
              aria-selected={tab === id}
              className={tab === id ? 'mdc-tab-btn mdc-tab-on' : 'mdc-tab-btn'}
              disabled={id !== 'data' && !sims.length}
              onClick={() => setTab(id)}
            >
              {name}
            </button>
          ))}
          {busy && (
            <span className="mdc-busy">
              <span className="mdc-spinner" /> {busy}
            </span>
          )}
        </div>
        <Timeline session={session} />
      </div>
      <div className="mdc-body">
        {/* the Data tab stays mounted: its jobs keep running while you look at results */}
        <div style={{ display: tab === 'data' ? undefined : 'none' }}>
          <DataTab session={session} server={server} />
        </div>
        {tab === 'contacts' && <ContactsTab session={session} />}
        {tab === 'stability' && <StabilityTab session={session} />}
        {tab === 'conformers' && <ConformersTab session={session} />}
        {tab === 'align' && <AlignTab session={session} />}
      </div>
    </div>
  );
}
