# StarMD: Mol* MD contacts

A Mol* plugin that takes GROMACS simulations (`.xtc` + `.tpr`), runs `gmx trjconv`, `gmx mdmat` and `gmx xpm2ps` on
them, and shows beside the 3D trajectory: an interactive residue contact map, a stability evaluation of every replica,
and conformers. It follows the "mining protein conformations" pipeline (steps 5 and 6: does a state hold, do the
replicas agree, what tells states apart).

## Build

Requirements: macOS or Linux, **Node.js ≥ 22.13** (Node 24 LTS recommended), `curl`, `tar`, `rsync`, and GROMACS
(any recent `gmx`; see [GROMACS](#gromacs) to build one locally). Foldseek is optional.

```bash
git clone https://github.com/jindrichlne/StarMD.git
cd StarMD
./setup.sh
```

`setup.sh` does, all inside the repository folder:

1. downloads Mol* at the pinned commit (`MOLSTAR_COMMIT` in `setup.sh`) into `molstar/`, around the files this
   repository tracks there (the MD contacts app and extension, plus `tsconfig.json`, `pnpm-lock.yaml` and
   `scripts/workspace/inventory.json`, which register them in the Mol* workspace);
2. installs the pnpm version Mol* pins into `.tools/`;
3. runs `pnpm install --frozen-lockfile` and builds the app into `molstar/apps/md-contacts/build/`.

Tool caches (npm, pnpm store, XDG) go to `.tools/` via `scripts/env.sh`; nothing is installed globally. After pulling
changes, rebuild with `REBUILD=1 ./run-md-contacts.sh`, or rerun `./setup.sh` if dependencies changed.

## Run

```bash
./run-md-contacts.sh            # simulations under ./MD, then open http://localhost:8765/
./run-md-contacts.sh /path/to/simulations --port 9000 --gmx /path/to/gmx
```

Everything is kept in this folder: tool caches in `.tools/`, GROMACS outputs in `md-contacts-work/` (safe to delete; it
is rebuilt on demand). The simulation files are only read.

`MD/` holds the inputs of the SLC35A2 example (two states, membrane systems, see `MD/README.md`); the trajectories
themselves (~12 GB), the minimisation and the replica folders are not in the repository: recreate them with
`MD/scripts/build_state.sh`, run them with `MD/run_all.sh` (or `MD/submit_lumi.sh` on LUMI), or point the server at
your own simulations.

1. **Data**: the server lists every trajectory under the folder, grouped by system (`state1/rep1`, `rep2`, ... are
   replicas of one group). Production runs are pre-selected, equilibration stages hidden behind a checkbox. "Prepare &
   load selected" runs GROMACS and loads the results; 100 ns of a 100 000-atom system takes about 10 s to extract and
   40 s for `gmx mdmat`, which continues in the background while you already look at the trajectory.
2. **Contacts**: the Cα distance map, a different quantity in each triangle (this frame, mean, std, time in contact,
   frame − mean, frame − start, `gmx mdmat` mean, disagreement between replicas, difference to another state). It
   follows the frame slider and playback. Hover highlights the residue pair in 3D; a click measures it in 3D and plots
   its distance over time for all replicas (click the plot to jump to a frame).
   Below the map, the **sequence** is shaded by the contacts of each residue in the frame shown (or by how many it
   gained or lost since the start); click a residue to pin it and see its partners along the sequence (shade: share of
   time in contact, outline: touching now), then click a partner to follow that pair.
3. **Stability**: per replica a verdict (holds / drifting / departed / too short) from the superposition-free dRMSD
   to the start and the share of starting contacts kept, with RMSD, radius of gyration and per-residue RMSF; per group
   whether the replicas agree; with two or more states loaded, which state every frame is closest to.
   A checkbox colours the 3D structure by RMSF (pale: rigid, dark orange: flexible; one scale for all replicas).
4. **Conformers**: frames are fingerprinted by their distance map, clustered by dRMSD (average linkage) and assigned;
   you get the number of conformers the merge heights support (one, when the motion is continuous), when each replica
   visits which, a 2D map of the frames, the medoid frames, and for any two conformers the difference map and the
   ranked switch pairs. Run it on one state or on all loaded simulations together.

   "Show all conformers superposed in 3D" draws the most typical frame of every conformer together, each in its
   conformer's colour, superposed on the first over the aligned Cα atoms (the simulations are hidden meanwhile; moving
   the frame slider goes back).
   The **clustering tree** is drawn as an icicle plot: every block is a cluster, as wide as its share of the frames,
   hanging from the height at which it joins its sibling. One tall split means separate states, a staircase of slivers
   one continuous motion. Click a block for its most typical frame; click two blocks to compare them (difference map,
   switch pairs) instead of two conformers.
5. **Align**: compare the frame shown with a frame of another replica or state. The residues of the two are aligned
   (as they are when the numbering matches, else by a global sequence alignment, so mutants, homologues and renumbered
   copies pair up), the reference frame is superposed on the frame shown over the aligned Cα atoms and drawn in grey
   in 3D, and the aligned sequences are shaded by how far each residue pair is apart after the fit. The numbers
   follow the frame slider; the reference frame has its own slider (and "medoid" buttons after a conformer run).

**3Di structural alphabet** (needs Foldseek, see below): every frame is encoded as one 3Di letter per residue, so a
trajectory becomes an alignment with frames as rows.

- Contacts tab, sequence panel: "3Di letters of this frame" (dark where a residue does not have its usual letter) and
  "3Di variability over time" (column entropy).
- Stability tab: share of residues whose letter changed since the start, column entropy per residue, and a
  residues-by-time map of where the letter is not the usual one.
- Conformers tab: "Compare frames by 3Di letters" as an alternative fingerprint, and "Local switch residues": the
  residues whose letter distribution differs between two conformers or tree blocks.

3Di sees local environment. Thermal motion alone changes about a quarter of the letters between two frames, and a
rigid-body change moves few of them: on the SLC35A2 example the 3Di fingerprint does not separate lumen-open from
cytosol-open (the distance fingerprint does), while the local-switch table still names the residues that differ.

**Smooth playback**: "Smooth" beside the play button draws 1 to 15 linearly interpolated frames between the saved
ones (made on demand, nothing extra is held in memory). Playback gets slower and smoother; every analysis still uses
the saved frames only.

### Without the server, or with analyses you already have

Drop files on "Already prepared files" (Data tab): trajectories (`.xtc`, `.trr`, `.dcd`) with a structure of the same
atoms (`.pdb`, `.gro`; same file name, or one structure for all), and `gmx mdmat` maps (`.xpm`) for the simulation
shown in 3D. The server also offers "Load as is" for small trajectories that have a structure beside them. Prepared
trajectories should hold the protein only, whole and fitted, e.g.

```bash
echo Protein Protein | gmx trjconv -s md.tpr -f md.xtc -o prot.xtc -pbc mol -center -dt 100
echo Protein Protein | gmx trjconv -s md.tpr -f md.xtc -o prot.pdb -pbc mol -center -dump 0
echo System | gmx mdmat -s prot.pdb -f prot.xtc -mean mean.xpm -t 4     # -t: the default 1.5 nm hides large motions
gmx xpm2ps -f mean.xpm -o mean.eps
```

### Foldseek (optional, for 3Di)

The server uses `--foldseek`, `$FOLDSEEK`, `foldseek` on the `PATH`, or `.tools/foldseek/bin/foldseek`;
`scripts/get-foldseek.sh` fetches the official build into `.tools/foldseek`. Without it everything else works and the
3Di views say so. For every prepared trajectory the server writes backbone + Cβ frames with `gmx trjconv -sep`,
renames force-field residue names (HISE, CYX, ...) to the standard ones, runs `foldseek createdb` and keeps the 3Di
strings as `3di.txt` (one line per frame) in the cache; such a file can also be dropped on the Data tab.

### GROMACS

The server uses `--gmx`, `$GMX`, or `gmx` / `gmx_mpi` on the `PATH`. On a machine without GROMACS,
`scripts/build-gromacs.sh` builds one into `.tools/gromacs` (CPU only, enough for the analysis tools).

## Where the code is

| Path | What |
|---|---|
| `molstar/packages/extensions/md-contacts/src/analysis.ts` | distance maps, stability, replica agreement, conformers, switch pairs (no DOM, unit tested) |
| `.../src/xpm.ts` | GROMACS XPM and XVG readers |
| `.../src/threedi.ts` | 3Di trajectories: column entropy, change since the start, local switch residues (unit tested) |
| `.../src/align.ts` | residue alignment of two simulations and superposition of two frames (unit tested) |
| `.../src/interpolate.ts` | the trajectory with in-between frames for smooth playback |
| `.../src/secondary-structure.ts` | one consensus secondary structure per trajectory (stops cartoon flicker) |
| `.../src/theme.ts` | the RMSF colour theme (also listed in Mol*'s own colour menu as "MD flexibility (RMSF)") |
| `.../src/session.ts` | the Mol* side: loading replicas as trajectories, frames, highlights, measurements |
| `.../src/ui/` | the panel (React): Data, Contacts, Stability, Conformers |
| `molstar/apps/md-contacts/` | the page, and `server/server.mjs` (no dependencies): scanning, running `gmx`, caching |

```bash
source scripts/env.sh && cd molstar
npx vitest run packages/extensions/md-contacts      # tests
npx tsc -b apps/md-contacts                         # types
node scripts/esbuild/app.mjs apps/md-contacts --prd # build (REBUILD=1 ./run-md-contacts.sh does it too)
```

## Limits

- The cartoon uses one secondary structure for the whole trajectory (DSSP on nine frames spread over it, each residue
  taking its most frequent assignment), so helices do not flicker; it therefore does not show a helix that really
  unfolds part-way through. Structure files with their own HELIX/SHEET records keep those.

- Distances are between Cα atoms (residues without a Cα, such as ACE/NME caps, are left out); `gmx mdmat` gives the
  smallest heavy-atom distance and is shown as its own layer.
- The stability verdict is a heuristic with fixed thresholds (listed in the Stability tab): read the plots.
- Conformer shares are shares of simulated frames, not populations.
- Trajectories are thinned (default one frame per 100 ps) and held in browser memory: about 60 MB per 1000 frames of
  a 5000-atom protein.
