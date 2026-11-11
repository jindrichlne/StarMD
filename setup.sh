#!/bin/bash
# One-time setup: completes molstar/ with the pinned upstream Mol* commit, installs pnpm and the workspace
# dependencies, and builds the MD contacts app. Everything stays inside this folder (.tools/ holds the caches).
#   ./setup.sh
# Needs Node.js >= 22.13 (Node 24 LTS recommended), curl, tar and rsync. Safe to run again.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Mol* master on 2026-11-11; the files this repository tracks under molstar/ are written against it
MOLSTAR_COMMIT=dabbf29f952880a034a9637e28c9844e3ae73d9f
PNPM_VERSION=12.15.0

node -e 'const [a,b]=process.versions.node.split(".").map(Number); if (a<22||(a===22&&b<13)) { console.error("Node.js >= 22.13 needed, found "+process.version); process.exit(1) }'
source "$HERE/scripts/env.sh"
mkdir -p "$HERE/.tools"

# 1. Mol* sources: upstream files are added around the ones tracked here, never over them
if [[ ! -f "$HERE/molstar/package.json" ]]; then
  echo "Downloading Mol* $MOLSTAR_COMMIT"
  TMP="$HERE/.tools/molstar-download"
  rm -rf "$TMP" && mkdir -p "$TMP"
  curl -fsSL "https://codeload.github.com/molstar/molstar/tar.gz/$MOLSTAR_COMMIT" | tar xz -C "$TMP" --strip-components 1
  rsync -a --ignore-existing "$TMP/" "$HERE/molstar/"
  rm -rf "$TMP"
fi

# 2. pnpm (the version Mol* pins), kept in .tools/node_modules
if [[ "$(pnpm --version 2>/dev/null)" != "$PNPM_VERSION" ]]; then
  npm install --prefix "$HERE/.tools" --no-audit --no-fund "pnpm@$PNPM_VERSION"
fi

# 3. dependencies and the app
cd "$HERE/molstar"
pnpm install --frozen-lockfile
node scripts/esbuild/app.mjs apps/md-contacts --prd
echo "Done. Start with ./run-md-contacts.sh and open http://localhost:8765/"
