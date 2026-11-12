#!/bin/bash
# Starts the Mol* MD contacts app: http://localhost:8765/
#   ./run-md-contacts.sh [simulations folder] [extra server flags, e.g. --port 9000 --gmx /path/to/gmx]
# Everything it writes stays in this folder (md-contacts-work/ is the GROMACS output cache, deleted when the app tab is closed or the
# server stops unless --keep-cache is given).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/MD"
if [[ $# -gt 0 && -d "$1" ]]; then ROOT="$1"; shift; fi
source "$HERE/scripts/env.sh"
APP="$HERE/molstar/apps/md-contacts"
if [[ ! -f "$APP/build/index.html" || -n "${REBUILD:-}" ]]; then
  (cd "$HERE/molstar" && node scripts/esbuild/app.mjs apps/md-contacts --prd)
fi
exec node "$APP/server/server.mjs" --root "$ROOT" --work "$HERE/md-contacts-work" "$@"
