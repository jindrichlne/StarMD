#!/bin/bash
# Optional: puts the official Foldseek build into this folder (.tools/foldseek) so the MD contacts server can encode
# trajectories in the 3Di structural alphabet. Nothing is installed outside this folder.
#   scripts/get-foldseek.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
case "$(uname -s)-$(uname -m)" in
  Darwin-*) BUILD=osx-universal ;;
  Linux-x86_64) BUILD=linux-avx2 ;;
  Linux-aarch64) BUILD=linux-arm64 ;;
  *) echo "No Foldseek build known for $(uname -s) $(uname -m): see https://github.com/steineggerlab/foldseek"; exit 1 ;;
esac
mkdir -p "$HERE/.tools" && cd "$HERE/.tools"
curl -fL -o foldseek.tar.gz "https://mmseqs.com/foldseek/foldseek-$BUILD.tar.gz"
tar xzf foldseek.tar.gz && rm foldseek.tar.gz
"$HERE/.tools/foldseek/bin/foldseek" version
echo "Foldseek is in $HERE/.tools/foldseek; restart ./run-md-contacts.sh"
