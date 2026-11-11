#!/bin/bash
# Optional: builds GROMACS inside this folder (.tools/gromacs) for the MD contacts app, for machines without one.
# Nothing is installed outside this folder. Needs cmake and a C++17 compiler; takes 10-30 minutes.
#   scripts/build-gromacs.sh [version, default 2026.3]
# The server finds .tools/gromacs/bin/gmx by itself when no other gmx is on the PATH (or pass --gmx).
set -euo pipefail
VERSION="${1:-2026.3}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
PREFIX="$HERE/.tools/gromacs"
BUILD="$HERE/.tools/gromacs-build"
mkdir -p "$BUILD" && cd "$BUILD"
TARBALL="gromacs-$VERSION.tar.gz"
[[ -f $TARBALL ]] || curl -fL -o "$TARBALL" "https://ftp.gromacs.org/gromacs/$TARBALL"
[[ -d gromacs-$VERSION ]] || tar xzf "$TARBALL"
mkdir -p "gromacs-$VERSION/build" && cd "gromacs-$VERSION/build"
# analysis tools only need a plain CPU build: own FFTW, no MPI, no GPU
cmake .. -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_BUILD_TYPE=Release \
  -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=OFF -DGMX_GPU=OFF -DGMX_OPENMP=OFF -DBUILD_TESTING=OFF -DGMX_INSTALL_NBLIB_API=OFF
cmake --build . -j "$(getconf _NPROCESSORS_ONLN)"
cmake --install .
echo "GROMACS $VERSION installed in $PREFIX; run ./run-md-contacts.sh"
