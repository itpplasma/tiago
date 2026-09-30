#!/usr/bin/env bash
# Build STELLOPT's xdiagno (LIBSTELL + DIAGNO) at a pinned commit.
#
# Usage: ./build_xdiagno.sh            # -> _work/bin/xdiagno          (stock)
#        ./build_xdiagno.sh --patched  # additionally _work/bin/xdiagno_patched
#
# --patched applies patches/*.patch, fixes for DIAGNO bugs found while
# benchmarking (each tracked in a Tiago issue with a draft upstream report).
#
# Ubuntu/Debian prerequisites (not installed by this script):
#   gfortran libopenmpi-dev openmpi-bin libscalapack-openmpi-dev
#   libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev libopenblas-dev
set -euo pipefail

STELLOPT_REPO=${STELLOPT_REPO:-https://github.com/PrincetonUniversity/STELLOPT.git}
STELLOPT_COMMIT=${STELLOPT_COMMIT:-2f181f0d4e71f28afe076967d606b2c020c9db33}
STELLOPT_MACHINE=${STELLOPT_MACHINE:-ubuntu}   # selects SHARE/make_<machine>.inc

here=$(cd "$(dirname "$0")" && pwd)
work="$here/_work"
src="$work/STELLOPT"
mkdir -p "$work/bin" "$work/home/bin"

if [ ! -d "$src/.git" ]; then
    git clone --filter=blob:none "$STELLOPT_REPO" "$src"
fi
git -C "$src" fetch --quiet origin "$STELLOPT_COMMIT" || true
git -C "$src" checkout --quiet "$STELLOPT_COMMIT"
git -C "$src" checkout --quiet -- .

build() {
    # STELLOPT's makefiles install into $HOME/bin; redirect HOME to keep it local.
    (cd "$src" && HOME="$work/home" MACHINE="$STELLOPT_MACHINE" STELLOPT_PATH="$src" \
        ./build_all -j"$(nproc)" LIBSTELL DIAGNO)
}

build
cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno"
echo "built $work/bin/xdiagno"

if [ "${1:-}" = "--patched" ]; then
    for patch in "$here"/patches/*.patch; do
        git -C "$src" apply "$patch"
        echo "applied $(basename "$patch")"
    done
    build
    cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno_patched"
    git -C "$src" checkout --quiet -- .
    echo "built $work/bin/xdiagno_patched"
fi
