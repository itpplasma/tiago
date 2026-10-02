#!/usr/bin/env bash
# Build STELLOPT's xdiagno (LIBSTELL + DIAGNO), optionally VMEC2000, at a pinned commit.
#
# Usage: ./build_xdiagno.sh [--patched] [--vmec]
#   always      -> _work/bin/xdiagno          (stock)
#   --patched   -> also _work/bin/xdiagno_patched
#   --vmec      -> also _work/bin/xvmec2000   (for the `equilibrium` suite)
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
    # Override STELLOPT's installation directory through make, keeping every
    # benchmark artifact local without changing the process's home directory.
    (cd "$src"
     export MACHINE="$STELLOPT_MACHINE" STELLOPT_PATH="$src"
     export FLAG_CALLED_FROM_BUILD_ALL=true
     make MYHOME="$work/home/bin" clean_release
     for code in LIBSTELL "$@"; do
         make -C "$code" MYHOME="$work/home/bin" \
             -j"${STELLOPT_JOBS:-$(nproc)}" clean_release
     done)
}

patched=0
vmec=0
for arg in "$@"; do
    case "$arg" in
        --patched) patched=1 ;;
        --vmec) vmec=1 ;;
        *) echo "unknown option: $arg" >&2; exit 1 ;;
    esac
done

build DIAGNO
cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno"
echo "built $work/bin/xdiagno"

if [ "$vmec" = 1 ]; then
    build VMEC2000
    cp "$src/VMEC2000/Release/xvmec2000" "$work/bin/xvmec2000"
    echo "built $work/bin/xvmec2000"
fi

if [ "$patched" = 1 ]; then
    for patch in "$here"/patches/*.patch; do
        git -C "$src" apply "$patch"
        echo "applied $(basename "$patch")"
    done
    build DIAGNO
    cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno_patched"
    git -C "$src" checkout --quiet -- .
    echo "built $work/bin/xdiagno_patched"
fi
