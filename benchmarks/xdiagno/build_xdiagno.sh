#!/usr/bin/env bash
# Build STELLOPT's xdiagno (LIBSTELL + DIAGNO) at a pinned commit.
#
# Usage: ./build_xdiagno.sh            # -> _work/bin/xdiagno
#        ./build_xdiagno.sh --vac-nfp  # additionally _work/bin/xdiagno_vacnfp
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
mkdir -p "$work/bin"

if [ ! -d "$src/.git" ]; then
    git clone --filter=blob:none "$STELLOPT_REPO" "$src"
fi
git -C "$src" fetch --quiet origin "$STELLOPT_COMMIT" || true
git -C "$src" checkout --quiet "$STELLOPT_COMMIT"

build() {
    # STELLOPT's makefiles install into $HOME/bin; redirect HOME to keep it local.
    (cd "$src" && HOME="$work/home" MACHINE="$STELLOPT_MACHINE" STELLOPT_PATH="$src" \
        ./build_all -j"$(nproc)" LIBSTELL DIAGNO)
}

mkdir -p "$work/home/bin"
git -C "$src" checkout --quiet -- DIAGNO/Sources/diagno.f90
build
cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno"
echo "built $work/bin/xdiagno"

if [ "${1:-}" = "--vac-nfp" ]; then
    # Stock DIAGNO leaves nfp=0 in -vac mode, so iflflg=1 loops evaluate to NaN.
    # This test-only patch takes nfp from the coil file's 'periods' line.
    git -C "$src" apply "$here/diagno_vac_nfp.patch"
    build
    cp "$src/DIAGNO/Release/xdiagno" "$work/bin/xdiagno_vacnfp"
    git -C "$src" checkout --quiet -- DIAGNO/Sources/diagno.f90
    echo "built $work/bin/xdiagno_vacnfp"
fi
