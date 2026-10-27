#!/bin/sh
set -eu
source_dir=$1
fortsym_build=$2
compiler=$3
if test -z "$fortsym_build"; then
    echo 'Set TIAGO_FORTSYM_BUILD to run FortSym regeneration'
    exit 77
fi
work=$(mktemp -d "${TMPDIR:-/tmp}/tiago-ece-codegen.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cmake -S "$source_dir/gen" -B "$work/build" -G Ninja \
    -DFORTSYM_BUILD="$fortsym_build" -DCMAKE_Fortran_COMPILER="$compiler"
cmake --build "$work/build" --target gen_ece
mkdir -p "$work/output/src/generated"
(cd "$work/output" && "$work/build/gen_ece")
for kernel in ece_hamiltonian_o ece_hamiltonian_x ece_hamiltonian_unmagnetized ece_test_medium; do
    diff -u "$source_dir/src/generated/$kernel.f90" "$work/output/src/generated/$kernel.f90"
done
