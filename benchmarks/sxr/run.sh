#!/usr/bin/env bash
# TIAGO soft X-ray chords against STELLOPT's chisq_sxr on the same VMEC
# equilibrium, chords and profiles (README.md).
#
# Usage: ./run.sh <TIAGO build directory>
#   builds STELLOPT at a pinned commit in its CI container (docker), runs one
#   STELLOPT evaluation of the LI383 low-resolution case with 36 SXR chords,
#   evaluates the same chords with TIAGO on STELLOPT's wout, and writes
#   _work/results.md, _work/results.csv and _work/sxr-stellopt.svg (matplotlib).
set -euo pipefail

STELLOPT_REPO=${STELLOPT_REPO:-https://github.com/PrincetonUniversity/STELLOPT.git}
STELLOPT_COMMIT=${STELLOPT_COMMIT:-d05a4a11231b191544016f24be1b2ead61a12c8f}
IMAGE=${STELLOPT_IMAGE:-lazerson/stellopt-compile@sha256:dd37050e6ebb1585cabce764d458db7edd8c4e02d91af3f22c468935d9e4f162}
INPUT_URL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/input.li383_low_res

build=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
work="$here/_work"
src="$work/STELLOPT"
run="$work/run"
mkdir -p "$work" "$run"

cmake --build "$build" --target benchmark_sxr

if [ ! -d "$src/.git" ]; then
    git clone --filter=blob:none "$STELLOPT_REPO" "$src"
fi
git -C "$src" fetch --quiet origin "$STELLOPT_COMMIT" || true
git -C "$src" checkout --quiet "$STELLOPT_COMMIT"

in_container() {
    docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e MACHINE=debian \
        -e STELLOPT_PATH=/stellopt -e MYHOME=/stellopt/bin \
        -v "$src":/stellopt -v "$run":/run_dir -w "$2" "$IMAGE" sh -c "$1"
}

if [ ! -x "$src/bin/xstelloptv2" ]; then
    mkdir -p "$src/bin"
    in_container "./build_all -o release -j ${STELLOPT_JOBS:-1}" /stellopt \
        > "$work/stellopt-build.log" 2>&1
fi

curl -fsSL "$INPUT_URL" -o "$run/input.li383_low_res"
"$build/benchmark_sxr" stellopt "$run/sxr.nml"
{
    sed '/^ *\/ *$/q' "$run/input.li383_low_res"
    printf "&OPTIMUM\n  NFUNC_MAX = 1\n  EQUIL_TYPE = 'VMEC2000'\n  OPT_TYPE = 'ONE_ITER'\n"
    printf "  FTOL = 1.0E-6\n  XTOL = 1.0E-6\n  GTOL = 1.0E-30\n  EPSFCN = 1.0E-4\n"
    printf "  MODE = 1\n  FACTOR = 100.0\n"
    cat "$run/sxr.nml"
    printf "/\n&END\n"
} > "$run/input.sxr"
in_container "/stellopt/bin/xstelloptv2 input.sxr" /run_dir > "$work/stellopt-run.log" 2>&1

wout=$(ls "$run"/wout_sxr*.nc | head -1)
"$build/benchmark_sxr" tiago "$wout" "$work/tiago.csv" > "$work/tiago-run.log"
MPLBACKEND=Agg python3 "$here/compare.py" "$run/stellopt.sxr" "$work/tiago.csv" "$work" \
    "$work/sxr-stellopt.svg"
{
    echo "STELLOPT $STELLOPT_COMMIT, image $IMAGE"
    echo "TIAGO $(git -C "$here" rev-parse HEAD), wout $(basename "$wout")"
} >> "$work/results.md"
cat "$work/results.md"
