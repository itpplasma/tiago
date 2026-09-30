#!/usr/bin/env bash
# Download the public coil sets, equilibria and VMEC inputs used by the benchmark
# (pinned commits).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
data="$here/_work/data"
mkdir -p "$data"

STELLOPT_RAW=https://raw.githubusercontent.com/PrincetonUniversity/STELLOPT/2f181f0d4e71f28afe076967d606b2c020c9db33
SIMSOPT_RAW=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347

fetch() {
    [ -s "$data/$2" ] && { echo "cached $2"; return; }
    curl -fsSL "$1" -o "$data/$2"
    echo "fetched $2"
}
fetch "$STELLOPT_RAW/BENCHMARKS/FIELDLINES_TEST/coils.NCSX"       coils.NCSX
fetch "$STELLOPT_RAW/BENCHMARKS/FIELDLINES_TEST/coils.NCSX_nfp1"  coils.NCSX_nfp1
fetch "$SIMSOPT_RAW/tests/test_files/coils.M16N08"                coils.M16N08
fetch "$STELLOPT_RAW/BENCHMARKS/DIAGNO_TEST/wout_ncsx.nc"        wout_ncsx.nc
fetch "$STELLOPT_RAW/BENCHMARKS/DIAGNO_TEST/input.ncsx"          input.ncsx
fetch "$SIMSOPT_RAW/tests/test_files/input.li383_low_res"      input.li383_low_res
