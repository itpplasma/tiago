#!/bin/sh
set -eu
cli=$1
coils=$2
geometry=$3
out=$4
mkdir -p "$out"
"$cli" --coils "$coils" --mse "$geometry" --output-dir "$out" \
    --mse-out angles.csv --response-out response.csv
# Analytic square-coil ratio at channel 1: atan(1/2), J=(0.4,-0.2).
awk -F, '
function abs(x) {return x < 0 ? -x : x}
NR == 1 {if ($0 != "label,gamma_rad") exit 1}
NR == 2 {if ($1 != "NOMINAL_01" || abs($2-atan2(1,2)) > 1e-13) exit 1}
END {if (NR != 3) exit 1}
' "$out/angles.csv"
awk -F, '
function abs(x) {return x < 0 ? -x : x}
$1 == "mse" && $2 == "NOMINAL_01" {
    if ($3 == 1) {if (abs($4-0.4) > 1e-13) exit 1; first=1}
    if ($3 == 2) {if (abs($4+0.2) > 1e-13) exit 1; second=1}
}
END {if (!first || !second) exit 1}
' "$out/response.csv"
printf '1\nBAD 0 0 0\n' > "$out/bad.txt"
if "$cli" --coils "$coils" --mse "$out/bad.txt" --output-dir "$out" \
    > "$out/error.log" 2>&1; then
    echo 'malformed MSE input accepted' >&2
    exit 1
fi
