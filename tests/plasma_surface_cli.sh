#!/bin/sh
# Loop integrals of the sampled-surface casing field against closed forms
# (see plasma_surface_fixture.f90); the tilted surface must be rejected.
set -eu
cli=$1
fixture=$2
out=$3
mkdir -p "$out"
"$fixture" "$out"
set -- $(cat "$out/expected.txt")
for case in poloidal:$1 toroidal_in:$2 toroidal_out:$3; do
    name=${case%%:*}
    "$cli" --plasma-surface "$out/surface.txt" --bprobes "$out/$name.diagno" \
        --output-dir "$out" --bprobe-out "$name.csv" > "$out/$name.log" 2>&1
    awk -F, -v want="${case#*:}" -v name="$name" '
    function abs(x) {return x < 0 ? -x : x}
    NR > 1 {sum += $2}
    END {
        printf "%s: %.16e expected %.16e\n", name, sum, want
        if (abs(sum - want) > 1e-9 * 12.566370614359172) exit 1
    }' "$out/$name.csv"
done
grep -q 'source = sampled flux surface' "$out/plasma_model.txt"
mkdir -p "$out/tilted"
"$fixture" "$out/tilted" tilted
if "$cli" --plasma-surface "$out/tilted/surface.txt" --bprobes "$out/poloidal.diagno" \
    --output-dir "$out/tilted" > "$out/tilted/error.log" 2>&1; then
    echo 'surface with a normal field accepted' >&2
    exit 1
fi
grep -q 'not a flux surface' "$out/tilted/error.log"
