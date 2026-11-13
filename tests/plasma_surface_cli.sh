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

for mode in nonfinite degenerate; do
    mkdir -p "$out/$mode"
    "$fixture" "$out/$mode" "$mode"
    if "$cli" --plasma-surface "$out/$mode/surface.txt" --bprobes "$out/poloidal.diagno" \
        --output-dir "$out/$mode" > "$out/$mode/error.log" 2>&1; then
        echo "invalid $mode surface accepted" >&2
        exit 1
    fi
    grep -q "$mode plasma surface" "$out/$mode/error.log"
done
mkdir -p "$out/bad-header"
printf 'missing dimensions\n' > "$out/bad-header/surface.txt"
if "$cli" --plasma-surface "$out/bad-header/surface.txt" --bprobes "$out/poloidal.diagno" \
    --output-dir "$out/bad-header" > "$out/bad-header/error.log" 2>&1; then
    echo 'malformed surface dimensions accepted' >&2
    exit 1
fi
grep -q 'malformed plasma surface grid size' "$out/bad-header/error.log"

mkdir -p "$out/bad-flux"
if "$cli" --plasma-surface "$out/surface.txt" --bprobes "$out/poloidal.diagno" \
    --plasma-linked-flux NaN --output-dir "$out/bad-flux" > "$out/bad-flux/error.log" 2>&1; then
    echo 'nonfinite linked flux accepted' >&2
    exit 1
fi
grep -q 'nonfinite linked flux' "$out/bad-flux/error.log"

for mode in null extra; do
    mkdir -p "$out/$mode"
    "$fixture" "$out/$mode" "$mode"
    if "$cli" --plasma-surface "$out/$mode/surface.txt" --bprobes "$out/poloidal.diagno" \
        --output-dir "$out/$mode" > "$out/$mode/error.log" 2>&1; then
        echo "malformed $mode surface row accepted" >&2
        exit 1
    fi
    grep -q 'missing or malformed row' "$out/$mode/error.log"
done
for mode in extra-header trailing-text; do
    mkdir -p "$out/$mode"
    if [ "$mode" = extra-header ]; then
        awk 'NR == 1 {$0 = $0 " 1"} {print}' "$out/surface.txt" > "$out/$mode/surface.txt"
    else
        cat "$out/surface.txt" > "$out/$mode/surface.txt"
        printf 'not a surface row\n' >> "$out/$mode/surface.txt"
    fi
    if "$cli" --plasma-surface "$out/$mode/surface.txt" --bprobes "$out/poloidal.diagno" \
        --output-dir "$out/$mode" > "$out/$mode/error.log" 2>&1; then
        echo "malformed $mode surface record accepted" >&2
        exit 1
    fi
done
