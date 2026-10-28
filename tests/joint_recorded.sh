#!/bin/sh
set -eu
binary=$1
input=$2
sensors=$3
work=$4
mkdir -p "$work"
awk -v R0=1.42 -v rs=0.85 -v nphi=4 -v nth=4 -v nprobe=12 -v dir="$work" -f "$sensors"
python3 "$(dirname "$0")/mock_archive.py" start "$work"
trap 'python3 "$(dirname "$0")/mock_archive.py" stop "$work"' EXIT
"$binary" "$input" "$work" "$work/work" "$work" "${5:-0}"
