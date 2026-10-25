#!/usr/bin/env bash
# Short tiago_infer run on LI383: MAP, Laplace, PSIS (50 draws) and NUTS
# (20 + 20 transitions) must finish without divergences, with a finite
# evidence and k-hat, and write the CSV files the plotting tools read.
# Usage: infer_smoke.sh <tiago_infer> <VMEC input> <sensors.awk> <workdir>
set -euo pipefail
exe=$1 input=$2 awkf=$3 work=$4
mkdir -p "$work"
cat > "$work/in.nml" <<NML
&infer
  vmec_input = '$input', sensors_awk = '$awkf', output_dir = '$work/out'
  seed = 3, n_warmup = 20, n_samples = 20, max_depth = 5, psis_draws = 50,
  posterior_solves = 3, laplace_draws = 50
/
NML
"$exe" "$work/in.nml" > "$work/log.txt" 2>&1 || { tail -20 "$work/log.txt"; exit 1; }
for f in summary.txt parameters.csv profiles.csv nuts_draws.csv iota.csv surfaces_truth.csv \
    surfaces_draws.csv truth_profiles.csv; do
    [[ -s "$work/out/$f" ]] || { echo "missing $f"; exit 1; }
done
cat "$work/out/summary.txt"
awk -F' = ' '$1 == "nuts_divergent" && $2 != 0 { bad = 1 }
    $1 == "psis_khat" && ($2 + 0 > 0.7 || $2 + 0 < -10) { bad = 1 }
    $1 == "evidence_proper" && $2 != "T" { bad = 1 }
    END { exit bad }' "$work/out/summary.txt"
