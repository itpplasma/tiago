#!/usr/bin/env bash
# Exercise the production namelist path with the pinned public CTH-like case.
set -euo pipefail
tiago=$1 fixture=$2 sensors=$3 work=$4
mkdir -p "$work"
awk -v R0=0.78 -v rs=0.5 -v nphi=3 -v nth=3 -v nprobe=12 -v dir="$work" -f "$sensors"
cat > "$work/recon.nml" <<NML
&reconstruction
  vmec_input = '$fixture/input.json'
  coils = '$fixture/coils.coils', coil_current_scale = 96.0, -36.0
  output_dir = '$work/out'
  flux = '$work/flux.diagno', segrog = '$work/seg.diagno', bprobes = '$work/probes.diagno'
  parameters = 'EXTCUR(1)', 'PHIEDGE', 'CURTOR', 'PRES_SCALE'
  prior_values = 4700.0, -0.035, 43229.08092460368, 432.29080924603676
  prior_sigma = 100.0, 0.001, 2000.0, 100.0
  truth_values = 4740.0, -0.0346, 44029.08092460368, 472.29080924603676
  synthesize_measurements = .true., add_noise = .false.
  sigma_relative = 0.003, sigma_flux = 1.0e-6, sigma_segrog = 1.0e-6, sigma_bprobe = 1.0e-6
  vmec_ftol = 1.0e-18
  plasma_nphi = 32, plasma_ntheta = 32, max_iterations = 12
/
NML
"$tiago" "$work/recon.nml" > "$work/log.txt" 2>&1 || { tail -20 "$work/log.txt"; exit 1; }
grep 'equilibrium derivatives: re-solved central finite differences' "$work/out/summary.txt"
awk -F'"' 'NR > 1 { split($3, v, ","); pull = (v[3] - v[5]) / v[4]
    printf "%-12s pull %10.2e\n", $2, pull; if (pull > 2 || pull < -2) bad = 1 }
    END { exit bad }' "$work/out/parameters.csv"
