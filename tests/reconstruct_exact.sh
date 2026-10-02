#!/usr/bin/env bash
# Noise-free synthetic reconstruction of LI383 must recover the true
# parameters (from start values 10-30 % off) to 0.05 of their uncertainty:
# the fit stops once the predicted chi^2 decrease is below the solver noise
# floor of 1e-3, i.e. within about sqrt(1e-3) = 0.03 sigma of the optimum.
# Usage: reconstruct_exact.sh <tiago_reconstruct> <VMEC input> <sensors.awk> <workdir> [conservative]
set -euo pipefail
tiago=$1 input=$2 sensors=$3 work=$4
conservative=.false.
if [[ ${5:-} == conservative ]]; then conservative=.true.; fi
mkdir -p "$work"
awk -v R0=1.42 -v rs=0.85 -v nphi=4 -v nth=4 -v nprobe=12 -v dir="$work" -f "$sensors"
cat > "$work/recon.nml" <<NML
&reconstruction
  vmec_input = '$input'
  output_dir = '$work/out'
  flux = '$work/flux.diagno', segrog = '$work/seg.diagno', bprobes = '$work/probes.diagno'
  parameters = 'PRES_SCALE', 'CURTOR', 'AC(1)', 'PHIEDGE'
  truth_values = 1.0, -174250.0, 1436035.6, 0.514386
  start_values = 0.7, -150000.0, 1200000.0, 0.50
  synthesize_measurements = .true., add_noise = .false.
  sigma_relative = 0.01, sigma_flux = 1.0e-4, sigma_segrog = 1.0e-8, sigma_bprobe = 1.0e-7
  plasma_nphi = 32, plasma_ntheta = 32, max_iterations = 20
  plasma_conservative = $conservative
/
NML
"$tiago" "$work/recon.nml" > "$work/log.txt" 2>&1 || { tail -20 "$work/log.txt"; exit 1; }
grep -A8 "chi^2 start" "$work/log.txt"
awk -F'"' 'NR > 1 { split($3, v, ","); pull = (v[3] - v[5]) / v[4]
    printf "%-12s pull %10.2e\n", $2, pull; if (pull > 0.05 || pull < -0.05) bad = 1 }
    END { exit bad }' "$work/out/parameters.csv"
