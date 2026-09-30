#!/usr/bin/env bash
# Synthetic equilibrium reconstructions of LI383 (VMEC++ + Tiago), see README.md.
# Usage: ./run.sh [build dir with tiago_reconstruct (default ../../build-vmecpp)]
#        [number of noise seeds for the pull study (default 20)]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
build=$(cd "${1:-$here/../../build-vmecpp}" && pwd)
seeds=${2:-20}
tiago="$build/tiago_reconstruct"
plot_pulls="$build/tiago_plot_pulls"
work="$here/_work"
results="$here/results"
jobs=${JOBS:-$(nproc)}
[ -x "$tiago" ] || { echo "$tiago not found (configure with -DTIAGO_ENABLE_RECONSTRUCTION=ON)"; exit 1; }
mkdir -p "$work" "$results"

SIMSOPT_RAW=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347
[ -s "$work/input.li383_low_res" ] ||
    curl -fsSL "$SIMSOPT_RAW/tests/test_files/input.li383_low_res" -o "$work/input.li383_low_res"
awk -v R0=1.42 -v rs=0.85 -v nphi=6 -v nth=4 -v nprobe=20 -v dir="$work" -f "$here/sensors.awk"

# case <name> <extra namelist lines>
case_nml() {
    local name=$1; shift
    mkdir -p "$work/$name"
    cat > "$work/$name/recon.nml" <<EOF
&reconstruction
  vmec_input = '$work/input.li383_low_res'
  output_dir = '$work/$name'
  flux = '$work/flux.diagno'
  segrog = '$work/seg.diagno'
  bprobes = '$work/probes.diagno'
  synthesize_measurements = .true.
  sigma_relative = 0.01, sigma_flux = 1.0e-4, sigma_segrog = 1.0e-8, sigma_bprobe = 1.0e-7
  samples = 4, plasma_nphi = 64, plasma_ntheta = 64
  max_iterations = 20
$(printf '  %s\n' "$@")
/
EOF
}
profile=(
  "parameters = 'PRES_SCALE', 'CURTOR', 'AC(1)', 'PHIEDGE'"
  "truth_values = 1.0, -174250.0, 1436035.6, 0.514386"
  "start_values = 0.7, -150000.0, 1200000.0, 0.50")
shape=(
  "parameters = 'PRES_SCALE', 'CURTOR', 'AC(1)', 'PHIEDGE', 'RBC(0,1)', 'ZBS(0,1)', 'RBC(1,1)', 'ZBS(1,1)'"
  "truth_values = 1.0, -174250.0, 1436035.6, 0.514386, 0.27073, 0.46465, -0.135, 0.16516"
  "start_values = 0.7, -150000.0, 1200000.0, 0.50, 0.26, 0.48, -0.14, 0.16")

run() {
    local name=$1 start
    start=$(date +%s.%N)
    OMP_NUM_THREADS=1 "$tiago" "$work/$name/recon.nml" > "$work/$name/log.txt" 2>&1 ||
        { tail "$work/$name/log.txt"; return 1; }
    echo "$(date +%s.%N) - $start" | bc | awk '{ printf "%.1f s wall (1 thread)\n", $1 }' \
        > "$work/$name/time.txt"
}

case_nml exact "${profile[@]}" "add_noise = .false."
case_nml noisy "${profile[@]}" "seed = 1"
case_nml shape "${shape[@]}" "seed = 1"
for s in $(seq 1 "$seeds"); do case_nml "pull_$s" "${profile[@]}" "seed = $s"; done
export -f run; export tiago work
printf '%s\n' exact noisy shape $(for s in $(seq 1 "$seeds"); do echo "pull_$s"; done) |
    xargs -P "$jobs" -I{} bash -c 'run {}'

for name in exact noisy shape; do
    mkdir -p "$results/$name"
    cp "$work/$name"/{summary.txt,parameters.csv,time.txt} "$results/$name/"
done
cp "$work/shape"/*.png "$results/shape/"   # figures of the most complete case
# pulls (fitted - truth) / sigma of every parameter over the noise seeds
{
    echo "seed,parameter,fitted,sigma,truth,pull"
    for s in $(seq 1 "$seeds"); do
        # the quoted name may contain commas: the numbers are the last six fields
        awk -F, -v s="$s" 'NR > 1 { name = $1; for (i = 2; i <= NF - 6; i++) name = name "," $i
            printf "%d,%s,%s,%s,%s,%.6f\n", s, name, $(NF-4), $(NF-3), $(NF-2), ($(NF-4) - $(NF-2)) / $(NF-3) }' \
            "$work/pull_$s/parameters.csv"
    done
} > "$results/pulls.csv"
awk -F'"' 'NR > 1 { split($3, v, ","); p = $2; n[p]++; m[p] += v[5]; q[p] += v[5] * v[5] }
    END { printf "%-12s %6s %8s %8s\n", "parameter", "seeds", "mean", "std";
          for (p in n) printf "%-12s %6d %8.3f %8.3f\n", p, n[p], m[p] / n[p],
              sqrt((q[p] - m[p] * m[p] / n[p]) / (n[p] - 1)) }' "$results/pulls.csv" | tee "$results/pulls.txt"
[ -x "$plot_pulls" ] && "$plot_pulls" "$results/pulls.csv" "$results/pulls.png"
echo "results in $results"
