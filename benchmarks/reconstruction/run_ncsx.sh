#!/usr/bin/env bash
# Full-resolution NCSX reconstruction with coil currents and free-boundary
# consistency (VMEC++ + Tiago), and the same equilibrium with VMEC2000 for
# comparison when it is built. See README.md.
# Usage: ./run_ncsx.sh [build dir with tiago_reconstruct (default ../../build-vmecpp)]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
build=$(cd "${1:-$here/../../build-vmecpp}" && pwd)
tiago="$build/tiago_reconstruct"
vacuum="$build/tiago_vacuum_cli"
xvmec="$here/../xdiagno/_work/bin/xvmec2000"
data="$here/../xdiagno/_work/data"
work="$here/_work/ncsx"
results="$here/results/ncsx"
[ -x "$tiago" ] || { echo "$tiago not found (configure with -DTIAGO_ENABLE_RECONSTRUCTION=ON)"; exit 1; }
"$here/../xdiagno/fetch_data.sh" > /dev/null   # pinned coils.NCSX, input.ncsx, wout_ncsx.nc
mkdir -p "$work" "$results"

# Fixed-boundary version of STELLOPT's NCSX input (MPOL 11, NTOR 6, NS up to 99);
# tiago_reconstruct takes the boundary from the free-boundary wout_ncsx.nc.
sed -e 's/LFREEB = T/LFREEB = F/' -e '/MGRID_FILE/d' "$data/input.ncsx" > "$work/input.ncsx_fixed"
awk -v R0=1.44 -v rs=0.55 -v nphi=6 -v nth=4 -v nprobe=20 -v dir="$work" -f "$here/sensors.awk"
# consistency points on the magnetic axis guess (RAXIS, ZAXIS of input.ncsx)
awk -v nfp=3 -v npts=12 \
    -v raxis="1.49569454253276 0.105806400912320 0.00721255454715878 -0.000387402652289249 -0.000202425864534069 -0.000162602353744308 -0.00000889569831063077" \
    -v zaxis="0 -0.0519492027001782 -0.00318814224021375 0.000226199929262002 0.000128803681387330 0.00000111266150452637 0.0000113732703961869" \
    -f "$here/axis_points.awk" > "$work/axis.txt"

# run <name> <command...>: wall time and the largest peak resident memory of
# the process and its descendants (from /proc; mpirun runs the solver as a child)
descendants() { local p; for p in $(pgrep -P "$1"); do echo "$p"; descendants "$p"; done; }
run() {
    local name=$1; shift
    local start peak=0 pid p hwm
    start=$(date +%s.%N)
    "$@" > "$work/$name.log" 2>&1 &
    pid=$!
    while kill -0 "$pid" 2> /dev/null; do
        for p in $pid $(descendants "$pid"); do
            hwm=$(awk '/VmHWM/ { print $2 }' "/proc/$p/status" 2> /dev/null || true)
            [ "${hwm:-0}" -gt "$peak" ] && peak=$hwm
        done
        sleep 1
    done
    wait "$pid"
    echo "$(date +%s.%N) - $start" | bc |
        awk -v p="$peak" '{ printf "%.1f s wall, %.0f MB peak\n", $1, p / 1024 }' > "$work/$name.time"
    cat "$work/$name.time"
}

nml() {
    local name=$1; shift
    cat > "$work/$name.nml" <<EOF
&reconstruction
  vmec_input = '$work/input.ncsx_fixed'
  vmec_boundary_wout = '$data/wout_ncsx.nc'
  output_dir = '$work/$name'
  coils = '$data/coils.NCSX'
  coil_extcur = '$work/input.ncsx_fixed'
  flux = '$work/flux.diagno', segrog = '$work/seg.diagno', bprobes = '$work/probes.diagno'
  consistency_points = '$work/axis.txt', consistency_sigma = 1.0e-3
  parameters = 'EXTCUR(1)', 'EXTCUR(2)', 'EXTCUR(3)', 'PHIEDGE', 'CURTOR', 'PRES_SCALE'
  truth_values = 652271.941985300, 651868.569367400, 537743.588647300, 0.497070205837336, -178606.25, 1.0
  synthesize_measurements = .true.
  sigma_relative = 0.01, sigma_flux = 1.0e-4, sigma_segrog = 1.0e-8, sigma_bprobe = 1.0e-7
  plasma_nphi = 64, plasma_ntheta = 64, vmec_ftol = 1.0e-12
$(printf '  %s\n' "$@")
/
EOF
}
# consistency floor: noise-free data at the truth, one LM step from the truth
nml truth "start_values = 652271.941985300, 651868.569367400, 537743.588647300, 0.497070205837336, -178606.25, 1.0" \
    "add_noise = .false., max_iterations = 1"
# the reconstruction: 3-20 % off, with noise
nml fit "start_values = 632000.0, 670000.0, 522000.0, 0.477, -160000.0, 0.8" \
    "add_noise = .true., seed = 1, max_iterations = 12"
run fit "$tiago" "$work/fit.nml"
run truth "$tiago" "$work/truth.nml"
cp "$work/fit"/{summary.txt,parameters.csv} "$work/fit"/*.png "$results/"
cp "$work/fit.time" "$results/time.txt"
grep "chi^2 start" "$work/truth/summary.txt" > "$results/consistency_floor.txt"

# VMEC2000 (Fortran VMEC, STELLOPT's solver) on the same fixed-boundary input,
# and the plasma signals of both equilibria on the same sensors
if [ -x "$xvmec" ]; then
    mkdir -p "$work/vmec2000"
    cp "$work/input.ncsx_fixed" "$work/vmec2000/input.ncsxf"
    (cd "$work/vmec2000" && OMP_NUM_THREADS=1 run vmec2000 mpirun --allow-run-as-root -np 1 "$xvmec" ncsxf)
    cat > "$work/vmecpp.nml" <<EOF
&reconstruction
  vmec_input = '$work/input.ncsx_fixed'
  output_dir = '$work/vmecpp'
  flux = '$work/flux.diagno', segrog = '$work/seg.diagno', bprobes = '$work/probes.diagno'
  parameters = 'PHIEDGE', truth_values = 0.497070205837336, start_values = 0.497070205837336
  synthesize_measurements = .true., add_noise = .false., vmec_ftol = 1.0e-12, max_iterations = 1
/
EOF
    OMP_NUM_THREADS=1 "$tiago" "$work/vmecpp.nml" > "$work/vmecpp.log" 2>&1
    for code in vmecpp vmec2000; do
        wout="$work/vmecpp/wout_truth.nc"
        [ "$code" = vmec2000 ] && wout="$work/vmec2000/wout_ncsxf.nc"
        "$vacuum" --plasma-wout "$wout" --flux "$work/flux.diagno" --segrog "$work/seg.diagno" \
            --bprobes "$work/probes.diagno" --output-dir "$work/signals_$code" > /dev/null 2>&1
    done
    {
        cat "$work/vmec2000.time" | sed 's/^/VMEC2000 solve (1 core): /'
        for k in flux segrog bprobes; do
            paste -d, "$work/signals_vmecpp/tiago_$k.csv" "$work/signals_vmec2000/tiago_$k.csv" |
                awk -F, -v k="$k" 'NR > 1 { d = $2 - $4; d = d < 0 ? -d : d; b = $4 < 0 ? -$4 : $4
                    if (d > md) md = d; if (b > mb) mb = b; n++ }
                    END { printf "%-8s %3d signals: max |VMEC++ - VMEC2000| / max |signal| = %.1e\n", k, n, md / mb }'
        done
    } | tee "$results/vmec2000.txt"
fi
echo "results in $results"
