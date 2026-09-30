# TIAGO: Toolkit for Inference and Analysis of Generalized Observables

Tiago is a modern Fortran toolkit for magnetic diagnostic studies. It reads
DIAGNO diagnostic files and STELLOPT coil files, and evaluates flux loops and
segmented Rogowski probes from coil currents (Biot–Savart) and, optionally, the
plasma currents of a VMEC equilibrium. Results are benchmarked against
STELLOPT's `xdiagno`.

## Highlights
- **Drop-in DIAGNO ingestion** – Load STELLOPT-compatible flux loops and
  segmented Rogowski files (including per-point effective areas) unchanged.
- **Vacuum solver + CLI** – Batch evaluation with `tiago_vacuum_cli` using
  libneo's Biot-Savart solver; comprehensive CMake test suite.
- **Cross-code validation** – `benchmarks/xdiagno/` builds STELLOPT's `xdiagno`
  at a pinned commit and compares accuracy and run time on identical inputs.

## Quick start
```bash
# Configure + build
cmake -S . -B build
cmake --build build
# optional: -DTIAGO_LIBNEO_TAG=<commit> to use another libneo, and
# -DCMAKE_Fortran_FLAGS="-march=native" for a machine-specific build
# (about 2x faster Biot-Savart kernels through AVX vectorisation)

# or the convenience wrappers
make              # configures + builds
make test         # rebuilds and runs ctest --output-on-failure
```
Pass extra cache entries through `CMAKE_ARGS`, e.g.
`make CMAKE_ARGS="-DTIAGO_LIBNEO_TAG=<commit>"`.

## Vacuum CLI usage
```
./build/tiago_vacuum_cli \
    tests/data/coil_sample.neo \
    tests/data/fluxloop_sample.diagno \
    tests/data/segrog_sample.diagno \
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8 \
    --nfp 3 --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv
```
This command emits `tiago_flux.csv` and `tiago_segrog.csv` in the output
directory (names via `--flux-out`, `--segrog-out`). `--samples` sets the points
per segment (midpoint rule, as DIAGNO's `int_type='midpoint'`); `--gauss` uses
Gauss-Legendre points instead, which converge far faster (8 Gauss points are
typically more accurate than thousands of midpoint samples). `--nfp` matches the field-period geometry,
while `--flux-turns` / `--segrog-turns` accept text files of `label scale`
pairs so Tiago mirrors DIAGNO's namelist-based `flux_turns` and
`segrog_turns` arrays.

## Magnetic probes and response matrices
`--bprobes file` evaluates DIAGNO magnetic probes (`x y z theta_inc phi_inc
eff_area` per row, angles in degrees; `--rphiz` for `R phi z`), written to
`tiago_bprobes.csv` as `eff_area * B . n`. `--response-out file` writes every
signal per unit EXTCUR of each coil group (`kind,label,group,value`, as
DIAGNO's `-mutual`), so for other coil currents the vacuum signals are
`sum_g M_g EXTCUR_g` without re-running Biot–Savart; this makes fitting coil
currents a linear least-squares problem.

## Plasma response
`--plasma-wout wout.nc` adds the field of the plasma currents of a VMEC
equilibrium to every flux loop and segmented Rogowski. The VMEC boundary is a
flux surface, so outside it the plasma acts like the sheet current
`mu0 K = n x B` on the boundary (virtual casing); Tiago evaluates its vector
potential and field directly on a `--plasma-nphi` (per field period) by
`--plasma-ntheta` grid, default 64 x 64. A warning is printed when a sensor is
closer to the boundary than two grid spacings. As in DIAGNO, loops that link
the plasma poloidally must be flagged `idia = 1`, which adds the plasma
toroidal flux `phiedge`. Pass `""` as coil file for plasma-only signals.

## Benchmark against STELLOPT xdiagno
A manual, reproducible accuracy and performance comparison with the reference
code lives in [`benchmarks/xdiagno/`](benchmarks/xdiagno/README.md). It is not
part of `ctest`.

## Directory layout
```
cmake/                     # Toolchain helpers
CMakeLists.txt             # FetchContent bridge + targets
benchmarks/xdiagno/        # Manual benchmark vs STELLOPT xdiagno
src/                       # Fortran diagnostics, solver, CLI
tests/                     # Sample inputs + regression drivers
```
