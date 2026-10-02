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
  exact straight-filament Biot–Savart kernels (Hanson & Hirshman); coil files
  are read with libneo.
- **Equilibrium reconstruction** – `tiago_reconstruct` fits VMEC input
  parameters (profiles, total current and flux, boundary shape, coil currents)
  to magnetic measurements with Levenberg–Marquardt. Equilibria come from
  VMEC++ through an ISO C adapter, and the Jacobian is exact: VMEC++'s implicit
  adjoint for the equilibrium, reverse mode in Tiago for the diagnostics. No
  Python is involved at any stage.
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
The standard `LIBNEO_BRANCH` CMake or environment variable overrides the
selected libneo ref for release validation; a nonempty CMake value takes
precedence over the environment. Removing it restores `TIAGO_LIBNEO_TAG`.

## Vacuum CLI usage
```
./build/tiago_vacuum_cli \
    --coils tests/data/coil_sample.neo \
    --flux tests/data/fluxloop_sample.diagno \
    --segrog tests/data/segrog_sample.diagno \
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8
```
Each of `--flux`, `--segrog` and `--bprobes` is optional (at least one is
needed). This command emits `tiago_flux.csv` and `tiago_segrog.csv` in the output
directory (names via `--flux-out`, `--segrog-out`). `--samples` sets the points
per segment (midpoint rule, as DIAGNO's `int_type='midpoint'`); `--gauss` uses
Gauss-Legendre points instead, which converge far faster (8 Gauss points are
typically more accurate than thousands of midpoint samples). `--nfp` matches the field-period geometry,
while `--flux-turns` / `--segrog-turns` accept text files of `label scale`
pairs so Tiago mirrors DIAGNO's namelist-based `flux_turns` and
`segrog_turns` arrays. The NCSX cases in `tests/cases/` come with such files;
their coil sets are fetched by `benchmarks/xdiagno/fetch_data.sh`:
```
benchmarks/xdiagno/fetch_data.sh
./build/tiago_vacuum_cli \
    --coils benchmarks/xdiagno/_work/data/coils.NCSX \
    --flux tests/cases/ncsx_nfp3/fluxloop.diagno \
    --segrog tests/cases/ncsx_nfp3/segrog.diagno \
    --output-dir build/ncsx --seg-area 3.40e-4 --gauss --samples 8 --nfp 3 \
    --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv
```

`--coil-extcur file` sets the coil-group currents from a VMEC `&INDATA` input
(`EXTCUR(i) = v`, `EXTCUR = a, b, ...`, slices `EXTCUR(i:j) = ...`, repeat
counts `n*v`) or a plain list of numbers. As in DIAGNO, each group's currents
are normalised to its first point and scaled by `EXTCUR(g)`; groups the file
leaves out are switched off (`EXTCUR = 0`, with a warning). Without
`--coil-extcur` the currents of the coil file are used as they are.

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
toroidal flux `phiedge`. Leave out `--coils` for plasma-only signals.

`--plasma-response-out file` writes the derivative of every signal's plasma
part with respect to the VMEC boundary field coefficients (`bsupumnc`,
`bsupvmnc` extrapolated to `s = 1`, per mode `m, n`). At fixed boundary shape
the plasma part is linear in them, so the Jacobian is exact and
`sum(value * coefficient)` reproduces it; together with `--response-out` for
the coil currents this is the linear part of an equilibrium reconstruction.
`--plasma-shape-response-out file` writes the derivatives with respect to the
boundary geometry coefficients (`rmnc`, `zmns`, and `rmns`, `zmnc` for
`lasym`) at fixed field coefficients, computed in reverse mode through the
sheet current and the Biot–Savart kernel.

Derivatives with respect to equilibrium parameters (`PHIEDGE`, `CURTOR`,
`PRES_SCALE`, ...) need the response of the equilibrium itself. That response
comes from `tiago_reconstruct` (below). The `equilibrium` suite of the xdiagno
benchmark checks the chain rule against central finite differences over VMEC
runs.

## Equilibrium reconstruction
`tiago_reconstruct input.nml` fits the parameters of a VMEC input to measured
flux loops, segmented Rogowskis and B-probes. It minimises

    chi^2 = sum_i ((S_i(x) - m_i) / sigma_i)^2 + consistency + priors

with Levenberg–Marquardt (in Tiago, on LAPACK), then reports the linearised covariance
`(J^T J)^-1`. Each evaluation of `S(x)` does three things. VMEC++ solves the
fixed-boundary equilibrium. Tiago builds the boundary sheet current from the
edge field, and evaluates the plasma signals as described above. Coil
contributions are added through the response matrix.

**Method.** The Jacobian `J = dS/dx` is exact. It costs one factorisation per
iterate, plus two cheap solves per residual row. It is built in three steps:

1. **Signals from the edge field.** `dS/dy` comes from the plasma field and
   shape responses. Here `y` is the edge field (`B^u`, `B^v` extrapolated to
   `s = 1`) together with the boundary `rmnc`, `zmns`.
2. **Edge field from the solver state.** `src/equilibrium/vmec_edge.f90`
   computes `y` from VMEC++'s spectral state in Fortran, and runs it in
   reverse mode. It reproduces VMEC++'s wout to 2e-14.
3. **Solver state from the parameters.** The adapter
   (`src/equilibrium/vmecpp_adapter.cc`) applies VMEC++'s implicit adjoint. It
   uses the exact force Jacobian generated by Enzyme, with the m = 1 gauge
   pinned. The profile and `PHIEDGE` chains are analytic, in Fortran.

**Scaling.** VMEC's forces on one flux surface depend only on the geometry of
that surface and its two neighbours. The one exception is the boundary, which
enters every surface through the spectral-condensation constraint. The
interior force Jacobian is therefore block tridiagonal in the radial surfaces,
and the adapter uses this in three ways:

- It extracts the blocks with `3 x (unknowns per surface)` forward-mode
  Hessian-vector products, probing every third surface at once.
- It probes the boundary columns one at a time. They give the coupling of the
  adjoint to the prescribed entries directly, so no per-row reverse product is
  needed.
- It factorises the blocks with a block LU. Time and memory grow as `ns m^3`
  and `ns m^2`, where `ns` is the number of surfaces and `m` the number of
  unknowns per surface. A dense LU grows as `(ns m)^3` and `(ns m)^2`.

Probing and the per-row passes run on thread-local copies of the VMEC++ model,
and VMEC++ solves use its own radial threads. Full-resolution NCSX (MPOL 11,
NTOR 6, NS 99, about 45k unknowns) would need 16 GB for a dense LU. With the
block LU, a Jacobian takes about 28 s on 4 cores. A 6-parameter reconstruction
finishes in under 3 minutes and needs 1.5 GB in total.

Three further savings:

- Solves after the first hot-restart from the last equilibrium.
- At a fixed boundary the plasma signals are exactly `R y_B`, with the field
  response `R` cached per boundary shape. Profile and coil fits therefore
  evaluate the sheet once.
- The fit stops when the predicted chi^2 decrease falls below the solver's
  noise floor, instead of probing it with further solves.

`TIAGO_VMECPP_VERIFY=1` checks every extracted block row and coupling column
against an isolated Hessian-vector product; `ctest` runs this on LI383.

**Accuracy.** The Jacobian agrees with central differences of re-solved
equilibria to about 1e-4 on LI383. At NCSX mode numbers the agreement is only
0.1–2.5 %, and VMEC++'s own `vmecpp.autodiff` reproduces the same gap (#27,
upstream). It does not bias the fit, which converges to the minimum of the
exactly evaluated chi^2.

**Parameters.** `PRES_SCALE`, `CURTOR`, `PHIEDGE`, `AM(k)`, `AI(k)`, `AC(k)`
(power-series profiles), `RBC(n,m)`, `ZBS(n,m)` and coil currents
`EXTCUR(g)`. The VMEC input must be fixed-boundary and stellarator-symmetric.
Free-boundary consistency is imposed instead through `consistency_points`:
virtual probes inside the plasma, where the coil field plus the sheet field
must vanish.

**Building.** Reconstruction is optional, because VMEC++'s exact force
Jacobian needs Clang with the matching Enzyme plugin:
```bash
sudo apt-get install clang-20 llvm-20-dev libclang-20-dev lld-20 libomp-20-dev \
    libzstd-dev liblapack-dev ninja-build
git clone --depth 1 --branch v0.0.264 https://github.com/EnzymeAD/Enzyme.git
cmake -G Ninja -S Enzyme/enzyme -B enzyme-build -DCMAKE_BUILD_TYPE=Release \
    -DLLVM_DIR=$(llvm-config-20 --cmakedir) \
    -DClang_DIR=$(llvm-config-20 --prefix)/lib/cmake/clang
ninja -C enzyme-build ClangEnzyme-20
CXX=clang++-20 cmake -S . -B build-vmecpp -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DTIAGO_ENABLE_RECONSTRUCTION=ON \
    -DTIAGO_ENZYME_PLUGIN=$PWD/enzyme-build/Enzyme/ClangEnzyme-20.so
cmake --build build-vmecpp
```
CMake fetches the following at pinned commits:
- VMEC++, built as a C++ library without its Python module;
- fortplot, for the figures.

VMEC++ also downloads HDF5 and netCDF-C. If those archives are not reachable,
point `-DTIAGO_VMECPP_HDF5_URL=` and `-DTIAGO_VMECPP_NETCDF_URL=` at local
copies.

**Input** (namelist `&reconstruction`; only `vmec_input` and `parameters` are
required):

| variable | meaning |
|---|---|
| `vmec_input` | VMEC `INDATA` file (converted with VMEC++'s `indata2json`) |
| `vmec_boundary_wout` | take the fixed boundary from this wout's last surface (e.g. a free-boundary equilibrium) |
| `parameters` | names as above, e.g. `'PRES_SCALE', 'CURTOR', 'RBC(0,1)', 'EXTCUR(2)'` |
| `start_values`, `parameter_scale` | start (default: input value) and scale per parameter |
| `prior_values`, `prior_sigma` | Gaussian priors (`prior_sigma > 0` enables one) |
| `flux`, `segrog`, `bprobes`, `seg_area`, `rphiz` | diagnostic files, as in `tiago_vacuum_cli` |
| `coils`, `coil_extcur` | coil file and currents, needed for `EXTCUR` parameters and consistency points |
| `measurements` | CSV `kind,label,value,sigma` (kind `flux`, `segrog` or `bprobe`) |
| `consistency_points`, `consistency_sigma` | rows `x y z` inside the plasma, and the tolerance on the field there [T] |
| `synthesize_measurements`, `truth_values`, `add_noise`, `seed` | synthetic data from the model at `truth_values` |
| `sigma_relative`, `sigma_flux`, `sigma_segrog`, `sigma_bprobe` | synthetic uncertainty `sigma_relative \|S\| + sigma_<kind>` |
| `plasma_nphi`, `plasma_ntheta`, `samples`, `gauss` | sheet grid and segment quadrature |
| `vmec_ftol`, `vmec_niter` | VMEC++ convergence (default `1e-14`, `20000`) |
| `max_iterations`, `check_jacobian`, `fd_step` | LM iterations; compare `J` with finite differences first |
| `output_dir` | results directory |

**Output** (in `output_dir`):
- `summary.txt`, `parameters.csv` (with start, fitted, sigma and truth),
  `covariance.csv`, `signals.csv`, `history.csv`;
- `wout_initial.nc`, `wout_fit.nc` and `wout_truth.nc`;
- figures: `signals_<kind>.png`, `residuals.png`, `convergence.png`,
  `profiles_{pressure,iota,current}.png` (with a pressure uncertainty band)
  and `boundary.png`.

[`benchmarks/reconstruction/`](benchmarks/reconstruction/README.md) holds the
synthetic reconstructions, with results, timings and figures:
- LI383, including a 20-seed pull study;
- full-resolution NCSX with coil currents and consistency points.

It also compares VMEC++ with VMEC2000 and gives cost estimates against
STELLOPT.

## Benchmark against STELLOPT xdiagno
A manual, reproducible accuracy and performance comparison with the reference
code lives in [`benchmarks/xdiagno/`](benchmarks/xdiagno/README.md). It is not
part of `ctest`.

## Directory layout
```
cmake/                     # Toolchain helpers
CMakeLists.txt             # FetchContent bridge + targets
benchmarks/xdiagno/        # Manual benchmark vs STELLOPT xdiagno
benchmarks/reconstruction/ # Synthetic reconstructions (tiago_reconstruct)
src/                       # Fortran diagnostics, solver, CLI
src/equilibrium/           # VMEC++ C adapter, Fortran bindings, edge field
src/reconstruct/           # Residuals, Jacobian, covariance, figures
tests/                     # Sample inputs + regression drivers
```
