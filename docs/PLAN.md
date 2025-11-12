# Tiago Implementation Plan

## References

- **STELLOPT / DIAGNO**
  - `STELLOPT/DIAGNO/Sources/diagno_flux.f90`, loop file format description in
    lines 9–34 and flux accumulation in lines 336–352.
  - `STELLOPT/DIAGNO/Sources/diagno_rogowski_new.f90`, segmented Rogowski I/O
    and signal assembly in lines 8–35 and 320–366.
  - `STELLOPT/DIAGNO/Sources/diagno_init_vmec.f90`, virtual-casing/volume
    integral bootstrapping in lines 71–197.
- **libneo**
  - `libneo/src/field/biotsavart_field.f90` lines 10–66 (`biotsavart_field_t`
    interface wrapping the Biot–Savart kernels).
  - `libneo/src/field/biotsavart.f90` lines 1–135 (vector potential and B-field
    evaluation on polygonal coils).
  - `libneo/src/spline_vmec_data.f90` lines 1–120 (existing VMEC spline + field
    synthesis routines we can reuse for Stage 4).

## Stage 0 – Scaffolding (done)
- Repository created under `itpplasma/tiago` (GitHub, private) with minimal
  CMake project fetching libneo and placeholder docs.
- No Fortran sources or tests checked in yet.

## Stage 1 – Diagnostic ingestion (geometry only)
Deliverables:
1. Diagnostic parsers in `src/diagnostics/flux_loops.f90` and
   `src/diagnostics/segmented_rogowski.f90` mirroring the file formats encoded
   in `diagno_flux.f90:L9-L34` and `diagno_rogowski_new.f90:L8-L35`.
2. Metadata registry (YAML/JSON) identifying diagnostic families, effective
   areas, and reference orientations.
3. CLI subcommand `tiago diag lint <file>` that validates coil files before
   they reach solvers.

Notes:
- We will keep the on-disk layouts identical to DIAGNO so existing pipelines can
  reuse their specification files.
- Parsers will be pure (return derived types) so they can feed either libneo or
  legacy workflows.

## Stage 2 – Vacuum Biot–Savart solver
Deliverables:
1. `src/solver/vacuum_forward.f90` that takes coil geometry + diagnostics and
   evaluates flux/voltage using libneo's `biotsavart_field_t` interface
   (see `biotsavart_field.f90:10-66`).
2. Pluggable numerical quadrature that operates per-diagnostic to avoid
   DIAGNO's global `int_step` cost (`diagno_flux.f90:336-352`).
3. `ctest` target `tiago_vacuum_smoke` running on the downloaded public coil
   file (see `docs/TESTING.md`).

Notes:
- Treat each diagnostic as an embarrassingly parallel job, but keep an OpenMP
  toggle identical to libneo's `ENABLE_OPENMP` option to match cluster setups.
- Separate accumulation of coil contributions from geometry traversal so we can
  reuse the same code path when VMEC plasma currents are introduced.

## Stage 3 – Cross-code validation harness
Deliverables:
1. `scripts/run_xdiagno.py` that shells out to `xdiagno -vac` with the same
   diagnostic files (reference STELLOPT program entry in
   `STELLOPT/DIAGNO/Sources/diagno.f90:60-188`).
2. `ctest` target `tiago_vs_xdiagno` that compares Tiago flux/voltage arrays
   against DIAGNO within configurable tolerances and produces a diff report in
   `tests/output/`.
3. CI recipe (GitHub Actions runner) that executes both tests on every PR.

Notes:
- The harness should treat the DIAGNO binary path as an environment variable
  (`TIAGO_XDIAGNO`) so cluster installs can reuse system builds.
- Comparison logic lives in a small Python helper to keep Fortran free from
  file-format concerns.

## Stage 4 – VMEC / virtual-casing support
Deliverables:
1. Wrapper module `src/solver/vmec_forward.f90` that calls libneo's VMEC spline
   routines (`spline_vmec_data.f90:1-120`) and exposes a `field_t` compatible
   with the Stage 2 diagnostic integrators.
2. Optional fallback to STELLOPT's virtual-casing kernels via ISO_C_BINDING for
   edge cases (the Fortran-to-Fortran interface targets
   `diagno_init_vmec.f90:71-197`).
3. Extended tests that run both Tiago and `xdiagno` with plasma response on a
   truncated W7-X `wout` file (kept small to satisfy runtime constraints).

Notes:
- Keep VC setup cached: once a `wout` is parsed, store the synthesized surface
  and share it across diagnostics; avoid recomputing per run.
- Use the same metadata structure from Stage 1 to bind diagnostics to VMEC
  surfaces so CLI UX stays uniform.

## Stage 5 – Performance + packaging polish
Deliverables:
1. Documented API + CLI reference.
2. Optional HDF5 output containing both Tiago and DIAGNO signals so downstream
   tools (SURROBIER, MaxEnt notebooks) can consume a single artifact.
3. Benchmarks on W7-X trimmed coil sets demonstrating speed-up relative to
   DIAGNO (targeting >10× for vacuum-only workflows).

