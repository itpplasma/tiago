# Toolkit for Inference and Analysis of Generalized Observables (Tiago)

Tiago is a modern Fortran toolkit for vacuum magnetic diagnostic studies. It
parses legacy DIAGNO coil descriptions, evaluates flux loops and segmented
Rogowski probes via libneo's Biot–Savart solver, and cross-validates every run
against STELLOPT's `xdiagno` binary.

## Highlights
- **Drop-in DIAGNO ingestion** – Stage 1 modules load STELLOPT-compatible flux
  loops, segmented Rogowski files, and registry metadata without altering the
  on-disk format.
- **Vacuum solver + CLI** – Stage 2 adds `tiago_vacuum_cli` for batch
  evaluation and `tiago_vacuum_smoke` for deterministic regression coverage.
- **Cross-code proof** – Stage 3 provides `scripts/run_xdiagno.py`, which now
  prepares control/input files automatically, runs both Tiago and `xdiagno`,
  writes CSV summaries, and emits comparison PNGs in `build/tests/output/` for
  visual inspection.

## Quick start
```bash
# Configure + build
cmake -S . -B build \
      -DTIAGO_LIBNEO_TAG=main \
      -DTIAGO_LIBNEO_GIT=git@github.com:itpplasma/libneo.git
cmake --build build

# or the convenience wrappers
make              # configures + builds
make test         # rebuilds and runs ctest --output-on-failure
```
Pass extra cache entries through `CMAKE_ARGS`, e.g.
`make CMAKE_ARGS="-DTIAGO_LIBNEO_TAG=dev-feature"`.

## Vacuum CLI usage
```
./build/tiago_vacuum_cli \
    tests/data/coil_sample.neo \
    tests/data/fluxloop_sample.diagno \
    tests/data/segrog_sample.diagno \
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8
```
This command emits `tiago_flux.csv` and `tiago_segrog.csv` in the output
directory. Override sample metadata via
`--flux-out`, `--segrog-out`, `--kind`, or the registry file referenced by
`docs/diagnostics/registry.json`.

## Cross-code validation & visual artifacts
```
python3 scripts/run_xdiagno.py \
    --coil tests/data/coil_sample.neo \
    --flux tests/data/fluxloop_sample.diagno \
    --segrog tests/data/segrog_sample.diagno \
    --output build/tests/output \
    --tiago-bin ./build/tiago_vacuum_cli \
    --seg-area 3.40e-4
```
When `TIAGO_XDIAGNO` is unset, the harness looks for `xdiagno` on `PATH`. The
script writes:
- `tiago_flux.csv`, `tiago_segrog.csv`
- `diagno_flux.csv`, `diagno_segrog.csv`
- `flux_diff.csv`, `segrog_diff.csv`
- `diagnostics.png` – a single figure with absolute traces (y-axis pinned to
  zero), relative-error subplot, and Tiago vs. `xdiagno` runtimes. Everything
  stays inside `build/tests/output/`.

`ctest` target `tiago_vs_xdiagno` wraps this flow so every test run produces the
PNG evidence automatically.

## Directory layout
```
cmake/                     # Toolchain helpers
CMakeLists.txt             # FetchContent bridge + targets
scripts/run_xdiagno.py     # Cross-code harness + PNG generator
src/                       # Fortran diagnostics, solver, CLI
tests/                     # Sample inputs + regression drivers
docs/diagnostics/          # Diagnostic registry metadata
docs/USER_GUIDE.md         # In-depth user documentation
```

## Further reading
See `docs/USER_GUIDE.md` for detailed CLI options, sample data notes, test
artifacts, and troubleshooting tips.
