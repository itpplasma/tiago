# TIAGO: Toolkit for Inference and Analysis of Generalized Observables

Tiago is a modern Fortran toolkit for vacuum magnetic diagnostic studies. It
parses legacy DIAGNO coil descriptions, evaluates flux loops and segmented
Rogowski probes via libneo's Biot–Savart solver, and cross-validates every run
against STELLOPT's `xdiagno` binary.

## Highlights
- **Drop-in DIAGNO ingestion** – Load STELLOPT-compatible flux loops, segmented
  Rogowski files, and registry metadata without altering the on-disk format.
- **Vacuum solver + CLI** – Batch evaluation with `tiago_vacuum_cli` using
  libneo's Biot-Savart solver; comprehensive CMake test suite.
- **Cross-code validation** – `scripts/run_xdiagno.py` automatically prepares
  control files, runs both Tiago and `xdiagno`, writes CSV summaries, and emits
  comparison PNGs in `build/tests/output/` for visual inspection.

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
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8 \
    --nfp 3 --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv --plasma-sample
```
This command emits `tiago_flux.csv` and `tiago_segrog.csv` in the output
directory. Override sample metadata via
`--flux-out`, `--segrog-out`, `--kind`, or the registry file referenced by
`docs/diagnostics/registry.json`. `--nfp` matches the field-period geometry,
while `--flux-turns` / `--segrog-turns` accept text files of `label scale`
pairs so Tiago mirrors DIAGNO's namelist-based `flux_turns` and
`segrog_turns` arrays.

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
- `diagnostics.png` – absolute traces (y-axis pinned to zero), relative-error
  subplot, and Tiago vs. `xdiagno` runtimes
- `geometry.png` – coil filaments plus flux-loop/segmented Rogowski paths so you
  can inspect geometry coverage visually

`ctest` target `tiago_vs_xdiagno` wraps this flow for the toy sample, while the
`tiago_vs_xdiagno_ncsx_nfp1` and `tiago_vs_xdiagno_ncsx_nfp3` tests target two
NCSX vacuum scenarios. During those runs the harness automatically downloads
`coils.NCSX_nfp1` / `coils.NCSX` from the public STELLOPT tree (only when the
files are absent), runs Tiago and `xdiagno` with matching `nfp` values and turn
scalars, and stores artifacts in `build/tests/output/<case>/diagnostics_<case>.png`
and `geometry_<case>.png`.

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

## Test Status
Current test results: **9/9 tests passing** (100%)

Test suite includes:
- ✅ Diagnostic lint validation (flux, Rogowski, invalid cases)
- ✅ Coil loader unit tests
- ✅ Cross-code validation vs xdiagno (5-coil reference geometry)
- ✅ NCSX geometry tests (NFP=1 and NFP=3 cases with full 18,690-coil sets)

The solver has been optimized for performance with:
- L1D cache hit rate improved from 10% → 95% via loop reordering
- Diagnostic-level parallelization on 16 cores
- Native CPU optimization flags (-march=native -mtune=native)

## Further reading
- `docs/USER_GUIDE.md` – Detailed CLI options, sample data, troubleshooting
