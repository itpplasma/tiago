# Tiago User Guide

This document covers everyday workflows for the Toolkit for Inference and
Analysis of Generalized Observables (Tiago).

## 1. Build & test
```
# configure + build
cmake -S . -B build \
      -DTIAGO_LIBNEO_TAG=main \
      -DTIAGO_LIBNEO_GIT=git@github.com:itpplasma/libneo.git
cmake --build build

# convenience targets
make              # wraps the two commands above
make test         # rebuilds and runs ctest --output-on-failure
```
`ctest` drives the following suites:
- `tiago_diag_lint_flux` / `tiago_diag_lint_segrog` / `tiago_diag_lint_invalid`
  ensure DIAGNO-format inputs parse cleanly.
- `tiago_vacuum_smoke` exercises the libneo Biot–Savart solver on the shipped
  sample coil/detector trio.
- `tiago_vs_xdiagno` runs Tiago and STELLOPT's `xdiagno`, compares results, and
  emits CSV + PNG diagnostics under `build/tests/output/`.

## 2. Vacuum solver CLI
`tiago_vacuum_cli` bridges libneo with DIAGNO ingestion:
```
./build/tiago_vacuum_cli COILS.fl DIAGNO.flux DIAGNO.segrog \
    --output-dir build/tests/output \
    --flux-out my_flux.csv \
    --segrog-out my_segrog.csv \
    --seg-area 3.40e-4 \
    --samples 8
```
Key flags:
- `--output-dir`: destination for CSV artifacts (defaults to `.`).
- `--flux-out` / `--segrog-out`: override filenames.
- `--samples`: quadrature sub-sampling per coil segment.
- `--seg-area`: effective area in m² for segmented Rogowski traces when the
  raw DIAGNO file omits metadata.

The CLI reuses `docs/diagnostics/registry.json`. Set
`TIAGO_DIAG_METADATA=/path/to/registry.json` to avoid passing
`--metadata` repeatedly to the lint command.

## 3. Cross-code harness (`scripts/run_xdiagno.py`)
The harness orchestrates Tiago and STELLOPT's `xdiagno`:
```
python3 scripts/run_xdiagno.py \
    --coil tests/data/coil_sample.neo \
    --flux tests/data/fluxloop_sample.diagno \
    --segrog tests/data/segrog_sample.diagno \
    --output build/tests/output \
    --tiago-bin ./build/tiago_vacuum_cli \
    --seg-area 3.40e-4
```
Behavior:
1. Runs Tiago's CLI to produce `tiago_flux.csv` / `tiago_segrog.csv`.
2. Auto-generates `diagno.control`, `input.` (VMEC stub), DIAGNO-style coils,
   and segmented Rogowski files inside `build/tests/output/`.
3. Invokes `xdiagno` (from `TIAGO_XDIAGNO` or `PATH`).
4. Converts `diagno_flux.*` and `diagno_seg.*` into CSVs.
5. Writes diff reports + PNG overlays:
   - `flux_diff.csv`
   - `segrog_diff.csv`
   - `diagnostics.png` (absolute traces + relative-error subplot + runtime info)

The PNGs stay under `build/tests/output/` for artifact-safe CI collection. They
plot Tiago (teal) vs. DIAGNO (orange) traces so discrepancies are visible at a
glance.

## 4. Sample data
- `tests/data/coil_sample.neo`: toy rectangular filament, used for smoke tests
  and harness conversions.
- `tests/data/fluxloop_sample.diagno`: classic DIAGNO flux loop with seven
  points and metadata flags.
- `tests/data/segrog_sample.diagno`: segmented Rogowski example; Tiago’s
  harness populates effective areas automatically via `--seg-area`.
- `docs/diagnostics/registry.json`: metadata for W7-X-esque diagnostics.

## 5. Troubleshooting
| Symptom | Resolution |
| --- | --- |
| `xdiagno not set; skipping cross-code comparison` | Ensure `xdiagno` exists on `PATH` or export `TIAGO_XDIAGNO=/abs/path/to/xdiagno`. |
| `DIAGNO COULD NOT OPEN A FILE` | Run the harness from the repository root so the generated `diagno.control`, `input.`, and `coils.tiago` remain visible to `xdiagno`. |
| `run_xdiagno failed: ... returned non-zero` | Inspect `flux_diff.csv` / `segrog_diff.csv` for the offending labels, and use the PNG overlays to confirm magnitude/sign of the mismatch. |

## 6. Artifact locations
After `ctest` (or `make test`) the following files live in `build/tests/output/`:
- `tiago_flux.csv`, `tiago_segrog.csv`
- `diagno_flux.csv`, `diagno_segrog.csv`
- `flux_diff.csv`, `segrog_diff.csv`
- `diagnostics.png`
- Harness control files (`diagno.control`, `input.`, `coils.tiago`, `segrog.diagno`, `fluxloop.diagno`).

These artifacts provide both numerical and visual evidence that Tiago and
DIAGNO agree within the configured tolerance (default `1e-3`). Adjust
`--tolerance` if your diagnostics require tighter or looser comparisons.
