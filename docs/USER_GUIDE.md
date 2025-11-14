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
  emits CSV + PNG diagnostics under `build/tests/output/sample/`.
- `tiago_vs_xdiagno_ncsx_nfp1` / `_ncsx_nfp3` cover the two NCSX benchmark
  cases, automatically downloading the coil files (if needed), applying
  diagnostic turn scaling, and publishing artifacts under
  `build/tests/output/<case>/`.

## 2. Vacuum solver CLI
`tiago_vacuum_cli` bridges libneo with DIAGNO ingestion:
```
./build/tiago_vacuum_cli COILS.fl DIAGNO.flux DIAGNO.segrog \
    --output-dir build/tests/output \
    --flux-out my_flux.csv \
    --segrog-out my_segrog.csv \
    --seg-area 3.40e-4 \
    --samples 8 \
    --nfp 3 \
    --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv \
    --plasma-sample
```
Key flags:
- `--output-dir`: destination for CSV artifacts (defaults to `.`).
- `--flux-out` / `--segrog-out`: override filenames.
- `--samples`: quadrature sub-sampling per coil segment.
- `--seg-area`: effective area in m² for segmented Rogowski traces when the
  raw DIAGNO file omits metadata.
- `--nfp`: field periods for replicated diagnostics (matches DIAGNO's
  `nfp_diagno`).
- `--flux-turns` / `--segrog-turns`: optional text files (one `label value`
  pair per line, `#` comments allowed) that apply DIAGNO's turn scaling to the
  Tiago outputs before CSV emission.
- `--plasma-wout`: VMEC equilibrium (wout) file. When present, segmented
  Rogowski diagnostics include plasma-response `B_external` sampled from the
  VMEC surface via the virtual-casing solver. Flux loops remain vacuum-only
  because no vector potential is available from the plasma solver.
- `--plasma-sample`: download (if needed) and reuse the lightweight Simsopt
  VMEC reference bundled with the regression tests. This flag is ignored when
  `--plasma-wout` is provided explicitly.
- `--plasma-nphi` / `--plasma-ntheta`: resolution of the VMEC surface grid fed
  to virtual casing (default 16×16). Increase these for higher-accuracy plasma
  response at the cost of setup time.

The CLI reuses `docs/diagnostics/registry.json`. Set
`TIAGO_DIAG_METADATA=/path/to/registry.json` to avoid passing
`--metadata` repeatedly to the lint command.

## 3. Cross-code harness (`scripts/run_xdiagno.py`)
The harness orchestrates Tiago and STELLOPT's `xdiagno`:
```
python3 scripts/run_xdiagno.py \
    --coil tests/cases/ncsx_nfp3/coils.NCSX \
    --coil-url https://raw.githubusercontent.com/PrincetonUniversity/STELLOPT/develop/BENCHMARKS/FIELDLINES_TEST/coils.NCSX \
    --flux tests/cases/ncsx_nfp3/fluxloop.diagno \
    --segrog tests/cases/ncsx_nfp3/segrog.diagno \
    --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv \
    --output build/tests/output/ncsx_nfp3 \
    --label ncsx_nfp3 \
    --tiago-bin ./build/tiago_vacuum_cli \
    --seg-area 3.40e-4 \
    --samples 14 \
    --nfp 3
```
Behavior:
1. Ensures the requested coil file exists (downloading via `--coil-url` if
   necessary) and runs Tiago's CLI to produce `tiago_flux.csv` / `tiago_segrog.csv`.
2. Auto-generates `diagno.control` (with matching `flux_turns`,
   `segrog_turns`, and `nfp`), `input.` (VMEC stub), DIAGNO-style coils,
   and segmented Rogowski files inside `build/tests/output/`.
3. Invokes `xdiagno` (from `TIAGO_XDIAGNO` or `PATH`).
4. Converts `diagno_flux.*` and `diagno_seg.*` into CSVs.
5. Writes diff reports + PNG overlays:
- `flux_diff.csv`
- `segrog_diff.csv`
- `diagnostics_<label>.png` (absolute traces + relative-error subplot + runtime info)
- `geometry_<label>.png` (3D plot of coils and diagnostic curves)

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
After `ctest` (or `make test`) each case stores artifacts under
`build/tests/output/<label>/`:
- `tiago_flux.csv`, `tiago_segrog.csv`
- `diagno_flux.csv`, `diagno_segrog.csv`
- `flux_diff.csv`, `segrog_diff.csv`
- `diagnostics_<label>.png`
- `geometry_<label>.png`
- Harness control files (`diagno.control`, `input.`, `coils.tiago`, `segrog.diagno`, `fluxloop.diagno`).

These artifacts provide both numerical and visual evidence that Tiago and
DIAGNO agree within the configured tolerance (default `1e-3`). Adjust
`--tolerance` if your diagnostics require tighter or looser comparisons.
