# Testing & Validation Strategy

All tests use `ctest` and run out of the CMake build tree.

## Data dependencies

| Dataset | Purpose | Source (planned) |
| --- | --- | --- |
| `tiago-testdata/w7x_trimmed_coils.json` | Small W7-X-like coil/diagnostic bundle for vacuum regression tests | Will live in a dedicated public repo (`https://github.com/itpplasma/tiago-testdata`) so the harness can download it with `cmake -P cmake/DownloadTestData.cmake` |
| `eval_2000/wout_vmec_aux.nc` | MaxEnt VMEC snapshot used for parity checks with `xdiagno` | Already published in `/home/ert/data/PAPERS/MAXENT/Koeberl_MaxEnt_2023_data/data/eval_2000`; we will mirror a truncated version into `tiago-testdata` to keep CI runtimes low |

The download helper will
1. Create `build/test-data/`.
2. Fetch each asset via `file(DOWNLOAD ...)`.
3. Verify SHA256 sums before handing files to tests.

## Test matrix

| Test | Stage | Driver | Notes |
| --- | --- | --- | --- |
| `tiago_vacuum_smoke` | Stage 2 | Fortran unit test | Loads the trimmed coil file, solves in vacuum, asserts sign conventions |
| `tiago_vs_xdiagno` | Stage 3 | Python helper | Runs Tiago and `xdiagno -vac` on the same inputs; diff `diagno_flux.*` and `diagno_seg.*` outputs |
| `tiago_vmec_surface` | Stage 4 | Fortran integration test | Exercises the VMEC spline-backed field class on the truncated MaxEnt `wout` |

`tiago_vs_xdiagno` requires the user to expose the legacy binary via
`TIAGO_XDIAGNO=/path/to/xdiagno`.  The script will skip the test (marking it as
`CTEST_CUSTOM_TESTS_IGNORE`) when the variable is unset so developers without
STELLOPT builds can still run the suite.

All tests produce machine-readable CSV summaries that the CI workflow stores as
artifacts.  This mirrors the MaxEnt diagnostic comparison done in
`/home/ert/data/W7X/DIAGNO/scripts/postprocess_diagnostics.py`.
