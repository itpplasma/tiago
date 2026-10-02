# TIAGO comparison measurements, October 2, 2026

Fresh public STELLOPT `2f181f0d4e71f28afe076967d606b2c020c9db33` and TIAGO
executions used GNU Fortran 16.2.1 on a shared Ryzen 9 5950X Linux workstation.
The JSON files retain every measured case; provenance includes exact input and
binary hashes. Wall times include startup, input loading and competing
workstation activity. Each row is one observation, so speed ratios describe
these runs rather than a statistical confidence interval.

## Accuracy and reference controls

`reference-default/` contains all six benchmark suites. Its reference build
uses STELLOPT's normal `-O2 -march=native -fcheck=all,no-array-temps` settings.
TIAGO uses its portable release build. The four documented reference patches
are retained only in designated controls; stock-reference defects have upstream
reports linked from the repository's issue tracker.

All 12 diagnostic-format controls agree. At 16 midpoint samples, realistic
NCSX and M16N08 geometry comparisons have maximum relative signal differences
of 2.21e-7 and 1.52e-7, respectively, with no missing or nonfinite signals.
The magnetic-probe maximum difference is 5.90e-7. The 350-entry coil response
matrix has median difference 2.34e-8 and maximum difference 4.26e-4, which must
remain visible alongside the median. The response-matrix signal reconstruction
agrees within 1.26e-14.

The plasma comparison uses different quadrature methods: DIAGNO's adaptive
virtual casing and TIAGO's fixed grids. At the finer TIAGO grid, maximum flux
and probe differences are approximately 4.80e-4 and 4.35e-4. Their wall-time
ratio is therefore not an equal-accuracy performance claim. Eighty boundary
field finite differences agree within 1.12e-11. VMEC2000 parameter-chain checks
agree within 3.88e-8, but these checks do not validate the VMEC++ implicit
adjoint. The default DIAGNO-compatible sheet retains spurious Ampere-current
derivatives; issue #25 tracks the alternative current formulation.

## Matched optimization settings

`reference-optimized/` reruns repository and realistic geometry cases with
portable `-O3 -fno-math-errno` on both implementations and reference runtime
checks disabled. Additional reference compatibility/BLAS flags are recorded
in its provenance. The optimized reference retains the same accuracy results.

| case, 16 samples | serial reference / TIAGO | four-process/thread reference / TIAGO |
| --- | ---: | ---: |
| NCSX repository, one field period | 2.09 | 5.31 |
| NCSX repository, three field periods | 3.39 | 6.94 |
| NCSX realistic geometry | 0.76 | 1.15 |
| M16N08 realistic geometry | 1.15 | 1.55 |

A ratio below one means TIAGO took longer. Small toy cases are dominated by
reference startup cost and cannot establish a kernel speed advantage. These
measurements do not support superiority on every workload.

Reproduce with `build_xdiagno.sh`, `fetch_data.sh`, and `bench.py`, following
the parent README. Save the result and declared build flags with
`snapshot_results.py`. Use a platform-specific STELLOPT make configuration
for the optimized comparison; the reference's normal runtime-checked build
should remain a separate correctness control.
