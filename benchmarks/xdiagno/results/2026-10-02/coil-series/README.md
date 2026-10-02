The coil potential kernel uses a five-term odd series for `abs(eps) <= 0.01`,
where `eps = L/(Ri+Rf)`. The relative truncation bound is
`eps^10/(11*(1-eps^2)) <= 9.1e-22`; larger ratios keep the existing log formula.
This avoids a log call for sufficiently short/distant segments and avoids the
loss of precision in `(1+eps)/(1-eps)` for distant points. The field kernel and
public mutable coil arrays keep their existing behavior.

Two rotated trials compared the fixed inactive-jumper baseline `3161170` and
the series candidate, with GNU 16.2.1, portable `-O3` and `-fno-math-errno` on
both kernels. M16N08 geometry and samples were identical. CPU affinity was set
before `fo`, with `OMP_PROC_BIND` and `OMP_PLACES` removed: four distinct worker
IDs each retained all four permitted CPUs. `results.json` records times,
paired ratios, worker affinity, and diagnostic errors; `provenance.json`
records source/data hashes, compiler, flags, and machine details.

| Measurement | One worker, paired speedup | Four workers, paired speedup |
| --- | --- | --- |
| Potential kernel | 1.537 / 1.448 | 1.467 / 1.454 |
| Full vacuum CLI, 16 samples | 1.162 / 1.330 | 1.494 / 1.687 |

The kernel compared 4096 potential points and 1024 field points in one binary,
with baseline/alternative order rotated. The field is bit-identical; the
potential change is at most `2.96e-15` relative to the largest baseline value.
The complete CLI includes `fo exec` launch/build checks, I/O and both diagnostic
integrations. All 277 CLI values are finite, none are missing, and maximum
relative DIAGNO error remains `1.53e-7`. DIAGNO timings here have different
reference affinity and support numerical comparison only. The existing
matched reference benchmark is separate.

These are two observations per setting on a shared host, not a universal
speed claim. Earlier exploratory four-thread timings that passed
`OMP_PROC_BIND=close` into `fo` were discarded: its child inherited one CPU.

`tests/test_coil_potential.f90` compares with the independent straight-wire
asinh integral at three length scales, signed currents, repeated nodes,
public-state mutations, `rho/L` through `1e12`, and both sides of `eps=0.01`.
The distant-point regression fails with the fixed baseline and passes with
the series. `series-certificate.json` records 80-digit controls for 21 far,
scaled and threshold cases: maximum candidate relative error `4.76e-15`,
versus `8.89e-5` for the old far-field calculation. Independent analytic field
and central potential-derivative checks also passed at the threshold/scales.

To repeat the CLI pairs, build both checkouts using `fo`, prepare the existing
`benchmarks/xdiagno/bench.py` geometry inputs and portable DIAGNO reference,
and run `repeat_cli.py` with the two checkout paths, data root and local CPU
masks. Its command uses `fo exec`, keeps exactly two trial orders, and records
raw affinities. Full pre-commit `OMP_NUM_THREADS=4 fo` passed all configured
checks with reconstruction disabled; root verifies the reconstruction build.
