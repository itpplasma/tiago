# Dispersion interferometer timing

Build `benchmark_interferometer` in an OpenMP-enabled Release build, then run
with `OMP_NUM_THREADS=1` and `4` (and `OPENBLAS_NUM_THREADS=1`). The benchmark
reports wall seconds per prediction, Jacobian-vector product (JVP) and
transposed Jacobian-vector product (VJP) for 1024 synthetic chords,
196608 quadrature points, 100 calls per observation and three observations
after warm-up. Mode numbers 1, 2 and 3 mean prediction, JVP and VJP.

The checksum column checks unchanged results; parallel calibration-adjoint
reductions may differ by roundoff. The workload includes both quadrature
rules and the recorded response, excludes provider sampling and initialization,
and does not measure a complete reconstruction. Retain host load, compiler and
all repetitions when interpreting the results.
