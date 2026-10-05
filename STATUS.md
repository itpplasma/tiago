# TIAGO delivery status

Updated October 2, 2026. Development stopped at the user's request; completed
work is committed and pushed to `main`. The tested production snapshot is
[`e14b45a`](https://github.com/itpplasma/tiago/commit/e14b45a0b0877163175b6180e01dfe47069576da).
This status file adds the handoff record without changing production code.

## Completed work

The original open PR and fetched branch histories are integrated into main.
[PR #2](https://github.com/itpplasma/tiago/pull/2) is merged. The integration
includes explicit raw-covariant and conservative projected plasma-boundary
models, their field/shape responses, reconstruction cache invalidation, and
mixed LLVM/GNU OpenMP runtime ordering repairs. DIAGNO-compatible contravariant
evaluation remains the default. A covariant input alone does not guarantee
conservation; the opt-in projection changes the modeled boundary field.

Inactive coil segments are skipped before singular logarithms. A bounded
five-term potential series improves distant-point accuracy and avoids a
logarithm for `abs(eps) <= 0.01`; larger arguments retain the existing formula.
Independent straight-wire integrals, scaling, current mutation, threshold and
derivative controls cover the change. The distant-point certificate reports
maximum relative error `4.76e-15`, compared with `8.89e-5` for the old formula.

Fresh benchmark runners and provenance retain all six DIAGNO comparison suites,
matched portable optimization settings, input/binary hashes, and controlled
coil-series repetitions. Reference defects were reported upstream:
[STELLOPT #506](https://github.com/PrincetonUniversity/STELLOPT/issues/506),
[#507](https://github.com/PrincetonUniversity/STELLOPT/issues/507),
[#508](https://github.com/PrincetonUniversity/STELLOPT/issues/508),
[#509](https://github.com/PrincetonUniversity/STELLOPT/issues/509), and
[VMEC++ #954](https://github.com/proximafusion/vmecpp/issues/954).
The corresponding TIAGO reporting tasks are closed.

## Verification and measurements

The final complete `fo` pipeline passes all 31 configured tests with
reconstruction enabled, LLVM/Clang 22, Enzyme 0.0.293, four OpenMP threads and
one BLAS thread. Earlier complete LLVM 20 and 22 serial/parallel pipelines
also pass. Introduced array-temporary diagnostics were removed; pre-existing
vendored-source parser hints and warnings are separate from these corrections.
No hosted TIAGO workflow existed at the audit, so local results are not a
claim of GitHub Actions execution.

The [DIAGNO measurement record](benchmarks/xdiagno/results/2026-10-02/README.md)
reports agreement for all twelve diagnostic-format controls. Realistic NCSX
and M16 geometry signal differences are at most `2.21e-7` and `1.52e-7`.
The 350-entry response matrix has median difference `2.34e-8` and maximum
`4.26e-4`; the maximum is retained rather than hidden by the median.
The plasma comparison uses different quadrature accuracy and does not support
an equal-accuracy speed ratio.

Matched portable `-O3` whole-process reference/TIAGO ratios range from 0.76
for serial realistic NCSX geometry to 6.94 for a four-worker repository case.
A ratio below one means TIAGO was slower. The
[coil-series record](benchmarks/xdiagno/results/2026-10-02/coil-series/README.md)
reports two paired observations per setting: potential-kernel speedups around
1.45–1.54 and complete vacuum-CLI speedups of 1.16–1.33 with one worker and
1.49–1.69 with four. All 277 CLI values are finite; maximum relative DIAGNO
error remains `1.53e-7`. These shared-host observations establish bounded
improvements, not superiority over every competitor or workload.

## Remaining work

[Issue #25](https://github.com/itpplasma/tiago/issues/25) remains open for the
default-policy and scientific-model decision. Explicit covariant/conservative
modes are implemented and their finite-difference responses are validated,
but the original universal raw-covariance conservation premise fails on the
pinned LI383 input. The measured Fourier curl is `0.001337677 T m`, with a
current-ripple bound of about `250.747 A`. The conservative projection closes
the modeled boundary one-form and preserves constant circulation coefficients;
it does not prove a globally force-balanced conservative physical field.
The [integration comment](https://github.com/itpplasma/tiago/issues/25#issuecomment-5956424158)
records this resolution and its limits.

The standalone VMEC parameter-chain benchmark does not by itself validate the
VMEC++ implicit adjoint. Reconstruction tests compare the complete TIAGO chain,
including that adjoint, with fully re-solved finite differences for all seven
parameter kinds on pinned LI383 input; this does not establish general
correctness. Equal-accuracy plasma performance comparisons, broader convergence
studies, and statistically repeated competitor measurements remain future work.
Upstream reports also remain subject to their maintainers' decisions.
No claim of complete scientific certification or universal superiority is made.

To reproduce the reconstruction pipeline, use the compatible LLVM/Enzyme
installation described in [README.md](README.md):

```sh
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 FO_JOBS=8 \
  FO_CMAKE_ARGS='-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DTIAGO_ENABLE_RECONSTRUCTION=ON -DTIAGO_ENZYME_PLUGIN=/usr/lib/ClangEnzyme-22.so' \
  fo
```

The benchmark READMEs provide their separate pinned reference builds and
reproduction commands. Further development is deferred after this handoff.


## Architecture decision, October 5, 2026

Documentation now records the intended long-term boundary with KIN6D. This
does not change the delivered production code or invalidate the October 2
verification record above.

- TIAGO owns diagnostic observation models, reconstruction/inference, priors
  and nuisance/calibration/measurement UQ.
- KIN6D owns general forward plasma physics, stationary/evolution solves,
  differentiability, numerical certification and model-reduction error.
- The existing VMEC++ exact-adjoint reconstruction remains a first-class
  provider/reference; it is not replaced merely because KIN6D is planned.
- A provider-neutral seam is introduced only when the first KIN6D consumer
  requires it, avoiding speculative abstraction churn.
- Several distinct equilibria remain several physical solutions. TIAGO may fit
  each locally and optionally report a finite mixture of Laplace/Gaussian
  approximations; it need not implement a general branch framework.
- The present covariance is a local Gaussian approximation. Correlated
  measurement covariance and nuisance/calibration parameters are higher
  priority than general posterior sampling.

See DESIGN.md and ROADMAP.md for the planned architecture.
