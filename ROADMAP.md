# TIAGO roadmap

Status: staged architecture plan, 2026-10-05. The October 2 VMEC++/magnetic
delivery remains the implemented baseline. This roadmap does not imply that
later items already exist.

## Rules

- Preserve the working VMEC++/DIAGNO-compatible path as a permanent reference.
- Extract abstractions only when a concrete second provider/diagnostic requires
  them.
- Keep exact/adjoint derivatives and finite-difference/pull-study oracles.
- Prefer matrix-free JVP/VJP composition over explicit full sensitivity
  matrices for large KIN6D states.
- Keep inference UQ distinct from forward numerical/model error.
- Do not introduce Python as a production dependency merely for inference.

## T0 — retain the delivered magnetic reconstruction

Keep current vacuum/plasma diagnostics, response matrices, VMEC++ implicit
adjoint, LM solver, local covariance, xdiagno comparison and LI383/NCSX
benchmarks green. Resolve scientific defects in that path on their own merits;
do not block maintenance on future KIN6D work.

## T1 — first provider seam, only with a real KIN6D consumer

Use the first native KIN6D Grad--Shafranov reconstruction as the trigger.
Factor the smallest interface that both VMEC++ and KIN6D need: solve/update,
provider state/observable view, parameter metadata, JVP/VJP/adjoint actions,
failure status and provenance.

Do not first design a universal plugin hierarchy for every future plasma model.
The existing VMEC++ implementation may remain internally specialized behind the
new seam.

Acceptance: the current VMEC benchmark is behaviorally unchanged and one
KIN6D Grad--Shafranov case reconstructs through the same TIAGO objective.

## T2 — proper local UQ

Retain Laplace/Gauss--Newton as the default, but generalize the residual model
from independent scalar sigmas to:

- correlated measurement covariance/factorization;
- shared diagnostic gains/offsets;
- sensor geometry/alignment nuisance parameters where differentiable;
- coil-current/calibration nuisance parameters;
- correlated or structured priors.

Add synthetic pull/coverage tests for each new uncertainty mechanism.

## T3 — integrated KIN6D equilibrium reconstruction

After native Grad--Shafranov, extend the provider contract only as needed for:

1. free-boundary axisymmetric MHD with explicit coil/vacuum data;
2. imported EQDSK checked/re-solved by KIN6D;
3. general 3-D MHD stationary states;
4. nested-surface VMEC/DESC/GVEC states as retained comparators;
5. later kinetic-MHD/DK/GK/FK stationary states once relevant diagnostics and
   derivatives exist.

TIAGO never assumes nested flux surfaces merely because one provider uses them.

## T4 — several equilibria, finite local mixtures

When a scientific case supplies several distinct equilibria, run independent
local reconstructions for each. Report each MAP/covariance and fit diagnostics
separately.

Optionally compute consistent local evidence/Laplace weights and expose a
finite Gaussian mixture over the local parameter posteriors. Do not average the
physical equilibria. Do not add a general branch-discovery engine to TIAGO.

## T5 — diagnostic breadth

Add diagnostics only with concrete experimental consumers. Magnetic diagnostics
remain the baseline. Future examples may include interferometry/polarimetry,
profile diagnostics, radiation/spectroscopy, ECE/CXRS, neutron/fusion-product
and fast-particle diagnostics, with their instrument/atomic/radiative-transfer
uncertainties owned explicitly.

Reuse KIN6D physical forward primitives when appropriate while keeping
instrument/calibration semantics in TIAGO.

## T6 — advanced posterior methods only on evidence

Do not schedule HMC/SMC/MCMC simply for completeness. First demonstrate a case
where the within-equilibrium posterior is materially non-Gaussian or the finite
mixture of local Laplace approximations fails a coverage/predictive check.
Then choose the cheapest gradient-aware method that addresses that failure.

## Permanent outputs

Every reconstruction records:
- measurement data and covariance model;
- nuisance/calibration parameterization;
- priors;
- provider identity/revision/configuration;
- fitted state/equilibrium identifier;
- optimizer history/status;
- Jacobian/adjoint verification status;
- local covariance/UQ method;
- any mixture weights and their approximation;
- forward numerical/model error components without silently combining them
  with posterior uncertainty.

Chris&AI
