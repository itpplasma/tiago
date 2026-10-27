# Joint equilibrium and profile inference

Status: in development on branch `inference/core-20261026`. The VMEC++
magnetic reconstruction (`tiago_reconstruct`, README and STATUS.md) is
unchanged and remains the benchmarked reference path. This document is the
contract for the inference core and for the profile diagnostics built on it.

## Generative model

The latent parameters are defined on the input side and never inverted from
outputs:

- kinetic profiles ne, Te, Ti, Zeff and a fast-particle pressure nuisance
  p_fast as functions of the provider's flux label s;
- the toroidal current (CURTOR and the AC coefficients of VMEC's current
  profile), the boundary and PHIEDGE parameters of the existing chain;
- coil-group currents (EXTCUR) with Gaussian priors from their measurements;
- diagnostic alignment nuisances (owned by each diagnostic) and calibration
  gains/offsets (owned by the likelihood).

For given latents the provider solves the equilibrium with the derived
pressure p = ne Te + (n_D + n_Z) Ti + p_fast (main ions D and one impurity of
charge Z_imp fixed by Zeff), all signals are predicted in real space, and the
likelihood compares them with the data. There is no outer fixed point between
profile fitting and equilibrium fitting.

## Profile-diagnostic interface

Module `tiago_profile_diagnostic` (`src/inference/profile_diagnostic.f90`,
library `tiago_inference`, built without VMEC++) defines the contract. A
diagnostic extends `profile_diagnostic_t` and implements:

| procedure | meaning |
|---|---|
| `n_signals()` | number of predicted signals |
| `n_nuisance()` | number of its own alignment parameters eta |
| `label(i)` | name of signal i |
| `sample_points(eta, x, dx_deta)` | real-space points x(3, n_points) as (R [m], phi [rad], Z [m]) and dx/deta (3, n_points, n_nuisance) |
| `predict(sample, eta, y)` | signals from the sampled quantities |
| `jvp(sample, eta, q_dot, eta_dot, y_dot)` | tangent at fixed points |
| `vjp(sample, eta, y_bar, q_bar, eta_bar)` | adjoint at fixed points |

`plasma_sample_t` holds, for each point, the quantities
`q(N_QUANTITIES, n)` with rows `Q_S` (normalized toroidal flux s, extrapolated
above 1 outside the plasma), `Q_B` (|B| in T), `Q_NE` (m^-3), `Q_TE` and
`Q_TI` (eV) and `Q_ZEFF`, their spatial derivatives `dq_dx(3, N_QUANTITIES, n)`
with respect to (R, phi, Z), and `inside(n)`. Outside the last closed surface
the kinetic profiles continue as f(1) exp(-(s - 1) / sol_decay) (Zeff
constant, sol_decay = 0.05 by default), so predictions remain continuous and
differentiable when a point crosses the boundary during sampling.

The core evaluates the sample, calls the diagnostic and composes the chain
rule: the diagnostic's q_bar goes through the profiles to the profile latents
and through s and |B| to the provider's implicit adjoint; the motion of the
points adds `sum(q_bar * dq_dx . dx_deta)` to eta_bar. A diagnostic therefore
needs only local derivatives and is testable without an equilibrium:
`check_diagnostic_derivatives` compares `jvp` with central differences of
`predict` and checks the adjoint identity <y_bar, J v> = <J^T y_bar, v>.

Rules for diagnostic authors:

- Instrument semantics (scattering volumes, spectral or antenna response,
  lines of sight and their quadrature, cold resonance location) live in the
  diagnostic. Line integrals choose their own quadrature points; a resonance
  search samples |B| along the line and differentiates the located position
  through q_bar(Q_B).
- The measurement covariance, shared calibration gains/offsets and the
  model-error term belong to the likelihood (`tiago_likelihood`), not to the
  diagnostic.
- Private machine geometry is read at run time; it never enters the repository.

`tiago_point_probe` is the reference implementation (ideal local ne/Te/Ti
values with a rigid (dR, dZ) alignment shift) used by the synthetic tests.

## Provider states and point location

Points are located in a provider state through `tiago_real_space_state`
(`real_space_state_t`: labels u = (rho, theta, phi) to Cartesian x and |B|,
parameter JVP/VJP at fixed labels, Newton location and implicit-function
derivatives; introduced with the synthetic-diagnostics work). There is one
location path. `tiago_state_sampler` fills the rows Q_S and Q_B of a
`plasma_sample_t` and maps their cotangents to the state parameters.

`tiago_vmecpp_state` implements the interface for a solved VMEC++
equilibrium. Its parameters are the solver's geometry coefficients (the
cotangent the implicit adjoint accepts) followed by the |B| spectrum bmnc,
each a C1 Hermite interpolant in rho. The joint model passes the geometry
cotangent to `tiago_equilibrium%jacobian` (`coef_cot`). The bmnc cotangent
has no adjoint to the VMEC++ inputs yet, so |B| is exact in position and
labels but its change with the equilibrium inputs at a fixed label is not
propagated; `joint_model_t%b_cotangent_dropped` records the size of the
dropped cotangent. This matters for diagnostics that read |B| (ECE).

## Likelihood and engines

`tiago_joint_model` stacks the magnetic rows of `tiago_reconstruction`
(whitened with a full covariance plus a model-error variance), the
free-boundary consistency rows and parameter priors, the profile diagnostics
L^-1 (g y - d) with their own full covariance, model-error term, alignment
and gain nuisances, and N(0, I) rows for the latent priors. The kinetic
pressure enters VMEC++ as a line-segment profile with knots at the half-grid
surfaces, so the solver's pressure cotangents apply directly.

`tiago_inference_engines` provides the MAP (Levenberg-Marquardt; VMEC++
hot-restarts from the previous equilibrium), the Laplace covariance and log
evidence, NUTS from fortnum with the Laplace covariance as dense metric, and
Pareto-smoothed importance sampling from the Laplace Gaussian (k-hat).
`tiago_infer` runs the synthetic study (truth drawn from the prior, synthetic
noise, MAP, Laplace, PSIS, NUTS) and writes CSV files for
`tools/inference_plots.py` and `tools/inference_coverage.py`.

## Synthetic coverage study (LI383, fixed boundary)

sCluster array 3369752 ran 50 seeds of `tiago_infer` at 22c28a4 (seeds
101-150; Clang 22.1.8 / Enzyme 0.0.293, 4 threads each, about 10 min per
seed). Each seed draws the truth from the prior (CURTOR, PHIEDGE, 7 B-spline
latents each for ne, Te, Ti, a (dR, dZ) probe alignment: 25 latents), adds
noise to 54 magnetic signals and 24 profile probes (5 %), and runs MAP,
Laplace, PSIS (400 draws) and NUTS (300 warm-up, 400 draws). Fractions of
truths inside the central intervals (nominal 0.683/0.954 for +-1/2 sd,
0.68/0.95 for quantile intervals):

| | 68 % | 95 % |
|---|---|---|
| latents, Laplace | 0.671 | 0.941 |
| latents, NUTS | 0.664 | 0.930 |
| ne, Te, Ti, p at 41 s values, Laplace | 0.669 | 0.938 |
| same, NUTS | 0.673 | 0.933 |
| iota (40 re-solved draws) | 0.565 | 0.904 |

Values within a seed are correlated, so the binomial intervals of
`tools/inference_coverage.py` understate the uncertainty. NUTS had no
divergences, mean acceptance 0.80-0.86, tree depth 3 and E-BFMI 0.9-1.3;
PSIS k-hat had median 0.42 and exceeded 0.7 in one seed. The iota bands come
from only 40 solves per seed and iota is dominated by the single CURTOR
latent, so its coverage is nearly all-or-nothing per seed; whether the
under-coverage is sampling noise or a bias has not been established.

## References and reference implementations

- Profile prior kernel: the non-stationary Gibbs kernel with a tanh
  core/edge length scale (Gibbs 1997; Paciorek and Schervish 2004, NIPS 16),
  as used in Minerva and in Köberl's surroBIER. TIAGO
  applies it to the log-profile values directly instead of to the log of the
  gradient, so plateaus (flattening) remain in the support of the prior.
- surroBIER (R. Köberl, MPCDF, MIT) and vmecnn (A. Merlo, MPCDF; Merlo et al.
  2023, JINST 18 P11012) are surrogate-based comparisons; TIAGO runs the
  VMEC++ solver and its implicit adjoint inside the inference instead.

Chris&AI
