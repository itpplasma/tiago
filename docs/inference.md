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
the kinetic quantities are ne = Te = Ti = 0 and Zeff = 1.

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

## References and reference implementations

- Profile prior kernel: the non-stationary Gibbs kernel with a tanh
  core/edge length scale used in Minerva and in Köberl's surroBIER
  (`src/surrobier/prior/GaussianProcessPrior.py`, `minerva_kernel`). TIAGO
  applies it to the log-profile values directly instead of to the log of the
  gradient, so plateaus (flattening) remain in the support of the prior.
- surroBIER (R. Köberl, MPCDF, MIT) and vmecnn (A. Merlo, MPCDF; Merlo et al.
  2023, JINST 18 P11012) are surrogate-based comparisons; TIAGO runs the
  VMEC++ solver and its implicit adjoint inside the inference instead.

Chris&AI
