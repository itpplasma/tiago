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

Recorded-signal extensions, versioned sources, reductions and calibration priors
are specified in [recorded-signals.md](recorded-signals.md). They use the same
state vector and derivative contract.

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
cotangent to `tiago_equilibrium%jacobian` (`coef_cot`) and the bmnc
cotangent as `bmnc_cot`. VMEC++ has no adjoint of its output stage, so Tiago
recomputes |B| on the half grid from the geometry coefficients, phi', iota and
the enclosed current as VMEC's bcovar does (|B|^2 = B^u B_u + B^v B_v,
`vmec_edge%modb`, equal to the wout bmnc to round-off) and differentiates that
by hand (`vmec_edge%modb_vjp`). The geometry part of its cotangent goes
through the VMEC++ implicit adjoint, the direct iota and current parts join
the profile cotangents, and PHIEDGE uses the same scaling argument as the
edge field. |B| at a fixed label therefore has complete derivatives with
respect to the equilibrium inputs, which diagnostics that read |B| (ECE)
need.

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

sCluster array 3370173 ran 50 seeds of `tiago_infer` at c2328c8 (seeds
101-150; Clang 22.1.8 / Enzyme 0.0.293, 4 threads each, 10 to 24 min per
seed). Each seed draws the truth from the prior (CURTOR, PHIEDGE, 7 B-spline
latents each for ne, Te, Ti, a (dR, dZ) probe alignment: 25 latents), adds
noise to 54 magnetic signals and 24 profile probes (5 %), and runs MAP,
Laplace, PSIS (400 draws), NUTS (300 warm-up, 400 draws) and 400 re-solved
draws for iota (flux surfaces are written for the first 40). Fractions of
truths inside the central intervals (nominal 0.683/0.954 for +-1/2 sd,
0.68/0.95 for quantile intervals; brackets are 95 % seed-bootstrap
intervals):

| | 68 % | 95 % |
|---|---|---|
| latents, Laplace | 0.671 [0.643, 0.699] | 0.938 [0.920, 0.954] |
| latents, NUTS | 0.670 [0.643, 0.694] | 0.929 [0.911, 0.946] |
| ne, Te, Ti, p at 41 s values, Laplace | 0.665 [0.636, 0.692] | 0.937 [0.920, 0.952] |
| same, NUTS | 0.675 [0.646, 0.703] | 0.940 [0.923, 0.955] |
| iota, 400 re-solved draws | 0.724 [0.641, 0.804] | 0.949 [0.908, 0.983] |

The empirical quantile band of N independent draws contains a calibrated
truth with probability 0.676/0.946 for N = 400 (Monte Carlo), so the
expected values for the draw-based rows are 0.676/0.945 (NUTS) and
0.676/0.945 (iota); all intervals contain them. Values within a seed are
correlated, so `tools/inference_coverage.py` resamples whole seeds, and
binomial intervals over all counted values understate the uncertainty, most
for iota, whose values at all s follow CURTOR. Per latent (one trial per
seed, Wilson 95 % intervals of roughly +-0.13) the Laplace coverage lies
between 0.54 and 0.82 at 68 % and between 0.84 and 1.00 at 95 %. The lowest
is the alignment latent align_1 (0.54 [0.40, 0.67], 0.84 [0.71, 0.92]),
which is the one outlier expected among 25 latents at this precision;
CURTOR has 0.64 [0.50, 0.76] and 0.92 [0.81, 0.97].
`tools/inference_coverage.py` writes this table (`coverage_parameters.txt`)
and the per-latent calibration figure. NUTS had no divergences, mean
acceptance 0.80-0.85, tree depth 3 and E-BFMI 0.9-1.4; PSIS k-hat had
median 0.43 and exceeded 0.7 in one seed.

The magnetics constrain CURTOR strongly: the median Laplace sd is 0.0185 of
the prior sd (about 370 A of the 2e4 A prior), and the joint Fisher
information from the magnetic rows equals that of `tiago_reconstruction`.
Consequently the iota posterior sd is small (median over seeds 0.002 at the
sixth and 0.005 at the third-last of the 15 half-grid surfaces, against about
0.02 before the fix) and the remaining profile uncertainty is that of the
kinetic latents: the posterior sd is 0.17 of the prior sd for the leading
B-spline latent of each of ne, Te, Ti, but 0.5 to 1.0 for the others, so the
profile diagnostics, not the magnetics, limit ne, Te, Ti at mid radius and
the edge, and PHIEDGE keeps 0.60 of its prior sd. The synthetic magnetics
must be whitened with the sigmas that `synthesize_magnetics` sets, which the
joint model now does (a regression test compares the joint magnetic Fisher
information with the reconstruction rows). `tools/inference_plots.py`
writes, for one seed, flux surfaces at three toroidal angles with 68/95 %
bands of the displacement from the truth, profiles (p, ne, Te, Ti, iota,
enclosed current) with bands, and the corner plot.

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
