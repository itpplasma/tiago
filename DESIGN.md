# TIAGO design

Status: architecture direction, 2026-10-05. The implemented production path is
the VMEC++ magnetic reconstruction described in README.md and STATUS.md.
Provider-neutral KIN6D integration is planned, not implemented.

## Mission and ownership

TIAGO owns inference from experimental/diagnostic data. Its responsibilities
are:

- diagnostic/instrument forward models and data ingestion;
- measurement covariance and systematic/calibration nuisance models;
- parameter priors and constraints;
- MAP/least-squares reconstruction;
- local posterior/Laplace covariance and validation;
- comparison of several candidate physical states/equilibria;
- optional design-of-diagnostic/experiment information built from the same
  sensitivities.

TIAGO does not own the general plasma equations. Grad--Shafranov, general 3-D
MHD, resistive/XMHD, kinetic-MHD and DK/GK/FK stationary states belong to their
forward provider. KIN6D is the planned general provider.

## Current VMEC++ provider

The current reconstruction chain is valuable and remains supported:

    VMEC++ equilibrium
        -> TIAGO edge/boundary representation
        -> TIAGO magnetic diagnostics
        -> exact diagnostic reverse mode
        -> VMEC++ implicit adjoint
        -> Levenberg--Marquardt
        -> local covariance

This path is independently benchmarked against STELLOPT/DIAGNO and finite
differences. It must not be rewritten merely to satisfy a speculative generic
interface.

The first provider abstraction should be extracted only when a real KIN6D
consumer exists. Preserve the existing VMEC path as a reference backend and
test it through the same behavioral contract.

## Forward-provider seam

The stable seam should be discovered from the first two concrete providers,
not designed exhaustively in advance. Conceptually TIAGO needs:

1. named/scaled physical parameters or parameter actions;
2. solve/update of a selected stationary/evolution state;
3. the physical state or observable information required by the diagnostic;
4. matrix-free tangent/JVP actions when forward sensitivities are useful;
5. matrix-free VJP/implicit-adjoint actions for many-parameter inference;
6. provider provenance and failure/status information.

TIAGO should not differentiate Newton/Krylov iterations or depend on a
provider's internal state layout when equation-level derivatives are
available.

The diagnostic layer may itself be differentiated in reverse and composed with
the provider adjoint. Physics-specific emission/response primitives may be
supplied by KIN6D while TIAGO owns the instrument/calibration/likelihood layer.

## Reconstruction objective

For data d, physical parameters theta and diagnostic nuisance parameters eta,
a typical problem is

\[
d = H(U(\theta),\eta)+\epsilon,
\qquad
R(U,\theta)=0.
\]

TIAGO owns H's instrument semantics and the probability/error model for
epsilon. The provider owns R and the state solve.

The current implementation minimizes a whitened least-squares objective with
Gaussian priors. Preserve this fast deterministic path as the default.

## Uncertainty quantification

### Current level: local Laplace/Gauss--Newton

For a whitened residual r and Jacobian J, the current covariance is

\[
C_{\rm post}\simeq (J^T J)^{-1}.
\]

The synthetic pull study is the acceptance oracle for claims about this
covariance, not merely visual agreement of one fit.

### Next level: correlated/systematic data errors

Generalize whitening to a declared covariance or factor without forcing users
to diagonal independent sigmas. Shared gain, offset, alignment, coil-current,
timing or other calibration uncertainties should normally enter as nuisance
parameters or correlated covariance structure.

Keep diagnostic-model uncertainty identifiable separately when it is not
already represented by nuisance parameters.

### Several equilibria

If KIN6D or another provider supplies distinct equilibria A and B, treat them as
two forward states. Fit each separately:

\[
p(\theta\mid d,A)\approx N(\hat\theta_A,C_A),
\qquad
p(\theta\mid d,B)\approx N(\hat\theta_B,C_B).
\]

When useful, attach explicit evidence/prior-derived weights and report a finite
mixture of local Laplace approximations. This may be called a Gaussian mixture,
but the equilibria themselves remain distinct and are not averaged.

TIAGO is not responsible for discovering every stationary solution. Providers,
continuation, different initial states and external scientific hypotheses
supply candidate equilibria.

### Advanced sampling

HMC/SMC/general MCMC is deferred. It becomes justified only when a concrete
posterior is materially non-Gaussian inside one solution basin or a finite
mixture of local Laplace approximations demonstrably fails.

## Orthogonal error ledger

Do not collapse conceptually different errors into one sigma. Report separately:

1. measurement/statistical uncertainty;
2. calibration/nuisance posterior uncertainty;
3. reconstruction uncertainty conditional on a selected forward model/state;
4. diagnostic-forward-model discrepancy where applicable;
5. provider numerical/discretization certification (KIN6D when used);
6. provider physical model/reduction error (KIN6D when used);
7. provider stochastic sampling uncertainty;
8. existence of several distinct candidate equilibria.

A deterministic interval and a posterior standard deviation are not combined by
root-sum-square without a declared joint probabilistic model.

## Data and reproducibility

Keep Fortran as the production implementation; no Python runtime is required.
Retain exact input, provider revision/configuration, diagnostic definitions,
measurement covariance/nuisance model, parameterization, optimizer settings and
UQ method in reconstruction provenance.

When KIN6D provides a self-describing NetCDF-4 result, prefer consuming its
normalized state/provenance contract rather than re-parsing source-code-specific
files. libneo remains the shared interchange/convention library where
appropriate.

Chris&AI
