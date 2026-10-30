# Two-fidelity electron-temperature prior

Use `tiago_multifidelity_profile::ar1_profile_t` when Thomson scattering
constrains the physical electron temperature and electron cyclotron emission
(ECE) constrains a calibrated lower-fidelity temperature. Both instruments
remain recorded-signal forward models in the joint likelihood. This prior
implements the autoregressive relation

$$T_{\rm high}(s)=C\,T_{\rm low}(s)+\delta(s),$$

where $s$ is normalized toroidal flux, $C$ is the ECE temperature calibration
scale and $\delta$ is an independent Gaussian-process (GP) shape discrepancy.
The high-fidelity temperature enters kinetic pressure and Thomson; the low
level enters ECE's medium, opacity and radiometer prediction. The convention
for the calibration record is **high = C times low**; invert a record using
the opposite convention before constructing this prior. Receiver voltage
calibrations remain in the recorded ECE model and describe a different
quantity. Do not add a second temperature-scale or slot gain nuisance.

The covariance model follows [Kennedy and O'Hagan
(2000)](https://doi.org/10.1093/biomet/87.1.1). Generic conditioning, normalized
log marginal likelihood and its gradients are in
[fortnum](https://github.com/lazy-fortran/fortnum/blob/main/docs/multifidelity.md).
Its Schur-complement recursion retains uncertainty at the low level for noisy,
non-nested observations. The independent-level simplification of
[Le Gratiet and Garnier (2014)](https://doi.org/10.1615/Int.J.UncertaintyQuantification.2014006914)
applies to noiseless nested designs; no such simplification is assumed here.

## Setup and parameters

Initialize `kinetic_profiles_t` first. Supply its `greville` nodes to
`prior%init(nodes, mean_low, calibration, log_mean, log_sd, hyper_source,
constraint_source, status)`, then call `profiles%set_multifidelity_te(prior)`
before `joint_model_t%setup`. `mean_low` contains physical temperature spline
coefficients in electronvolts (eV), not logarithms. A spline's coefficients
are associated with the Greville nodes but generally differ from its values
at those nodes. The calibration is a scalar `calibration_prior_t` with a
nonempty `record%record_id`; its mean and factor define the normal prior on
$C$. Device-specific records stay in the device adapter.

`log_mean` and `log_sd` specify independent normal priors on log low-level
amplitude, log low-level length, log discrepancy amplitude and log discrepancy
length, in that order. Amplitudes are standard deviations in eV; lengths are
in $s$. Both kernels are squared exponential. Record the physical rationale
for these priors in `hyper_source`; there are no default fitted lengths.
Every amplitude, length and unfixed calibration scale is inferred jointly
with the profiles and equilibrium. The posterior uses independent standard
normal coordinates, with parameter order: low GP coefficients, discrepancy
coefficients, two or four hyperprior coordinates, and the calibration
coordinate when unfixed. `hyper_index()` gives the offset before the first
hyperprior coordinate; `scale_index()` gives the calibration index, or zero
when fixed. Add the profile's `offset` when indexing the full profile vector.

The discrepancy lies in the Euclidean orthogonal complement of the constant
vector and `mean_low`, plus optional columns `protected`. With orthonormal
basis $Q$ of that complement and discrepancy kernel $K_\delta$, it is
$\delta=Q\,\operatorname{chol}(Q^T K_\delta Q)z_\delta$. These are constraints
on the spline coefficients. They protect a fixed local calibration tangent,
not every possible temperature shape, and do not prove global
identifiability. Describe the chosen protected space in `constraint_source`.
No numerical jitter is inserted; singular covariance factors return failure.

This is a Gaussian prior in physical-temperature coefficients. Unlike the
existing log-temperature prior, it has negative temperatures in its support;
it is suitable when temperature is well separated from zero relative to its
prior uncertainty. Recorded instruments retain their physical-domain checks.
The existing exponential scrape-off continuation applies beyond $s=1$.

## Joint likelihood and local uncertainty

ECE slots automatically use the low level when the prior is active. Other
slots use the high level unless `low_fidelity_te` is explicitly set, which
also permits synthetic point-probe tests. The pressure always uses the
physical high level. Jacobians include calibration, kernel and pressure
sensitivities and the recorded ECE moving-ray chain.

All GP and hyperprior coordinates contribute normalized standard-normal prior
rows. In these non-centered coordinates no additional kernel determinant
belongs in the joint density: the Cholesky transformation generates the
conditional GP measure. The MAP is defined in these coordinates. For an
active AR(1) prior, `laplace` uses the observed posterior Hessian, obtained by
central differences of the analytic score, including residual curvature
omitted by Gauss--Newton. Its evidence includes the observation and prior
normalization and the observed Hessian determinant. It describes one local
mode; multiple equilibrium modes still require separate fits and evidence
weights.

`discrepancy=.false., fixed_scale=.true.` removes discrepancy and calibration
coordinates, leaving a single GP seen at the low level and at $C$ times that
level. Its two kernel hyperparameters remain posterior parameters.

The synthetic CTest `tiago_multifidelity_profile` checks calibration and shape
recovery with local Laplace intervals, repeated independent measurement noise,
profile and hyperparameter derivatives, and the single-GP limit. Its CSV
contains the calibration, one discrepancy coefficient, their standard errors
and normalized log evidence for each trial. This coverage study conditions on
one planted truth and does not establish coverage over arbitrary hyperpriors
or equilibrium modes. `tiago_ece_recorded` checks the composed ECE derivatives;
`tiago_joint_model`, in reconstruction builds, checks the pressure/equilibrium
chain, MAP and evidence with both fidelities.
