# Joint inference from recorded signals

Use one `tiago_joint_model::joint_model_t` for a common equilibrium and kinetic
profile state. Interferometer, recorded Thomson and electron cyclotron emission
(ECE) slots share its density
and electron-temperature B-splines; magnetic fluxes constrain the equilibrium,
whose pressure is computed from the same kinetic profiles. The equilibrium
adjoint carries geometry, magnetic-field magnitude and contravariant field
sensitivities back to equilibrium and profile parameters. ECE uses the composed
moving-ray tape, including launch, refraction, cell lengths, opacity, walls and
receiver response. The diagnostic kernels and calibration leaves are reused.

A device adapter supplies logical signal names, units, versioned response
metadata, geometry and prior distributions. It implements `signal_adapter_t`
and keeps archive addresses and device configuration outside this repository.
Read bounded windows through `online_source_t`, optionally reduce them online,
or read `tiago-signal-1` offline files through `netcdf_source_t`. Feed the
resulting samples and full measurement covariance into diagnostic slots.
This joint likelihood conditions on that fixed observed covariance; it does
not differentiate the ECE receiver's predicted thermal variance.
Avoid retaining complete raw archives. Check source versions, units, response
validity intervals and the reduction covariance before inference.

Initialize every recorded diagnostic with its explicit `calibration_prior_t`.
The joint standard-normal nuisance coordinates are mapped through `eta_factor`;
its default diagonal contains each slot's `eta_sigma`. Off-diagonal blocks can
share a calibration or correlate repeated observations. `campaign_factor`
constructs a shared initial offset and Brownian drift increments at strictly
increasing times. Insert its factor in the corresponding nuisance rows; these
rows are in diagnostic nuisance coordinates, before each instrument's physical
calibration factor. A common magnetic integrator has gain, offset and linear
drift priors in `mag_calibration`, observation times in `mag_time` (seconds),
and `mag_time_origin`. Its rows must have a common physical signal unit.
Repeated observations here share a kinetic state. Evolving profiles require
separate state blocks; a drift prior alone does not supply profile dynamics.

The existing `energy_consistency_t` can add the diamagnetic/kinetic-energy
constraint in joules. Supply its physical nuisance means in `eta_mean` and
their prior factor in `eta_factor`, together with fixed real-space volume
quadrature. Use it within its documented diamagnetic approximation. Account
for shared measurements before combining an energy estimate with its source
flux signals; treating both as independent observations double-counts data.

Set a slot's `student_nu` or `mag_student_nu` to positive degrees of freedom
for independent Student likelihood rows after covariance whitening; zero
selects Gaussian rows. Prior rows remain Gaussian. Add an instrument's proposed
discrepancy columns in `discrepancy_basis` and its protected profile/calibration
sensitivity columns in `protected_jacobian` before `setup`. The covariance-metric
projection removes their overlap, including dependent protected columns.
Freeze that projection during a posterior evaluation. This implements a local
constraint motivated by [Kennedy and O'Hagan (2001)](https://doi.org/10.1111/1467-9868.00294);
it does not establish global nonlinear identifiability. Protect sensitivities
across the intended parameter region and validate on independent truths.

For two-fidelity electron temperature, install the
[AR(1) profile prior](multifidelity-profile.md) before joint setup. Thomson and
kinetic pressure use the physical high level; recorded ECE uses the low level,
linked by its temperature calibration scale and a constrained shape GP. Kernel
hyperparameters are posterior parameters, with sourced lognormal priors. Local
Laplace evidence uses the observed Hessian in standard-normal coordinates.

For a cheap proposal model, call `train_surrogate` with paired exact design
sites, fixed input lengths, correction amplitude, autoregression coefficient
and nugget standard deviation. Its inexpensive level is the joint tangent
model, corrected by the [Kennedy and O'Hagan (2000)](https://doi.org/10.1093/biomet/87.1.1)
Gaussian process in [fortnum](https://github.com/lazy-fortran/fortnum).
`use_surrogate` enables the conditional mean only
inside the componentwise training box; prior rows remain exact. Retrain after
changing data, calibration factors, discrepancy or covariance. Validate on
held-out exact evaluations. Disable the surrogate for exact final likelihoods,
No-U-Turn sampling targets or `psis_validate` importance correction. Emulator variance is
conditional uncertainty, not automatically added observation noise, and its
mean alone must not be reported as the exact posterior.

Validate a new device module first with synthetic responses: reproduce every
instrument's analytic limits and published references, check complete joint
finite differences, and require offline/HTTP posterior parity. Inject density
scale, temperature/gain and integrator errors. Repeat noisy experiments, report
pointwise coverage with binomial uncertainty, compare fixed-calibration and
uncoupled fits, and plot truth and posterior density and temperature with units.
A converged local fit or a positive covariance is insufficient: check residuals,
competing equilibrium/radiation modes and exact likelihood corrections. The
public `tiago_joint_recorded` test uses the pinned public low-resolution VMEC
case already used by reconstruction tests; all signals and geometry are
synthetic. It checks the complete derivative and source parity. Extended
coverage and ablation studies use the same fixture and must report the tested
parameter/prior range and the limitations of a Laplace approximation.
