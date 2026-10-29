# Recorded electron cyclotron emission

`tiago_ece_recorded::record_ece` composes the existing equilibrium/profile
medium, cold rays, perpendicular weakly relativistic opacity, thermal transport
and calibrated radiometer. It returns voltage in volts (V) and thermal noise
variance in V squared. No device configuration or measurements are built in.
This is a restricted near-perpendicular electron cyclotron emission (ECE)
model; check the domain below before using it for inference.

## Configuration and evaluation

Create one `ece_channel_t` for each frequency/antenna ray. Its frequency is
in Hz, origin and direction are Cartesian, `length` is the positive
Hamiltonian integration duration in metres, and `nscan` is the number of
transport cells. `length` is not arc length. Select ordinary or extraordinary
polarization with `ECE_O_MODE` or `ECE_X_MODE` and harmonic 1 or 2.
A direction is normalized at launch and the index is placed on the cold
mode's dispersion surface, including its medium dependence.

Construct `ece_response_channel_t` objects using the existing
[receiver quadrature contract](ece-transfer.md). Their ray indices refer to
these `ece_channel_t` entries; the frequency at each response node must match
the referenced rays. Several frequencies and antenna nodes can feed one
recorded channel, and rays can be shared.

```fortran
use tiago_ece_radiometer, only: ece_t
use tiago_ece_recorded, only: recorded_ece_t, record_ece
! diagnostic is ece_t; tape is recorded_ece_t.
call diagnostic%init(rays, response, receiver_prior, wall_prior, code)
! Check code == 0. medium is the existing state_medium_t.
call record_ece(diagnostic, medium, eta, tape, status)
! Check status%code == 0 and tape%valid before using any tape output.
! tape%voltage(:), tape%variance(:)
call tape%jvp(parameter_dot, eta_dot, voltage_dot, variance_dot, status)
call tape%vjp(voltage_bar, variance_bar, parameter_bar, eta_bar, status)
```

The medium and its profile latents follow the
[equilibrium-medium interface](ece-geometry.md). The combined parameter vector
is the provider's output coefficients followed by the existing kinetic-profile
latents. Inputs are copied into the tape; a subsequent provider update does
not change a previous evaluation or its derivatives. JVP means Jacobian-vector
product; VJP returns its transpose. Adjoint outputs are overwritten.

`eta` contains the standard-normal receiver calibration latents first, followed
by standard-normal wall latents. The correlated affine Gaussian maps are the
existing `calibration_prior_t`; the two blocks are independent. Receiver
physical parameters are `(gain [V/eV], offset [V], noise temperature [eV])`
for each recorded channel. Wall parameters are `(power reflectivity [1],
brightness temperature [eV])` for each ray. Every one has a supplied prior;
there is no implicit unity-gain calibration. Physical domains remain those of
the transfer and receiver components. Response shapes, launch geometry,
frequency and integration duration are supplied fixed inputs; their uncertain
calibration would require additional parameters and priors.

Pass observations and covariance from the generic
[signal sources and reductions](recorded-signals.md). Do not add the computed
thermal variance to a covariance that already contains it. A likelihood that
uses parameter-dependent variance must include its normalization and the
variance derivative, as the posterior test does.

## Sampling and the profile contract

Each ray is traced in contiguous intervals of fixed Hamiltonian duration.
Transport uses the distance between interval endpoints and the plasma at their
Cartesian midpoint. Refine `nscan` to resolve curvature, emission layers and
source gradients; adaptive ray accuracy alone does not resolve the transport
integral. Samples remain ordered from the observer into the plasma.

The composed tape includes launch-index derivatives, every ray interval,
midpoint displacement, profile sampling, vector-field magnitude, cell length,
opacity, reciprocal walls, response, calibration and thermal noise. Ray
sensitivities retain the existing accepted-step schedule convention;
adaptive decisions, duration and terminal events are inactive.
Independent rays use OpenMP in forward, tangent and reverse evaluations.
No opacity or transport physics is duplicated.

`tiago_ece_radiometer::ece_t` now implements the recorded
`profile_diagnostic_t` contract on these prepared cells, replacing the old
cold-resonance interpolation. `tape%sampled(prepared,sample)` exports this
sampled diagnostic for contract checks or fixed-path applications. Its
`predict`, `jvp` and `vjp` use the existing sampled rows `Q_B`, `Q_NE` and
`Q_TE`, with receiver and wall priors in `eta`. `record` also returns variance;
`record_jvp` and `record_vjp` additionally expose cell-length sensitivities.
`sample_points` returns cylindrical `(R,phi,Z)` points; calibration latents
leave those points unchanged.

For changing-provider inference, call `record_ece` at every evaluation and
use its composed tape. Exporting a prepared sampled diagnostic alone freezes
its paths and omits their provider dependence. The existing `joint_model_t`
slot does not automatically dispatch this composition. The supplied
`plasma_sample_t` from the tape uses the norm of the vector field, consistently
with ray tracing and opacity; it does not substitute the separately interpolated
VMEC scalar magnitude spectrum. Provider-output cotangents still need the
provider's own implicit adjoint to infer equilibrium input parameters. That
solved-equilibrium chain is not verified here.

## Physical domain and angle bound

All [opacity guards](ece-opacity.md) and [cold fast-branch guards](ece-rays.md)
remain active. Every sampled cell must lie in the provider's valid chart and
in the near-resonance opacity domain. An unsupported exterior field is not
replaced by a guessed vacuum field. There is no cutoff tunnelling, mode
conversion or oblique opacity.

Let `theta=Te/(510998.95069 eV)` and let `mu` be the cosine of the angle to
the magnetic field. The adapter requires

```text
abs(mu) <= 1e-4 * sqrt(theta)
```

for both the midpoint wave-index vector and the physical ray segment. It checks
all cells, a stronger requirement than checking only the emitting layer.
At Te = 1000 eV this allows about 4.42e-6 radians from perpendicular.
Larger angles fail with an explicit perpendicular-opacity error. This is a
conservative small-Doppler domain, not a general oblique model.

The first harmonic of the extraordinary mode additionally requires the
**tenuous** regime `r <= 0.2*theta`, where `r=(omega_p/omega_c)**2` and
`omega_p`, `omega_c` are electron plasma and cyclotron angular frequencies.
The adapter rejects larger density, including the intermediate regime;
it does not enable the separate finite-density opacity branch. A launch
preflight reports this restriction before tracing. Propagation restrictions
may still prevent access to an otherwise valid local emission layer.

The angle control uses the Maxwellian weakly relativistic resonant-shell
limit from [Bornatici et al. (1983), sections 2.3 and 3.1](https://doi.org/10.1088/0029-5515/23/9/005).
With velocity `u=v/sqrt(2*Te/m_e)`, resonance coordinate
`z=(1-n*omega_c/omega)/theta`, and `a=N_parallel/sqrt(2*theta)`, the resonance
surface is `u_perp**2+(u_parallel-a)**2=a**2-z`.
Here `n` is the harmonic and `N_parallel` the index parallel to the field.
The test integrates the Maxwellian and perpendicular-velocity weighting on
this displaced sphere for tenuous first/second-harmonic extraordinary modes.
Across `-10 <= z <= -0.1` at the angle bound its line-core difference from the
perpendicular result is below 1.2e-7 relative. This is a leading tenuous-X
limit check; it does not certify ordinary-mode finite-angle errors, arbitrary
density, the zero-absorption edge at z=0, or remote wings.

## Verification and limits

`tiago_ece_recorded` checks an optically thick isothermal synthetic VMEC chart
through the complete chain, vacuum wall radiation, mesh refinement, the sampled
contract, combined provider/profile/calibration finite differences, voltage and
variance transpose identities, and angle/density rejection. The isothermal
limit has `T_rad=Te` at the resonance and throughout the emitting region.
For nonuniform temperatures at finite linewidth, the radiation temperature
is a transport-weighted value and need not equal Te at the cold resonance.

`tiago_ece_posterior` performs 128 repeated independent noise realizations
with six synthetic frequency channels, a known common gain error and a shifted
magnetic resonance. It jointly fits one temperature-profile direction, the
common gain and the field shift with standard-normal priors. Each fit starts
at the maximum for the averaged data, which is found from the prior mean, and
iterates to the same convergence criterion as a cold start. The likelihood
retains parameter-dependent thermal noise and its normalization. Intervals
use the local Fisher/Laplace covariance. This is conditional synthetic recovery
in a three-dimensional subspace, not coverage of every B-spline coefficient,
a solved-equilibrium posterior, a continuum-discretization certificate or
experimental validation. The CSV contains every fit and interval width.

`tiago_ece_recorded_codegen` regenerates and byte-compares the six FortSym
kernels and all derivative rows when `TIAGO_FORTSYM_BUILD` is configured.
`benchmark_ece_recorded` measures forward, tangent and reverse throughput
on 64 independent synthetic rays with 48 cells, including tape construction
in the forward timing. Existing component tests retain the published slab,
cold, vacuum, Kirchhoff and transfer checks.
