# Thermal transfer and radiometer response

`tiago_ece_transfer` and `tiago_ece_receiver`, in `tiago_inference`, compose
thermal radiation transfer on supplied rays with a calibrated total-power
receiver. They accept runtime absorption, paths and instrument response.
The [recorded ECE composition](ece-recorded.md) connects these primitives
to the equilibrium medium and perpendicular relativistic opacity. This page
describes transfer and receiver inputs separately.

## Transfer on a supplied ray

Call `thermal_transfer(te,alpha,length,rho,wall_te,trace,status)`, with one
entry per cell, ordered from the observer into the plasma. `te` is electron
temperature in electronvolts (eV), `alpha` is absorption in inverse metres,
and `length` is physical arc length in metres. Cells have constant absorption
and temperature. The thermal source function equals `te` in the
Rayleigh--Jeans limit; the emissivity in brightness-temperature units is
`alpha*te` in eV/m. The transfer equation and its constant-source solution
are given in Rybicki and Lightman, *Radiative Processes in Astrophysics*
(1979), section 1.4, Eq. (1.30).

The generated cell map is `Tout = Tin + (1-exp(-alpha*length))*(te-Tin)`.
The observer-directed pass processes the cells in reverse order. This
retains foreground absorption of emission from deeper cells. Small optical
depths use the tenth-order alternating series for `1-exp(-tau)` below
`tau=1e-3`; its value remainder is below `2.6e-41` at the switch, before
floating-point rounding. Values and derivatives come from the same generated
expression. A spatially varying source still requires cell refinement.

Check zero `status` and `trace%valid` before using `trace%temperature` (eV)
or `trace%optical_depth` (dimensionless). Negative or nonfinite inputs,
inconsistent shapes, and optical-depth overflow fail. Zero absorption,
zero-length cells and empty rays are allowed. A new failed evaluation
invalidates the tape. `thermal_transfer_batch` takes arrays `(cell,ray)`,
one wall pair per ray and one status per ray, and uses OpenMP over rays.

## Optional reciprocal-wall model

`rho` is power reflectivity, `0 <= rho < 1`; `wall_te` is wall brightness
temperature in eV, common to both ends. Both are explicit inputs. A caller
inferring them must map them through a calibration prior and keep them in
their physical domain. `rho=0,wall_te=0` gives no incoming radiation.

This is a two-stream model: each wall returns radiation along the same
reciprocal ray, with no mode conversion or diffuse redistribution. If `En`
and `Ef` are the zero-input emissions towards the near and far ends, and
`P=exp(-sum(alpha*length))`, the boundary equations are

```text
Tnear = En + P*((1-rho)*wall_te + rho*Tfar)
Tfar  = Ef + P*((1-rho)*wall_te + rho*Tnear).
```

The generated solution sums all reflections. Equal wall and plasma
temperatures are an equilibrium, including near the transparent, reflective
limit. A perfectly closed cavity is excluded. A cold-ray turning point is
part of the supplied path, not an extra wall. Real wall geometry, multiple
different paths, diffuse scattering and nonreciprocal propagation need a
different boundary model.

## Frequency, antenna response and calibration

Construct `ece_response_channel_t` entries with positive frequency nodes in
Hz, nonnegative frequency quadrature `weights` in Hz, dimensionless power
`response`, nonnegative antenna quadrature `antenna`, and integration
`duration` in seconds. `rays(angle,frequency)` indexes the supplied ray
temperatures. The antenna weights include pattern and solid-angle quadrature;
they are normalized on initialization. Frequencies identify the rays; the
receiver does not evaluate an opacity or trace a path at those frequencies.
The normalized frequency response averages brightness before calibration.
Uncertain response shapes, pointing and durations are not inferred here.

Call `receiver%init(channels,n_rays,prior,status)` with the existing
`calibration_prior_t`. Its three physical parameters per channel are, in
order: gain (V/eV), offset (V), and receiver noise temperature (eV). They
use the existing correlated affine Gaussian-prior map from standard-normal
latents. Gain must be positive and receiver temperature nonnegative at
every evaluation. The offset includes the receiver's mean background.
`evaluate(temperature,eta,voltage,variance,status)` returns V and V squared;
both outputs must be ignored on failure. `n_nuisance()` gives the latent count.

For antenna-averaged brightness `Tf` at each frequency and power response
`R`, the ideal single-polarization total-power variance is

```text
variance = gain**2 * integral R**2*(Tf+receiver_te)**2 df
           / (duration * (integral R df)**2).
```

For constant brightness this reduces to the
[ideal radiometer equation, Eq. (3.154)](https://www.cv.nrao.edu/~sransom/web/Ch3.html),
with effective bandwidth `(integral R df)**2/integral R**2 df` in Hz,
returned by `effective_bandwidth()`. The full quadrature retains frequency
variation of the noise temperature. This is one receiving spatial mode;
antenna nodes integrate its brightness, not independent receivers. Flicker
noise, digitization, switching penalties and inter-channel covariance are
excluded. Do not add this variance to observation covariance that already
contains it. A likelihood using the parameter-dependent variance must also
include its normalization and derivative.

## Derivatives and verification scope

The transfer tape's `jvp` maps perturbations of every temperature, absorption,
length, reflectivity and wall temperature to radiation temperature. Its `vjp`
returns the transpose in the same order. The receiver's `jvp` and `vjp`
differentiate both voltage and variance with respect to ray temperatures
and calibration latents, including shared ray indices and correlated priors.
These names mean Jacobian-vector and transposed Jacobian-vector products.
All adjoint outputs are overwritten. Channel loops use OpenMP.

[Perpendicular thermal opacity](ece-opacity.md) can now supply the absorption
and its temperature, density, field-magnitude and frequency derivatives.
The caller must supply path derivatives. In particular,
length sensitivity alone does not account for profiles sampled on a moving
ray. The [recorded adapter](ece-recorded.md) supplies these moving-ray derivatives
and the `profile_diagnostic_t` implementation.

`tiago_ece_transfer` checks vacuum, thin emission, thick thermal equilibrium,
nonuniform-source refinement, independent repeated wall reflections,
frequency/antenna averages, the radiometer equation, every derivative column,
transpose identities, calibration domains and per-ray errors. The optional
`tiago_ece_transfer_codegen` regenerates all six FortSym kernels and compares
them with committed output. `benchmark_ece_transfer` measures independent
256-cell rays including allocation and derivative-tape construction.
`benchmark_ece_receiver` measures 256 channels with four frequency and four
antenna nodes each, including calibration-prior transformation and noise.

These are transfer/receiver verifications with synthetic prescribed opacity.
Published optical-depth checks belong to the [opacity tests](ece-opacity.md);
shifted-resonance recovery belongs to the [recorded tests](ece-recorded.md).
General oblique relativistic opacity remains unavailable.
