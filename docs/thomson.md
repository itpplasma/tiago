# Thomson scattering to recorded pulse integrals

`tiago_thomson_polychromator::recorded_thomson_t` predicts photoelectrons and
voltage integrals from the existing real-space kinetic-profile sample contract.
It extends `recorded_diagnostic_t`, reuses `thomson_t` for spatial quadrature and
alignment, and uses `calibration_prior_t` for correlated calibration parameters.
It can therefore be composed with the existing kinetic-profile latents, field
providers and generic signal sources without a device-specific adapter.

## Spectrum and polarization

`tiago_thomson_spectrum::thomson_spectrum_t` implements the fully relativistic,
incoherent Maxwell--Jüttner spectrum for electron temperatures of 10--20000 eV.
Incident light and the detected polarization are perpendicular to the scattering
plane. The Selden spectrum includes the Naito depolarization integral; the latter
uses Gauss--Laguerre quadrature (24 nodes by default). Exponentially scaled
modified Bessel functions from [fortnum](https://github.com/itpplasma/fortnum)
avoid underflow in the cold limit.

After `init`, `cross_section(Te,lambda0,theta,wavelength,s,dsdt)` returns the
photon differential cross-section divided by the squared classical electron
radius, in inverse metres, and its temperature derivative. Temperatures are in
electronvolts, wavelengths in metres and scattering angles in radians. This is
the **photon** spectrum: conversion from the energy spectrum includes the
scattered-to-incident wavelength ratio. Its wavelength integral retains the
relativistic angular cross-section correction; it need not equal one.

`normalized` has the same arguments and an optional output `area`. It returns a
unit-area photon probability density, its quotient-rule temperature derivative,
and the unnormalized area. Its composite logarithmic wavelength rule is separate
from the channel integration rule. Use the unnormalized cross-section for
absolute signal prediction. Zero scattering angle is excluded. Normalization
is tested at 0.1, pi/2 and pi radians, not arbitrarily close to zero angle.

The physical basis is [Selden (1980)](https://doi.org/10.1016/0375-9601(80)90276-5),
[Naito et al. (1993)](https://doi.org/10.1063/1.860593), and the equivalent
velocity-space expression in
[Mirnov and Den Hartog (2017), Eq. (60)](https://doi.org/10.1088/1361-6587/aa61e2).
Collective scattering, drift, arbitrary polarization, multiple scattering,
finite collection-angle averaging and finite-transit-time corrections are not
included.

## Channels, pulses and priors

Construct `polychromator_channel_t` entries with a volume index, laser wavelength
and pulse energy (joules), scattering angle, collection solid angle (steradians),
and wavelength nodes with quadrature weights (metres), filter transmission and
detector quantum efficiency (both dimensionless). Entries can share a volume and
represent different filters or pulses. Spatially varying density and temperature
are integrated before forming the signal. The volume length is in metres.

Call `init(volumes,channels,prior,status,aligned,gas)`; zero status means success.
The optional `aligned` flag has the existing `thomson_t` meaning: true fixes the
alignment, false exposes its two displacement parameters. The physical optical
calibration vector has `3+2*n_channels` entries in this order:

| Entries | Parameter | Units |
|---|---|---|
| 1 | Absolute calibration scale | 1 |
| 2 | Laser out-coupling scale | 1 |
| 3 | Common optical transmission | 1 |
| 4, 6, ... | Channel spectral-response scale | 1 |
| 5, 7, ... | Channel integral gain | V s / photoelectron |

Every entry is mapped from a standard Gaussian latent through the supplied W0
prior, including correlations. Alignment parameters precede these latents.
The map is affine; callers must keep the physical parameters positive. Shared
scale factors are degenerate without informative calibration records or priors.
Measured pulse energies and supplied response curves are inputs; uncertainties
in their absolute scale belong in the out-coupling and response priors. Uncertain
response **shape**, wavelength, geometry and pulse-to-pulse correlated energy
errors require an adapter or an extension beyond these scalar scales.

`predict`, `jvp` and `vjp` implement the existing diagnostic contract in V s.
`photoelectrons`, `photon_jvp` and `photon_vjp` provide the corresponding count
map. JVP and VJP mean Jacobian-vector and vector-Jacobian products. The density
and temperature rows and every calibration chain are differentiated; positional
alignment chains through the existing provider sampling contract. Independent
channels use OpenMP. A vacuum sample with both density and temperature zero
contributes zero; that boundary has zero derivative by convention.

## Absolute Rayleigh calibration record

Supplying a `gas` prior to `init` selects an elastic Rayleigh calibration
measurement with the same optics and gains. Its three physical parameters are
pressure (Pa), temperature (K), and dynamic polarizability volume (m^3).
Their three Gaussian latents follow the optical latents. `elastic_response`
is the filter transmission times quantum efficiency at the laser wavelength;
it is supplied explicitly because a notch cannot safely be inferred from sparse
quadrature nodes. The common optical and channel response scales still apply.

This is the polarized classical dipole law for an isotropic gas, with the ideal
gas density. It supports an argon-like monatomic calibration, not nitrogen's
anisotropic Rayleigh or rotational Raman spectrum. Polarizability is a parameter
with a prior, not a built-in gas constant. The argon reference test uses
[Sneep and Ubachs (2005), Table 2](https://doi.org/10.1016/j.jqsrt.2004.07.025).
Share optical latents between plasma and calibration records in the inference
problem; the gas record does not identify all optical factors separately.

## Trace reduction and noise

`tiago_thomson_reduction::pulse_integral` constructs an existing
`linear_reduction_t` from time stamps (seconds), inclusive pulse and baseline
endpoint indices and optional extra variance (V^2 s^2). It uses trapezoidal
weights and subtracts a time-averaged baseline. Gates cannot share an interval.
The existing reduction propagates full input covariance, derivatives and
provenance from a `signal_window_t`; output units are V s. Windows may come from
any W0 source and can be reduced before the next window is read.

`thomson_noise(photons,bg,gain,digitvar,excess,bgscale,variance,row)` returns the
pulse variance and its six partial derivatives in argument order. `bg` is the
expected in-gate background count, `bgscale` is pulse duration divided by an
independent baseline duration, `excess >= 1` is the detector excess-noise factor,
and `digitvar` is the already-integrated digitiser variance. The row gives JVP
by a dot product and VJP by multiplication with the output seed. This model
assumes stationary Poisson background and independent baseline counts. Use the
trace's full covariance for correlated electronics or shaped pulses. Add only
noise contributions absent from that covariance to avoid double counting.

## Verification

`tiago_thomson` checks the cold Gaussian width, independent velocity-space
reference values, unit spectral area, ideal-filter temperature ratios, vacuum,
nonuniform profiles, every optical calibration column, photon and voltage
products, gas calibration, trace covariance and noise derivatives. A synthetic
joint plasma/Rayleigh inversion checks biased spectral and absolute scales
against a local Gaussian posterior. One synthetic realization does not establish
frequentist coverage, and no experimental validation is claimed.

The [FortSym](https://github.com/itpplasma/fortsym) generator and all pointwise
kernel derivatives are committed together. See [generation instructions](../gen/README.md)
for the optional regeneration CTest check. Ordinary builds need no generator.
