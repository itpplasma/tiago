# Soft X-ray cameras

`tiago_soft_xray::soft_xray_t` predicts the signals of soft X-ray (SXR)
pinhole cameras: line integrals of a filtered emissivity along the chord each
detector sees through its pinhole. It implements the
[profile-diagnostic contract](profile-diagnostics.md), so it runs in the
[joint inference](recorded-joint-inference.md), in
[`tiago_synth`](synthetic-forward.md) and on any provider state, including
point tables with island plateaus.

## Signal

Channel `c` records

```text
y_c = G_c / (4 pi) int eps dl        [W]
```

along the straight chord that starts at the pinhole and runs along
detector -> pinhole. `eps` [W m^-3] is the emissivity after the camera's filter
and detector absorption, taken as isotropic and zero outside the last closed
flux surface. `G_c` [m^2 sr] is the etendue. Without an etendue the code uses
`G_c = 4 pi`, so `y_c = int eps dl` [W m^-2].

`pinhole_camera_t` describes one camera: pinhole centre, area and normal,
detector centres, areas and normal, chord length from the pinhole and filter.
`camera_chords` returns its chords and etendues
`G = A_d A_a cos(theta_d) cos(theta_a) / d^2` (detector area `A_d`, pinhole
area `A_a`, their distance `d`, angles between the chord and the normals),
the usual approximation for pinhole and detector small compared with `d`
([Ingesson et al., Fusion Sci. Technol. 53, 528 (2008)](https://doi.org/10.13182/FST53-528)).
`set_cameras(cameras, npanel)` sets up all cameras of one diagnostic; the
chords get `npanel` panels of four Gauss points each.

Only the central chord of each channel is integrated. The finite width of the
viewing cone, detector nonlinearity, time response and reflections are not
modelled.

## Emissivity models

`model = SXR_PROFILE_MODEL` (default): `eps = exp(u(s))` with `u` a cubic
B-spline in the flux label `s`, the parametrisation and Gibbs/Matérn
prior of the kinetic profiles (`tiago_kinetic_profiles`).
`set_emissivity_basis(breakpoints)` chooses the basis,
`set_emissivity_prior(mean, sigma_f, l_core, l_edge, s0, sw[, kernel])` makes
the profile latent with prior mean values at the Greville abscissae
(`emissivity%greville`), and `set_emissivity_table(s, eps, ierr)` fixes it to a
least-squares fit of a table. The latent coordinates are standard normal and
are this diagnostic's nuisance parameters `eta`; in a joint model they need no
further prior setup. This model makes no assumption about the radiation
mechanism. Use it when line and recombination radiation of impurities
contribute to the filtered signal.

`model = SXR_BREMSSTRAHLUNG_MODEL`: free-free emission from the shared ne, Te
and Zeff,

```text
eps = C ne^2 Zeff Te^(-1/2) int t(E) exp(-E / Te) dE,    C = 1.69e-38 W m^3 eV^(-1/2),
```

with ne in m^-3 and Te, E in eV. `C` is the total bremsstrahlung coefficient of
the NRL Plasma Formulary (J. D. Huba, Naval Research Laboratory, 2019,
radiation section), which assumes a Gaunt factor of one; without a filter (`t = 1`) the integral
gives `C ne^2 Zeff Te^(1/2)`. The filter transmission `t(E)` (`sxr_filter_t`,
one per camera) is zero below its first energy node, linear between nodes and
constant above the last one, and should include the detector's absorption
efficiency. The integral is evaluated in closed form per interval, with a
series for intervals narrower than `0.01 Te`. Line and recombination radiation
are not included, and a Gaunt factor different from one only rescales the
signal, which the camera gain absorbs.

In a field with islands the flux label is the provider's: a
[point table](synthetic-forward.md#equilibria-on-a-grid-and-islands) with
island plateaus, or one whose `rho` comes from
[field-line labels](field-line-labels.md), makes the emissivity flat across
each island in both models.

## Calibration and joint inference

Detector responsivity and amplifier gain are calibration parameters of the
recorded signal. Wrap the diagnostic in `calibrated_diagnostic_t` with one
log-gain group per camera:

```fortran
call sxr%set_cameras(cameras, 100)
call sxr%set_emissivity_basis([0.0_dp, 0.3_dp, 0.6_dp, 1.0_dp])
call sxr%set_emissivity_prior(mean, 0.5_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
call prior%init(log(gain_mean), gain_covariance, status)   ! one entry per camera
call recorded%init(sxr, prior, status, groups=sxr%camera, log_gains=.true.)
allocate (model%diags(k)%d, source=recorded)
```

The prior comes from the cameras' calibration record. With the emissivity-profile
model, a common factor of all gains is not separable from the emissivity scale;
the data fix the gain ratios between cameras. With the bremsstrahlung model a
common gain trades against the density and temperature scale, which Thomson
scattering, interferometry and their calibrations determine.

The JVP and VJP with respect to the sampled quantities (s, ne, Te, Zeff) and the
emissivity latents are exact, so the joint model carries the derivatives
through the equilibrium adjoint like the other line integrals. Emission stops
at the last closed flux surface; the derivative does not include the motion
of that boundary across a quadrature point, which matters only when the
emissivity there is not small.

## Synthetic signals

`tiago_synth` evaluates an `sxr` group of the geometry file (see
[synthetic forward driver](synthetic-forward.md)): `detector(channel, xyz)`,
`aperture(channel, xyz)`, `length(channel)` [m] from the aperture, optional
`etendue(channel)` [m^2 sr], `camera(channel)` (1, 2, ...), and per camera a
filter `filter_energy(filter, node)` [eV] and `filter_transmission(filter, node)`
whose first index is the camera number. If the profile table has
`sxr_emissivity(s)` [W m^-3], the emissivity-profile model uses it; otherwise the
bremsstrahlung model uses the table's ne, Te and Zeff.

## Verification

- `tiago_soft_xray`: filter response against Simpson sums over photon energy
  (1e-10) for Te from 30 eV to 20 keV, the step filter in closed form, and the
  series and closed-form kernels where they meet; the etendue by hand; on
  concentric circles a Gaussian emissivity profile (closed form with `erf`,
  1e-12), unfiltered bremsstrahlung (closed form, 1e-8) and filtered
  bremsstrahlung (Simpson sums, 1e-8); contract derivatives (JVP within 7e-8,
  adjoint 1e-14), the chain through the provider parameters against central
  differences (3e-8 on LI383) and its transposition identity, and per-camera
  log gains. `tiago_soft_xray_vmec` repeats the derivative checks on LI383.
- `tiago_point_table`: an SXR chord through an island with an emissivity
  plateau agrees with the hand integral to 1e-4, within the grid's
  interpolation error.
- `tiago_synthetic_forward`: the `sxr` group against Simpson sums (1e-10).
- `tiago_joint_soft_xray`: two calibrated cameras, one per model, in the joint
  model of the recorded-signal test: complete directional derivative through
  re-solved VMEC++ equilibria against central differences (4e-8 relative),
  adjoint identity (2e-17), and a local fit that recovers the injected gain ratio
  of the bremsstrahlung cameras.
- [STELLOPT benchmark](../benchmarks/sxr/README.md): STELLOPT's own SXR signal
  (`chisq_sxr`, emissivity `Zeff ne^2 Te^(1/2)`) and TIAGO agree to 9e-5 of the
  largest signal on 36 chords of the LI383 equilibrium.

![TIAGO and STELLOPT soft X-ray signals of 36 chords on LI383 and their
difference](figures/sxr-stellopt.svg)

DIAGNO, the STELLOPT magnetic diagnostic code that TIAGO's magnetics are
benchmarked against, has no emissivity or soft X-ray model; the reference here
is STELLOPT's optimizer target.
