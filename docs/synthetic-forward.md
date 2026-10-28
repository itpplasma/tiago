# Synthetic forward driver

`tiago_synth <input.nml>` evaluates every configured diagnostic once, from a
VMEC wout, tabulated kinetic profiles and real diagnostic geometry, and writes
one netCDF-4 file with one group per diagnostic. It adds no physics: profile
instruments are the existing `interferometer_t`, `thomson_t` and `xics_t` on
the VMEC state (`vmec_state_t`), and magnetics are the sensor-set evaluators
with the virtual-casing plasma model.

```text
&synth
  wout     = 'wout.nc'          ! required
  profiles = 'profiles.nc'      ! ne, te, ti (and zeff) on s; omitted: profiles stay at their defaults
  geometry = 'geometry.nc'      ! one group per profile instrument
  sensors  = 'machine.sensors'  ! magnetic sensor set (docs/sensor-sets.md)
  coils    = ''                 ! coil file; empty: plasma contribution only
  extcur   = ''                 ! coil currents for the coil file
  output   = 'synthetic.nc'
  plasma_sheet = 'contravariant' ! or covariant, conservative (docs/magnetics.md)
  npanel = 200, thomson_points = 4, ece_points = 4000, profile_knots = 12
/
```

## Inputs

All lengths are in metres, Cartesian (x, y, z). Dimension order is the one
written by netCDF4-python or `ncgen`; `label(n, strlen)` is an optional
character variable in every group.

| group | variables |
|---|---|
| `interferometer` | `start(chord, xyz)`, `end(chord, xyz)` |
| `xics` | `start(chord, xyz)`, `end(chord, xyz)` |
| `thomson` | `center(volume, xyz)`, `direction(volume, xyz)` (laser), `length(volume)` |
| `ece` | `frequency(channel)` [Hz], `harmonic(channel)`, `origin(channel, xyz)`, `direction(channel, xyz)`, `length(channel)` |

The profile table has one dimension `s` (normalized toroidal flux in [0, 1])
and the variables `ne` [m^-3], `te` [eV], `ti` [eV] and optionally `zeff`,
all positive (1 < zeff < 6). Each profile is fixed to the least-squares fit of
the cubic B-spline of `tiago_kinetic_profiles` (`profile_knots` intervals) to
the logarithm of the table (logit for Zeff), so a table needs at least
`profile_knots + 3` points. A profile whose logarithm is a cubic in s is
reproduced exactly. Beyond s = 1 the profiles decay as described in
`docs/profile-diagnostics.md`; profiles given in r_eff are converted by the
caller with s = (r_eff / a)^2, a = `Aminor_p` of the wout.

Missing groups are skipped; malformed groups, NaN values, non-positive
profiles and missing variables stop the run with a message naming the file,
group and variable. Nothing is filled in.

## Output

Global attributes: `tiago_revision` (`git describe` of the build),
`aminor_p`, `plasma_sheet`, and for every input `<name>_path` and
`<name>_sha256`. Every group has a `label` variable and these signals:

| group | variables |
|---|---|
| `interferometer` | `line_density` [m^-2] = int ne dl along the chord |
| `thomson` | `ne` [m^-3] (volume average), `te` [eV] (density weighted), `s`, `r_eff` [m], `inside` at the volume centre |
| `xics` | `brightness` (relative, unit scale), `ti` [eV] (emissivity-weighted along the chord) |
| `ece` | `te` [eV] at the first cold resonance n f_ce(abs(B)) = f inside the plasma along the line of sight, with its `s`, `r_eff`, `R`, `phi`, `Z`; NaN if there is none |
| `magnetics` | `signal` (flux loops [Wb], Rogowski coils and B-probes [T m^2]), `kind` (sensor kind code) |

The ECE value is the electron temperature at the cold resonance; it neglects
relativistic downshift and optical depth. The radiation temperature of the
recorded-ECE model (`docs/ece-recorded.md`) needs a near-perpendicular line of
sight within its opacity angle bound and calibrated receiver priors, and is not
part of this driver.

## Example

`examples/synthetic/` holds a public example on the low-resolution LI383
equilibrium of the Simsopt test files (the wout of the
`tiago_download_wout_li383` test fixture):

```sh
cd examples/synthetic
ncgen -4 -o geometry.nc geometry.cdl
ncgen -4 -o profiles.nc profiles.cdl
cp <build>/tests/cases/wout_li383_low_res_reference.nc .
<build>/tiago_synth li383.nml
```

## Verification

`tiago_synthetic_forward` checks the driver end to end, through the geometry
and profile files and the output groups: on concentric circles
(s = r^2, abs(B) = B0 R0 / R) the chord integrals of ne and of the XICS
emissivity and its Ti moment agree with independent Simpson quadratures to
1e-11, Thomson points and a finite volume reproduce the tabulated profiles to
1e-11, and the ECE cold resonance lies at R = B0 R0 / B_res; on the LI383 wout,
Thomson points placed at X(u) for known labels u give s = u_1^2 and the
tabulated ne and Te to 1e-9. It also checks the SHA-256 implementation
against the FIPS 180-4 examples. `tiago_synth_example` runs the example above.
