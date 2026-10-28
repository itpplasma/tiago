# Perpendicular thermal electron-cyclotron opacity

`tiago_ece_opacity` evaluates weakly relativistic Maxwellian absorption and
emission for the first and second harmonics of the ordinary (O) and
extraordinary (X) modes. It supplies the absorption to the existing
[thermal transfer and receiver components](ece-transfer.md); prescribed
absorption remains available through `thermal_transfer`.

These are perpendicular, near-resonance, weak-damping approximations. They
have no parallel wave-vector argument and must not be used as an oblique
Doppler-broadened model. The [recorded composition](ece-recorded.md) joins them to the
[cold-ray solver](ece-rays.md) and the profile-diagnostic adapter, with an
explicit small-angle guard.

## Local interface and units

```fortran
call ece_opacity(te, ne, bfield, frequency, mode, harmonic, &
    alpha, emissivity, index2, jacobian, status, x1_regime)
```

Inputs are electron temperature `te` in eV, electron density `ne` in inverse
cubic metres, positive magnetic-field magnitude `bfield` in tesla, and
`frequency` in Hz. `mode` is `ECE_O_MODE` or `ECE_X_MODE`; `harmonic` is 1 or 2.
Outputs are absorption `alpha` in inverse metres, emissivity in
W m^-3 sr^-1 Hz^-1, and the squared real propagation index `index2`.
`jacobian(3,4)` contains their derivatives in that output order with respect
to `(te,ne,bfield,frequency)`. Only status zero makes the outputs usable.
`ECE_OPACITY_INPUT` denotes malformed input; `ECE_OPACITY_DOMAIN` denotes an
excluded physical regime or a nonfinite intermediate result.

Frequency is cycles per second, not angular frequency. The emissivity is
per single polarization and per Hz. In the Rayleigh–Jeans limit the vacuum
blackbody intensity for that polarization is `e*te*frequency**2/c**2`, where
`e` is the elementary charge and `c` the vacuum speed of light. The generated
emissivity obeys Kirchhoff's law:

```text
emissivity / alpha = index2 * e*te*frequency**2/c**2.
```

This emission approximation assumes photon energy is much less than thermal
energy; it is not a Planck-law model at quantum temperatures.

The batch interface accepts one value of each local input, mode and harmonic
per sample, returns `(3,4,n)` derivative rows and one status per sample, and
parallelizes the sample loop with OpenMP. The optional X1 regime is shared
across the batch. Physical constants are fixed published SI constants;
instrument gain, offset, wall reflectivity, bandpass and noise parameters
remain explicit inputs of the existing transfer/receiver components.

## Approximations and domain

The formulas follow [Bornatici et al., Nuclear Fusion 23, 1153 (1983)](https://doi.org/10.1088/0029-5515/23/9/005):

| Mode | Formula | Propagation approximation |
| --- | --- | --- |
| O1 | Eqs. (3.1.25), (3.1.40a) | Hermitian thermal index and weak damping |
| O2 | Eq. (3.1.41), harmonic 2 | Leading cold index near the second harmonic |
| X1, tenuous | Eqs. (3.1.15a), (3.1.28) | Unit index |
| X1, finite density | Eqs. (3.1.30)–(3.1.32) | Finite-Larmor-radius correction and leading index |
| X2 | Eqs. (3.1.17)–(3.1.20), (3.1.35)–(3.1.38) | Hermitian thermal root and real polarization |

Let `theta=te/(m_e*c**2/e)`, with `m_e` the electron mass, and let
`r=(omega_p/omega_c)**2`, with electron plasma and cyclotron angular
frequencies `omega_p` and `omega_c`. The real resonance coordinate is
`z=(1-harmonic*omega_c/(2*pi*frequency))/theta`. Fortnum evaluates the
upper-lip Dnestrovskii functions of orders 5/2 and 7/2, including their
real parts on both sides of resonance. Their imaginary parts vanish for
`z>=0`; the propagation index remains available there.

The interface requires `1e-12 <= theta <= 0.01` and relative detuning
`abs(theta*z)<=0.1`. O1 requires `r<1`, O2 `r<4`, and X2 `r<2`.
The default `ECE_X1_FINITE` requires `20*theta <= r < 2`;
`ECE_X1_TENUOUS` requires `r <= 0.2*theta`. The intermediate density regime
is rejected rather than interpolated between incompatible asymptotic models.
The X1 finite-density approximation therefore excludes vacuum; select the
tenuous regime for its vacuum limit.

Nonpositive propagation indices and strong damping are rejected. The
weak-damping checks require `alpha*c/(4*pi*frequency*sqrt(index2))<=0.1` and,
for O1/X2, `abs(r*Im(F_7/2))<=0.1`. These conservative guards delimit the
implemented approximations; they are not an accuracy certificate for every
accepted point. Cutoff tunnelling, mode conversion, oblique propagation,
non-Maxwellian distributions and general complex dispersion are excluded.

## Composition with transport

```fortran
call thermal_ece_transfer(te, ne, bfield, frequency, mode, harmonic, length, &
    rho, wall_te, trace, status, x1_regime)
```

Cell arrays run observer first, as in `thermal_transfer`. `length` is in
metres, `rho` is the supplied wall power reflectivity, and `wall_te` is the
wall radiation temperature in eV. The result contains the existing transfer
radiation temperature and optical depth. Its `jvp` takes perturbations of
`te(:),ne(:),bfield(:),frequency,length(:),rho,wall_te`; `vjp` returns the
transpose in that order. The methods overwrite their outputs. The local
absorption derivative rows feed the existing transfer tape. No transport
physics or receiver calibration is duplicated.

Transport evolves radiation temperature, proportional to intensity divided
by the squared index. Its local thermal source is therefore `te`, consistently
with the emissivity above. This composition holds the sampled path fixed.
The caller must propagate ray motion, displaced profile samples and path-length
changes when connecting a geometry-dependent diagnostic. The field input is
a magnitude, not a vector-field provider. Receiver and calibration derivatives
remain available on the existing receiver tape.

## Verification and reproducibility

`test_ece_opacity` checks O1, O2, tenuous X1 and X2 slab optical depths against
the perpendicular tenuous limits in Bornatici's Table IV, p.1192. With
`L_B=B/abs(dB/ds)` at resonance and `theta<<1`, these are respectively
`(pi/2)*r*theta*omega_c*L_B/c`, `pi*r*theta**2*omega_c*L_B/c`,
`(pi/2)*r*omega_c*L_B/c`, and `pi*r*theta*omega_c*L_B/c`.
Refining temperature toward zero verifies the asymptotic limit rather than
requiring a finite-temperature model to equal its limiting formula.

The same test checks vacuum, the cold off-resonance limit, Kirchhoff's law,
every local derivative column on both resonance wings, coupled transport
finite differences and its transpose identity. The existing transfer tests
cover optically thin/thick limits, wall reflections and receiver calibration.
All fixtures are synthetic or evaluations of published equations.

The X1 finite-Larmor-radius correction is checked against 48 independent
90-digit principal-value quadratures of the defining functions, produced by
`gen/reference_ece_opacity.py`. There is a numerical inconsistency between
Bornatici's Eq. (3.1.30) and Table III: at `z=-3,r=1.9`, the equation gives
`0.24956268465884969`, while the table prints `0.2909`. Independent defining
integrals and the recurrence simplification give the former value. The test
uses the defining equation, not that inconsistent table entry; no agreement
with Table III is claimed.

`gen/gen_ece_opacity.f90` generates all ten local physics kernels and their
partial rows with FortSym. Configure `TIAGO_FORTSYM_BUILD` to include
`tiago_ece_opacity_codegen`, which regenerates and byte-compares each file.
`benchmark_ece_opacity` measures the OpenMP batch including emissivity,
propagation index, all derivative rows and status checks. These verifications
do not establish finite-angle accuracy or experimental signal recovery.
