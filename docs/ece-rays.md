# Cold electron-cyclotron rays

`tiago_ece_ray` supplies a device-independent cold-plasma ray primitive with
forward and reverse sensitivities. It does not yet predict a radiometer signal.
The existing `tiago_ece_radiometer::ece_t` remains the cold-resonance profile
diagnostic. Standalone [thermal transfer and receiver components](ece-transfer.md) now
supply prescribed-opacity transport, reciprocal walls, frequency/antenna
averaging and calibrated voltage/noise. Relativistic emission/absorption,
their coupling to the rays and profile sampler, and posterior validation
are still needed for a complete electron cyclotron emission (ECE) diagnostic.

## Medium and coordinates

A caller implements `cold_medium_t%n_params()` and `%sample`. The count is a
pure function. At a Cartesian point `x(3)` in metres, `sample` supplies:

| Output | Shape | Meaning |
| --- | --- | --- |
| `q` | `(4)` | Cartesian magnetic field components in tesla, then electron density in m^-3 |
| `dq` | `(3,4)` | Spatial first derivatives of `q`, per metre |
| `d2q` | `(3,3,4)` | Spatial second derivatives of `q`, per metre squared |
| `qp` | `(4,n_params)` | Parameter derivatives of `q` at fixed `x` |
| `dqp` | `(3,4,n_params)` | Parameter derivatives of `dq` at fixed `x` |

The medium is immutable and its callback must be safe for concurrent calls.
Every output must be initialized, including zeros. Parameter derivatives use
the caller's physical units. These are field values and derivatives, not
coefficients of a second equilibrium or profile representation.

There is no VMEC adapter for this primitive yet. The existing profile sampler
supplies the field magnitude and its first derivatives; it does not supply the
field vector or second spatial derivatives required here. It would be incorrect
to identify the magnetic field direction with the ray direction, or to use the
exterior virtual-casing field inside the plasma. The existing VMEC flux map and
kinetic profile parameterization should supply a future adapter's quantities.

## Hamiltonian and launch

The collisionless, electron-only Appleton-Hartree relation is used, with fixed
ions. Let `N = c k / omega` be the dimensionless refractive-index vector,
where `c` is the speed of light, `k` the wavevector and
`omega = 2 pi frequency`. The electron plasma and cyclotron angular
frequencies are `omega_pe` and `omega_ce`. Define `X = omega_pe^2 / omega^2`,
`Y = omega_ce / omega` and `mu = N.B/(|N| |B|)`. The Hamiltonian is
`H = (|N|^2 - n_mode^2(X,Y,mu))/2`. O and X denote the ordinary and
extraordinary branches, respectively; at perpendicular propagation their
squared indices reduce to `1-X` and `1-X(1-X)/(1-X-Y^2)`.
See Stix, *Waves in Plasmas* (1992), chapter 1, and the cold-ray limit of
[Marushchenko et al., Computer Physics Communications 185, 165 (2014)](https://doi.org/10.1016/j.cpc.2013.09.002).

`trace_cold_ray` takes the observing frequency in Hz, the mode (`ECE_O_MODE`
or `ECE_X_MODE`), an initial six-vector `(x,N)` and a positive Hamiltonian
integration duration. Launch on the dispersion surface: for a unit direction
`d`, evaluate `cold_hamiltonian(d,B,ne,f,mode,...)`, then set
`N = sqrt(1-2*H)*d` when the radicand is positive. The launch guard requires
`abs(H) <= 1e-7`; accepted states are also checked for Hamiltonian drift.

Hamilton's equations are `dx/dt = H_N`, `dN/dt = -H_x`, where subscripts
denote partial derivatives. Here `t` is a
Hamiltonian parameter with units of metres, not elapsed time or arc length.
A later transport calculation must use `ds = |H_N| dt`. The adaptive
Cash-Karp integrator comes from fortnum; no numerical integrator is duplicated.
Turning and reflection arise from these equations without a discontinuous
specular reset. A normal cutoff is regular in the unmagnetized limit; an
oblique O-mode slab turns before its normal-incidence cutoff.

The supported rays remain on the fast electromagnetic branch `|N| <= 1`
(with numerical allowance `max(1e-8,20*rtol)`). A vacuum ray has `|N|=1`.
The intended magnetized domain is below the ordinary density cutoff (`X<1`),
away from resonances and branch degeneracies. Singular stages, nonfinite
medium data, negative density, off-shell launches and slow branches are errors.
The unmagnetized kernel provides the regular `B=0` limit, but derivatives that
turn on a magnetic field at a field null are not defined for the individual
magnetoionic branches. Exact zero-index magnetized stages are rejected.
There is no tunnelling, mode conversion, ion response or collision model.

`trace_cold_rays` evaluates independent rays with OpenMP. It returns a status
per ray and the first error as the aggregate status. An empty batch succeeds.
Check status and `ray%valid` before using the trace. Sensitivities also require
`ray%derivatives_valid`; a field null with active magnetic-field parameters or
spatial derivatives is rejected for differentiation.

## Derivatives and verification

`cold_ray_jvp` propagates tangents of the initial six-vector, medium parameters
and frequency to the final six-vector. `cold_ray_vjp` applies its transpose.
The medium must remain at the nominal parameters used to compute the trace.
The generated Hamiltonian gradient and Hessian are composed with the medium
jets. Fortnum differentiates the accepted Runge-Kutta schedule with that
schedule held fixed, including the ray path. Duration, tolerances, mode and
step acceptance decisions are inactive. These are endpoint derivatives;
there is no terminal-event sensitivity or transport derivative in this API.

`tiago_ece_ray` checks the perpendicular dispersion limits, straight vacuum
rays, an analytic linear slab with a turning point and reflected path, the
normal unmagnetized cutoff, a cylindrical parabolic-density oscillator,
oblique O/X rays against central finite differences, the adjoint identity and
invalid launches. The synthetic magnetic field is divergence-free and spatially
varying. `tiago_ece_codegen` regenerates the kernels and compares them with the
committed output when `TIAGO_FORTSYM_BUILD` is set; otherwise it is skipped.
The independent-ray benchmark is `benchmark_ece_ray`.

These tests verify the cold-ray primitive. They do not validate relativistic
optical depths, radiation temperatures, the full VMEC inference chain or
experimental measurements.
