# Recorded magnetic signals

Link `tiago_inference` for the magnetic response and recorded integrator model;
link `tiago_signals` for drift reduction. Geometry and calibration records are
supplied at run time. All examples and tests are synthetic.

## Forward model and parameter order

`tiago_magnetic_response::magnetic_response_t%init(solver,sensors,rule)` reuses
`tiago_sensor_signals` for flux loops, diamagnetic loops, segmented Rogowski
coils and B-probes, including turns and sensor-set combinations. It implements
`profile_diagnostic_t`. There are no additional real-space sample points:
vector potential, vector field and linked-flux evaluation remain in the existing
magnetic solver, rather than being represented as scalar kinetic samples.

At fixed coil and boundary geometry the magnetic signals are exactly affine
in coil currents, boundary-field coefficients and linked toroidal flux. The
adapter's `eta` contains **physical increments about the supplied solver** in
that order: one per coil group, one per VMEC boundary-field coefficient, then
one linked-flux increment in webers. The field-coefficient order is the existing
`mode_column_name` order. External field sources are held fixed and supply no
field-coefficient columns. The JVP (Jacobian-vector product) and VJP (its
transpose action) use the existing exact response matrices. The zero increment
returns the solver's signals, including its existing linked flux.

Reinitialize after changing geometry. This block does not approximate shape
changes by an affine response; use the existing shape-response/provider chain
when fitting geometry. An external source without derivative support cannot
be varied through this block. Do not introduce an independent coefficient
prior if these coefficients are already outputs of an equilibrium provider.

`tiago_integrator_diagnostic::integrator_diagnostic_t%init(model,prior,time,origin,status)`
wraps this response, or another profile diagnostic, as a
`recorded_diagnostic_t`. `time` has one timestamp in seconds per output;
`origin` is the time reference in the same convention. The map is

```text
recorded(i) = gain * model(i) + offset + drift * (time(i) - origin).
```

The `calibration_prior_t` has three correlated physical parameters in the
order gain, offset, drift. Their units are recorded/model units, recorded units,
and recorded units per second. Standard-normal calibration latents follow the
wrapped model's `eta`; all sampled-field, alignment and model-parameter actions
are composed exactly. Use separate wrappers for independently calibrated
channel groups. Do not share a calibration across incompatible units.
The wrapper copies its model: refresh that copy after a provider/geometry update.

This is the signal **after integration**, with the configured integrator gain
and polarity. It does not differentiate the signal to a pickup voltage or model
electronic bandwidth, eddy currents or saturation. Instrument-specific transfer
functions can be separate recorded diagnostics.

## Baseline drift and covariance

`tiago_integrated_drift::integrated_drift_t%init(input,baseline,keep,origin,id,status)`
constructs a `reduction_step_t` for a bounded integrated-signal window.
`baseline` and `keep` are disjoint, unique, one-based sample indices; `keep`
is increasing. At least two baseline times must differ. The physical signal
must be zero on the baseline; any physical baseline must first be predicted and
accounted for consistently in both observations and predictions.

The baseline offset and linear drift are fitted by generalized least squares
using their full covariance. Applying the step subtracts that fitted line and
returns only the `keep` samples, with their original timestamps. With `W` the
resulting linear operator, the covariance is exactly `W C W^T`, including the
cross-correlations between the baseline and kept samples. Adding only the
fitted-drift variance would double-count common offset noise and omit cross
terms. This is linear covariance propagation as in
[JCGM 100:2008, section 5.2](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf).
`baseline_covariance()` returns the offset/drift fit covariance at `origin`.

The baseline covariance and two-parameter normal matrix must be positive
definite. The operator, times and units are fixed at initialization; changed
times or units are rejected. A new input covariance is propagated with the
same fixed operator; reinitialize to refit the statistically optimal weights.
Use an `id` version that identifies the baseline selection, kept indices and
reference time. The output retains the source record and parent checksum.

Dense window covariance and reduction weights are intentional here: select a
bounded window or reduce incoming chunks first. This is not a full-record
streaming regression. Apply the same reduction to observations and predictions,
compose its JVP/VJP, and use the returned covariance once in the likelihood.
A fitted baseline removes two calibration directions; do not also infer them
from that reduced trace without an independent calibration constraint.

## Diamagnetic and kinetic energy consistency

`tiago_energy_consistency::energy_consistency_t` is a profile diagnostic whose
single output is `Wdia - Wkin` in joules. Supply real-space `points(3,n)`,
fixed volume quadrature `volume_weights(n)` in cubic metres and the impurity
charge `zimp`, as in the existing two-ion kinetic-profile closure. Its `eta`
is `(delta_flux [Wb], B0 [T], R0 [m], paramagnetic_flux [Wb])`.
These are physical quantities, not standardized calibration latents; propagate
their provider/calibration derivatives and priors through the caller.

```text
Wdia = 3 pi R0 B0 (paramagnetic_flux - delta_flux) / mu0
Wkin = integral (3/2) e [ne Te + (n_D + n_Z) Ti] dV.
```

`delta_flux` is the compensated, vacuum-subtracted, turn-normalized change
in flux, negative for diamagnetism with positive `B0`. The paramagnetic flux
correction is explicit and has to come from the chosen equilibrium model.
The relation is the circular-cylinder, low-beta, large-aspect-ratio limit of
[Pustovitov, NIFS annual report 2010, p. 64, equation (1)](https://www.nifs.ac.jp/report/annrep10/pdf/064.pdf),
with isotropic thermal energy and toroidal length `2 pi R0`. It is a consistency
check with declared model error, not a universal three-dimensional energy
reconstruction. Finite-aspect-ratio effects, boundary motion, anisotropy and
fast particles require a provider-specific correction/model. The kinetic
pressure and its derivatives reuse TIAGO's generated `kinetic_pressure`.
Temperatures are in electronvolts, densities in inverse cubic metres; `e` is the
exact elementary charge. Vacuum permeability `mu0` uses the
[published CODATA 2018 value](https://physics.nist.gov/cuu/pdf/wall_2018.pdf)
`1.25663706212e-6 H/m`.

Only the fixed quadrature is differentiated here. A varying integration volume
requires its geometry derivative outside this diagnostic. Compare against zero
with an explicit consistency uncertainty, retaining measurement, calibration
and model error separately. Do not count the magnetic measurement twice by
adding an independent energy datum made from that same measurement.

## Verification and generation

`tiago_recorded_magnetics` checks the circular loop in a uniform field,
Rogowski and B-probe limits, the linked-flux term, vacuum current cancellation,
correlated calibration recovery, coefficient responses and finite differences.
The pressure-cylinder control independently uses
`Bz(r)^2 + 2 mu0 p(r) = B0^2`, `p(r)=p0(1-r^2/a^2)` and checks the low-beta
energy relation. `tiago_integrated_drift` checks drift removal, independent
endpoint covariance, correlated baseline propagation, derivatives and invalid
baseline selections. These tests validate no experimental calibration.

All new scalar maps and derivative rows are generated with FortSym by
`gen/gen_magnetics.f90`. Configure `-DTIAGO_FORTSYM_BUILD=/path/to/fortsym/build`
to register `tiago_magnetic_codegen`, which regenerates into an isolated build
directory and compares all three files with the checked-in sources. The
standalone generator build follows [gen/README.md](../gen/README.md).
`benchmark_magnetics` measures the recorded calibration and energy prediction,
JVP and VJP on one million synthetic samples, with OpenMP channel/sample loops.

Chris&AI
