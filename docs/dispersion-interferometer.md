# Dispersion interferometer recorded signals

`tiago_dispersion_interferometer::dispersion_interferometer_t` predicts
interleaved lock-in signals `(I_1,Q_1,I_2,Q_2,...)` from the existing
line-integral diagnostic and the existing [recorded-diagnostic calibration prior](recorded-signals.md).
It links with `tiago_inference`; no device geometry or archive is built in.

For electron column `N = integral ne dl` in m^-2 and fundamental vacuum
wavelength `lambda` in metres, the adopted positive-electron phase is

```text
phase = (3/2) re lambda N - (4 pi/lambda) neutral_scale D
        + env_gain E + phase0
I = offset_I + amplitude_I cos(phase)
Q = offset_Q + amplitude_Q sin(phase)
```

Here `re = 2.8179403205e-15 m` is the
[CODATA 2022 classical electron radius](https://physics.nist.gov/cgi-bin/cuu/Value?re),
`D = integral [n_neutral(omega)-n_neutral(2 omega)] dl` is the neutral
index-difference column in metres, and `E` is a supplied environmental proxy
with a declared runtime unit. `omega` is the fundamental angular frequency.
The positive-electron convention is `phi(2 omega)-2 phi(omega)`, where `phi`
is propagation phase. It agrees with the positive density coefficient in
Bagryansky et al., *Review of Scientific Instruments* 77, 053501 (2006),
[Eq. (3)](https://doi.org/10.1063/1.2202922).
Common vacuum path changes cancel exactly before code emission; neutral
dispersion does not cancel. The common-path principle is also described by
[Knauer et al., EPS 2016, P4.017](https://scipub.euro-fusion.org/wp-content/uploads/WPS1CP16_16154_submitted.pdf).
No machine layout or calibration is taken from these references.

This is the dilute, cold, unmagnetized straight-ray phase approximation.
It omits refraction, absorption, finite-temperature/magnetization corrections,
optical misalignment, electronic filtering and the detector modulation waveform.
The I/Q response describes demodulated quadratures with separate gains and
baselines, rather than an instrument-specific harmonic-demodulation algorithm.
Neutral and environmental inputs are observations/configuration supplied at
runtime; this module does not infer a gas composition or atmospheric index law.
A reference-baseline subtraction must use differences of these same inputs.

## Calibration and the profile chain

`init(chords, prior, neutral, environment, status)` takes existing `chord_t`
objects, one `calibration_prior_t` with eight physical parameters, and
one `neutral` and `environment` value per chord. Zero status means success.
One object shares a calibration across its chords. Use separate objects for
separately calibrated groups. Repeat a chord with different runtime inputs to
represent a sequence of measurements at a fixed plasma state.

The physical parameter order is fixed; the public `DI_*` constants name it:

| Index | Constant | Physical parameter and units |
| --- | --- | --- |
| 1 | `DI_WAVELENGTH` | fundamental wavelength, m |
| 2 | `DI_AI` | I amplitude, recorded signal units |
| 3 | `DI_AQ` | Q amplitude, recorded signal units |
| 4 | `DI_OI` | I offset, recorded signal units |
| 5 | `DI_OQ` | Q offset, recorded signal units |
| 6 | `DI_PHASE0` | phase offset, rad |
| 7 | `DI_NEUTRAL` | neutral-column correction scale, dimensionless |
| 8 | `DI_ENV` | environmental correction, rad per proxy unit |

All eight standard-normal calibration latents are `eta`. Correlated physical
priors use `prior%mean + matmul(prior%factor,eta)`, as the recorded-signal contract specifies.
Set the joint-model slot's `eta_sigma=1` and leave its separate `gain_sigma=0`.
Wavelength and both amplitudes must remain positive; invalid evaluations stop
with an error. Use Gaussian priors narrow enough for the intended local fit.
The universal physical constant `re` is fixed, rather than an instrument
calibration. No wavelength, gain or correction coefficient is hidden in the
operator. The wavelength/density scale degeneracy is constrained by the
wavelength prior, and a constant environmental proxy is degenerate with
`phase0`; varying reference measurements are needed to separate them.

The operator delegates sample locations and column integration to
`tiago_line_integrals::interferometer_t`. Its JVP (Jacobian-vector product)
and VJP (transposed Jacobian-vector product) compose with the existing kinetic
B-spline profiles, Gaussian-process priors and provider/VMEC flux map. It
integrates the supplied `Q_NE` row, including the core's differentiable
outside-plasma continuation, without imposing another boundary mask.
Each chord describes one physical passage. Supply separate paths and combine
their columns before the phase map for a multipass optical arrangement; this
object does not implicitly double the density column.

## Quadrature control

For each configured chord, `init` retains the requested composite
Gauss-Legendre rule and a second rule with twice as many panels. The predicted
I/Q uses the finer rule. `sample_points` returns both sets to the existing
sampler; the coarse points contribute only to the convergence check.

Call `quadrature_check(sample,atol,rtol,columns,error,status)` before accepting
a prediction. It returns the fine columns and `abs(fine-coarse)` in m^-2,
with success only if every difference is at most `atol+rtol*abs(columns)`.
Increase `chords%npanel`, reinitialize and resample if this fails. The
comparison is an empirical convergence indicator, not a rigorous bound:
unresolved narrow structures can be missed by both rules. Check refinements
against the spatial scales of the profile and equilibrium. Derivatives are
exact for the fixed accepted discretization; changing quadrature nodes is not
differentiated. Choose tolerances before optimization and retain sufficient
resolution throughout it. Channel sums and point adjoints use OpenMP.

## Reduction and noise

`tiago_iq_phase_reduction::iq_phase_reduction_t` implements the existing
`reduction_step_t`. Input samples are interleaved I/Q pairs with identical
time stamps within each pair and nondecreasing time across pairs.
Set `id` to a versioned reduction identity, `amplitude(2)`, `offset(2)`, and
`anchor` to a known phase within pi of the first true phase. The calibrated
principal phase is `atan2((Q-offset_Q)/amplitude_Q,
(I-offset_I)/amplitude_I)`. Each phase is lifted to the fringe closest to
the preceding unwrapped phase. Set the next chunk's anchor to the last output
of the preceding chunk. No raw trace need be retained in full.

Correct unwrapping requires adjacent true phases to differ by less than pi,
with adequate noise margin. An absolute fringe reference is required;
a single quadrature pair cannot identify the fringe number. Missing fringes,
large gaps, clipping and near-zero phasors need independent quality control.
Phasors whose normalized radius is at most `min_radius` are rejected.

The reduction returns radians and the first-order conditional covariance
`J C_IQ transpose(J) + noise_variance * identity`, including correlated input
noise. This is a local high-signal-to-noise approximation, not a circular
likelihood near zero amplitude or fringe ambiguity. Calibration is held fixed
in this reduction contract; retain calibration uncertainty in the recorded
I/Q likelihood when inferring it, rather than adding an independent copy to
phase noise. Output provenance records the parent checksum and reduction
parameters. JVP/VJP are exact within a fixed fringe assignment and exclude the
discrete decision at an unwrapping boundary.

`tiago_dispersion_column::dispersion_column` inverts the unwrapped total
phase, using the same wavelength and correction inputs. Like the phase, I/Q
and atan2 kernels, it exposes generated partial-derivative rows. Its output
is an electron column, rather than a profile reconstructed from one chord.

## Verification and limits

CTest covers slab columns and cylinder profiles `ne=n0(1-r^2)^alpha`, with
radius normalized to one and impact parameter `b`. The closed-form column
is `n0 (1-b^2)^(alpha+1/2) sqrt(pi) Gamma(alpha+1)/Gamma(alpha+3/2)`, from
[DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1). Tests use alpha 0.5, 1, 1.3 and
2 at nine impact parameters, the vacuum/neutral/dilute cold-plasma limits,
finite differences of generated rows, calibration derivatives, full B-spline
and circular/VMEC state chains, covariance propagation and a multi-fringe
I/Q-to-column round trip, including a chunk boundary.

A 400-realization synthetic prior-predictive study fits density scale and all
eight calibration latents, and checks marginal Laplace coverage, pull bias
and contraction of density and I/Q offsets. This tests a locally identifiable,
high-signal-to-noise regime with a known density shape. It does not validate
experimental measurements, arbitrary profile identifiability, strong-noise
unwrapping, multimodal posterior coverage or a complete VMEC equilibrium
re-solve in that coverage study. The existing VMEC state derivative test
perturbs its geometry coefficients, not a new implicit equilibrium solve.
