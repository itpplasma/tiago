# Coherent harmonic demodulation of a detector waveform

`tiago_harmonic_lockin` converts a detector waveform and a calibrated sinusoidal
reference into in-phase and quadrature (I/Q) samples for TIAGO's existing
`iq_phase_reduction_t`. It is a generic coherent harmonic lock-in, with one
complete modulation period per output. It does not reproduce instrument firmware,
recover a clock from an arbitrary reference, or estimate calibration parameters.

The physical starting point is Eq. (1) of [Brunner et al., *Review of Scientific
Instruments* **93**, 023506 (2022)](https://doi.org/10.1063/5.0070041), available as
an [institutional full text](https://publikationen.bibliothek.kit.edu/1000143338/148149761).
Their Eq. (2) identifies the first and second harmonics for lock-in demodulation;
Eq. (3) describes a different, sectioned-integration method. Section II A discusses
filtering in the lock-in implementation. This implementation specifies its own
coherent finite-window filter below, rather than importing unspecified filter
coefficients or a phase-locked loop from that experiment.

## Signal, normalization and sign

Write the calibrated detector sample as

\[
s(t)=D+A\cos[\rho\sin\theta(t)+\phi],\qquad
\theta(t)=2\pi f(t-t_0-\tau)+\alpha.
\]

Here `D` is a constant detector background within a period, `A` is the fringe
amplitude in calibrated detector units, `rho` is the modulation depth in radians,
`phi` is the measured phase in radians, `f` is the modulation frequency in hertz,
`t0` is the reference epoch, `tau` is the detector delay relative to the reference
in seconds, and `alpha` is the reference phase at the epoch in radians. Positive
`tau` shifts the modulation at the detector later. Detector calibration is
`s=(detector-detector_offset)/detector_gain`.

The Jacobi–Anger expansion of the cosine gives the second cosine-harmonic
coefficient `2 A J2(rho) cos(phi)` and the first sine-harmonic coefficient
`-2 A J1(rho) sin(phi)`. `Jn` denotes the Bessel function of the first kind of
integer order `n`. An average over a complete period therefore gives

\[
 I={\langle(s-\bar s)\cos2\theta\rangle\over J_2(\rho)}
   =A\cos\phi,\qquad
 Q={\langle(s-\bar s)\sin\theta\rangle\over-J_1(\rho)}
   =A\sin\phi.
\]

The minus sign in the quadrature denominator is essential. The factor two in the
Fourier coefficients cancels the mean square of a unit sine or cosine. This
convention makes the downstream phase `atan2(Q,I)`. In particular, it does not
use the zeroth-order Bessel/Struve normalization of Eq. (3). The finite-sample
implementation replaces each average with an equally weighted sum over `N`
samples and divides by `N J2` or `-N J1`.

There are `N` samples per period with integer nanosecond cadence `cadence_ns`.
Write `dt = cadence_ns * 1e-9` seconds, so that `f N dt = 1`. Aliasing from higher
detector harmonics is a finite-grid error. Initialization checks the discrete
responses of both `cos(rho sin(theta))` and `sin(rho sin(theta))` against their analytical
Bessel responses, including leakage into the opposite quadrature. All four
absolute discrepancies must be at most `alias_tolerance`. Either Bessel
normalizer at or below `minimum_normalizer` in absolute value is rejected.
These user-supplied tolerances define admissible grids and noise amplification;
they are not hardware defaults.

The implementation subtracts each period's detector mean algebraically:
`sum(s*w)-sum(s)*sum(w)/N`. This rejects constant backgrounds even when floating
point reference weights have a small nonzero sum. It cannot reject a background
or amplitude variation synchronous with the selected harmonics. The identities
`I=A cos(phi)` and `Q=A sin(phi)` assume that the background `D`, amplitude `A`
and phase `phi` are constant within each period.

## Reference and calibration contract

The supplied reference must satisfy

\[
 {r_k-r_{\rm offset}\over r_{\rm gain}}
 =p\sin(2\pi k/N+\alpha),\quad k=0,1,\ldots,
\]

within `reference_tolerance`, where `p` is `reference_polarity`, either `+1` or
`-1`, and the reference gain is positive. The polarity describes inversion of
the observed reference channel; it does not change the detector mixing phase
`theta`. Delay is applied to the detector mixing weights, not to this reference
validation. The caller supplies the calibrated epoch, phase, polarity and delay.
A noisy square reference or a drifting clock
requires a separately specified reference reconstruction stage and is rejected
by this API when it violates this contract.

The final calibration is componentwise
`iq=(harmonic_iq-iq_offset)/iq_gain`. Both gains must be nonzero. This permits
explicit measured offset, gain and sign corrections; it does not fit a rotated
ellipse. `valid_from_ns` and `valid_to_ns` bound every accepted input timestamp.
All calibration values, tolerances, frequency, depth, cadence and period length
must be supplied explicitly. A completed output is marked invalid when either
component is nonfinite or its radius does not exceed `minimum_radius`.

## Configuration and calls

`harmonic_lockin_config_t` contains the following fields. Integer nanosecond
values use `int64` and real values use `real64`, the 64-bit kinds from
`iso_fortran_env`. All real values must be finite. The timestamp validity interval
is inclusive, must be ordered, and must contain `epoch_ns`.

| Field | Meaning, units and bounds |
| --- | --- |
| `cadence_ns` | Positive sample spacing in integer nanoseconds; `N * cadence_ns` must fit in a 64-bit integer. |
| `samples_per_period` | Integer sample count `N`, at least 8. |
| `frequency_hz` | Positive modulation frequency `f` in hertz; `abs(f * N * dt - 1)` must be at most 64 machine epsilons of `real64`. |
| `epoch_ns` | Reference epoch `t0` in integer nanoseconds; the first input timestamp. |
| `delay_ns` | Signed detector delay `tau` in integer nanoseconds. |
| `valid_from_ns`, `valid_to_ns` | Inclusive bounds on input timestamps in integer nanoseconds. |
| `modulation_depth` | Positive modulation depth `rho` in radians. |
| `detector_offset` | Offset in raw detector units, removed before harmonic extraction. |
| `detector_gain` | Nonzero gain in raw detector units per calibrated detector unit; either sign is allowed. |
| `reference_offset`, `reference_gain` | Offset and positive peak amplitude in raw reference units. |
| `reference_phase` | Modulation phase `alpha` at the epoch in radians, before reference-channel polarity. |
| `reference_polarity` | Reference-channel polarity `p`, either `+1` or `-1`. |
| `reference_tolerance` | Positive, dimensionless absolute error allowed after reference offset/gain calibration. |
| `alias_tolerance` | Positive, dimensionless absolute discrepancy in each harmonic-response check, before Bessel normalization. |
| `minimum_normalizer` | Positive, dimensionless cutoff for `abs(J1(rho))` and `abs(J2(rho))`, before multiplication by `N`. |
| `iq_offset(2)` | Offsets in harmonic I/Q units, ordered as I then Q. |
| `iq_gain(2)` | Nonzero component gains in harmonic I/Q units per final output unit, ordered as I then Q. |
| `minimum_radius` | Nonnegative radius threshold in final I/Q units, after `iq_offset` and `iq_gain`. |

Use an instance `lockin` of `harmonic_lockin_t`, a `config` of
`harmonic_lockin_config_t` and an `output` of `harmonic_lockin_output_t`:

```fortran
call lockin%initialize(config, status)
call lockin%consume(time_ns, detector, reference, output, status)
n_buffered = lockin%buffered_samples()
```

`status` and `n_buffered` are default integers. Both subroutines return
`status=0` on success and `status=1` on rejection. `consume` requires equal-length
one-dimensional input arrays: `time_ns` has kind `int64`, while `detector` and
`reference` have kind `real64`, from `iso_fortran_env`. `buffered_samples()`
returns the number of samples retained in the current partial period.

Each `consume` call replaces `output`. On rejection its allocatable arrays are
unallocated. On success with `M` completed periods, the arrays have shapes
`time_ns(M)`, `iq(2,M)` and `valid(M)`, with `iq(1,:)=I` and `iq(2,:)=Q`.
Successful calls with no completed period allocate zero-sized output arrays.
An invalid completed I/Q sample is still emitted with `valid=.false.` and does
not cause the chunk to be rejected.

## Time, filtering and streaming

The filter is a coherent rectangular window, followed by decimation by `N`.
In the following formula `dt` is the cadence converted to seconds and `Omega`
is baseband angular frequency in radians per second. The averaging response is

\[
 H(\Omega)=e^{-i\Omega(N-1)dt/2}
 {\sin(N\Omega dt/2)\over N\sin(\Omega dt/2)}.
\]

The expression has its continuous limit at a zero denominator. Its group delay
is `(N-1)dt/2` seconds, and its output cadence is `N dt` seconds. The output
`time_ns` records the last input sample of the completed window. The physical
window center is this timestamp minus the group delay, which may be a half
nanosecond. There is no additional bandpass or infinite-impulse-response filter.
This window attenuates higher frequencies but does not guarantee elimination of
all baseband aliasing for an arbitrarily varying phase. `H` describes the
averaging stage, not the full transfer from a varying physical phase to
`atan2(Q,I)`.

`initialize` starts a new stream. `consume` accepts consecutive chunks and emits
only completed periods; empty chunks are permitted. A partial period remains in
state and is never emitted by padding. Oscillator position and partial sums carry
across calls. The first timestamp must equal `epoch_ns`; every later timestamp
must equal the preceding timestamp plus `cadence_ns`. Gaps, overlaps, changed
cadence, out-of-validity samples, nonfinite inputs and a mismatched reference
reject the whole chunk before advancing state. A rejected initialization also
preserves the previous stream.

Timestamps remain 64-bit integers. Mixing phases use the bounded sample position
within the period, and delay is reduced modulo the integer period before real
conversion. No absolute epoch is converted to floating point. State stores three
real values per modulation sample (two mixer weights and one expected reference)
and a fixed number of accumulators; it never stores the run history. Output
storage is proportional to completed periods in the current caller chunk.

The output I/Q columns can be interleaved into a `signal_window_t`, with the
same timestamp for each pair, before applying `iq_phase_reduction_t`. That
reduction's noise/provenance contract must be supplied by the adapter. This
module does not infer waveform noise covariance. For phase unwrapping across
successive I/Q windows the adapter must carry the previous phase as the next
reduction's `anchor`; the demodulator itself does not unwrap phase.

## Verification and regeneration

The native CTest `tiago_harmonic_lockin` uses only synthetic detector samples
from the cosine model, independent of the generated mixer. It checks multiple
phases, amplitudes, modulation depths, backgrounds, calibration gains, polarity,
reference phase, delays, large integer epochs, three chunk partitions, rejected
chunk retries, singular normalizers, noncoherent grids, harmonic aliasing and a
long stream. `tiago_harmonic_codegen` checks reproducible FortSym generation
when `FORTSYM_BUILD` is supplied. Ordinary builds use the committed generated
kernel and need no computer algebra installation.
