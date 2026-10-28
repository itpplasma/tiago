# Periodic magnetic pickup voltage

`tiago_magnetic_pickup:magnetic_pickup_t` maps a single channel's signed linked
flux in webers (`Wb`) to continuous recorded voltage in volts (`V`). Use the
existing static magnetic response to supply linked flux, including winding
turns, orientation and any compensation combination. This operator does not
recompute sensor geometry. Electronics response is supplied per channel; no
device coefficients are built in.

The supported signal is the periodic, band-limited trigonometric interpolant
of an odd number `N >= 3` of uniformly spaced samples. It is not a streaming
filter or a model of start-up transients. The period is `N*dt`, where `dt` is
the sample interval in seconds. Even lengths are rejected because the Nyquist
bin has no unique real-valued derivative at the sample positions under an
arbitrary delay. Input outside the represented band can alias and must be
excluded by the caller's acquisition model.

For zero-based indices `j,k`, the forward discrete Fourier transform (DFT) is
`X[k] = sum_j flux[j]*exp(-2*pi*i*j*k/N)`. The inverse has the positive sign
and a factor `1/N`. Here `i` is the imaginary unit. Signed frequencies are
`f[k] = k/(N*dt)` for `k <= N/2` and `(k-N)/(N*dt)` otherwise. This matches the
[FFTW DFT convention](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html).
The voltage spectrum is

```
Y[k] = (-i*2*pi*f[k]) * gain * H[k] * exp(-i*2*pi*f[k]*delay) * X[k].
```

The leading factor is Faraday's minus time derivative. A positive delay means
the output at time `t` uses the input at `t-delay`. Delay is in seconds and is
bounded to half a record period. The response `H` is the full complex voltage
ratio, with units `V/V`, from ideal pickup voltage to the reported continuous
voltage; it can include coil departures from ideal sensitivity, cable,
preamplifier, rack and analog acquisition response. Any exposure or sample
aperture response must already be in `H`. `H[0]` must be real and its remaining
bins must have exact Hermitian conjugate pairs. Construct each negative bin
by conjugation of its positive partner. DC voltage is zero. The caller supplies
`max_frequency` in hertz; every nonzero signed grid frequency must be within
this calibrated band. No interpolation, extrapolation or truncation occurs.
A calibration available only above a nonzero lower cutoff cannot be used to
claim this full-grid model without additional low-frequency response evidence.

Call `init(utc_ns, reference_ns, prior, H, max_frequency, 'Wb', 'V', 'V/V',
status)`, then `predict(flux, z, voltage, status)`. Status zero indicates success;
any nonzero status rejects the operation. A failed initialization invalidates
the object. UTC means Coordinated Universal Time: timestamps and the reference
are nonnegative 64-bit integer nanoseconds since the Unix epoch. Intervals are
checked for exact positivity and uniformity before conversion to seconds.
Every timestamp minus the reference and the step must fit within `2^53`
nanoseconds, so integer differences are exactly representable in binary64.
This deliberately limits reference separation to about 104 days. No absolute
epoch timestamp is converted to floating point.

`prior` reuses `calibration_prior_t`, with two physical parameters in order:
dimensionless `gain` and `delay` in seconds. Its full two-parameter covariance
may correlate gain and delay. The standardized nuisance coordinates `z` map
through its existing Cholesky factor. Record source and record identifier are
required. The caller's record must identify the full chain represented by `H`,
its channel, stages, calibration convention and band; these descriptions are
retained in the existing provenance contract. The prior validity interval is
in seconds relative to `reference_ns` and must cover every sample. Changing
the reference requires changing this interval consistently. Finiteness,
units, dimensions, timing, response symmetry and validity are checked.

`jvp(flux,z,flux_dot,z_dot,voltage_dot,status)` propagates directional
perturbations; `vjp(flux,z,voltage_bar,flux_bar,z_bar,status)` is its transpose
under ordinary real sample-space dot products. The tabulated `H`, clock and
reference are fixed inputs to these derivatives. Uncertainty in `H` is not
represented by the gain/delay prior. Derivative products need linear storage
in `N` and use the shared FFT; no sample-by-sample covariance matrix is
formed. They can be composed with existing calibration and magnetic response
products. This voltage operator does not differentiate clipping, ADC
(analog-to-digital converter) quantization, resets, integration or saturation.
Those producer operations require their own recorded-signal contract.

Full-chain amplitude and phase calibration, rather than coil-only flatness,
is supported by Rahbarnia et al., *Plasma Physics and Controlled Fusion* 63
(2021) 015005, [doi:10.1088/1361-6587/abc395](https://doi.org/10.1088/1361-6587/abc395),
section 2.2 and figures 2–4. The authors' component calibration includes
preamplifier, long cables and rack electronics. The method motivates `H`; its
example curves are not coefficients for another channel or campaign. This
generic mechanism does not certify an experiment: dated channel responses,
producer timing and integration-state records remain necessary.
