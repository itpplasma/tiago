# Recorded signals and calibration

A device module supplies channel configuration and instrument geometry at run
time. TIAGO supplies source, reduction and recorded-diagnostic contracts, with
no archive or device configuration compiled into the library. Link
`tiago_signals` for ingestion and `tiago_inference` for diagnostics.
The offline reader requires the netCDF C library at configure and link time.

## Sources and adapters

`tiago_signal_source::signal_source_t%read(request, window, status)` accepts a
`signal_request_t`: logical `channel`, closed time interval `[first,last]` in
seconds, and optional allocatable `version`. An absent or empty version asks
for the source's current version. A successful `signal_window_t` contains
nondecreasing `time`, scalar `samples`, nonnegative independent `variance`,
`units`, optional additive correlated `covariance`, and `provenance`:
`source_id`, resolved immutable `version`, `request_key`, content `checksum`,
`record_id` and free-form `metadata`. Time zero is defined by the source and
recorded in metadata. Empty windows succeed with zero-length arrays.

`total_covariance()` returns `diag(variance) + covariance` (an absent correlated
term is zero). This allocates a dense matrix; reduce large traces before
calling it. Producers must supply a positive semidefinite total covariance.
Consumers must check status before using any output: `SIGNAL_OK=0`,
`SIGNAL_ERROR=1`; other nonzero source statuses are also failures.

`signal_adapter_t` is the interface for a separately implemented transport:

| Method | Contract |
| --- | --- |
| `configure(path,status)` | Load a runtime configuration file; keep private values outside TIAGO. |
| `map_request(request,backend,status)` | Map the logical request to `backend_request_t` (`endpoint`, `method`, `body`, canonical `request_key`). |
| `decode(request,bytes,final,window,status)` | Consume an `int8` byte chunk. Return `SIGNAL_MORE=2` when incomplete, or `SIGNAL_OK` with a complete window, or an error. `final` marks the last chunk. |

An adapter owns incremental decoder state; one instance handles one request
at a time. A transport feeds chunks, consumes successful windows immediately,
and never concatenates the entire response. It must reject an incomplete final
chunk, enforce the requested version and time interval, and record response
provenance. Header/version negotiation and HTTP transport belong to the online
source, which is not implemented here. A private device module (for example
W7-X or LHD) implements this interface and configures its own logical-channel
mapping; no private endpoint, stream name or measurement belongs in TIAGO.

`cached_source_t%init(upstream,capacity)` creates a bounded in-memory FIFO cache
around any source, including a reduced source. It reads before fetching and
keys exact channel/time/version requests. Entries retain all provenance and
return independent copies. Errors and inconsistent pinned versions or checksums
are not cached. An unpinned request is a session snapshot until `clear()`;
use pinned versions for reproducibility. The upstream target must outlive the
cache. Clear the cache after changing upstream configuration or reduction
parameters. Persistence and cross-process sharing are left to future sources.

## Offline layout: `tiago-signal-1`

`tiago_netcdf_source::netcdf_source_t%path` names one netCDF file containing
one logical scalar channel and one version. A device-side registry may choose
files for many channels; the same read contract applies to every file.
Both classic netCDF and netCDF-4 with the following root-level layout work.

| Item | Required content |
| --- | --- |
| Dimension `sample` | Number of samples. |
| `time(sample)` | Double, nondecreasing seconds. |
| `samples(sample)` | Double, in the declared units. |
| `variance(sample)` | Double, independent variance in units squared. |
| `covariance(sample,sample)` | Optional double, symmetric correlated covariance in units squared. |
| Global text attributes | `format="tiago-signal-1"`, `channel`, `source_id`, `version`, `record_id`, `metadata`, `units`, `checksum`. |
| Global double attributes | `valid_from`, `valid_to`, seconds; use finite extremes when unrestricted. |

Times must be globally sorted by the producer. The reader uses binary searches
and the standard [netCDF hyperslab API](https://docs.unidata.ucar.edu/netcdf-c/current/group__variables.html)
through ISO C bindings, so it reads only the requested data and covariance
submatrix. It validates the returned window, channel, version and layout.
The file must remain immutable during a read. A partial read cannot certify the
ordering or integrity of data outside the requested window.

The `checksum` attribute is computed by `window%seal()` before writing the full
file. `content_checksum` specifies a portable noncryptographic content checksum:
start lookup3 seed at zero, add the bytes of `units`, then uppercase 16-digit
IEEE binary64 hexadecimal strings for `valid_from`, `valid_to` and, per sample,
`time`, `samples`, `variance`. If covariance is present, add the text
`covariance`, then its entries in Fortran column order. Each addition uses
`fortio_checksum::lookup3_checksum` seeded by the preceding result. The result
is `lookup3-hex64:` followed by eight uppercase hexadecimal digits. This detects
accidental changes; it is not an authenticity signature. Whole-record reads
check the stored checksum; partial reads seal the returned content instead.

## Reduction and likelihood

`reduction_step_t` declares `apply(input,output,status)` and exact `jvp`/`vjp`
with respect to input samples at fixed times and fixed step parameters. Each
step owns its additional noise model, returns the observable covariance, and
records its identity and parent checksum. A nonlinear step must document any
linearized covariance approximation separately from its exact derivatives.
`reduced_source_t` composes a source and a step through non-owning pointers;
these targets must outlive the wrapper. Wrappers can be chained or cached.

The reference `linear_reduction_t` has `weights`, `noise_covariance`, `units`
and a versioned `id`. It returns `weights * samples`, covariance propagated as
`weights * total_covariance * transpose(weights) + noise_covariance`, and a
time tag equal to the final input sample. Supply integration weights (including
sample durations) for a pulse integral. Its covariance propagation is exact
for linear maps, including correlated inputs, following
[JCGM 100:2008, section 5.2](https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf).
No pulse extraction or phase-unwrapping physics is implied by this reference.

Pass the reduced `samples` as a diagnostic slot's `data` and
`total_covariance()` as its `covariance` in `tiago_joint_model`. This includes
reduction noise in the existing likelihood whitening; do not add it again as
model error. Reductions on predictions must use the same parameters as those
on observations, and compose their JVP and VJP by the chain rule.

## Calibration and recorded diagnostics

`calibration_prior_t%load(source,request,status)` reads a calibration record
through any signal source. All parameter means are stored as `samples` at one
record timestamp, their prior covariance as the window's total covariance,
and their order and units in `metadata`/`units`. `record_id` identifies the
calibration, and its validity interval must cover the requested interval.
The covariance must be positive definite; exact known constants should be
removed from the inferred parameter block. Keep the parameter order fixed
within a record version.

The prior stores the mean, covariance, Cholesky factor, complete source record
and validity. Standard-normal latents map to physical parameters as
`mean + factor * z`; `jvp` and `vjp` supply this linear transformation and its
transpose. Validity is checked when loading; callers must reload for an
observation outside that interval.

`recorded_diagnostic_t` extends `tiago_profile_diagnostic::profile_diagnostic_t`
and carries this prior. Its calibration latents follow alignment nuisances in
`eta`, so the existing sampled-field JVP/VJP and equilibrium/profile chain
continue to apply. The reference `calibrated_diagnostic_t%init(model,prior,status)`
composes any existing profile diagnostic with shared gain and offset, in that
order: each predicted recorded value is `gain * model_signal + offset`.
Signals sharing a calibration must have compatible units. Its scalar map and
partial derivatives are generated with [FortSym](https://github.com/itpplasma/fortsym).
Channel loops use OpenMP above 1024 signals.

In a joint-model slot, set `eta_sigma` for the appended calibration latents to
one, retain the intended alignment scales, and leave `gain_sigma=0` for a gain
already represented here. The existing standard-normal nuisance-prior rows
then give the full correlated physical calibration prior without a second
calibration factor. Neither the profile sampler nor the equilibrium solver is
duplicated. Specialized recorded instruments can implement different response
maps through the same abstract type.

Tests use synthetic ramps, correlated sums, calibration records and the
existing point probe. They check the offset-only zero-input limit, finite
differences, adjoint identities and point-motion composition. No archive,
device measurement, instrument spectrum or radiation-transport model is
validated by these interface tests.
