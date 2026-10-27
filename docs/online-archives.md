# Writing an archive adapter

Link `tiago_signals` and implement the existing
[`signal_adapter_t`](recorded-signals.md) contract in a separate device module.
Keep archive addresses, channel names, geometry, calibration records and data
in private runtime configuration. The public example uses only the loopback
server in `tests/mock_archive.py` and synthetic signals.

## REST/JSON configuration

`tiago_rest_adapter::rest_adapter_t%configure(path,status)` reads JSON. Each
entry maps a logical channel to a GET URL and paths understood by
[JSON-Fortran 9.3.1](https://github.com/jacobwilliams/json-fortran/tree/9.3.1):

```json
{
  "source_id": "synthetic",
  "channels": [{
    "logical": "temperature",
    "url": "http://127.0.0.1:8080/synthetic?first={first}&last={last}&version={version}",
    "time_path": "data.time",
    "value_path": "data.values",
    "variance_path": "data.variance",
    "version_path": "version",
    "version": "v1",
    "units": "eV",
    "scale": 1.0,
    "offset": 0.0,
    "time_scale": 1.0
  }],
  "policy": {"max_window": 60.0, "max_response_bytes": 8388608}
}
```

`{first}` and `{last}` are required and represent seconds in the source's time
origin, using round-trip decimal binary64 values. `{version}` and `{channel}`
are optional and URL escaped. A nonempty request version overrides the configured
version; the response must contain that version. An unpinned read resolves the
first response version and pins every subsequent window to it. Sources must
supply immutable revision identifiers rather than a moving alias for reproducible
reads. The configured version is the URL selector for an initially unpinned read.

Time/value arrays must have equal length, finite values, nondecreasing times,
and lie inside the requested closed interval. Empty arrays succeed.
`variance_path` is optional; omission explicitly means zero independent variance.
When supplied, it must contain nonnegative variances with the same length.
`time_scale` is positive and converts response timestamps to seconds. Values
are converted by `scale * value + offset`, and variances by `scale**2`.
These are fixed unit conversions, using the existing FortSym affine kernel;
inferred instrument gains and offsets belong in
[`calibration_prior_t`](recorded-signals.md), with their calibration priors.
No physical calibration should be hidden in the JSON scaling fields.

## Transport, pacing and memory

Assign the adapter and its policy to an `online_source_t`:

```fortran
online%adapter => adapter
online%policy = adapter%policy
call cache%init(online, 8)
call cache%read(request, observation, status)
```

Both `adapter` and `online` must have `target` and outlive their wrappers.
Use the existing `cached_source_t` as the caller-facing source: it checks the
cache before any network access, including before splitting a long request.
Its cache is a bounded session snapshot; clear it after changing adapter or
reduction configuration. The transport itself retains no cache of raw responses.
A source/adapter instance handles one read at a time. Independent instances can
run concurrently, sharing the same process-wide host scheduler.

The libcurl client uses ISO C bindings and a small C shim for variadic options,
POSIX clocks and mutexes. It does not execute external commands. TLS peer and
hostname verification are always enabled. HTTP/HTTPS are the only protocols;
redirects and credentials embedded in URLs are rejected. Authentication and
non-GET requests require a separate transport extension.

Each adapter can configure every field of `http_policy_t` in its JSON `policy`:

| Setting | Default | Meaning |
| --- | ---: | --- |
| `min_interval` | 0.5 s | Minimum spacing between request starts. |
| `timeout` | 30 s | Maximum transfer time, including connection. |
| `connect_timeout` | 10 s | Connection timeout. |
| `backoff`, `max_backoff` | 0.5, 30 s | Initial and maximum exponential base delay. |
| `jitter` | 0.25 | Uniform additional fraction of the base delay. |
| `max_retries` | 3 | Total retry budget across a complete split read. |
| `max_elapsed` | 120 s | Complete read deadline, including pacing/retries. |
| `max_window` | 60 s | Maximum seconds per backend request. |
| `max_response_bytes` | 8388608 | Hard transfer-byte cap per response. |
| `max_windows` | 10000 | Maximum backend windows in one read. |
| `max_output_samples` | 1000000 | Maximum retained output sample count. |
| `user_agent` | TIAGO/1 and repository URL | Identifying HTTP User-Agent. |

There is one transfer per hostname in flight, across ports and adapter instances
inside one process. Different processes have independent schedulers. The next
start respects both the previous and current caller's spacing. HTTP 429/5xx and
transient connection failures use exponential backoff with jitter. libcurl
parses both integer and HTTP-date Retry-After headers; the requested delay is
never shortened to the backoff cap. It also blocks other callers on that host.
If the delay exceeds the complete read deadline, the read fails instead of
retrying early. TLS/configuration failures and other 4xx responses are not retried.

Long closed intervals are partitioned at consecutive representable binary64
times, so a boundary sample is requested once. Each response is streamed into
the adapter in libcurl blocks. The default JSON decoder buffers at most the
configured response cap because JSON-Fortran parses a complete JSON document.
An oversized or incomplete response fails. Interrupted success bodies fail
without replay, since the W0 decoder has no rollback operation. Transient failures
before a success body can retry safely. An interrupted decoder is finalized
and its output discarded; custom decoders must reset their state on `final`.

Assign an existing `reduction_step_t` to `online%reduction` to process each
emitted window before retaining results; raw arrays are then released between
windows. Choose windows appropriate to the step's sample-count requirements.
Reductions must document any lost correlations between different windows.
Without reduction, `read` necessarily retains the requested raw samples up to
the output cap. Correlated within-window covariance is joined block diagonally;
cross-window correlated uncertainty needs a specialized source/reduction.
Dense covariance can consume quadratic memory, so reduce to small outputs.

## Binary and specialized formats

A device adapter can directly implement W0's incremental `decode` method. It
receives `int8` blocks, including NUL and high-bit bytes unchanged. Alternatively,
assign a `tiago_rest_adapter::payload_decoder_t` implementation to
`adapter%decoder`; this retains file-configured URL mapping but replaces JSON
parsing entirely. The decoder supplies units, validated windows, resolved version,
and W0 provenance/checksum. Return `SIGNAL_MORE` until enough bytes are available,
`SIGNAL_OK` for each complete emitted window, and an error for incomplete final
input. Reset internal state on final input, including an interrupted body.
Preserve sample ordering across emissions and implement reductions with the
existing exact JVP/VJP contract. No binary format is mandated by TIAGO.

## Synthetic example and verification

Build normally, then run `ctest --test-dir build -R 'tiago_(http|online|archive)'
--output-on-failure`. CTest starts an ephemeral loopback HTTP server, exercises
the clients and example, and shuts the server down as a cleanup fixture.
`examples/online_signals.f90` demonstrates online data, a three-sample mean,
the W0 cache, and comparison with the existing point-probe diagnostic.
The fixture writes its generated `adapter.json` in the build's `tests/mock_archive`.

Tests check analytic ramp values and variances, bounded windows and memory,
empty/malformed responses, binary bytes, version pinning, retry budgets,
Retry-After, connection failure, timeout and shared-host timing. Server-side
counters independently check cache hits and retry counts. Four concurrent
callers check serialization. Existing recorded-diagnostic tests check the
reused affine kernel and its derivatives. Optional
`-DTIAGO_FORTSYM_BUILD=/path/to/fortsym/build` enables a regeneration comparison.
These tests validate transport and synthetic ingestion, not an experimental
archive, instrument physics, TLS certificate deployment or archive authentication.
