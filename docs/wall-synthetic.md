# Wall-loss synthetic infrared signals

Link `tiago_wall` and use `tiago_wall_synthetic` for an explicit loss-event to
surface to raw infrared (IR) forward model. This input is independent of the
six-quantity profile-diagnostic interface. All numerical inputs are SI units,
except `kinetic_energy_eV`, which is explicitly in electronvolts (eV).
`tests/test_wall_synthetic.f90` is a complete synthetic example.

The first supported thermal model is a semi-infinite solid with constant thermal
properties, uniform initial temperature within each cell, prescribed deposited
surface flux, and no lateral conduction, active cooling, radiative or convective
boundary loss. The caller must establish that its material, thickness and time
range satisfy this approximation. Setting a regime string does not validate it
for an experimental component. Other thermal regimes require a separate validated
solver; its surface temperatures can be passed directly to `synthesize_ir`.
No machine geometry, calibration, or experimental material records are supplied.

## Event and normalization contract

Initialize `surface_deposition_t` with `init_deposition(surface, edge_s, area_m2,
mesh_id, version, status)`. Cell areas are square metres; bin edges are seconds.
Accumulate bounded `loss_event_block_t` blocks with `accumulate_losses`. Each
block supplies `time_s`, `kinetic_energy_eV`, `weight`, `deposited_fraction`,
one-based `facet`, `wall_hit`, `exit_code`, and the strings `source_id`, `version`,
`frame`, `position_kind`. Accepted position kinds are `guiding-centre` and
`physical-impact`. Only `wall_hit == 1` and `exit_code == 2` are deposited.
Numerical failures and non-wall records are excluded, even when their unused
energy/facet values are invalid. Selected malformed records reject the block
without changing the surface. Blocks must share source identity, normalization,
frame and position representation. A facet is a finite-area deposition cell;
this contract does not infer sub-facet spreading or gyro-orbit displacement.

`PARTICLE_COUNT` means weight in particles per event and absolute impact time.
`PARTICLE_RATE` means response-marker weight in particles per second and time
after birth; it requires `injection_history_t` with second-valued `edge_s` and
dimensionless piecewise-constant `amplitude`. The rate response is integrated
against the shifted source history exactly in each deposition bin. Count events
reject an injection history. Relative source weights alone are insufficient.

An eV is exactly 1.602176634e-19 joule. Deposited energy is particle count times
residual impact energy times `deposited_fraction`, which must be supplied in
[0,1]. `deposition_flux` divides deposited joules by cell area and bin duration
once, yielding watts per square metre. A source rate is integrated to particle
count before this division. Output bins must cover all selected count events or
the nonzero shifted injection support; partial coverage rejects the block.
Bins include their left edge, exclude their right edge, and include the final
right edge in the last bin. Total absorbed energy is conserved over these bins.
Rejected/outside events are not silently reclassified as a zero heat load.

`semi_infinite_material_t` requires `conductivity_W_mK`, `capacity_J_m3K`,
`initial_K` per cell and `version`, with `regime = SEMI_INFINITE_REGIME`.
`capacity_J_m3K` is volumetric heat capacity, not specific heat per unit mass.
`surface_temperature` evaluates the exact square-root response to each
piecewise-constant flux bin at the supplied times, including diffusion into the
solid after the flux stops. The deposition time grid represents a bin-averaged
load, so an instantaneous event is spread uniformly over its chosen bin.
Refine this grid and the exposure quadrature to establish temporal convergence.

## SIMPLE NetCDF input

`tiago_simple_wall_source:accumulate_simple_netcdf(path, surface, source_strength,
normalization, source_id, source_version, frame, deposited_fraction, status,
history, external_weight, block_size)` reads bounded blocks, with a default
maximum of 65536 markers. The last three arguments are optional.
`source_strength` is the total source particles (`PARTICLE_COUNT`) or total
source particles per second (`PARTICLE_RATE`), over all source markers.
`deposited_fraction` is required per facet. The all-marker relative-weight sum
normalizes each marker, including markers that are retained, thermalized or
numerically failed. The reader never normalizes over losses alone.

For SIMPLE schema 1.1 it reads explicit `kinetic_energy_eV` (`eV`),
`marker_weight` (`1`), `times_lost` (`s`), `wall_hit_triangle` (zero-based with
`start_index=0`), `wall_hit`, `orbit_exit_code`, `cartesian_frame` and
`stl_facet_count`. Only `supplied_relative_weights` are accepted automatically.
When sampling weights are unknown, `external_weight` is required. An explicit
physical source strength is always required. No default equal marker weight is
promoted into physical normalization.

Legacy files without explicit energy use the recorded `particle_energy_eV`
times the square of the fourth component of `zend` (normalized speed).
External source weights may override recorded weights; the caller supplies the
legacy Cartesian frame explicitly. No completed orbit calculation needs to be
repeated to use this adapter. Coordinates are not read or converted: deposition
uses the recorded facet identity and the caller's areas. The reader checks the
facet count, frame and index convention, but SIMPLE files do not provide a mesh
content digest. The caller must verify the exact archived wall mesh and retain
its content identity in `mesh_id`; matching names or counts alone is insufficient.
The entire read is transactional if a later block is invalid.

## Raw camera contract

`ir_camera_t` requires all calibration arrays; allocation with zero defaults is
not a calibration. `synthesize_ir(camera, temperature, result, status)` consumes
temperatures in kelvin with shape `(cell,quadrature_time)` and integrates spectral
Planck radiance plus reflected background for an opaque diffuse surface. Emissivity is wavelength/cell
resolved, and `background_W_m3sr` is incident reflected-background radiance in
watts per cubic metre per steradian. The emitted radiance is emissivity times
Planck radiance plus one minus emissivity times background radiance.

Supply positive `wavelength_m`, positive integration weights `wavelength_weight_m`
in metres, and `transmission` in [0,1] for the composite normalized spectral
energy response. The absolute radiometric conversion belongs in the gain.
The spectral energy response must include the detector's wavelength sensitivity
when that sensitivity is not flat; optical transmission alone is insufficient.
Supply sparse `link_pixel`, `link_cell`, `etendue_m2sr`, `link_visible` arrays,
sorted by pixel. Etendue is effective viewed area times solid angle, in square
metres steradian; it includes the projection and optical pixel mapping. Pixel
spreading and point-spread functions are represented by multiple links. Do not
include the same view factor again in gain or spectral transmission.

Supply `exposure_edge_s(2,exposure)` with ordered nonoverlapping exposures,
`quadrature_time_s`, `quadrature_weight_s` and `quadrature_exposure`. Each
positive quadrature weight is seconds, its time lies inside its exposure, and
weights sum to that exposure's duration. Spectral and temporal quadrature are
caller-controlled approximations, not an implicit integration tolerance.
The numerical Planck kernel accepts wavelength times temperature up to
100 metre kelvin; values beyond that bound require a different resolved
long-wavelength representation. Nonfinite outputs are rejected.

`integrated=.true.` produces gain times exposure-integrated radiant energy plus
offset; gain has output units per joule. `integrated=.false.` divides by exposure
duration and uses gain in output units per watt. Allowed `units` are `counts`
and `V`. Both modes require per-pixel positive `gain`, finite `offset`, explicit
upper `saturation`, and `noise_variance(pixel,exposure)` in squared output units.
The output is the deterministic expected raw signal; noise is represented by
the supplied variance, not by random draws. Only affine, unsaturated detector
response is implemented; nonlinear detector
calibration must be modeled separately. Geometry, calibration, source and thermal
versions are mandatory strings (`geometry_version`, `calibration_version`,
`source_version`, `thermal_version`); they should identify archived records,
including normalization, quadrature and material assumptions.

`wall_ir_result_t` contains one existing `signal_window_t` per pixel, the exposure
edges, a validity mask and flags. Missing/occluded views use `IR_UNMAPPED`;
saturated values use `IR_SATURATED`. Their stored finite detector values are
not trustworthy measurements. Use `ir_valid_window` to select valid samples
before applying existing `linear_reduction_t` reductions. It preserves variance
and selects the corresponding covariance submatrix. All-invalid pixels reject
extraction. Window timestamps are exposure midpoints; full edges remain in the
result and must be compared to matching measured exposures.

The supplied variance describes detector noise only. Reused or split Monte
Carlo markers, calibration nuisance covariance and model discrepancy do not
become independent detector noise automatically. Retain the source's sampling
and split-marker relationships in its archived provenance; supply justified
correlated contributions through the existing window covariance interface.
The forward model does not estimate orbit-sampling uncertainty or attribute a
measured wall signal uniquely to fast ions rather than other heat sources.

## Conditional derivatives and regeneration

`temperature_energy_jvp` and `temperature_energy_vjp` are the Jacobian-vector
product and its transpose for deposited joules with fixed bins, areas and
material. `ir_jvp` and `ir_vjp` differentiate temperature, pixel gain and offset
with fixed wavelength/exposure quadrature, optical links and visibility. They
reject invalid or saturated operating points. Facet assignment, hard time-bin
crossings, visibility and saturation transitions are nonsmooth. There is no
Jacobian of the stochastic orbit solver and no marker-by-pixel dense matrix.

The pointwise kernels and their partial derivatives are emitted by FortSym from
`gen/gen_wall.f90`; ordinary consumers build the checked-in generated sources.
From the TIAGO repository root, regenerate with:

```sh
cmake -S gen -B build-codegen -G Ninja -DFORTSYM_BUILD=/path/to/fortsym/build
cmake --build build-codegen --target gen_wall -j2
build-codegen/gen_wall
```

Configure the native TIAGO build with `TIAGO_FORTSYM_BUILD` to enable the exact
`tiago_wall_codegen` CTest regeneration comparison. Independent energy balance,
source normalization, thermal pulse, Planck, covariance and derivative oracles
are registered as `tiago_wall_synthetic` and `tiago_simple_wall_source`.
`benchmark_wall` measures the cell/pixel OpenMP path with synthetic data.

## Public physical basis and limits

[Cornelissen et al. (2022), section 2.1, equations 2–3](https://doi.org/10.1088/1361-6587/aca0bf)
give the restricted semi-infinite square-root thermal model. This implementation
uses its constant-property limit and linear superposition for flux bins; it
does not adopt the paper's material fits or experimental calibration.

[LeViness et al. (2024), section III](https://doi.org/10.1063/5.0214788) distinguish
source-rate weighting, entrance phase space and scintillator impact density from
camera response. [Jansen van Vuuren et al. (2024), sections 5–6](https://doi.org/10.1016/j.fusengdes.2024.114520)
likewise require optical and sensor response. No published design geometry is
used as an experimental calibration here. The
[FC-FILD measurement report](https://doi.org/10.1063/5.0223953) establishes a
different charge-collection observable. No FILD response is implemented by this
wall-IR API: guiding-centre position and energy do not specify physical
entrance velocity or gyrophase, so a calibrated incident-phase-space kernel or
an appropriate full-orbit entrance model is required for that extension.
