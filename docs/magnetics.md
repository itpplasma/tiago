# Magnetic diagnostics and `tiago_vacuum_cli`

`tiago_vacuum_cli` evaluates DIAGNO-compatible magnetic diagnostics from coil
currents and, optionally, the plasma currents of an equilibrium. Coil files are
read with libneo; the field of each straight filament uses the exact
Hanson–Hirshman Biot–Savart kernels. Run `tiago_vacuum_cli --help` for the full
option list.

## Coil signals

```bash
./build/tiago_vacuum_cli \
    --coils tests/data/coil_sample.neo \
    --flux tests/data/fluxloop_sample.diagno \
    --segrog tests/data/segrog_sample.diagno \
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8
```

At least one of `--flux`, `--segrog`, `--bprobes` and `--mse` is required. The
example writes `tiago_flux.csv` and `tiago_segrog.csv` (names set with
`--flux-out` and `--segrog-out`).

- `--samples n` sets the points per segment with the midpoint rule (DIAGNO's
  `int_type='midpoint'`). `--gauss` uses Gauss–Legendre points instead; eight
  Gauss points are typically more accurate than thousands of midpoint samples.
- `--nfp` sets the field-period geometry.
- `--flux-turns` and `--segrog-turns` read text files of `label scale` pairs,
  equivalent to DIAGNO's `flux_turns` and `segrog_turns` namelist arrays.

The NCSX cases in `tests/cases/` include such files; their coil sets are
downloaded by `benchmarks/xdiagno/fetch_data.sh`:

```bash
benchmarks/xdiagno/fetch_data.sh
./build/tiago_vacuum_cli \
    --coils benchmarks/xdiagno/_work/data/coils.NCSX \
    --flux tests/cases/ncsx_nfp3/fluxloop.diagno \
    --segrog tests/cases/ncsx_nfp3/segrog.diagno \
    --output-dir build/ncsx --seg-area 3.40e-4 --gauss --samples 8 --nfp 3 \
    --flux-turns tests/cases/ncsx_nfp3/flux_turns.csv \
    --segrog-turns tests/cases/ncsx_nfp3/segrog_turns.csv
```

### Coil currents

`--coil-extcur file` sets coil-group currents from a VMEC `&INDATA` input
(`EXTCUR(i) = v`, `EXTCUR = a, b, ...`, slices `EXTCUR(i:j) = ...`, repeat
counts `n*v`) or from a plain list of numbers. As in DIAGNO, each group's
currents are normalised to its first point and scaled by `EXTCUR(g)`; groups
missing from the file are switched off with a warning. Without
`--coil-extcur` the currents in the coil file are used unchanged.

### Magnetic probes and response matrices

`--bprobes file` evaluates DIAGNO magnetic probes (`x y z theta_inc phi_inc
eff_area` per row, angles in degrees; `--rphiz` for `R phi z` positions) and
writes `eff_area * B . n` to `tiago_bprobes.csv`.

`--response-out file` writes every signal per unit `EXTCUR` of each coil group
(`kind,label,group,value`, like DIAGNO's `-mutual`). Vacuum signals for other
currents are then `sum_g M_g EXTCUR_g` without re-running Biot–Savart, so
fitting coil currents is a linear least-squares problem.

### MSE

`--mse file` evaluates motional Stark effect pitch angles and writes
`tiago_mse.csv`; see [mse.md](mse.md).

## Plasma response (VMEC)

`--plasma-wout wout.nc` adds the field of the plasma currents of a VMEC
equilibrium to every flux loop and segmented Rogowski coil. Because the VMEC
boundary is a flux surface, the plasma acts outside it like the sheet current
`mu0 K = n x B` on the boundary (virtual casing). Tiago evaluates the vector
potential and field of this sheet on a grid of `--plasma-nphi` (per field
period) by `--plasma-ntheta` points, default 64 x 64, and warns when a sensor
is closer to the boundary than two grid spacings. As in DIAGNO, loops that link
the plasma poloidally must be flagged `idia = 1`, which adds the toroidal flux
`phiedge`. Omit `--coils` for plasma-only signals.

### Boundary field models

- `--plasma-contravariant` (default) follows DIAGNO.
- `--plasma-covariant` uses VMEC's covariant boundary field
  (`bsubumnc`, `bsubvmnc`) unchanged and reports the Fourier curl norm and a
  bound on the toroidal current ripple. The covariant field alone does not
  guarantee current conservation.
- `--plasma-conservative` additionally projects each Fourier pair
  `(B_u, B_v)` onto `n B_u + m B_v = 0` (with `n` including the field
  periodicity). This is the smallest correction in the Euclidean norm of the
  two coefficients and preserves the `m = n = 0` circulation coefficients; the
  resulting sheet is closed and conserves toroidal current. It changes the
  field model and can differ from DIAGNO at finite resolution.

`plasma_model.txt` in the output directory records the mode and its
diagnostics. Field and shape responses include the projection; covariant
response columns refer to the original `bsub*` coefficients.

### Response matrices

`--plasma-response-out file` writes the derivative of every signal's plasma
part with respect to the boundary field coefficients of the selected model
(`bsup*` or `bsub*`, extrapolated to `s = 1`, per mode `m, n`). At fixed
boundary shape the plasma part is linear in them, so the Jacobian is exact.
`--plasma-shape-response-out file` writes derivatives with respect to the
boundary geometry (`rmnc`, `zmns`, and `rmns`, `zmnc` for `lasym`) at fixed
field coefficients, computed in reverse mode through the sheet current and the
Biot–Savart kernel.

Derivatives with respect to equilibrium parameters (`PHIEDGE`, `CURTOR`,
`PRES_SCALE`, ...) need the response of the equilibrium itself; these come
from `tiago_reconstruct` ([reconstruction.md](reconstruction.md)).

## General plasma field sources

`tiago_plasma_field_source::plasma_field_source_t` is the library interface for
plasma fields. Its batch methods take Cartesian points in metres and return
the exterior plasma B in tesla and A in tesla metres. `plasma_support_t`
implements it with the VMEC virtual casing above, including its responses. A
solver keeps a copy of a source installed with
`solver%set_plasma_source(source)`.

### SPECTRE

With `-DTIAGO_ENABLE_SPECTRE=ON` (see [building.md](building.md)):

```bash
./build-spectre/tiago_vacuum_cli --plasma-spectre equilibrium.h5 \
    --bprobes probes.diagno --output-dir signals
```

This uses libneo's `spectre_field_t`, its coordinate inverse and Cartesian
transforms; the file must pass libneo's chart validity checks. The exterior
plasma field is computed with Tiago's virtual-casing kernels on the outer
plasma interface, so internal islands and interface currents are allowed.

- A free-boundary file uses the vacuum-side boundary field, including a
  physical boundary sheet current. A fixed-boundary file supplies only the
  plasma-side field, which assumes no additional current sheet on the
  boundary.
- The outer boundary must be a flux surface; a sampled relative normal field
  above `1e-6` is rejected.
- Sensors must lie outside the plasma and away from the boundary. Increase
  `--plasma-nphi` and `--plasma-ntheta` to check quadrature convergence.

The interior total field is available from the library:

```fortran
use tiago_spectre_source, only: spectre_source_t
! source is type(spectre_source_t); solver is type(vacuum_solver_t)
call source%init('equilibrium.h5', 64, 64)
call solver%set_plasma_source(source)
call solver%sample_bfield(points, bfield, interior=.true.)
```

Initialise the solver's coils first with `solver%init(coil_file)` (or `''`).
The interior call returns the SPECTRE total field without adding coils;
exterior diagnostics add the coil field once. Interior A, interior VMEC
fields and SPECTRE coefficient or shape responses are not implemented.
libneo cannot invert at the magnetic axis, so use off-axis interior points.

Loops with `idia = 1` need the signed global toroidal flux, which libneo does
not currently read from SPECTRE files. Pass it with
`--plasma-linked-flux <signed_Wb>` (or `linked_flux=` to `source%init`), using
the sign convention of the VMEC correction `phiedge * signgs`. Library callers
apply `finalize_flux_signals` with
`solver%plasma_flux(required=any(loops%idia == 1))`.

Tests cover analytic diagnostic integrals, identity with the VMEC path, SPECTRE
field samples, finite differences of the exterior A, and a nested-surface
cross-check against VMEC. The vacuum fixed-interface test case is not a
force-relaxed island equilibrium. Tolerances and limitations are listed in the
[fixture provenance](../tests/data/spectre/PROVENANCE.md).
