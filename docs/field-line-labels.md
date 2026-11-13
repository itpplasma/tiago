# Field-line labels

Kinetic profiles are constant on flux surfaces. In a field with islands the
flux surfaces are no longer the nested surfaces of an equilibrium code, so a
diagnostic point cannot be given the label s of a nested chart. The
field-line label replaces it: the field line through the point is followed
for many toroidal turns, and its label is the mean of the nested label
s = rho^2 of a reference chart over its returns to the starting plane.

- On a surface of the chart the label is s itself.
- On a flux surface deformed by a perturbation it is the mean of s over that
  surface.
- Inside an island all points of one island tube share a label, so a profile
  Te = f(label) is flat across the island; on the X side of an island the
  label changes faster with s than s itself.

The label is a model of the profile, not a measurement. It assumes that the
profile is a function of the flux surface alone and that parallel transport
equilibrates each surface or island tube; it carries no island-specific
transport.

## Library

`tiago_field_line_label` (`src/inference/field_line_label.f90`):

- `line_field_t`: abstract field, `eval(R, phi, Z, b, inside)` with
  b = (B_R, B_phi, B_Z).
- `grid_field_t`: a field tabulated on a cylindrical grid over one field
  period, tricubic splines periodic in phi (libneo `batch_interpolate`).
  `init` takes the three components, `load_mgrid(path, currents, ...)` reads
  a makegrid file through the netCDF C library and sums the coil groups
  with the given currents (in raw mode divided by `raw_coil_cur`). Lines
  closer than `margin` (default 2 cm) to the R or Z edge of the grid are lost.
- `trace_returns`: classical fourth-order Runge-Kutta in phi, with the
  generated direction dR/dphi = R B_R / B_phi, dZ/dphi = R B_Z / B_phi
  (`gen/gen_field_line.f90`).
- `field_line_labels(field, chart, R, Z, phi0, turns, steps, out)`: start
  points in the plane phi0; the chart is any `real_space_state_t`, for
  example `vmec_state_t`. `out` holds s at the start, the label, the standard
  deviation of s over the returns, the number of returns and whether the
  line was lost (left the grid or could not be located in the chart). The
  loop over start points is OpenMP-parallel.

The memory of a `grid_field_t` is that of the tricubic coefficients,
64 values per node and component: 2.4 GB for a 65 x 360 x 65 mgrid.

## Command line

```text
tiago_field_line_labels input.nml

&field_line_labels
  mgrid = 'mgrid.nc'
  currents = -4.95e6, -4.95e6, 2.589e6, 6.2535e5, -2.695e6, 0.0  ! one per group
  scan_group = 6                      ! group replaced by each scan value
  scan_currents = 0.0, 500.0, 1000.0  ! 0 for scan_group: a single case
  wout = 'wout.nc'                    ! reference chart
  radii = 'radii.txt'                 ! start radii [m] at Z = 0, one per line
  phi_deg = 18.0
  turns = 200, steps = 720
  output = 'labels.json'
/
```

The JSON output holds `R_m`, `s_own` (s of the start points), one entry per
scan current in `cases` (`current_A`, `label`, `label_spread`, `returns`,
`lost`; `null` where undefined), the TIAGO revision and the SHA-256 of the
inputs.

## Verification

`tiago_field_line_label` uses an analytic field on circles about
(R0, 0) whose lines follow the Hamiltonian
H = iota0 p + iota1 p^2 / 2 + eps p cos(2 theta - phi), p = r^2 / 2, so that
K = H - p / 2 is constant along each line and a 2/1 island sits at
iota = 1/2. The chart is the circular state, s = r^2 / a^2.

- Without the perturbation, the label equals s of the start point to 1e-9.
- For circulating lines, the label agrees within 1e-5 (4000 turns) with the
  mean of s over the invariant curve K = const, weighted by
  dxi / (dxi/dphi) and computed by quadrature without tracing. The
  difference falls as 1/turns (9e-5 at 250 turns, 2e-6 at 4000).
- Inside the island, two lines started at mirror points about the O point
  get one label within 2e-4, closer to the resonant s than their start.
- The field tabulated on a 61 x 72 x 61 grid reproduces the analytic field
  to 1.3e-6 relative to B_phi and its labels to 6e-8.
- An mgrid with two coil groups, read with currents (2, -0.5), reproduces
  the weighted sum at the grid nodes, also three field periods away.
