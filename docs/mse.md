# MSE channel geometry and pitch angles

The point diagnostic computes the **pi-component polarization angle** of the
motional Stark effect (MSE). It follows the beam-polarimetry principle of
F. M. Levinton et al., *Physical Review Letters* **63**, 2060–2063 (1989),
[doi:10.1103/PhysRevLett.63.2060](https://doi.org/10.1103/PhysRevLett.63.2060).
The vector derivation below specifies TIAGO's conventions; the coefficient
numbering and angle zero need not match another instrument's calibration.

## Coordinate and polarization conventions

Positions are Cartesian `(x,y,z)` in metres, with right-handed cylindrical
basis `e_R=(cos(phi),sin(phi),0)`, `e_phi=(-sin(phi),cos(phi),0)`,
`e_Z=(0,0,1)` and `phi=atan2(y,x)`. Magnetic fields are in tesla.
The normalized beam direction `v` points along the neutral velocity. Beam
speed cancels from the angle. The normalized sight direction `n` points from
the emitting volume **towards the observer**; reverse an observer-to-plasma
ray before using it. Both geometrical paths are infinite lines, not rays:
negative line parameters are allowed.

The reference vector `r` defines the image-plane zero axis

```
z0 = normalize(r - n (n . r))
h  = n x z0
```

The default reference is global `+Z`. A different reference supports an optical
roll or vertical viewing geometry. It must not be parallel to `n`.
The angle increases from `z0` towards `h`. For pi emission the electric-vector
polarization is along the projection of `E_motional = v x B` onto this plane:

```
N = h  . (v x B) = (h  x v) . B = a . B
D = z0 . (v x B) = (z0 x v) . B = c . B
gamma = atan2(N, D), reduced modulo pi into [-pi/2, pi/2)
```

The plasma electric field, finite beam width/divergence, line integration,
spectral mixing and optical calibration are omitted. A sigma-component
instrument measures the perpendicular polarization; convert by pi/2 modulo
pi and apply its measured optical calibration. These choices must be supplied
before comparing a synthetic angle with experimental data.

Expanding `B = B_R e_R + B_phi e_phi + B_Z e_Z` gives exactly

```
tan(gamma) = (A1 B_Z + A2 B_R + A3 B_phi) /
             (A4 B_Z + A5 B_R + A6 B_phi)
(A1,A2,A3) = (a . e_Z, a . e_R, a . e_phi)
(A4,A5,A6) = (c . e_Z, c . e_R, c . e_phi)
```

For horizontal beam `v = v_R e_R + v_phi e_phi`, radial viewing `n=e_R`
and vertical reference, `h=-e_phi`, so
`A=(v_R,0,0,0,-v_phi,v_R)` and
`tan(gamma)=v_R B_Z/(v_R B_phi-v_phi B_R)`.
At `B_R=0` and `v_R != 0` this is the usual `B_Z/B_phi` limit.
Dropping `B_R` for an arbitrary three-dimensional stellarator field is not
justified. Reversing `B` or `v` leaves the angle unchanged modulo pi;
reversing `n` changes its sign with this reference convention.

## Closest approach and input format

Given beam line `p+t v` and sightline `s+u n`, let `w=s-p`, `c=v.n` and
`Delta=|v x n|^2`. The least-distance solution is
`t=(v.w-c n.w)/Delta`, `u=(c v.w-n.w)/Delta`.
The sample position is the midpoint of these two closest points. Their
separation must not exceed the channel's `max_gap`. Lines with
`Delta <= 1e-12`, zero/nonfinite directions, an undefined reference projection
or a midpoint on the cylindrical axis (within roundoff) are rejected.
A zero field or vanishing projected motional electric field has no defined
polarization and returns an error, rather than an angle of zero.

The version-1 text format is a positive channel count followed by one complete
channel per line:

```
label px py pz vx vy vz sx sy sz nx ny nz rx ry rz max_gap
```

Directions and references need not be normalized. `max_gap` is nonnegative,
in metres. Labels are unique, contain at most 127 characters, and use only
letters, digits, `_`, `.` and `-`. Blank lines and `#`/`!` comments are ignored.
All 16 numeric columns are required and finite; extra rows/columns are errors.
See [the example](../examples/mse_channels.txt), whose first midpoint is
`(4,0,0)` and second is `(0,4,0.1)` m. **These are illustrative channels,
not measured LHD BL3 geometry.** Later BL3 data can share a beam line across
channels while supplying each sightline, reference and acceptable gap.

## APIs, field sources and responses

`tiago_mse` exports `read_mse_file`, `mse_geometry` (midpoint, six coefficients,
optional gap and Cartesian projection rows) and
`mse_angles(channels, field, angles, ierr, message, field_jacobian)`.
Here `field(3,nchannel)` contains Cartesian B **at those midpoints**;
`field_jacobian(3,nchannel)` is optional and contains `d gamma / d(Bx,By,Bz)`
in rad/T. This interface accepts an analytic field or an internal equilibrium
field supplied by the caller. Check `ierr == 0` before using any output.

`vacuum_solver_t%mse(channels, angles, ierr, message, coil_jacobian,
plasma_jacobian, shape_jacobian)` evaluates the existing coil kernels plus
the selected exterior plasma field source at the midpoints. Its optional Jacobians
have `(nchannel,nparameter)` shape. Coil columns follow EXTCUR group order;
plasma columns follow the existing `mode_column_name`/`shape_column_name` APIs.
Coil derivatives hold the plasma source fixed; they do not re-solve equilibrium.
The contraction is `d gamma / dp = (d gamma / dB) . (dB/dp)`, where
`d gamma / dB = (D a - N c)/(N^2+D^2)`.
These are **local nonlinear angle derivatives**, not linear response matrices
whose weighted sum reconstructs the angles. They are undefined across the
chosen modulo-pi branch cut; use wrapped angle differences for residuals.
The local kernels and derivatives are generated with FortSym; regeneration is
specified in [gen/README.md](../gen/README.md).

The virtual-casing source is physically valid **outside the plasma only**.
Inside, its sheet field cancels the external coil field rather than giving
the physical internal field (see `tiago_plasma_support`). Therefore the CLI's
`--plasma-wout` path cannot predict an internal BL3 MSE measurement. Use
`mse_angles` with a genuine internal field provider for that calculation.
No MSE measurements are wired into `tiago_reconstruct` yet; its current VMEC
boundary-field representation does not supply this internal field.

## CLI and verification

```
build/tiago_vacuum_cli --coils tests/data/mse_coils.coils \
    --mse examples/mse_channels.txt --output-dir build/mse \
    --response-out response.csv
ctest --test-dir build -R 'tiago_mse' --output-on-failure
```

`--mse-out` changes the default `tiago_mse.csv` (`label,gamma_rad`).
`--response-out`, `--plasma-response-out` and `--plasma-shape-response-out`
include `kind=mse` rows in rad/A, rad per boundary-field coefficient, and rad
per boundary-shape coefficient, respectively, at the supplied field state.
For example the two orthogonal square coils produce `gamma=atan(1/2)` at
channel 1, with EXTCUR derivatives `(0.4,-0.2)` rad/A.

The registered tests use `B_phi=2/R`, `B_Z=0.2` with hand-calculated angles at
three radii and toroidal orientations, the explicit electric vector and Stokes
parameters for nonplanar geometry, analytic square-coil ratios, and central
finite differences for Cartesian B, EXTCUR, VMEC boundary-field and shape
coefficients. The plasma tests cover contravariant, covariant and conservative
representations at an exterior point. They validate this field operator and
its derivatives, not an internal equilibrium or experimental calibration.
