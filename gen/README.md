# Regenerate the generated kernels

Ordinary TIAGO builds use the checked-in Fortran kernels without a CAS.
To regenerate from the symbolic cross products and the differentiated
`atan2` expression, use a built FortSym tree at commit
`7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6`:

```sh
cmake -S gen -B build-codegen -G Ninja -DFORTSYM_BUILD=/path/to/fortsym/build
cmake --build build-codegen -j4
build-codegen/gen_mse
build-codegen/gen_profile_kernels
git diff --exit-code -- src/generated
```

`gen_profile_kernels` (FortSym `48dda45a34c0da0db1c78cd4130b3a97375e3c8e`)
emits the cubic Hermite basis jets, the shifted-circle test state with its
label and parameter derivatives, the electron-cyclotron resonance field and
the coronal XICS line emissivity with its partial derivatives.
At that FortSym revision `gen_mse` reproduces the checked-in MSE kernels.

Run the executable from the TIAGO repository root. FortSym must have its
SymEngine backend and runtime dependencies available. No generated file should
be edited by hand. Validate with TIAGO's `tiago_mse` and
`tiago_profile_instruments` CTest tests after generation; symbolic generation
does not replace their independent physical and FD oracles.

# Regenerate the kinetic pressure kernel

`gen_kinetic` emits `src/generated/kinetic_pressure.f90` (pressure from ne,
Te, Ti, Zeff and p_fast with its partial derivatives). It was generated with
FortSym `48dda45a34c0da0db1c78cd4130b3a97375e3c8e`:

```sh
cmake -S gen -B build-codegen -G Ninja -DFORTSYM_BUILD=/path/to/fortsym/build
cmake --build build-codegen --target gen_kinetic
build-codegen/gen_kinetic
git diff --exit-code -- src/generated/kinetic_pressure.f90
```

`tiago_kinetic_profiles` checks the kernel against quasi-neutrality solved
independently.

# Regenerate the recorded-signal reference kernel

`gen_recorded` emits shared gain/offset signals and their partial derivatives.
It was run with FortSym commit
`21d157b394e775a999a8d9ed0be409eac1eb9574`:

```sh
cmake --build build-codegen --target gen_recorded
build-codegen/gen_recorded
git diff --exit-code -- src/generated/recorded_affine.f90
```

`tiago_recorded_diagnostic` independently checks analytic values, the zero-input
limit and derivatives through the existing profile-diagnostic contract.

# Regenerate the recorded Thomson kernels

`gen_thomson` uses FortSym `21d157b394e775a999a8d9ed0be409eac1eb9574` and emits
`ts_*.f90`: the relativistic photon spectrum, Naito depolarization integrand,
scaled Bessel derivative chain, photoelectron and voltage maps, Rayleigh
calibration and count-noise variance, with their required derivative rows.

```sh
cmake --build build-codegen --target gen_thomson
build-codegen/gen_thomson
git diff --exit-code -- src/generated/ts_*.f90
```

To run the exact regeneration comparison in CTest, configure the ordinary TIAGO
build with `-DTIAGO_FORTSYM_BUILD=/path/to/fortsym/build`, then build and run
`ctest --test-dir build -R thomson`. Without this optional path, the physical
and derivative tests still run against the committed kernels. The generator
uses FortSym's native expression differentiation and emitter; it needs no
symbolic simplification backend.

# Regenerate recorded magnetics and energy kernels

`gen_magnetics`, using FortSym `21d157b394e775a999a8d9ed0be409eac1eb9574`,
emits `integrator_signal.f90`, `diamagnetic_energy.f90` and
`thermal_energy_density.f90` under `src/generated/`. It derives every partial
row from the forward expression; the energy diagnostic reuses the existing
kinetic-pressure kernel rather than duplicating its ion closure.

```sh
cmake --build build-codegen --target gen_magnetics
build-codegen/gen_magnetics
git diff --exit-code -- src/generated
```

With `-DTIAGO_FORTSYM_BUILD=/path/to/fortsym/build` in the native TIAGO build,
`ctest -R tiago_magnetic_codegen` regenerates these files in the build tree and
checks exact equality. Ordinary builds consume the committed kernels.

# Regenerate the cold ECE ray kernels

`gen_ece` uses FortSym `21d157b394e775a999a8d9ed0be409eac1eb9574` to emit
Appleton-Hartree O/X Hamiltonians, their gradients and Hessians, the regular
unmagnetized limit and a synthetic divergence-free polynomial test medium.
The value, gradient and Hessian share generated subexpressions.

```sh
cmake --build build-codegen --target gen_ece
build-codegen/gen_ece
cmake -S . -B build -DTIAGO_FORTSYM_BUILD=/path/to/fortsym/build
ctest --test-dir build -R '^tiago_ece' --output-on-failure
```

The regeneration test runs in a temporary directory and compares all four
outputs without changing the checkout. Ordinary library builds need no
computer algebra system.
The physical scope and derivative domain are in [docs/ece-rays.md](../docs/ece-rays.md).

# Regenerate the dispersion-interferometer kernels

`gen_interferometer` emits phase, calibrated I/Q, calibrated atan2 and
phase-to-column maps with complete partial-derivative rows, using FortSym
`21d157b394e775a999a8d9ed0be409eac1eb9574`. It uses the same build as above:

```sh
cmake --build build-codegen --target gen_interferometer
build-codegen/gen_interferometer
git diff --exit-code -- src/generated/dispersion_*.f90
```

The optional native regeneration test runs the generator in an isolated build
directory and byte-compares all four files. Configure the normal TIAGO build
with `-DFORTSYM_BUILD=/path/to/fortsym/build`, build, then run
`ctest --test-dir build -R tiago_interferometer_generated --output-on-failure`.
FortSym and TIAGO must use compatible Fortran compilers for this option.
Ordinary builds retain no CAS dependency. Independent physical and derivative
oracles are described in [the instrument documentation](../docs/dispersion-interferometer.md).

# Regenerate thermal ECE transfer and receiver kernels

`gen_ece_transfer` uses FortSym
`21d157b394e775a999a8d9ed0be409eac1eb9574` and emits six kernels:
`ece_loss`, `ece_loss_small`, `ece_thermal_cell`, `ece_reflection_loss`,
`ece_two_stream_wall`, and `ece_radiometer_noise`, including every partial row.
The small-depth alternating polynomial and its derivatives share one source.

```sh
cmake --build build-codegen --target gen_ece_transfer
build-codegen/gen_ece_transfer
ctest --test-dir build -R tiago_ece_transfer --output-on-failure
```

Configure the ordinary build with `TIAGO_FORTSYM_BUILD` to include the exact
regeneration comparison. These are thermal transfer and receiver kernels;
they do not generate an ECE absorption coefficient. See
[the component contract and limits](../docs/ece-transfer.md).

# Regenerate the equilibrium ray-medium kernels

`gen_ece_geometry` uses FortSym
`21d157b394e775a999a8d9ed0be409eac1eb9574` for Fourier geometry through third
label derivatives, higher Hermite basis derivatives, the contravariant vector
field contraction, density jets and inverse-map spatial/parameter rows.

```sh
cmake --build build-codegen --target gen_ece_geometry
build-codegen/gen_ece_geometry
ctest --test-dir build -R tiago_ece_geometry --output-on-failure
```

With `TIAGO_FORTSYM_BUILD` set, `tiago_ece_geometry_codegen` regenerates all
nine outputs in an isolated build directory and compares the committed files.
The [medium contract](../docs/ece-geometry.md) records the knot, axis and
provider-adjoint limitations.
