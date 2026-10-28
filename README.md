# TIAGO: Toolkit for Inference and Analysis of Generalized Observables

Tiago is a modern Fortran toolkit for plasma diagnostic forward models,
equilibrium reconstruction and uncertainty quantification. It is aimed at
stellarator physicists who need DIAGNO-compatible synthetic diagnostics with
exact derivatives, and inference built on them.

Status: research code under active development; interfaces may change.

## Features

- Magnetic diagnostics: reads DIAGNO flux-loop, segmented Rogowski and
  magnetic-probe files and STELLOPT coil files, and evaluates them from coil
  currents with exact straight-filament Biot–Savart kernels. Coil response
  matrices make coil-current fits linear.
- Plasma response from a VMEC equilibrium by virtual casing, with exact
  derivatives with respect to boundary field and shape; optional SPECTRE
  equilibria through libneo. MSE pitch angles.
- Equilibrium reconstruction (`tiago_reconstruct`): Levenberg–Marquardt fit of
  VMEC++ input parameters to magnetic data, with an exact Jacobian from the
  VMEC++ implicit adjoint and reverse-mode diagnostics, and local covariance.
- Recorded-signal models with calibration and derivatives: magnetics, Thomson
  scattering, dispersion interferometry and electron cyclotron emission (cold
  rays, thermal opacity and transfer, receiver).
- Joint inference of equilibrium and kinetic profiles from several
  diagnostics, with archive access over HTTP/REST or from netCDF files.
- Validation against STELLOPT `xdiagno` and finite differences.

## Building

Requirements: CMake 3.24 or newer, Fortran and C compilers, netCDF-C,
libcurl, and Python 3 for the test suite. libneo, fortnum and json-fortran are
fetched at pinned commits.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

`make` and `make test` wrap the same steps. Reconstruction
(`-DTIAGO_ENABLE_RECONSTRUCTION=ON`) additionally needs Clang with the Enzyme
plugin and LAPACK. Build options, libneo revision selection and the Enzyme
toolchain are described in [docs/building.md](docs/building.md).

## Usage

Evaluate flux loops and segmented Rogowski coils from a coil file:

```bash
./build/tiago_vacuum_cli \
    --coils tests/data/coil_sample.neo \
    --flux tests/data/fluxloop_sample.diagno \
    --segrog tests/data/segrog_sample.diagno \
    --output-dir build/tests/output --seg-area 3.40e-4 --samples 8
```

This writes `tiago_flux.csv` and `tiago_segrog.csv`. Adding
`--plasma-wout wout.nc` includes the plasma contribution of a VMEC
equilibrium. `tiago_cli diag lint <file>` checks diagnostic files and sensor
sets, and `tiago_reconstruct input.nml` runs a reconstruction.

## Documentation

Magnetics and reconstruction:
- [Building](docs/building.md)
- [Magnetic diagnostics and plasma response](docs/magnetics.md)
- [Equilibrium reconstruction](docs/reconstruction.md)
- [MSE](docs/mse.md), [sensor sets](docs/sensor-sets.md)

Recorded signals and inference:
- [Recorded signals and calibration](docs/recorded-signals.md),
  [archive adapters](docs/online-archives.md)
- [Recorded magnetics](docs/recorded-magnetics.md),
  [Thomson scattering](docs/thomson.md),
  [dispersion interferometer](docs/dispersion-interferometer.md)
- ECE: [recorded ECE](docs/ece-recorded.md), [rays](docs/ece-rays.md),
  [equilibrium medium](docs/ece-geometry.md), [opacity](docs/ece-opacity.md),
  [transfer and radiometer](docs/ece-transfer.md)
- [Profile diagnostics](docs/profile-diagnostics.md),
  [joint inference](docs/inference.md),
  [joint inference from recorded signals](docs/recorded-joint-inference.md)

Benchmarks: [STELLOPT xdiagno](benchmarks/xdiagno/README.md),
[reconstruction](benchmarks/reconstruction/README.md).

Design and plans: [DESIGN.md](DESIGN.md), [ROADMAP.md](ROADMAP.md). Tiago owns
the inverse-problem layer (diagnostic models, calibration, priors, posterior
uncertainty); forward equilibrium and kinetic physics stay in their own codes.

## Related projects

- [libneo](https://github.com/itpplasma/libneo): coil, field and SPECTRE readers.
- [VMEC++](https://github.com/proximafusion/vmecpp): equilibria and implicit
  adjoint for reconstruction.
- [fortnum](https://github.com/lazy-fortran/fortnum): profile bases and samplers.
- [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT): reference
  `xdiagno` implementation.

## License

MIT, see [LICENSE](LICENSE).

Energetic-particle wall-loss events can be converted to surface loads and raw
infrared signals with the separate [`tiago_wall` API](docs/wall-synthetic.md).
It requires explicit physical source normalization, material assumptions and
camera calibration; the repository supplies synthetic fixtures only.
