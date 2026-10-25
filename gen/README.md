# Regenerate MSE local kernels

Ordinary TIAGO builds use the checked-in Fortran kernels without a CAS.
To regenerate from the symbolic cross products and the differentiated
`atan2` expression, use a built FortSym tree at commit
`7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6`:

```sh
cmake -S gen -B build-codegen -G Ninja -DFORTSYM_BUILD=/path/to/fortsym/build
cmake --build build-codegen -j4
build-codegen/gen_mse
git diff --exit-code -- src/generated/mse_rows.f90 src/generated/mse_local.f90
```

Run the executable from the TIAGO repository root. FortSym must have its
SymEngine backend and runtime dependencies available. No generated file should
be edited by hand. Validate with TIAGO's `tiago_mse` CTest tests after generation;
symbolic generation does not replace their independent physical and FD oracles.

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
