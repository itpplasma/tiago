# Tiago — Trimmed Input for Agile Geometric Observables

Tiago is a green-field forward-modeling tool for stellarator magnetic diagnostics.  It
reuses the modern Fortran infrastructure in [`libneo`](https://github.com/itpplasma/libneo)
while reproducing the feature set of STELLOPT's `xdiagno` utility.

The two short-term priorities are:

1. **Vacuum-mode drop-in for DIAGNO** – parse coil/diagnostic geometry, call
   libneo's Biot–Savart kernels, and emit flux-loop / segmented-Rogowski
   predictions fast enough for daily workflows.
2. **Cross-code validation harness** – compare Tiago outputs against
   `xdiagno -vac` on tiny, published coil files so we can track regressions.

> ⚠️ 12 Nov 2025 – Implementation deliberately stops at design documents.
> The Fortran sources, drivers, and tests will be added in staged PRs once the
> plan in [`docs/PLAN.md`](docs/PLAN.md) is approved.

## Repository layout

```
tiago/
├── CMakeLists.txt         # FetchContent bridge into libneo (no sources yet)
├── cmake/                 # Reserved for toolchain helpers
├── docs/                  # Project plan, testing strategy, references
├── scripts/               # Future helper scripts (empty for now)
├── tests/                 # ctest scaffolding & data download recipes
└── README.md              # This file
```

## Building (scaffolding only)

Tiago already bootstraps `libneo` via `FetchContent`.  The following commands
only verify that the dependency graph resolves:

```bash
cmake -S . -B build \
      -DTIAGO_LIBNEO_TAG=main \
      -DTIAGO_LIBNEO_GIT=git@github.com:itpplasma/libneo.git
cmake --build build
```

No Fortran target is compiled yet; `cmake --build` simply ensures libneo is
fetchable.  The placeholder `tiago_tests` target will be replaced with explicit
`ctest` entries in Stage 2 (see [`docs/PLAN.md`](docs/PLAN.md)).

## Documentation

- [`docs/PLAN.md`](docs/PLAN.md) — staged implementation roadmap with file- and
  line-level references into STELLOPT and libneo.
- [`docs/TESTING.md`](docs/TESTING.md) — describes the cross-code validation
  strategy, planned data sources, and how `xdiagno` will be invoked.

## Contributing

1. Discuss intended feature branches on the plasma/dev channel (to avoid
   duplicated work with the SURROBIER + STELLOPT teams).
2. Keep commits small and reference plan stages.
3. Add/update the system tests described in `docs/TESTING.md` before landing
   functionality.
