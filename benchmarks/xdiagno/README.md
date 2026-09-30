# Tiago vs. STELLOPT xdiagno benchmark

Compares `tiago_vacuum_cli` against the reference code it replaces,
[STELLOPT](https://github.com/PrincetonUniversity/STELLOPT)'s `xdiagno`, on
identical inputs: accuracy, performance, and DIAGNO file-format semantics.

This is a manual benchmark. It is **not** run by `ctest` or CI: building
STELLOPT takes several minutes and the full run takes about 15 minutes on 4 cores.

## Quick start

```bash
# 0. Build Tiago (see the top-level README), giving build/tiago_vacuum_cli
# 1. System packages for STELLOPT (Ubuntu/Debian)
sudo apt-get install gfortran libopenmpi-dev openmpi-bin libscalapack-openmpi-dev \
     libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev libopenblas-dev
# 2. Build xdiagno at a pinned STELLOPT commit      -> _work/bin/xdiagno
#    --patched also builds it with patches/*.patch  -> _work/bin/xdiagno_patched
./build_xdiagno.sh --patched
# 3. Download the public coil sets (pinned commits) -> _work/data/
./fetch_data.sh
# 4. Run everything, or pick suites: repo geometry semantics features plasma
python3 bench.py                 # add --quick for one sample count only
```

Everything is written below `_work/` (git-ignored). Results go to
`_work/results/results.md` (the tables below) and `results.json`.
Running as root under OpenMPI is allowed automatically.

## What is compared

Both codes get the same STELLOPT coil file, the same EXTCUR values, and the same
DIAGNO-format diagnostic files:
- EXTCUR is written as `input.`, reproducing the coil file's own group currents.
- Segmented Rogowskis carry a per-point `eff_area = 3.4e-4 / (npts-1)`.

xdiagno runs with `-vac -coil <file>`, `int_type='midpoint'`, and
`int_step` equal to Tiago's `--samples`.

| suite | contents |
|---|---|
| `repo` | the cases in `tests/`: toy square coil (3 current variants), NCSX nfp=1, NCSX nfp=3 |
| `geometry` | generated sensor sets on a torus between plasma and coils: 12–16 poloidal (diamagnetic) loops, 5 toroidal loops, 96–128 saddle loops, 96–128 segmented Rogowskis. Coil sets: NCSX (18.7k points, 10 groups) and M16N08 (33k points, 256 groups) |
| `semantics` | one small input per DIAGNO-format feature. Checks marked "DIFFER" are open Tiago issues |
| `features` | 40 magnetic probes and the per-coil-group response matrix (350 entries, `xdiagno -mutual`) on the NCSX coils |
| `plasma` | plasma-only signals of the NCSX equilibrium from STELLOPT's `DIAGNO_TEST` (`xdiagno -vmec`, adaptive virtual casing): 35 flux loops (diamagnetic loops with `idia=1`), 25 Rogowskis, and a closed loop checked against Ampère's law with the VMEC toroidal current. xdiagno needs about 6 minutes on 4 ranks here |

Metrics:
- **Tiago vs xdiagno:** relative difference per signal, normalised by
  `max(|xdiagno|, 1e-6 · max|xdiagno|)`.
- **Quadrature error:** each code's median deviation from its own run at 64
  samples per segment. This shows how many samples are actually needed.
- **Wall time:** the whole process, including coil loading. `1` is serial; `n`
  is 4 MPI ranks for xdiagno and 4 OpenMP threads for Tiago.

The geometry suite uses only `iflflg=0`, `idia=0` and closed loops (first
point repeated), so the semantic differences below don't contaminate the
accuracy numbers.

### Patched xdiagno and upstream bugs

Bugs found in DIAGNO while benchmarking are fixed by small patches in
`patches/`, one per bug, each tracked in a Tiago issue with a draft upstream
report (label `upstream`). `build_xdiagno.sh --patched` builds
`_work/bin/xdiagno_patched` with all of them. It is used only where stock
DIAGNO fails: the `ncsx_nfp3` repo case and the `iflflg_period` check.

| patch | DIAGNO bug |
|---|---|
| `diagno_vac_nfp.patch` | `-vac` leaves `nfp = 0`, so every `iflflg=1` flux loop is NaN |
| `diagno_segrog_without_coils.patch` | segmented Rogowskis crash without a coil file (`coil_group` not allocated) |
| `diagno_nextcur_exceeds_coil_groups.patch` | `-vmec` with a coil file crashes when the wout has more EXTCUR values than coil groups |
| `diagno_mut_file_per_diagnostic.patch` | naming any `*_mut_file` makes the B-probe and Mirnov routines read their (unnamed) matrix files and crash |

The `plasma` suite works around the last two with stock xdiagno by passing a
zero-current coil file with one group per EXTCUR value.

If you find a bug in a reference code, add a reproducible patch here, apply it
from `build_xdiagno.sh`, and open a tracking issue with a draft upstream report.

## Results

Run on 2026-09-30 on a 4-core x86_64 VM with gfortran 13.3 and STELLOPT
2f181f0. Tiago was built with the default flags (portable, no `-march=native`);
xdiagno was built with STELLOPT's `make_ubuntu.inc`, which uses `-O2 -march=native`
and `-fcheck=all`.

### Accuracy and performance

Wall time in seconds including coil loading; `n` = 4 MPI ranks (xdiagno) or OpenMP threads (Tiago). *Tiago vs xdiagno* is the relative difference at equal `int_step`/`--samples`; *quadrature error* is each code's median deviation from its own run at 64 samples per segment.

| case | coil pts | flux/seg | samples | xdiagno 1 | xdiagno n | Tiago 1 | Tiago n | Tiago vs xdiagno median / max | missing / non-finite | quad. err xdiagno / Tiago |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| repo_sample | 5 | 4/3 | 2 | 0.40 | 0.46 | 0.01 | 0.01 | 1.1e-08 / 3.5e-08 | 0 / 0 | 3.6e-03 / 3.6e-03 |
| repo_sample | 5 | 4/3 | 6 | 0.40 | 0.48 | 0.01 | 0.01 | 3.4e-08 / 5.0e-08 | 0 / 0 | 3.8e-04 / 3.8e-04 |
| repo_sample | 5 | 4/3 | 16 | 0.40 | 0.48 | 0.00 | 0.00 | 5.6e-09 / 1.0e-08 | 0 / 0 | 5.0e-05 / 5.0e-05 |
| repo_varying_current | 15 | 4/3 | 2 | 0.40 | 0.45 | 0.00 | 0.00 | 1.5e-08 / 2.8e-08 | 0 / 0 | 3.5e-03 / 3.5e-03 |
| repo_varying_current | 15 | 4/3 | 6 | 0.39 | 0.45 | 0.01 | 0.01 | 3.1e-08 / 4.9e-08 | 0 / 0 | 3.9e-04 / 3.9e-04 |
| repo_varying_current | 15 | 4/3 | 16 | 0.40 | 0.45 | 0.00 | 0.00 | 2.4e-09 / 9.2e-09 | 0 / 0 | 5.1e-05 / 5.1e-05 |
| repo_negative_current | 15 | 4/3 | 2 | 0.41 | 0.45 | 0.01 | 0.01 | 1.7e-08 / 4.9e-08 | 0 / 0 | 2.4e-03 / 2.4e-03 |
| repo_negative_current | 15 | 4/3 | 6 | 0.39 | 0.45 | 0.00 | 0.00 | 3.1e-08 / 1.8e-07 | 0 / 0 | 2.7e-04 / 2.7e-04 |
| repo_negative_current | 15 | 4/3 | 16 | 0.40 | 0.48 | 0.01 | 0.00 | 4.5e-09 / 2.4e-08 | 0 / 0 | 3.5e-05 / 3.5e-05 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 2 | 0.59 | 0.67 | 0.12 | 0.11 | 2.2e-08 / 8.7e-08 | 0 / 0 | 8.6e-03 / 8.6e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 6 | 0.62 | 0.67 | 0.13 | 0.12 | 3.6e-08 / 1.0e-07 | 0 / 0 | 1.4e-03 / 1.4e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 16 | 0.70 | 0.67 | 0.20 | 0.13 | 2.2e-09 / 2.8e-08 | 0 / 0 | 1.8e-04 / 1.8e-04 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 2 | 0.58 | 0.68 | 0.11 | 0.10 | 1.1e-08 / 1.0e-07 | 0 / 0 | 1.8e-02 / 1.8e-02 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 6 | 0.63 | 0.65 | 0.14 | 0.11 | 3.9e-08 / 1.9e-07 | 0 / 0 | 2.0e-03 / 2.0e-03 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 16 | 0.73 | 0.72 | 0.21 | 0.13 | 2.1e-09 / 5.0e-08 | 0 / 0 | 3.7e-04 / 3.7e-04 |
| geom_ncsx | 18690 | 113/96 | 2 | 1.63 | 0.88 | 0.98 | 0.33 | 7.3e-09 / 5.9e-07 | 0 / 0 | 3.6e-04 / 3.6e-04 |
| geom_ncsx | 18690 | 113/96 | 6 | 3.59 | 1.47 | 2.79 | 0.79 | 3.0e-08 / 6.1e-07 | 0 / 0 | 4.2e-05 / 4.2e-05 |
| geom_ncsx | 18690 | 113/96 | 16 | 9.24 | 2.94 | 7.43 | 2.15 | 2.7e-09 / 2.2e-07 | 0 / 0 | 5.7e-06 / 5.7e-06 |
| geom_m16n08 | 33024 | 149/128 | 2 | 3.54 | 1.64 | 2.23 | 0.75 | 2.0e-09 / 7.0e-07 | 0 / 0 | 3.7e-05 / 3.7e-05 |
| geom_m16n08 | 33024 | 149/128 | 6 | 9.04 | 3.09 | 6.44 | 1.78 | 2.6e-08 / 2.3e-07 | 0 / 0 | 4.1e-06 / 4.1e-06 |
| geom_m16n08 | 33024 | 149/128 | 16 | 23.10 | 7.02 | 16.74 | 4.44 | 7.0e-10 / 1.5e-07 | 0 / 0 | 5.5e-07 / 5.5e-07 |

### Quadrature: midpoint vs Gauss-Legendre (Tiago, `--gauss`)

Median deviation from the converged 64-sample midpoint result, and wall time on 4 threads, at equal points per segment.

| case | samples | midpoint error | midpoint time | Gauss error | Gauss time |
|---|---:|---:|---:|---:|---:|
| geom_ncsx | 2 | 3.6e-04 | 0.33 | 4.0e-06 | 0.33 |
| geom_ncsx | 6 | 4.2e-05 | 0.79 | 3.6e-07 | 0.81 |
| geom_ncsx | 16 | 5.7e-06 | 2.15 | 3.8e-07 | 2.06 |
| geom_m16n08 | 2 | 3.7e-05 | 0.75 | 4.8e-08 | 0.70 |
| geom_m16n08 | 6 | 4.1e-06 | 1.78 | 3.7e-08 | 1.78 |
| geom_m16n08 | 16 | 5.5e-07 | 4.44 | 3.7e-08 | 4.69 |

### Magnetic probes and response matrices (NCSX coils, vacuum)

| quantity | Tiago vs xdiagno median / max | missing / non-finite |
|---|---|---|
| 40 B-probes | 2.4e-08 / 5.9e-07 | 0 / 0 |
| response matrix (350 entries, xdiagno -mutual) | 2.3e-08 / 4.3e-04 | 0 / 0 |
| sum_g M_g EXTCUR_g vs Tiago signals | 2.2e-15 / 1.2e-14 | 0 / 0 |

### Plasma response (NCSX, plasma only)

35 flux loops, 25 segmented Rogowskis and 20 B-probes; xdiagno -vmec (adaptive virtual casing, tol 1e-6) on 4 ranks took 347 s. Tiago uses the VMEC boundary sheet current on a grid of `grid` x `grid` points per field period.

| grid | Tiago n threads | flux: median / max vs xdiagno | Rogowski: median / max | B-probe: median / max |
|---:|---:|---|---|---|
| 32 | 0.05 | 3.9e-06 / 4.6e-04 | 2.4e-06 / 7.7e-06 | 2.3e-05 / 3.7e-04 |
| 64 | 0.17 | 4.3e-06 / 4.8e-04 | 1.5e-06 / 7.6e-06 | 1.3e-05 / 4.4e-04 |

Ampere loop around the plasma (x eff_area): mu0 I_tor = 7.633075e-07, xdiagno 7.633317e-07, Tiago 7.633352e-07.

Boundary-field Jacobian (`--plasma-response-out`) vs a finite difference in `bsupvmnc(0,0)` over all 80 signals: median 4.8e-13, max 1.1e-11.

### DIAGNO-format semantics

| check | issue | what | xdiagno | Tiago | verdict |
|---|---|---|---|---|---|
| closed_loop | - | reference: closed square loop | flux:SQ=6.8325e-08 | flux:SQ=6.8325e-08 | agree |
| open_polygon | #7 | loop without repeated first point | flux:SQ_OPEN=6.8325e-08 | flux:SQ_OPEN=6.8325e-08 | agree |
| iflflg_period | #8 | iflflg=1 loop over one field period (nfp=3) | flux:TOR_PERIOD=6.6715e-07 | flux:TOR_PERIOD=6.6715e-07 | agree |
| idia_plus1 | #9 | idia=1 on a horizontal loop | flux:HORIZ=6.8325e-08 | flux:HORIZ=6.8325e-08 | agree |
| idia_minus | #9 | idia=-1: subtract flux of loop 1 | flux:SMALL=6.8325e-08<br>flux:BIG_MINUS_SMALL=2.0605e-07 | flux:SMALL=6.8325e-08<br>flux:BIG_MINUS_SMALL=2.0605e-07 | agree |
| zero_length_coil_segment | #10 | duplicated coil point | flux:SQ=6.8325e-08 | flux:SQ=6.8325e-08 | agree |
| extcur_zero | #11 | EXTCUR(2)=0 switches group B off | flux:SQ=3.0566e-08 | flux:SQ=3.0566e-08 | agree |
| extcur_array | #11 | EXTCUR = 1.0, 5.0 (namelist array) | flux:SQ=2.1936e-07 | flux:SQ=2.1936e-07 | agree |
| label_with_space | #12 | label 'Loop A/upper' | flux:Loop A/upper=6.8325e-08 | flux:Loop A/upper=6.8325e-08 | agree |
| segrog_point_area | #12 | per-point eff_area 1e-4 / 5e-4 | seg:SEG_NONUNIF=4.5539e-10 | seg:SEG_NONUNIF=4.5539e-10 | agree |

### Reading the results

- **Accuracy.** On identical inputs, Tiago reproduces xdiagno to about 1e-8
  median and 7e-7 worst case on 470 vacuum signals with up to 33k coil
  points, and the two codes have identical quadrature errors. All ten
  DIAGNO-format semantics checks agree.
- **Performance.** On realistic coil sets Tiago is 1.3–1.4× (serial) and
  1.7–1.9× (4 threads) faster than xdiagno with portable flags. Built with
  `-DCMAKE_Fortran_FLAGS=-march=native`, as xdiagno is, the M16N08 case at 6
  samples takes 3.4 s serial and 1.0 s on 4 threads, against xdiagno's 9.0 and
  3.1 s. The ~0.4 s of xdiagno on the tiny cases is MPI start-up.
- **Quadrature.** Midpoint (DIAGNO's rule) needs about 16 samples per segment
  for a 1e-5 median error. `--gauss` reaches the same with 2, at the same cost
  per point. The Gauss errors of about 4e-7 (NCSX) and 4e-8 (M16N08) are those
  of the 64-sample midpoint reference itself.
- **Response matrices.** The largest relative differences are round-off zeros
  (1e-15 vs 1e-22 at a matrix scale of 6e-6).
- **Plasma.** Tiago's sheet-current model and xdiagno's adaptive virtual casing
  agree to a few 1e-6 median; the diamagnetic loops (`idia=1`) are compared
  after a 0.5 Wb `phiedge` cancels. Tiago's Ampère loop agrees with μ0·I_tor
  to 4e-6 relative (limited by VMEC's half-mesh B), and it runs in 0.2 s
  instead of 6 minutes.
