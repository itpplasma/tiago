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
# 4. Run everything, or pick suites: repo geometry semantics plasma
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

The `plasma` suite works around the last two with stock xdiagno by passing a
zero-current coil file with one group per EXTCUR value.

If you find a bug in a reference code, add a reproducible patch here, apply it
from `build_xdiagno.sh`, and open a tracking issue with a draft upstream report.

## Results

Run on 2026-09-30 with a 4-core x86_64 VM, gfortran 13.3, STELLOPT 2f181f0,
Tiago d60c0fb. Tiago was built with `-O3 -march=native`, plasma off.
STELLOPT's `make_ubuntu.inc` release flags include `-fcheck=all`.

### Accuracy and performance

Wall time in seconds including coil loading; `n` = 4 MPI ranks (xdiagno) or OpenMP threads (Tiago). *Tiago vs xdiagno* is the relative difference at equal `int_step`/`--samples`; *quadrature error* is each code's median deviation from its own run at 64 samples per segment.

| case | coil pts | flux/seg | samples | xdiagno 1 | xdiagno n | Tiago 1 | Tiago n | Tiago vs xdiagno median / max | missing / non-finite | quad. err xdiagno / Tiago |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| repo_sample | 5 | 4/3 | 2 | 0.43 | 0.45 | 0.01 | 0.01 | 1.4e-08 / 5.0e-08 | 0 / 0 | 3.6e-03 / 3.6e-03 |
| repo_sample | 5 | 4/3 | 6 | 0.40 | 0.46 | 0.00 | 0.00 | 4.9e-08 / 6.5e-08 | 0 / 0 | 3.8e-04 / 3.8e-04 |
| repo_sample | 5 | 4/3 | 16 | 0.40 | 0.49 | 0.00 | 0.01 | 9.5e-09 / 2.4e-08 | 0 / 0 | 5.0e-05 / 5.0e-05 |
| repo_varying_current | 15 | 4/3 | 2 | 0.40 | 0.45 | 0.01 | 0.00 | 1.4e-08 / 4.2e-08 | 0 / 0 | 3.5e-03 / 3.5e-03 |
| repo_varying_current | 15 | 4/3 | 6 | 0.42 | 0.48 | 0.01 | 0.01 | 4.5e-08 / 6.4e-08 | 0 / 0 | 3.9e-04 / 3.9e-04 |
| repo_varying_current | 15 | 4/3 | 16 | 0.41 | 0.47 | 0.01 | 0.01 | 1.4e-08 / 2.4e-08 | 0 / 0 | 5.1e-05 / 5.1e-05 |
| repo_negative_current | 15 | 4/3 | 2 | 0.39 | 0.50 | 0.00 | 0.00 | 2.7e-08 / 6.3e-08 | 0 / 0 | 2.4e-03 / 2.4e-03 |
| repo_negative_current | 15 | 4/3 | 6 | 0.41 | 0.46 | 0.01 | 0.01 | 4.5e-08 / 1.9e-07 | 0 / 0 | 2.7e-04 / 2.7e-04 |
| repo_negative_current | 15 | 4/3 | 16 | 0.40 | 0.45 | 0.01 | 0.01 | 1.4e-08 / 3.8e-08 | 0 / 0 | 3.5e-05 / 3.5e-05 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 2 | 0.60 | 0.74 | 0.49 | 0.44 | 1.4e-08 / 9.9e-08 | 0 / 0 | 8.6e-03 / 8.6e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 6 | 0.64 | 0.73 | 0.57 | 0.50 | 4.4e-08 / 1.1e-07 | 0 / 0 | 1.4e-03 / 1.4e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 16 | 0.72 | 0.71 | 0.77 | 0.55 | 1.4e-08 / 2.6e-08 | 0 / 0 | 1.8e-04 / 1.8e-04 |
| repo_ncsx_nfp3 (xdiagno_vacnfp) | 18690 | 9/5 | 2 | 0.56 | 0.64 | 0.52 | 0.49 | 2.1e-08 / 1.2e+00 | 0 / 0 | 2.0e-02 / 2.0e-02 |
| repo_ncsx_nfp3 (xdiagno_vacnfp) | 18690 | 9/5 | 6 | 0.73 | 0.68 | 0.66 | 0.53 | 5.5e-08 / 1.3e+00 | 0 / 0 | 2.0e-03 / 1.6e-03 |
| repo_ncsx_nfp3 (xdiagno_vacnfp) | 18690 | 9/5 | 16 | 0.74 | 0.78 | 0.90 | 0.57 | 1.5e-08 / 1.3e+00 | 0 / 0 | 4.0e-04 / 2.4e-04 |
| geom_ncsx | 18690 | 113/96 | 2 | 1.68 | 0.91 | 3.09 | 1.79 | 1.5e-08 / 6.1e-07 | 0 / 0 | 3.6e-04 / 3.6e-04 |
| geom_ncsx | 18690 | 113/96 | 6 | 4.02 | 1.62 | 7.75 | 4.18 | 4.4e-08 / 6.0e-07 | 0 / 0 | 4.2e-05 / 4.2e-05 |
| geom_ncsx | 18690 | 113/96 | 16 | 9.10 | 2.92 | 19.00 | 9.47 | 1.4e-08 / 2.1e-07 | 0 / 0 | 5.7e-06 / 5.7e-06 |
| geom_m16n08 | 33024 | 149/128 | 2 | 3.54 | 1.60 | 7.75 | 4.68 | 1.2e-08 / 7.2e-07 | 0 / 0 | 3.7e-05 / 3.7e-05 |
| geom_m16n08 | 33024 | 149/128 | 6 | 8.92 | 3.30 | 18.38 | 9.49 | 4.0e-08 / 2.2e-07 | 0 / 0 | 4.1e-06 / 4.1e-06 |
| geom_m16n08 | 33024 | 149/128 | 16 | 23.19 | 7.05 | 44.52 | 22.61 | 1.4e-08 / 1.5e-07 | 0 / 0 | 5.5e-07 / 5.5e-07 |

### DIAGNO-format semantics

| check | issue | what | xdiagno | Tiago | verdict |
|---|---|---|---|---|---|
| closed_loop | - | reference: closed square loop | flux:SQ=6.8325e-08 | flux:SQ=6.8325e-08 | agree |
| open_polygon | #7 | loop without repeated first point | flux:SQ_OPEN=6.8325e-08 | flux:SQ_OPEN=5.1244e-08 | DIFFER |
| iflflg_period | #8 | iflflg=1 loop over one field period (nfp=3) | flux:TOR_PERIOD=6.6715e-07 | flux:TOR_PERIOD=4.3971e-07 | DIFFER |
| idia_plus1 | #9 | idia=1 on a horizontal loop | flux:HORIZ=6.8325e-08 | flux:HORIZ=1.3397e-10 | DIFFER |
| idia_minus | #9 | idia=-1: subtract flux of loop 1 | flux:SMALL=6.8325e-08<br>flux:BIG_MINUS_SMALL=2.0605e-07 | flux:SMALL=6.8325e-08<br>flux:BIG_MINUS_SMALL=1.3910e-09 | DIFFER |
| zero_length_coil_segment | #10 | duplicated coil point | flux:SQ=6.8325e-08 | flux:SQ=nan | DIFFER |
| extcur_zero | #11 | EXTCUR(2)=0 switches group B off | flux:SQ=3.0566e-08 | flux:SQ=6.8325e-08 | DIFFER |
| extcur_array | #11 | EXTCUR = 1.0, 5.0 (namelist array) | flux:SQ=2.1936e-07 | flux:SQ=6.8325e-08 | DIFFER |
| label_with_space | #12 | label 'Loop A/upper' | flux:Loop A/upper=6.8325e-08 | flux:Loop=6.8325e-08 | DIFFER |
| segrog_point_area | #12 | per-point eff_area 1e-4 / 5e-4 | seg:SEG_NONUNIF=4.5539e-10 | seg:SEG_NONUNIF=2.4649e-10 | DIFFER |

### Reading the results

- **Accuracy.** Wherever both codes use the same file semantics, Tiago
  reproduces xdiagno to about 1e-7 relative (max 7e-7), on 400+ signals with up
  to 33k coil points. At equal samples, both codes have identical quadrature
  error, so the Hanson–Hirshman kernel is correct.
- **`repo_ncsx_nfp3` max error ~1.3.** These are the two `iflflg=1` loops;
  Tiago implements the repeat flag differently (#8). They need
  `xdiagno_vacnfp`, because stock DIAGNO returns NaN for them (upstream bug
  #20, see `diagno_vac_nfp.patch`).
- **Performance.** Serial, Tiago is about 2× slower than xdiagno on realistic
  coil sets, and scales worse with threads (≈2× on 4 threads vs ≈3.3× for
  xdiagno with 4 MPI ranks). The causes and fixes are tracked in #19.
  For the tiny `repo_sample` cases, xdiagno's ~0.4 s is MPI start-up.
- **Quadrature.** Midpoint needs about 16 samples per segment for a 1e-5
  median error. xdiagno's `simpson`/`bode` rules reach ~1e-8 with 6.
- **Semantics.** Every `DIFFER` row is an open Tiago issue; the issue column
  links to it. The table is the acceptance test for those fixes.
