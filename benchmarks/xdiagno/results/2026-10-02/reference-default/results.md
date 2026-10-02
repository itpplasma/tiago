### Accuracy and performance

Wall time in seconds including coil loading; `n` = 4 MPI ranks (xdiagno) or OpenMP threads (Tiago). *Tiago vs xdiagno* is the relative difference at equal `int_step`/`--samples`; *quadrature error* is each code's median deviation from its own run at 64 samples per segment.

| case | coil pts | flux/seg | samples | xdiagno 1 | xdiagno n | Tiago 1 | Tiago n | Tiago vs xdiagno median / max | missing / non-finite | quad. err xdiagno / Tiago |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| repo_sample | 5 | 4/3 | 2 | 1.49 | 2.10 | 0.01 | 0.00 | 1.1e-08 / 3.5e-08 | 0 / 0 | 3.6e-03 / 3.6e-03 |
| repo_sample | 5 | 4/3 | 6 | 1.47 | 1.64 | 0.00 | 0.00 | 3.4e-08 / 5.0e-08 | 0 / 0 | 3.8e-04 / 3.8e-04 |
| repo_sample | 5 | 4/3 | 16 | 1.31 | 1.50 | 0.00 | 0.00 | 5.6e-09 / 1.0e-08 | 0 / 0 | 5.0e-05 / 5.0e-05 |
| repo_varying_current | 15 | 4/3 | 2 | 1.91 | 1.83 | 0.01 | 0.00 | 1.5e-08 / 2.8e-08 | 0 / 0 | 3.5e-03 / 3.5e-03 |
| repo_varying_current | 15 | 4/3 | 6 | 1.36 | 1.53 | 0.00 | 0.00 | 3.1e-08 / 4.9e-08 | 0 / 0 | 3.9e-04 / 3.9e-04 |
| repo_varying_current | 15 | 4/3 | 16 | 1.30 | 1.52 | 0.00 | 0.00 | 2.4e-09 / 9.2e-09 | 0 / 0 | 5.1e-05 / 5.1e-05 |
| repo_negative_current | 15 | 4/3 | 2 | 1.31 | 1.49 | 0.00 | 0.00 | 1.7e-08 / 4.9e-08 | 0 / 0 | 2.4e-03 / 2.4e-03 |
| repo_negative_current | 15 | 4/3 | 6 | 1.30 | 1.50 | 0.00 | 0.00 | 3.1e-08 / 1.8e-07 | 0 / 0 | 2.7e-04 / 2.7e-04 |
| repo_negative_current | 15 | 4/3 | 16 | 1.30 | 1.50 | 0.00 | 0.00 | 4.5e-09 / 2.4e-08 | 0 / 0 | 3.5e-05 / 3.5e-05 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 2 | 1.42 | 1.62 | 0.10 | 0.09 | 2.2e-08 / 8.7e-08 | 0 / 0 | 8.6e-03 / 8.6e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 6 | 1.45 | 1.61 | 0.12 | 0.10 | 3.6e-08 / 1.0e-07 | 0 / 0 | 1.4e-03 / 1.4e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 16 | 1.50 | 1.64 | 0.18 | 0.12 | 2.2e-09 / 2.8e-08 | 0 / 0 | 1.8e-04 / 1.8e-04 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 2 | 1.44 | 1.79 | 0.11 | 0.11 | 1.1e-08 / 1.9e-07 | 0 / 0 | 7.3e-03 / 7.3e-03 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 6 | 1.50 | 1.66 | 0.14 | 0.10 | 4.9e-08 / 4.1e-07 | 0 / 0 | 1.1e-03 / 1.1e-03 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 16 | 1.56 | 1.70 | 0.23 | 0.14 | 3.5e-09 / 1.6e-07 | 0 / 0 | 1.5e-04 / 1.5e-04 |
| geom_ncsx | 18690 | 113/96 | 2 | 2.80 | 1.93 | 0.87 | 0.31 | 7.3e-09 / 5.9e-07 | 0 / 0 | 3.6e-04 / 3.6e-04 |
| geom_ncsx | 18690 | 113/96 | 6 | 3.53 | 2.24 | 2.46 | 0.91 | 3.0e-08 / 6.1e-07 | 0 / 0 | 4.2e-05 / 4.2e-05 |
| geom_ncsx | 18690 | 113/96 | 16 | 10.53 | 3.28 | 6.57 | 1.93 | 2.7e-09 / 2.2e-07 | 0 / 0 | 5.7e-06 / 5.7e-06 |
| geom_m16n08 | 33024 | 149/128 | 2 | 3.23 | 2.17 | 1.80 | 0.58 | 2.0e-09 / 7.0e-07 | 0 / 0 | 3.7e-05 / 3.7e-05 |
| geom_m16n08 | 33024 | 149/128 | 6 | 6.28 | 3.63 | 6.91 | 2.68 | 2.6e-08 / 2.3e-07 | 0 / 0 | 4.1e-06 / 4.1e-06 |
| geom_m16n08 | 33024 | 149/128 | 16 | 22.41 | 8.11 | 21.74 | 6.43 | 7.0e-10 / 1.5e-07 | 0 / 0 | 5.5e-07 / 5.5e-07 |

### Quadrature: midpoint vs Gauss-Legendre (Tiago, `--gauss`)

Median deviation from the converged 64-sample midpoint result, and wall time on 4 threads, at equal points per segment.

| case | samples | midpoint error | midpoint time | Gauss error | Gauss time |
|---|---:|---:|---:|---:|---:|
| geom_ncsx | 2 | 3.6e-04 | 0.31 | 4.0e-06 | 0.37 |
| geom_ncsx | 6 | 4.2e-05 | 0.91 | 3.6e-07 | 1.02 |
| geom_ncsx | 16 | 5.7e-06 | 1.93 | 3.8e-07 | 1.87 |
| geom_m16n08 | 2 | 3.7e-05 | 0.58 | 4.8e-08 | 0.59 |
| geom_m16n08 | 6 | 4.1e-06 | 2.68 | 3.7e-08 | 2.71 |
| geom_m16n08 | 16 | 5.5e-07 | 6.43 | 3.7e-08 | 4.70 |

### Magnetic probes and response matrices (NCSX coils, vacuum)

| quantity | Tiago vs xdiagno median / max | missing / non-finite |
|---|---|---|
| 40 B-probes | 2.4e-08 / 5.9e-07 | 0 / 0 |
| response matrix (350 entries, xdiagno -mutual) | 2.3e-08 / 4.3e-04 | 0 / 0 |
| sum_g M_g EXTCUR_g vs Tiago signals | 2.2e-15 / 1.3e-14 | 0 / 0 |

### Plasma response (NCSX, plasma only)

35 flux loops, 25 segmented Rogowskis and 20 B-probes; xdiagno -vmec (adaptive virtual casing, tol 1e-6) on 4 ranks took 371 s. Tiago uses the VMEC boundary sheet current on a grid of `grid` x `grid` points per field period.

| grid | Tiago n threads | flux: median / max vs xdiagno | Rogowski: median / max | B-probe: median / max |
|---:|---:|---|---|---|
| 32 | 0.04 | 3.9e-06 / 4.6e-04 | 2.4e-06 / 7.7e-06 | 2.3e-05 / 3.7e-04 |
| 64 | 0.14 | 4.3e-06 / 4.8e-04 | 1.5e-06 / 7.6e-06 | 1.3e-05 / 4.4e-04 |

Ampere loop around the plasma (x eff_area): mu0 I_tor = 7.633075e-07, xdiagno 7.633317e-07, Tiago 7.633352e-07.

Boundary-field Jacobian (`--plasma-response-out`) vs a finite difference in `bsupvmnc(0,0)` over all 80 signals: median 4.8e-13, max 1.1e-11.

### Equilibrium-parameter derivatives (LI383 low resolution, plasma only)

80 signals; 13 fixed-boundary VMEC runs of 1.6 s each. Central differences with step h and h/2; *chain rule* is the boundary-field Jacobian (`--plasma-response-out`) times the finite difference of the boundary coefficients (plus d phiedge on idia = 1 loops); *dI/dp* is the change of the current enclosed by the AMPERE loop (signal / (mu0 eff_area)) next to that of VMEC's `ctor`.

| parameter | value | h | h vs h/2 median / max | chain rule median / max | dI/dp Tiago / VMEC ctor |
|---|---:|---:|---|---|---|
| PHIEDGE | 0.514386 | 5.14e-05 | 8.2e-08 / 4.5e-06 | 2.4e-10 / 3.9e-08 | -108.8 / -6.2238e-06 |
| CURTOR | -174250 | 1.74e+01 | 5.1e-09 / 2.8e-07 | 3.9e-11 / 1.1e-09 | 1.0156 / 1.0168 |
| PRES_SCALE | 1 | 1.00e-04 | 2.1e-08 / 1.5e-06 | 1.7e-10 / 1.4e-08 | -46.083 / 4.6566e-06 |

Radial resolution: enclosed current minus CURTOR = -174250 A (the converged value). VMEC's `ctor` extrapolates the covariant B_u to the boundary; Tiago's sheet uses the extrapolated contravariant B^u, B^v.

| ns | VMEC ctor - CURTOR [A] | Tiago Ampere loop - CURTOR [A] |
|---:|---:|---:|
| 16 | -2930.7 | -2875.0 |
| 32 | -592.5 | -747.7 |
| 64 | -120.9 | -353.4 |
| 128 | -26.1 | -298.9 |

The remaining difference between the two is angular truncation: the extrapolated contravariant field does not conserve the sheet current exactly (the current through a poloidal cross-section varies with phi), which vanishes with MPOL/NTOR:

| MPOL / NTOR (ns 64) | Tiago Ampere loop - VMEC ctor [A] |
|---|---:|
| 4 / 3 | -232.5 |
| 6 / 4 | +144.0 |
| 8 / 6 | +25.2 |

Relative sensitivity (p / S) dS/dp of some signals:

| signal | PHIEDGE | CURTOR | PRES_SCALE |
|---|---:|---:|---:|
| flux:DIA_00 | -2.945 | -2.623 | +3.284 |
| flux:TOR_02 | -0.093 | +0.824 | +0.134 |
| flux:SAD_00_00 | -0.270 | +0.557 | +0.356 |
| seg:SEG_00_00 | -0.096 | +0.820 | +0.138 |
| seg:AMPERE | +0.000 | +0.999 | +0.000 |

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
| extcur_slice | #11 | EXTCUR(1:2) = 2.0 5.0 (array slice) | flux:SQ=2.4993e-07 | flux:SQ=2.4993e-07 | agree |
| extcur_partial | #11 | only EXTCUR(2:) = 5.0 given | flux:SQ=1.8880e-07 | flux:SQ=1.8880e-07 | agree |
| label_with_space | #12 | label 'Loop A/upper' | flux:Loop A/upper=6.8325e-08 | flux:Loop A/upper=6.8325e-08 | agree |
| segrog_point_area | #12 | per-point eff_area 1e-4 / 5e-4 | seg:SEG_NONUNIF=4.5539e-10 | seg:SEG_NONUNIF=4.5539e-10 | agree |
