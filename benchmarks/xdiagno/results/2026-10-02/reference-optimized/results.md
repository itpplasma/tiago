### Accuracy and performance

Wall time in seconds including coil loading; `n` = 4 MPI ranks (xdiagno) or OpenMP threads (Tiago). *Tiago vs xdiagno* is the relative difference at equal `int_step`/`--samples`; *quadrature error* is each code's median deviation from its own run at 64 samples per segment.

| case | coil pts | flux/seg | samples | xdiagno 1 | xdiagno n | Tiago 1 | Tiago n | Tiago vs xdiagno median / max | missing / non-finite | quad. err xdiagno / Tiago |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| repo_sample | 5 | 4/3 | 2 | 1.43 | 1.87 | 0.01 | 0.02 | 1.1e-08 / 3.5e-08 | 0 / 0 | 3.6e-03 / 3.6e-03 |
| repo_sample | 5 | 4/3 | 6 | 1.47 | 1.74 | 0.01 | 0.01 | 3.4e-08 / 5.0e-08 | 0 / 0 | 3.8e-04 / 3.8e-04 |
| repo_sample | 5 | 4/3 | 16 | 1.61 | 2.74 | 0.02 | 0.04 | 5.6e-09 / 1.0e-08 | 0 / 0 | 5.0e-05 / 5.0e-05 |
| repo_varying_current | 15 | 4/3 | 2 | 1.42 | 1.87 | 0.01 | 0.03 | 1.5e-08 / 2.8e-08 | 0 / 0 | 3.5e-03 / 3.5e-03 |
| repo_varying_current | 15 | 4/3 | 6 | 1.44 | 1.93 | 0.01 | 0.03 | 3.1e-08 / 4.9e-08 | 0 / 0 | 3.9e-04 / 3.9e-04 |
| repo_varying_current | 15 | 4/3 | 16 | 1.58 | 1.97 | 0.03 | 0.01 | 2.4e-09 / 9.2e-09 | 0 / 0 | 5.1e-05 / 5.1e-05 |
| repo_negative_current | 15 | 4/3 | 2 | 1.49 | 1.99 | 0.01 | 0.01 | 1.7e-08 / 4.9e-08 | 0 / 0 | 2.4e-03 / 2.4e-03 |
| repo_negative_current | 15 | 4/3 | 6 | 1.48 | 1.86 | 0.01 | 0.01 | 3.1e-08 / 1.8e-07 | 0 / 0 | 2.7e-04 / 2.7e-04 |
| repo_negative_current | 15 | 4/3 | 16 | 1.47 | 1.83 | 0.01 | 0.05 | 4.5e-09 / 2.4e-08 | 0 / 0 | 3.5e-05 / 3.5e-05 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 2 | 2.30 | 2.58 | 0.23 | 0.20 | 2.2e-08 / 8.7e-08 | 0 / 0 | 8.6e-03 / 8.6e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 6 | 1.81 | 2.17 | 0.28 | 0.62 | 3.6e-08 / 1.0e-07 | 0 / 0 | 1.4e-03 / 1.4e-03 |
| repo_ncsx_nfp1 | 18690 | 8/4 | 16 | 2.46 | 3.91 | 1.18 | 0.74 | 2.2e-09 / 2.8e-08 | 0 / 0 | 1.8e-04 / 1.8e-04 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 2 | 2.43 | 3.08 | 0.42 | 0.31 | 1.1e-08 / 1.9e-07 | 0 / 0 | 7.3e-03 / 7.3e-03 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 6 | 2.73 | 3.44 | 0.45 | 0.38 | 4.9e-08 / 4.1e-07 | 0 / 0 | 1.1e-03 / 1.1e-03 |
| repo_ncsx_nfp3 (xdiagno_patched) | 18690 | 9/5 | 16 | 1.96 | 3.02 | 0.58 | 0.43 | 3.5e-09 / 1.6e-07 | 0 / 0 | 1.5e-04 / 1.5e-04 |
| geom_ncsx | 18690 | 113/96 | 2 | 3.77 | 4.05 | 1.96 | 1.53 | 7.3e-09 / 5.9e-07 | 0 / 0 | 3.6e-04 / 3.6e-04 |
| geom_ncsx | 18690 | 113/96 | 6 | 10.89 | 3.92 | 5.05 | 2.11 | 3.0e-08 / 6.1e-07 | 0 / 0 | 4.2e-05 / 4.2e-05 |
| geom_ncsx | 18690 | 113/96 | 16 | 13.25 | 7.23 | 17.46 | 6.30 | 2.7e-09 / 2.2e-07 | 0 / 0 | 5.7e-06 / 5.7e-06 |
| geom_m16n08 | 33024 | 149/128 | 2 | 4.86 | 2.50 | 3.48 | 1.39 | 2.0e-09 / 7.0e-07 | 0 / 0 | 3.7e-05 / 3.7e-05 |
| geom_m16n08 | 33024 | 149/128 | 6 | 14.28 | 5.00 | 8.37 | 1.93 | 2.6e-08 / 2.3e-07 | 0 / 0 | 4.1e-06 / 4.1e-06 |
| geom_m16n08 | 33024 | 149/128 | 16 | 18.15 | 5.82 | 15.85 | 3.74 | 7.0e-10 / 1.5e-07 | 0 / 0 | 5.5e-07 / 5.5e-07 |

### Quadrature: midpoint vs Gauss-Legendre (Tiago, `--gauss`)

Median deviation from the converged 64-sample midpoint result, and wall time on 4 threads, at equal points per segment.

| case | samples | midpoint error | midpoint time | Gauss error | Gauss time |
|---|---:|---:|---:|---:|---:|
| geom_ncsx | 2 | 3.6e-04 | 1.53 | 4.0e-06 | 1.60 |
| geom_ncsx | 6 | 4.2e-05 | 2.11 | 3.6e-07 | 1.24 |
| geom_ncsx | 16 | 5.7e-06 | 6.30 | 3.8e-07 | 4.48 |
| geom_m16n08 | 2 | 3.7e-05 | 1.39 | 4.8e-08 | 1.73 |
| geom_m16n08 | 6 | 4.1e-06 | 1.93 | 3.7e-08 | 1.83 |
| geom_m16n08 | 16 | 5.5e-07 | 3.74 | 3.7e-08 | 3.82 |
