# Bootstrap current on G-EQDSK flux surfaces

`tiago_bootstrap` evaluates the bootstrap current of Sauter, Angioni and Lin-Liu
([Phys. Plasmas 6, 2834 (1999)](https://doi.org/10.1063/1.873240), erratum
[Phys. Plasmas 9, 5140 (2002)](https://doi.org/10.1063/1.1517052)) on the
flux surfaces of a G-EQDSK equilibrium. Its rows `jdotb@<psi_N>` serve as an
edge-current prior of the free-boundary fit
([axisymmetric free boundary](axisymmetric-free-boundary.md)).

## Equations

TIAGO's convention is `B = F grad phi + grad psi x grad phi`, with `psi` per
radian and the plasma current along `+phi` (Sauter's `psi` is `-psi`). With
`j = F' B / mu0 + R p' e_phi` the flux-surface average of the equilibrium is

```text
<j.B> = F p' + F' <B^2> / mu0,                                             (1)
```

and the bootstrap current of electrons and singly charged ions is

```text
<j.B>_bs = F [L31 (p_e dln n_e/dpsi + p_i dln n_i/dpsi) + (L31 + L32) p_e dln T_e/dpsi
              + (L31 + alpha L34) p_i dln T_i/dpsi],                        (2)
```

with `L31`, `L32`, `L34`, `alpha` of Sauter's Eqs. 14-17 as functions of the
effective trapped fraction `f_t`, `Z_eff` and the collisionalities `nu_e*`,
`nu_i*` of Eqs. 18b, c (ion charge in `nu_i*` after the erratum, Coulomb
logarithms of Eqs. 18d, e), with `R = (R_max + R_min) / 2` and
`epsilon = (R_max - R_min) / (R_max + R_min)` of the surface.

Averages use the weight `dl / B_p`,

```text
<A> = oint A dl / B_p / oint dl / B_p,                                       (3)
```

on the surface traced along rays from the magnetic axis, `R = R_axis + rho cos(theta)`,
`Z = Z_axis + rho sin(theta)`, with 512 rays and the periodic trapezoidal rule
in `theta`. Each ray marches out in steps of 2 mm to the first crossing of the
label, since near the X-point a longer step can cross the separatrix leg into
the private flux region, where `psi_N < 1` again. The trapped fraction is the
bounce integral of Lin-Liu and Miller
([Phys. Plasmas 2, 1666 (1995)](https://doi.org/10.1063/1.871315), Eq. 1),

```text
f_t = 1 - 3/4 <B^2> int_0^{1/B_max} lambda dlambda / <sqrt(1 - lambda B)>,  (4)
```

after `lambda = (1 - y^2) / B_max`, with 48 Gauss-Legendre points in `y`.

## Usage

```text
tiago_bootstrap equilibrium.geqdsk kinetic.csv relative_sigma psi_N...
```

`kinetic.csv` holds `psi_N,n_e,T_e,T_i,n_i,Z_eff` [m^-3, eV] on a uniform
`psi_N` grid from 0 to 1. Standard output gets `label,value,sigma` rows
`jdotb@<psi_N>` with `sigma = relative_sigma |<j.B>_bs|`, standard error the
surface quantities, (1), (2), `F` and `<B^2>`. In the measurements of
`tiago_axisym_fit` such rows are predicted by (1) on the solved equilibrium.

## Check

`tests/bootstrap_reference.py` compares with TORAX 2.0.0 (Apache 2.0) on the
FreeGS `TestTokamakSensor` equilibrium. TORAX evaluates (4) on matplotlib
contours of `psi_N` and (2) with its own Sauter coefficients
(`tests/freegs/generate_bootstrap_reference.py`). At nine labels with
`nu_e*` from 0.03 to 100:

| quantity | largest relative difference |
|---|---|
| `f_t` | `5e-4` |
| `nu_e*` | `4e-3` (at `psi_N = 0.2`, from the interpolation of `q`) |
| `<j.B>_bs` | `1.5e-3` |
| `<j.B>` of (1), `tiago_bootstrap` and `jdotb` rows of `tiago_axisym_fit` | `5e-5` of the largest `abs(<j.B>)` |

The reference evaluates (1) with its own contour averages. Near the edge
`F p'` and `F' <B^2> / mu0` cancel, so the difference is taken relative to the
largest value.

## Limits

Only one main ion species of charge 1. The Ohmic part of `<j.B>` is not
included.
