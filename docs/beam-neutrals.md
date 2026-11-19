# Neutral beam densities

`tiago_beam_neutrals` (`src/inference/beam_neutrals.f90`) returns the density
of the full, half and third energy components of a neutral heating beam and of
its thermal halo at Cartesian points. Charge-exchange, FIDA and neutron models
share it.

## Model

Component j has the energy E/j per atom of mass m. With power P, power
fractions f_j, axial distance l from the beam origin and transverse offsets
(a, b) in the horizontal and vertical directions of the beam frame,

```text
n_j(x)  = Gamma_j / (v_j pi w_h w_v) exp(-a^2/w_h^2 - b^2/w_v^2) exp(-tau_j(l)),   (1)
Gamma_j = P f_j j / E,        v_j = sqrt(2 E / (j m)),                             (2)
w_k(l)^2 = w_k^2 (1 - l/F_k)^2 + (l alpha_k)^2,     k = h, v,                      (3)
tau_j(l) = int_0^l ne sigma_s(E/(j A), ne, Te, Zeff) ds   along the axis,           (4)
```

with 1/e half-widths w_k at the origin, focal distances F_k and 1/e
divergences alpha_k. The axis attenuation applies to the whole cross-section.
`sigma_s` is the analytic fit of Suzuki et al., Plasma Phys. Control. Fusion
40, 2097 (1998), eqs. (24), (26), (28), with tables 2a and 3a (H, D, T at
B = 5 T, low table below 100 keV/amu) and the impurity tables 2c and 3c
(C, N, O, 1 < Zeff < 5) or 2b and 3b (He, B). The fit holds for
10 keV/amu < E/(jA) < 10 MeV/amu, 1e18 < ne < 1e21 m^-3 and, below
100 keV/amu, E/(100 jA) < Te/keV < E/(2 jA). Outside this range in Te the fit
is extrapolated. `beam_check` rejects component energies below 9 keV/amu.

Eq. (4) is a trapezoid sum over `nodes` equally spaced axis points from the
origin (`beam_axis_nodes`), linear between nodes. The origin must lie outside
the plasma.

## Halo

Beam atoms charge-exchange with thermal ions, density
`n_i = ne (Z - Zeff)/(Z - 1)` for one impurity of charge Z, with the
H+ + H(1s) cross-section of Janev, Reiter and Samm, Report Jül-4105 (2003),
eq. (44) and table 9. Halo atoms are ionized by electrons at
`nu = ne <sigma v>_e(Te)` (Maxwellian average of eq. (14), table 3) and
undergo further charge exchange at `nu_cx = n_i sigma_cx(3 Ti/A) sqrt(6 Ti/m)`.
Their random walk spreads them by `L^2 = (Ti/m)/((nu_cx + nu) nu)` per
direction, so

```text
n_halo(x) = n_i S / nu  exp(-a^2/W_h^2 - b^2/W_v^2) / (pi W_h W_v),
W_k^2 = w_k^2 + 4 L^2,     S = sum_j Gamma_j exp(-tau_j(l)) sigma_cx(E/(jA)),        (5)
```

evaluated with the plasma at x. Proton-impact ionization and charge exchange
of halo atoms with impurity ions are neglected, which overestimates the halo.
Excited beam and halo states are not resolved; their effect on stopping is in
`sigma_s`.

## Interface

```fortran
type(neutral_beam_t) :: beam          ! origin, direction, energy, power, fractions, widths, ...
class(beam_plasma_t) :: plasma        ! extend with eval(x, ne, te, ti, zeff)
call beam_densities(beam, plasma, x, nb, nh, status)   ! nb(3, n), nh(n) [m^-3]
```

For inference the steps are separate: sample the plasma at
`beam_axis_nodes(beam)`, call `beam_attenuation(beam, ne, te, zeff, att)`,
then `beam_density(beam, att, x, nb)`. `beam_density_jvp` and
`beam_density_vjp` give the derivatives of `nb` with respect to (ne, Te, Zeff)
at the axis nodes. The halo has no derivatives yet.

## Verification

`tiago_beam_neutrals` checks

- the stopping fit against ASCOT5's independent implementation
  (itpplasma/ascot5 29e093f8, `src/suzuki.c`) for D beams into D with N7+ at
  seven points from 9.58 to 500 keV/amu, to 1e-9, and its partial derivatives
  against central differences;
- the attenuation in a sin^2 density slab against composite Gauss sums of the
  same integrand, and the unattenuated density by hand;
- the transverse flux integral of `n_j v_j` against `Gamma_j` with focus and
  divergence, to 2e-3;
- the JVP against central differences and the VJP transposition identity;
- the Maxwellian ionization rate against a sum over electron speed;
- the halo atoms per unit length against `n_i S / nu`.
