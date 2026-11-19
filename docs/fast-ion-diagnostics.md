# Fast-ion diagnostics: D-D neutrons and FIDA

Forward models of the D-D neutron rate and of fast-ion D-alpha (FIDA) spectra
from a gyrotropic fast-ion distribution and the sampled plasma profiles, with
exact derivatives for the profile-diagnostic contract
([inference](inference.md#profile-diagnostic-interface)).

## Fast-ion distribution

`tiago_fast_ion_distribution` holds f(s, E, p): fast ions per m^3, per eV
and per unit pitch p = v_par/v, on flux-label nodes s_k and energy-pitch
cells (E_i, p_j) of widths dE_i, dp_j. Between nodes f is linear in s, above
the last node it vanishes. The density in cell (i, j) is

$$ n_{ij}(s) = f(s, E_i, p_j)\, dE_i\, dp_j . \tag{1}$$

An orbit or Fokker-Planck code supplies f. The text format is documented in
the module header and read by `read_fast_ion_distribution`.

## D-D neutrons

`tiago_neutron_rate` (`neutron_t`) sums three emissivities [m^-3 s^-1]:

$$ \varepsilon_\mathrm{tn} = \tfrac12 n_D^2 \langle\sigma v\rangle_\mathrm{BH}(T_i), \tag{2}$$

$$ \varepsilon_\mathrm{bt} = n_D \sum_i F_i(s)\, \langle\sigma u\rangle(v_i, T_i),
\qquad F_i = \sum_j n_{ij}, \tag{3}$$

$$ \varepsilon_\mathrm{bb} = \tfrac12 \sum_{a,b} n_a n_b\,
\overline{\sigma(|\mathbf v_a-\mathbf v_b|)\,|\mathbf v_a-\mathbf v_b|}^{\,\psi} . \tag{4}$$

Symbols: n_D = n_e (Z_imp - Z_eff)/(Z_imp - 1) the deuterium density with one
impurity of charge Z_imp, <sigma v>_BH the Bosch-Hale reactivity, sigma the
Bosch-Hale D(d,n)3He cross-section at the centre-of-mass energy m_D u^2/4
([Bosch and Hale 1992](https://doi.org/10.1088/0029-5515/32/4/I07), Eqs. 8, 9,
12, Tables IV, VII). a and b run over cells, and the bar is the average over
the gyro-phase difference psi with

$$ |\mathbf v_a-\mathbf v_b|^2 = v_a^2+v_b^2-2v_av_b\,(p_ap_b+q_aq_b\cos\psi),
\qquad q=\sqrt{1-p^2}. \tag{5}$$

The rate on Maxwellian deuterons at rest with thermal speed w = sqrt(2 T_i/m_D)
is

$$ \langle\sigma u\rangle(v,T_i) = \frac{1}{\sqrt\pi\, v}\int_{-\infty}^{\infty}
g(v+wt)\, e^{-t^2}\, dt, \qquad g(u)=\mathrm{sgn}(u)\,u^2\sigma(|u|), \tag{6}$$

evaluated with 96 Gauss-Legendre nodes on |t| <= 7. Its T_i derivative follows
from

$$ \partial_w\langle\sigma u\rangle = \frac{1}{\sqrt\pi\, v}\int t\, g'(v+wt)\,
e^{-t^2}\,dt, \qquad g'(u)=2|u|\sigma+u^2\,\partial_u\sigma . \tag{7}$$

Because f is linear in s, Eq. (4) is a quadratic form of the node
distributions. Its node and neighbour terms are computed once at setup.

A channel returns scale_c sum_m w_m eps(x_m). `neutron_volume_points` gives
the axisymmetric total rate with w = 2 pi R dR dZ, and `neutron_chord_points`
the line-integrated emission per steradian with w = dl/(4 pi).
`tiago_neutrons` evaluates Eqs. (2)-(4) on flux-surface profiles and
integrates them over the volume.

Not modelled: target rotation, emission anisotropy and detector response.

## FIDA

`tiago_fida` (`fida_t`) predicts the spectral radiance
[photons s^-1 m^-2 sr^-1 nm^-1] of channel c in wavelength bin k:

$$ S_{ck} = \frac{\mathrm{scale}_c}{4\pi\, \Delta\lambda_k}\int_c dl
\sum_{ij} n_{ij}(s)\, \Big\langle Y\nu\,\chi_k(\lambda) \Big\rangle_\gamma , \tag{8}$$

$$ \nu = \sum_l n_l\,\sigma_\mathrm{cx}(|\mathbf v-\mathbf v_l|)\,|\mathbf v-\mathbf v_l|
+ n_h\langle\sigma_\mathrm{cx} u\rangle(v, T_i), \tag{9}$$

$$ \lambda = \lambda_0\,(1+\mathbf v\cdot\hat{\mathbf d}/c). \tag{10}$$

Symbols: v the fast-ion velocity at gyro-phase gamma about the field
direction, n_l and v_l the density and velocity of beam component l, n_h the
halo density, sigma_cx the H+ + H(1s) charge-exchange cross-section (Janev and
Smith 1993, Eq. 44, Table 9, coefficients as in FIDASIM's `p_cx_1_janev` and
`tiago_hydrogen_1s_cross_sections`, regenerated here with its derivative),
Y the photons per charge exchange, lambda_0 = 656.1 nm, d the unit sight
direction away from the detector and chi_k the indicator of bin k. The
gyro-phase average uses `ngyro` midpoints.

`fida_set_beams` evaluates the component and halo densities of each beam at
the sample points with the Suzuki-stopping model of
[beam neutrals](beam-neutrals.md), for a given plasma. `fida_set_neutrals`
accepts densities and velocities from another source. Neither is
differentiated with respect to the sampled profiles.

Weight functions: `fida_weights` returns W_ck(E_i, p_j) such that
S_ck = sum_ij W_ck n_ij for densities n_ij uniform along the sight line.

Not modelled: the collisional-radiative population of n = 3 (absorbed in Y),
Stark splitting, the instrument function and the backgrounds.

## Verification

| Check | Reference | Result |
|---|---|---|
| sigma at E_cm = 10, 30, 100, 400 keV | Bosch-Hale Table V | max rel. error 5.5e-5 |
| <sigma v>_BH at T = 1-20 keV | Bosch-Hale Table VIII (R-matrix integration) | 1.4e-4 |
| Eq. (3) for a Maxwellian population at 5 keV, divided by n_f n_D | Table VIII | 2.0e-3 |
| Eq. (4) for the same population, divided by n_f^2/2 | Table VIII | 2.3e-4 |
| Eq. (6) at T_i = 0.01 eV | sigma(E/2) v | 5.9e-7 |
| FIDA spectrum of one (E, p) cell, halo only | closed-form gyro-phase distribution, [Salewski et al. 2011](https://doi.org/10.1088/0029-5515/51/8/083014) | 1.4e-3 |
| FIDA spectrum of an isotropic monoenergetic population | flat box lambda_0 (1 +- v/c) | 6e-3, zero outside |
| sum W n against S, uniform distribution | Eq. (8) | 2e-16 |
| jvp, vjp of both diagnostics | central differences, transposition | 8e-10, 3e-16 |

The Bosch-Hale fits agree with their R-matrix tables to 0.3 %, which bounds
the Maxwellian checks. Tests: `tests/test_fast_ion_neutrons.f90`,
`tests/test_fida.f90`, `tests/neutrons_cli.cmake`. The pointwise cross-sections
and reactivity are generated with FortSym by `gen/gen_fast_ion.f90`.

## Open

- Comparison with FIDASIM on its synthetic tokamak case (spectra, weights,
  neutron rate).
- Derivatives of the beam and halo densities with respect to n_e, T_e, T_i
  and Z_eff in the FIDA chain.
- n = 3 population, Stark splitting and instrument function for absolute FIDA
  radiance.
