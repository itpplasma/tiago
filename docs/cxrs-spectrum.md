# CXRS spectra and line moments

`tiago_cxrs_spectrum` (`src/inference/cxrs_spectrum.f90`) predicts the
charge-exchange spectrum of a hydrogen-like impurity line, and its brightness,
apparent line-of-sight velocity and apparent temperature, along straight sight
lines that cross one or more [neutral beams](beam-neutrals.md). The default
line is N VII n = 9 -> 8 at 566.937 nm.

## Local emission

An impurity of charge Z, density n_Z, temperature T_Z and flow V captures an
electron from beam component j (velocity u_j, density n_j) or from the halo
(density n_h). With the effective emission rate q(E) of the line, tabulated
against the collision energy per amu E and supplied by the user, and e the unit
vector along the sight line away from the observer,

```text
eps_j F_j(v_l) = n_Z n_j / (4 pi) int f_Z(v) q(E(|v - u_j|)) delta(v_l - v.e) d^3v.   (1)
```

The weight q(|v - u_j|) shifts and narrows the distribution of emitting ions
(cross-section effect, von Hellermann et al., Plasma Phys. Control. Fusion 37,
71 (1995)). Writing v = V + w with w Maxwellian of per-direction variance
s^2 = T_Z/m_Z, e1 along u_j - V and e = a e1 + b e2,

```text
M0 = <q>,   M1 = a <q w1>,   M2 = a^2 <q w1^2> + (b^2/2) <q w_perp^2>,               (2)
```

integrated by 16-point Gauss-Hermite quadrature in w1 and 10-point
Gauss-Laguerre quadrature in w_perp^2/(2 s^2). The local mean velocity is
V.e + M1/M0 and the local variance M2/M0 - (M1/M0)^2. Halo atoms are thermal at
the ion temperature and flow, so their weight is q at the mean collision
energy E = 1.5 T_i (1/A_D + 1/A_Z) and has no velocity dependence. Gyro-motion
during the lifetime of the excited state is not modelled.

## Rest line

The rest line is the zero-field fine-structure pattern of the n -> n'
transition: Dirac energies with the reduced mass, statistically populated
upper j levels, and component intensities proportional to the hydrogenic line
strengths

```text
S(nlj -> n'l'j') = (2j+1)(2j'+1) {l j 1/2; j' l' 1}^2 max(l, l') R(nl, n'l')^2,      (3)
```

with the radial integral R computed by Gauss-Laguerre quadrature. The pattern
is centred on its intensity-weighted mean and convolved with the normal Zeeman
triplet: sigma components at +-lambda^2 mu_B |B|/(h c) carry the fraction
(1 + cos^2 theta)/2 of the intensity, theta the angle between B and the sight
line. For N VII 9 -> 8 the fine structure adds a variance equivalent to 183 eV
and the Zeeman triplet at 2.5 T perpendicular to the sight line 29 eV. A
statistical upper-level population overestimates the fine-structure width when
capture favours low l.

## Channel moments and spectra

Sight points carry composite Gauss-Legendre rules. A channel returns

```text
B      = sum_i w_i eps_i,                         brightness [photons m^-2 s^-1 sr^-1]
v_app  = sum_i w_i eps_i <v_l>_i / B,
T_app  = m_Z c^2 var(lambda) / lambda_0^2,       var(lambda) = Doppler + fine structure + Zeeman + instrument,
```

with `temperature_doppler` the Doppler part alone, and the emission-weighted
flow `velocity_local`, temperature `temperature_local` and position, so that
cross-section and line-integration effects can be separated. `spectrum`
returns the spectral radiance on a wavelength grid as a sum of Gaussians, one
per sight point, fine-structure component and Zeeman component, convolved
with a Gaussian instrument function of the channel's FWHM. A single-Gaussian
fit of that spectrum approximates the moments.

The plasma is supplied by extending `cxrs_plasma_t` (`eval` for ne, Te, Ti,
Zeff and `impurity` for n_Z, T_Z, flow and field). Derivatives with respect to
the profiles are not yet provided.

## Verification

`tiago_cxrs_spectrum` checks the quadrature rules on monomials, the hydrogen
2p-1s radial integral 128 sqrt(6)/243 a0, two Wigner 6j values, the hydrogen
n = 2 fine structure against its (Z alpha)^4 term, the 2:1 Lyman-alpha
intensities, eq. (2) against closed-form Gaussian moments for q linear in E,
the brightness of a Gaussian beam crossed by a chord in closed form, apparent
velocity and temperature equal to V.e and T_Z for constant q, and the moments
of the synthetic spectrum against the channel moments with fine structure,
Zeeman effect, halo and instrument function.
