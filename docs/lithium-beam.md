# Lithium-beam emission

`tiago_lithium_beam` predicts the line emission of a neutral alkali beam,
such as the Li I 2s-2p line of lithium-beam emission spectroscopy, from the
electron density and temperature along the beam.

## Model

The collisional-radiative model of [Brix et al. (2001)](https://scipub.euro-fusion.org/wp-content/uploads/2014/11/EFDC010248.pdf)
(who cite Schweinzer et al., Plasma Phys. Control. Fusion 34, 1173 (1992),
for it) evolves the level populations \(N_i\) per
injected atom along the beam coordinate \(z\),

\[
 v\,\frac{dN_i}{dz} = \sum_j M_{ij} N_j ,
 \qquad
 M = n_e\,E(T_e) + n_i\,I + A ,
\]

with beam speed \(v\), electron and ion rate coefficients
\(\langle\sigma v\rangle\) in \(E\) and \(I\) (gains off the diagonal,
excitation, de-excitation and ionisation losses on it) and spontaneous rates
\(A\). The ion density is \(n_i = r\,n_e\). The emission of the line
\(u\to l\) per metre and per injected atom is

\[
 \epsilon(z) = A_{lu} N_u(z)/v .
\]

A calibrated signal is \(\epsilon\) times the beam current over \(e\) and the
optical transmission of the channel.

`beam_atomic_data_t` holds the rate tables: electron rates on ascending
temperatures, interpolated linearly in \(\ln T_e\), ion rates for the beam
energy and the Einstein coefficients. TIAGO ships no atomic data. The caller
supplies the tables, for example from Schweinzer et al., At. Data Nucl. Data
Tables 72, 239 (1999), the database cited by Brix et al. (2001).

`lithium_beam_emission` integrates on the caller's nodes \(z_k\) with
Crank-Nicolson steps whose rate matrix uses the node-averaged \(n_e\) and
\(T_e\) of each cell. It returns the populations, \(\epsilon(z_k)\) and the
exact derivatives of the discrete \(\epsilon(z_k)\) in \(n_e(z_j)\) and
\(T_e(z_j)\), from the tangent recursion

\[
 (1-hM_k)\,\dot N_{k+1} = (1+hM_k)\,\dot N_k + h\,\dot M_k\,(N_k+N_{k+1}),
 \qquad h = \frac{z_{k+1}-z_k}{2v}.
\]

## Verification

`tiago_lithium_beam` checks

- a two-level beam in a uniform plasma against the closed-form \(2\times2\)
  matrix exponential: relative error \(3.6\times10^{-6}\) at 100 steps and
  \(9.1\times10^{-7}\) at 200, the second-order ratio 4,
- conservation of the total population to \(10^{-13}\) without ionisation,
- density and temperature derivatives on pedestal profiles against central
  differences, to \(1.5\times10^{-7}\) relative.

The rates in the tests are synthetic. The model has no fine-structure,
charge-exchange to impurities or beam divergence, and no observation optics.
