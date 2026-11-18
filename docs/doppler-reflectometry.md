# Doppler reflectometry

A Doppler reflectometer measures the angular frequency shift of scattered
microwaves. The local operator takes the signed Bragg wavenumber \(k\) in
\(\mathrm{m^{-1}}\), its unit direction \(\hat{\boldsymbol k}\), a radial electric
field \(E_r\) in \(\mathrm{V\,m^{-1}}\), an outward unit normal \(\hat{\boldsymbol n}\)
and the magnetic field \(\boldsymbol B\) in tesla. The scattering layer and
wavenumber come from ray tracing. A nuisance velocity \(v_\mathrm{ph}\) describes
the propagation of the scattering fluctuations along \(\hat{\boldsymbol k}\).

\[
 B_\mathrm{proj} = \hat{\boldsymbol k}\cdot(\hat{\boldsymbol n}\times\boldsymbol B),
 \qquad
 \omega_D = k \left(E_r\,\frac{B_\mathrm{proj}}{|\boldsymbol B|^2}
                      + v_\mathrm{ph}\right).
\]

`tiago_doppler_reflectometry::doppler_reflectometry` returns \(\omega_D\) in rad/s
and its gradient in \(E_r, v_\mathrm{ph}, k, |\boldsymbol B|^2, B_\mathrm{proj}\).
Pass \(2\pi f_D\) for a measured shift \(f_D\) in Hz. The field squared must be
positive. The predicted shift has zero electric-field sensitivity when either
\(k\) or \(B_\mathrm{proj}\) is zero.

An electric-field fit needs a constraint on the phase velocity and the
scattering-layer position. To compare a predicted shift with an observation,
propagate the ray and magnetic covariance through the returned gradient,
including covariance between \(k\), \(B_\mathrm{proj}\) and \(|\boldsymbol B|^2\).
The caller supplies those quantities and the receiver frequency convention.
The local scattering approximation follows [Conway et al. (2006)](https://doi.org/10.1088/0029-5515/46/9/S15).
