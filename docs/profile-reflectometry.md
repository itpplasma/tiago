# Profile reflectometry

`tiago_profile_reflectometry` predicts the phase and group delay of an
ordinary-mode (O-mode) reflectometer at normal incidence and inverts group
delays to cutoff positions.

## Forward model

The caller implements `sight_density_t`: the electron density \(n_e(s)\) in
\(\mathrm{m^{-3}}\), its slope and its parameter derivatives along the sight
line, with \(s\) in metres from the plasma edge, where \(n_e(0)=0\). With
\(X = n_e/n_c\), \(n_c = \varepsilon_0 m_e (2\pi f)^2/e^2\) and
\(N = \sqrt{1-X}\), the cutoff \(s_c\) solves \(n_e(s_c) = n_c\), and

\[
 \varphi = \frac{4\pi f}{c}\left(d + \int_0^{s_c} N\,ds\right) - \frac{\pi}{2},
 \qquad
 \tau = \frac{1}{2\pi}\frac{d\varphi}{df} = \frac{2}{c}\left(d + \int_0^{s_c}\frac{ds}{N}\right),
\]

with \(d\) the vacuum distance from antenna to edge. The phase integral is
Morales et al. (2017), Eq. (5), for the round trip. The \(-\pi/2\) is the
reflection phase of the WKB field of
[Gusakov and Surkov (2004)](https://arxiv.org/abs/physics/0410106), Sec. 2. The second form of
\(\tau\) uses \(\partial(fN)/\partial f = 1/N\). The substitution
\(s = s_c - u^2\) removes the \(1/\sqrt{\cdot}\) singularity, and
Gauss-Legendre nodes in \(u\) (128 by default) evaluate both integrals.

`omode_reflection` returns \(s_c\), \(\varphi\), \(\tau\) and their
derivatives in the profile parameters,

\[
 \partial_p s_c = -\frac{\partial_p n_e}{n_e'}\Big|_{s_c},\quad
 \partial_p\varphi = \frac{4\pi f}{c}\int_0^{s_c}\frac{\partial N}{\partial n_e}\,\partial_p n_e\,ds,
\]

and \(\partial_p\tau\) from the \(u\) form, whose integrand
\(\partial_{n_e}N^{-1}\,(\partial_p n_e + n_e'\,\partial_p s_c)\) stays finite at
the cutoff. The density must rise monotonically to the cutoff. A swept
(FMCW) reflectometer records the beat frequency
\(f_b = (\tau + \tau_\mathrm{sys})\,df/dt\), so the sweep law \(df/dt\) and the
system delay \(\tau_\mathrm{sys}\) of the waveguides are calibration inputs.

## Inversion

`invert_omode_group_delay` reconstructs cutoff positions from plasma group
delays at ascending frequencies by the recursion of Bottollier-Curtet as
written by [Morales et al. (2017)](https://doi.org/10.1063/1.4979513), Eqs.
(6) and (7): phase increments \(\pi(\tau_{n-1}+\tau_n)(f_n-f_{n-1})\),
trapezoids over the reconstructed steps with
\(N(f_n, s_i) = \sqrt{1-f_i^2/f_n^2}\) and the weight 2/3 in the last step. The
first cutoff \(s_1\) is an input, and the first phase assumes a linear layer
below it. Errors in \(s_1\) propagate to all later positions.

## Verification

`tiago_profile_reflectometry` checks

- a linear layer, where \(s_c\), \(\varphi\) and \(\tau\) have closed forms, to
  \(7\times10^{-11}\) relative (rounding of \(1-X\) next to the cutoff),
- \(\tau\) against the central difference of \(\varphi\) in frequency on a
  pedestal profile, to \(10^{-10}\) relative,
- parameter derivatives of \(s_c\), \(\varphi\) and \(\tau\) against central
  differences, to \(10^{-5}\) relative,
- the inversion of 180 synthetic group delays from 2 to 46.75 GHz on the
  pedestal profile, to 45 µm.

The model has no magnetic field (X-mode), oblique incidence, density
fluctuations or full-wave corrections.
