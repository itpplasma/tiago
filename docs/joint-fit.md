# Joint fit of measured profiles: `tiago_joint_fit`

`tiago_joint_fit input.nml` fits a fixed-boundary VMEC++ equilibrium and the
kinetic profiles n_e, T_e, T_i to measured point values and interferometer
line integrals read from files. It uses the [joint model](inference.md): the
equilibrium pressure is the kinetic pressure of the fitted profiles, every
point is located in the solved equilibrium, and each instrument group has its
own alignment shift and calibration. It requires a build with
`-DTIAGO_ENABLE_RECONSTRUCTION=ON`.

The posterior is the local Laplace approximation at the maximum a posteriori
(MAP) point. Profile bands and posterior draws come from that Gaussian. The
command does not run magnetic diagnostics; boundary, rotational transform or
current profile come from the VMEC input, and only the parameters named in
`names` are fitted.

```text
&joint_fit
  vmec_input = 'vmec.json'          ! stellarator-symmetric, kinetic-pressure compatible
  output_dir = 'fit'
  probes_file = 'probes.csv'        ! group,quantity,R,phi,Z,value,sigma
  chords_file = 'chords.csv'        ! x0,y0,z0,dx,dy,dz,length,value,sigma (optional)
  prior_table = 'prior.csv'         ! s,ne,Te,Ti,Zeff: prior means, Zeff fixed
  names = 'PHIEDGE'
  prior = 1.0
  prior_sigma = 0.02
  breakpoints = 0, 0.3, 0.6, 0.85, 0.95, 1.0
  z_imp = 6.0
  sigma_f = 1.0, 1.0, 1.0           ! prior sd of log ne, log Te, log Ti
  l_core = 0.3, l_edge = 0.05, s0 = 0.85, sw = 0.05
  align_sigma = 0.01, 0.0           ! per group: rigid (dR, dZ) shift prior [m]
  cal_sigma = 0.05, 0.0             ! per group: log calibration factor prior
  model_rel = 0.05, 0.05            ! per group: relative model-error sd
  student_nu = 4.0, 0.0             ! per group: Student likelihood, 0 = Gaussian
  chord_cal_sigma = 0.03
  chord_model_rel = 0.01
  laplace_draws = 400
  posterior_solves = 40             ! draws re-solved and written as wout
/
```

Quantities in `probes_file`: 1 n_e [m^-3], 2 T_e [eV], 3 T_i [eV]. R, Z in m,
phi in rad. Groups are numbered from one; each group is one
[point-probe](inference.md#profile-diagnostic-interface) diagnostic. Line
integrals are in m^-2 along straight chords in Cartesian coordinates.
The VMEC input must satisfy the kinetic-pressure conditions of
[reconstruction.md](reconstruction.md#kinetic-and-additional-pressure).

Outputs in `output_dir`:

| File | Content |
|---|---|
| `parameters.csv` | MAP and Laplace sd per latent, physical units for equilibrium parameters |
| `covariance.csv` | full Laplace covariance of all latents |
| `predictions.csv` | data, prediction (calibration applied) and sigma per signal |
| `profiles.csv` | MAP and 2.5, 16, 50, 84, 97.5 % quantiles of n_e, T_e, T_i, p on 101 values of s |
| `draws/latents.csv`, `draws/profiles.csv` | Laplace draws, their profiles |
| `draws/wout_NNNN.nc`, `wout_map.nc` | re-solved equilibria of the first `posterior_solves` draws and of the MAP |

`predict_only = .true.` writes `predictions.csv` at the prior mean without
fitting, for synthetic data.

The `tiago_joint_fit` test predicts point values and one chord on the public
VMEC++ Solov'ev tokamak from known profiles, fits them from a prior whose
mean is 40 % too high and requires the MAP profiles to recover the known
ones within three Laplace standard deviations or 3 %.
