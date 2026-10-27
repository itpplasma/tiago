program test_interferometer_profiles
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use tiago_dispersion_interferometer, only: dispersion_interferometer_t
   use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_NE, Q_S
   use tiago_sight_geometry, only: chord_t
   use tiago_calibration, only: calibration_prior_t
   use tiago_analytic_state, only: circular_state_t
   use tiago_real_space_state, only: located_t
   use tiago_state_sampler, only: sample_state, sample_state_vjp
   use tiago_kinetic_profiles, only: kinetic_profiles_t, N_PROFILES, P_NE
   implicit none
   type(dispersion_interferometer_t) :: model
   type(calibration_prior_t) :: prior
   type(circular_state_t) :: state
   type(kinetic_profiles_t) :: profiles
   type(plasma_sample_t) :: sample
   type(located_t), allocatable :: located(:)
   real(dp), allocatable :: x(:, :), dx(:, :, :), z(:), zdot(:), zbar(:)
   real(dp), allocatable :: f(:, :), df(:, :), fbar(:, :), sbar(:), qbar(:, :)
   real(dp) :: mean(8), cov(8, 8), eta(8), edot(8), ebar(8), ybar(2), yp(2), ym(2)
   real(dp) :: p0(4), pdot(4), pbar(4), rhs, lhs, h, error
   type(chord_t) :: chord(1)
   integer :: i, status, k
   mean = [1.0e-5_dp, 1.1_dp, 0.9_dp, 0.17_dp, -0.23_dp, 0.3_dp, 1.0_dp, 0.2_dp]
   cov = 0
   do i = 1, 8
      cov(i, i) = (0.01_dp*mean(i))**2
   end do
   call prior%init(mean, cov, status)
   if (status /= 0) error stop 'profile chain prior'
   chord(1) = chord_t('synthetic', [3.3_dp, 0.0_dp, -0.5_dp], &
                      [0.0_dp, 0.0_dp, 1.0_dp], 1.0_dp, 8, 4)
   call model%init(chord, prior, [2.0e-8_dp], [0.4_dp], status)
   if (status /= 0) error stop 'profile chain init'
   call profiles%init([0.0_dp, 0.4_dp, 1.0_dp])
   call profiles%set_profile(P_NE, &
                 [(5.0e19_dp*(1 - 0.7_dp*profiles%greville(i)), i=1, profiles%ncoef)], &
                             .true., 0.3_dp, 0.8_dp, 0.2_dp, 0.6_dp, 0.3_dp)
   allocate (z(profiles%n_latent), zdot(profiles%n_latent), zbar(profiles%n_latent))
   z = 0.1_dp
   zdot = [(cos(real(i, dp)), i=1, size(z))]
   eta = 0.1_dp
   edot = [(sin(real(i, dp)), i=1, 8)]
   call state%get_params(p0)
   pdot = [0.2_dp, -0.1_dp, 0.05_dp, 0.3_dp]
   call model%sample_points(eta, x, dx)
   call sample_state(state, x, sample, located)
   allocate (f(N_PROFILES, sample%n), df(N_PROFILES, sample%n), &
             fbar(N_PROFILES, sample%n), sbar(sample%n), qbar(N_QUANTITIES, sample%n))
   call profiles%evaluate(z, sample%q(Q_S, :), f, df)
   sample%q(Q_NE, :) = f(P_NE, :)
   ybar = [0.6_dp, -0.2_dp]
   call model%vjp(sample, eta, ybar, qbar, ebar)
   fbar = 0
   fbar(P_NE, :) = qbar(Q_NE, :)
   zbar = 0
   call profiles%vjp(z, sample%q(Q_S, :), fbar, zbar, sbar)
   qbar(Q_S, :) = qbar(Q_S, :) + sbar
   pbar = 0
   call sample_state_vjp(state, located, qbar, pbar)
   rhs = dot_product(pbar, pdot) + dot_product(zbar, zdot) + dot_product(ebar, edot)
   do k = 4, 6
      h = 10.0_dp**(-k)
      call state%set_params(p0 + h*pdot)
      call predict(z + h*zdot, eta + h*edot, yp)
      call state%set_params(p0 - h*pdot)
      call predict(z - h*zdot, eta - h*edot, ym)
      lhs = dot_product(ybar, (yp - ym)/(2*h))
      error = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
      print '(A,2ES14.6)', 'GP/profile/state/calibration VJP FD step,error: ', h, error
      if (error > 1.0e-6_dp) error stop 'interferometer profile chain'
   end do
contains
   subroutine predict(latent, calibration, y)
      real(dp), intent(in) :: latent(:), calibration(:)
      real(dp), intent(out) :: y(:)
      type(plasma_sample_t) :: smp
      type(located_t), allocatable :: loc(:)
      real(dp), allocatable :: values(:, :), ds(:, :)
      call sample_state(state, x, smp, loc)
      allocate (values(N_PROFILES, smp%n), ds(N_PROFILES, smp%n))
      call profiles%evaluate(latent, smp%q(Q_S, :), values, ds)
      smp%q(Q_NE, :) = values(P_NE, :)
      call model%predict(smp, calibration, y)
   end subroutine
end program
