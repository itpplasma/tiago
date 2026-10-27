program test_dispersion_interferometer
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use tiago_dispersion_interferometer, only: dispersion_interferometer_t
   use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_NE, &
                                       check_diagnostic_derivatives
   use tiago_sight_geometry, only: chord_t
   use tiago_calibration, only: calibration_prior_t
   use tiago_dispersion_phase, only: dispersion_phase
   use tiago_dispersion_column, only: dispersion_column
   use tiago_dispersion_iq, only: dispersion_iq
   use tiago_dispersion_atan, only: dispersion_atan
   use tiago_iq_phase_reduction, only: iq_phase_reduction_t
   use tiago_signal_source, only: signal_window_t
   implicit none
   real(dp), parameter :: pi = acos(-1.0_dp), n0 = 5.0e19_dp
   type(calibration_prior_t) :: prior
   real(dp) :: mean(8), covariance(8, 8), sd(8)
   integer :: i, status
   mean = [1.0e-5_dp, 1.1_dp, 0.9_dp, 0.17_dp, -0.23_dp, 0.3_dp, &
           1.0_dp, 0.2_dp]
   sd = [1.0e-8_dp, 0.03_dp, 0.03_dp, 0.1_dp, 0.1_dp, 0.1_dp, &
         0.1_dp, 0.1_dp]
   covariance = 0.0_dp
   do i = 1, 8
      covariance(i, i) = sd(i)**2
   end do
   covariance(4, 5) = 0.002_dp
   covariance(5, 4) = covariance(4, 5)
   call prior%init(mean, covariance, status)
   call check(status == 0, 'correlated calibration prior')
   call integrals()
   call phase_limits()
   call inverse_derivatives()
   call roundtrip()
   call check_invalid()
   print '(A)', 'dispersion interferometer passed'
contains
   subroutine check(ok, text)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: text
      if (.not. ok) then
         print '(A)', 'FAIL '//text
         error stop 1
      end if
   end subroutine

   subroutine integrals()
      type(dispersion_interferometer_t) :: model
      type(chord_t) :: chord(1)
      type(plasma_sample_t) :: smp
      real(dp), allocatable :: x(:, :), dx(:, :, :)
      real(dp) :: b, alpha, h, exact, columns(1), err(1), maxerr, jerr, aerr
      integer :: a, j
      real(dp), parameter :: alphas(4) = [0.5_dp, 1.0_dp, 1.3_dp, 2.0_dp]
      maxerr = 0.0_dp
      do a = 1, 4
         alpha = alphas(a)
         do j = 0, 8
            b = 0.1_dp*j
            h = sqrt(1.0_dp - b*b)
            chord(1) = chord_t('synthetic', [3.0_dp + b, 0.0_dp, -h], &
                               [0.0_dp, 0.0_dp, 1.0_dp], 2*h, 2048, 4)
            call model%init(chord, prior, [0.0_dp], [0.0_dp], status)
            call check(status == 0, 'cylinder init')
            call model%sample_points(spread(0.0_dp, 1, 8), x, dx)
            smp%n = size(x, 2)
            if (allocated(smp%q)) deallocate (smp%q)
            allocate (smp%q(N_QUANTITIES, smp%n))
            smp%q = 0.0_dp
            smp%q(Q_NE, :) = n0*max(0.0_dp, 1 - b*b - x(3, :)**2)**alpha
            ! Beta integral, DLMF 5.12.1, t = z/sqrt(1-b^2).
            exact = n0*(1 - b*b)**(alpha + 0.5_dp)*sqrt(pi) &
                    *gamma(alpha + 1)/gamma(alpha + 1.5_dp)
            call model%quadrature_check(smp, 0.0_dp, 1.0e-7_dp, &
                                        columns, err, status)
            print '(A,5ES24.16)', 'cylinder,', alpha, b, columns(1), &
               exact, abs(columns(1)/exact - 1)
            call check(status == 0, 'cylinder convergence indicator')
            maxerr = max(maxerr, abs(columns(1)/exact - 1))
            call check(abs(columns(1)/exact - 1) < 1.0e-7_dp, &
                       'cylinder beta integral')
         end do
      end do
      print '(A,ES12.4)', 'cylinder max relative error: ', maxerr
      smp%q(Q_NE, :) = n0
      call model%quadrature_check(smp, 0.0_dp, 1.0e-12_dp, columns, err, status)
      call check(abs(columns(1)/(n0*2*h) - 1) < 1.0e-12_dp, 'slab integral')
      print '(A,ES12.4)', 'slab relative error: ', abs(columns(1)/(n0*2*h) - 1)
      call check_diagnostic_derivatives(model, smp, spread(0.1_dp, 1, 8), &
                                        1.0e-4_dp, jerr, aerr)
      print '(A,2ES12.4)', 'I/Q JVP FD and adjoint errors: ', jerr, aerr
      call check(jerr < 1.0e-7_dp .and. aerr < 1.0e-12_dp, 'I/Q derivatives')
      ! A coarse rule must report an unresolved non-polynomial profile.
      chord(1)%npanel = 1
      call model%init(chord, prior, [0.0_dp], [0.0_dp], status)
      call model%sample_points(spread(0.0_dp, 1, 8), x, dx)
      smp%q = spread(n0*sqrt(max(0.0_dp, 1 - b*b - x(3, :)**2)), &
                     1, N_QUANTITIES)
      call model%quadrature_check(smp, 0.0_dp, 1.0e-10_dp, columns, err, status)
      call check(status /= 0, 'unresolved quadrature rejected')
   end subroutine

   subroutine phase_limits()
      real(dp) :: phase(1), jp(1, 8), other(1), scratch(1, 8), exact, ratio
      real(dp) :: arg(8), plus(1), minus(1), h, maxerr, cold, xplasma
      integer :: k
      arg = [1.0e20_dp, 1.0e-5_dp, 1.0e-8_dp, 0.9_dp, &
             0.4_dp, 0.2_dp, 0.3_dp, 0.1_dp]
      call phase_at(arg, phase, jp)
      maxerr = 0.0_dp
      do k = 1, 8
         h = 1.0e-5_dp*max(abs(arg(k)), 1.0e-20_dp)
         arg(k) = arg(k) + h
         call phase_at(arg, plus, scratch)
         arg(k) = arg(k) - 2*h
         call phase_at(arg, minus, scratch)
         arg(k) = arg(k) + h
         ratio = abs((plus(1) - minus(1))/(2*h) - jp(1, k)) &
                 /max(1.0_dp, abs(jp(1, k)))
         maxerr = max(maxerr, ratio)
      end do
      call check(maxerr < 1.0e-7_dp, 'phase partials FD')
      print '(A,ES12.4)', 'phase partials scaled FD error: ', maxerr
      arg(8) = 1.0e6_dp
      call phase_at(arg, other, scratch)
      call check(other(1) == phase(1), 'common-path vibration cancellation')
      call check(jp(1, 8) == 0, 'common-path zero derivative')
      arg(1) = 0
      arg(3) = 0
      arg(5) = 0
      arg(7) = 0
      call phase_at(arg, phase, jp)
      call check(phase(1) == 0, 'vacuum phase')
      ! Bagryansky et al., RSI 77, 053501 (2006), Eq. (3),
      ! DOI 10.1063/1.2202922: positive electron phase convention.
      arg(1) = 1.0e20_dp
      call phase_at(arg, phase, jp)
      exact = 1.5_dp*2.8179403205e-15_dp*arg(2)*arg(1)
      call check(abs(phase(1)/exact - 1) < 1.0e-14_dp, 'published phase coefficient')
      ! Independent cold-plasma index difference, length = 1 m.
      xplasma = arg(1)*2.8179403205e-15_dp*arg(2)**2/pi
      cold = 4*pi/arg(2)*(sqrt(1 - xplasma/4) - sqrt(1 - xplasma))
      ratio = abs(phase(1)/cold - 1)
      print '(A,ES12.4)', 'dilute cold-plasma index relative difference: ', ratio
      call check(ratio < 1.0e-5_dp, 'dilute cold-plasma limit')
      arg(1) = 0
      arg(3) = 2.0e-8_dp
      call phase_at(arg, phase, jp)
      exact = -4*pi*arg(3)*arg(4)/arg(2)
      call check(abs(phase(1)/exact - 1) < 1.0e-14_dp, 'neutral dispersion')
   end subroutine

   subroutine phase_at(a, value, jac)
      real(dp), intent(in) :: a(8)
      real(dp), intent(out) :: value(1), jac(1, 8)
      call dispersion_phase(a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), &
                            value, jac)
   end subroutine

   subroutine inverse_derivatives()
      real(dp) :: args(7), saved(7), value(1), jac(1, 7), scratch(1, 7)
      real(dp) :: ajac(1, 6), ascratch(1, 6), plus(1), minus(1), h, error
      integer :: k
      args = [4.0_dp, 1.0e-5_dp, 2.0e-8_dp, 0.9_dp, 0.4_dp, 0.2_dp, 0.3_dp]
      saved = args
      call column_at(args, value, jac)
      error = 0
      do k = 1, 7
         h = 1.0e-5_dp*abs(saved(k))
         args = saved
         args(k) = saved(k) + h
         call column_at(args, plus, scratch)
         args(k) = saved(k) - h
         call column_at(args, minus, scratch)
         error = max(error, abs((plus(1) - minus(1))/(2*h) - jac(1, k)) &
                     /max(1.0_dp, abs(jac(1, k))))
      end do
      print '(A,ES12.4)', 'column inverse partial FD relative error: ', error
      call check(error < 1.0e-7_dp, 'column inverse derivatives')
      args(:6) = [0.7_dp, 0.5_dp, 1.1_dp, 0.9_dp, 0.17_dp, -0.23_dp]
      saved = args
      call atan_at(args(:6), value, ajac)
      error = 0
      do k = 1, 6
         h = 1.0e-5_dp*abs(saved(k))
         args = saved
         args(k) = saved(k) + h
         call atan_at(args(:6), plus, ascratch)
         args(k) = saved(k) - h
         call atan_at(args(:6), minus, ascratch)
         error = max(error, abs((plus(1) - minus(1))/(2*h) - ajac(1, k)))
      end do
      print '(A,ES12.4)', 'atan inverse partial FD absolute error: ', error
      call check(error < 1.0e-8_dp, 'atan inverse calibration derivatives')
   end subroutine

   subroutine column_at(a, value, jac)
      real(dp), intent(in) :: a(7)
      real(dp), intent(out) :: value(1), jac(1, 7)
      call dispersion_column(a(1), a(2), a(3), a(4), a(5), a(6), a(7), value, jac)
   end subroutine

   subroutine atan_at(a, value, jac)
      real(dp), intent(in) :: a(6)
      real(dp), intent(out) :: value(1), jac(1, 6)
      call dispersion_atan(a(1), a(2), a(3), a(4), a(5), a(6), value, jac)
   end subroutine

   subroutine make_window(w, n)
      type(signal_window_t), intent(out) :: w
      integer, intent(in) :: n
      integer :: k
      allocate (w%time(2*n), w%samples(2*n), w%variance(2*n))
      w%time = [(real((k - 1)/2, dp), k=1, 2*n)]
      w%samples = 0
      w%variance = 1.0e-6_dp
      w%units = 'V'
      w%provenance%source_id = 'synthetic'
      w%provenance%version = '1'
      w%provenance%request_key = 'synthetic-iq'
      w%provenance%record_id = 'synthetic'
      w%provenance%metadata = 'interleaved I,Q'
      call w%seal()
   end subroutine

   subroutine roundtrip()
      integer, parameter :: n = 201
      type(signal_window_t) :: raw, reduced, perturbed, plus, minus
      type(iq_phase_reduction_t) :: step
      real(dp) :: phase(n), value(1), jp(1, 8), ji(2, 5), iq(2), col(1), jc(1, 7)
      real(dp) :: v(2*n), vb(2*n), w(n), jv(n), h, maxerr, cov(n, n), row(2*n)
      integer :: k, j
      call make_window(raw, n)
      step%id = 'synthetic-phase-v1'
      step%amplitude = mean(2:3)
      step%offset = mean(4:5)
      allocate (raw%covariance(2*n, 2*n))
      raw%covariance = 2.0e-7_dp
      do k = 1, n
         call dispersion_phase(1.0e20_dp*(k - 1)/10, mean(1), 1.0e-8_dp, &
                               mean(7), 0.2_dp, mean(8), mean(6), 0.0_dp, value, jp)
         phase(k) = value(1)
         call dispersion_iq(phase(k), mean(2), mean(3), mean(4), mean(5), iq, ji)
         raw%samples(2*k - 1:2*k) = iq
      end do
      call raw%seal()
      step%anchor = phase(1)
      call step%apply(raw, reduced, status)
      call check(status == 0, 'phase reduction apply')
      maxerr = maxval(abs(reduced%samples - phase))
      print '(A,ES12.4)', 'multi-fringe phase roundtrip absolute error [rad]: ', maxerr
      call check(maxerr < 1.0e-12_dp, 'unwrapped phase roundtrip')
      do k = 1, n, 10
         print '(A,3ES24.16)', 'phase,', real(k - 1, dp), phase(k), &
            reduced%samples(k)
      end do
      maxerr = 0
      do k = 1, n
         call dispersion_column(reduced%samples(k), mean(1), 1.0e-8_dp, &
                                mean(7), 0.2_dp, mean(8), mean(6), col, jc)
         maxerr = max(maxerr, abs(col(1) - 1.0e20_dp*(k - 1)/10)/1.0e20_dp)
      end do
      print '(A,ES12.4)', 'column roundtrip error / 1e20 m^-2: ', maxerr
      call check(maxerr < 1.0e-12_dp, 'column roundtrip')
      v = [(sin(0.7_dp*k), k=1, 2*n)]
      w = [(cos(0.3_dp*k), k=1, n)]
      call step%jvp(raw, v, jv)
      call step%vjp(raw, w, vb)
      call check(abs(dot_product(w, jv) - dot_product(vb, v)) < 1.0e-12_dp, &
                 'phase reduction adjoint')
      h = 1.0e-6_dp
      perturbed = raw
      perturbed%samples = raw%samples + h*v
      call perturbed%seal()
      call step%apply(perturbed, plus, status)
      call check(status == 0, 'phase plus perturbation')
      perturbed%samples = raw%samples - h*v
      call perturbed%seal()
      call step%apply(perturbed, minus, status)
      call check(status == 0, 'phase minus perturbation')
      maxerr = maxval(abs(jv - (plus%samples - minus%samples)/(2*h)))
      print '(A,ES12.4)', 'phase reduction FD absolute error: ', maxerr
      call check(maxerr < 2.0e-8_dp, 'phase reduction FD')
      do j = 1, n
         w = 0
         w(j) = 1
         call step%vjp(raw, w, row)
         vb = matmul(raw%total_covariance(), row)
         call step%jvp(raw, vb, cov(:, j))
      end do
      maxerr = maxval(abs(reduced%total_covariance() - cov))
      print '(A,ES12.4)', 'phase covariance J C J^T absolute error [rad^2]: ', maxerr
      call check(maxerr < 1.0e-18_dp, 'phase covariance')
      ! A streaming caller advances the anchor between bounded chunks.
      perturbed = raw
      perturbed%samples = raw%samples(201:)
      perturbed%time = raw%time(201:)
      perturbed%variance = raw%variance(201:)
      perturbed%covariance = raw%covariance(201:, 201:)
      call perturbed%seal()
      step%anchor = phase(100)
      call step%apply(perturbed, plus, status)
      call check(status == 0, 'phase streaming chunk')
      call check(maxval(abs(plus%samples - phase(101:))) < 1.0e-12_dp, &
                 'phase chunk continuity')
   end subroutine

   subroutine check_invalid()
      type(signal_window_t) :: raw, output
      type(iq_phase_reduction_t) :: step
      type(dispersion_interferometer_t) :: model
      type(chord_t) :: chord(1)
      call make_window(raw, 2)
      step%id = 'synthetic'
      call step%apply(raw, output, status)
      call check(status /= 0, 'zero phasor rejected')
      raw%samples = 1
      raw%time(2) = 0.5_dp
      call raw%seal()
      call step%apply(raw, output, status)
      call check(status /= 0, 'unpaired times rejected')
      chord(1)%label = 'invalid'
      chord(1)%direction = 0
      call model%init(chord, prior, [0.0_dp], [0.0_dp], status)
      call check(status /= 0, 'zero chord direction rejected')
   end subroutine
end program
