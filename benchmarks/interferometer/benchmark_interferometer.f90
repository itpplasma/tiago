program benchmark_interferometer
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   use omp_lib, only: omp_get_max_threads
   use tiago_dispersion_interferometer, only: dispersion_interferometer_t
   use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_NE
   use tiago_sight_geometry, only: chord_t
   use tiago_calibration, only: calibration_prior_t
   implicit none
   integer, parameter :: nc = 1024, reps = 100
   type(dispersion_interferometer_t) :: model
   type(calibration_prior_t) :: prior
   type(chord_t) :: chords(nc)
   type(plasma_sample_t) :: sample
   real(dp), allocatable :: x(:, :), dx(:, :, :), qdot(:, :), qbar(:, :)
   real(dp) :: mean(8), cov(8, 8), y(2*nc), dy(2*nc), ybar(2*nc), ebar(8)
   real(dp) :: eta(8), edot(8), elapsed, checksum
   integer(int64) :: start, finish, rate
   integer :: i, k, status, mode
   mean = [1.0e-5_dp, 1.1_dp, 0.9_dp, 0.17_dp, -0.23_dp, 0.3_dp, 1.0_dp, 0.2_dp]
   cov = 0
   do i = 1, 8
      cov(i, i) = (0.01_dp*mean(i))**2
   end do
   call prior%init(mean, cov, status)
   if (status /= 0) error stop 'benchmark prior'
   do i = 1, nc
      chords(i) = chord_t('synthetic', [3.2_dp, 0.0_dp, -1.0_dp], &
                          [0.0_dp, 0.0_dp, 1.0_dp], 2.0_dp, 16, 4)
   end do
   call model%init(chords, prior, spread(0.0_dp, 1, nc), &
                   spread(0.0_dp, 1, nc), status)
   if (status /= 0) error stop 'benchmark model'
   eta = 0
   edot = 0.1_dp
   call model%sample_points(eta, x, dx)
   sample%n = size(x, 2)
   allocate (sample%q(N_QUANTITIES, sample%n), qdot(N_QUANTITIES, sample%n), &
             qbar(N_QUANTITIES, sample%n))
   sample%q = 0
   sample%q(Q_NE, :) = 5.0e19_dp*exp(-x(3, :)**2)
   qdot = 0
   qdot(Q_NE, :) = 0.1_dp*sample%q(Q_NE, :)
   ybar = [(cos(0.3_dp*i), i=1, 2*nc)]
   call system_clock(count_rate=rate)
   print '(A,I0,A,I0,A,I0)', 'threads=', omp_get_max_threads(), &
      ',chords=', nc, ',quadrature_points=', sample%n
   do mode = 1, 3
      checksum = 0
      do k = 0, 3
         call system_clock(start)
         do i = 1, reps
            select case (mode)
            case (1)
               call model%predict(sample, eta, y)
            case (2)
               call model%jvp(sample, eta, qdot, edot, dy)
            case (3)
               call model%vjp(sample, eta, ybar, qbar, ebar)
            end select
         end do
         call system_clock(finish)
         elapsed = real(finish - start, dp)/rate/reps
         select case (mode)
         case (1)
            checksum = sum(y)
         case (2)
            checksum = sum(dy)
         case (3)
            checksum = sum(qbar) + sum(ebar)
         end select
         if (k > 0) print '(I0,A,I0,A,ES14.6,A,ES24.16)', &
            mode, ',', k, ',', elapsed, ',', checksum
      end do
   end do
end program
