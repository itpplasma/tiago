program test_recorded_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_TE, Q_B, Q_S, &
                                        check_diagnostic_derivatives
    use tiago_point_probe, only: point_probe_t
    use tiago_recorded_diagnostic, only: calibrated_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_analytic_state, only: circular_state_t
    use tiago_real_space_state, only: located_t
    use tiago_state_sampler, only: sample_state, sample_state_vjp
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_TE, N_PROFILES
    implicit none
    type(point_probe_t) :: probe
    type(calibrated_diagnostic_t) :: diag
    type(calibration_prior_t) :: prior
    type(plasma_sample_t) :: sample
    real(dp) :: eta(4), y(2), err, adj, motion_error, analytic_error, zero_error
    real(dp), allocatable :: x(:, :), dx(:, :, :)
    integer :: status, fails
    fails = 0
    probe%x0 = reshape([2.0_dp, 0.0_dp, 0.0_dp, 2.1_dp, 0.0_dp, 0.1_dp], [3, 2])
    probe%quantity = [Q_TE, Q_B]
    probe%aligned = .false.
    call prior%init([2.0_dp, 0.3_dp], &
                    reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2]), status)
    if (status /= 0) error stop 'prior setup'
    call diag%init(probe, prior, status)
    if (status /= 0) error stop 'diagnostic setup'
    call check(diag%n_signals() == 2, 'signal count')
    call check(diag%n_nuisance() == 4, 'alignment and calibration in state')
    call check(diag%n_calibration() == 2, 'explicit calibration count')
    eta = [0.01_dp, -0.02_dp, 0.0_dp, 0.0_dp]
    call diag%sample_points(eta, x, dx)
    call check(maxval(abs(dx(:, :, 3:))) == 0.0_dp, 'calibration fixes geometry')
    call check(abs(x(1, 1) - 2.01_dp) < 1.0e-14_dp, 'delegated alignment')
    sample%n = 2
    sample%x = x
    sample%inside = [.true., .true.]
    allocate (sample%q(N_QUANTITIES, 2), sample%dq_dx(3, N_QUANTITIES, 2))
    sample%q = 0.0_dp
    sample%q(Q_TE, 1) = 100.0_dp
    sample%q(Q_B, 2) = 3.0_dp
    sample%dq_dx = 0.0_dp
    call diag%predict(sample, eta, y)
    analytic_error = maxval(abs(y - [200.3_dp, 6.3_dp]))
    call check(analytic_error < 1.0e-12_dp, &
               'analytic gain and offset')
    eta(3:) = [0.4_dp, -0.2_dp]
    call check_diagnostic_derivatives(diag, sample, eta, 1.0e-5_dp, err, adj)
    call check(err < 1.0e-8_dp, 'all sampled quantities and calibration FD')
    call check(adj < 1.0e-12_dp, 'complete adjoint identity')
    call check_motion()
    call check_provider()
    ! Vacuum/zero-input limit still records the offset.
    sample%q = 0.0_dp
    eta(3:) = 0.0_dp
    call diag%predict(sample, eta, y)
    zero_error = maxval(abs(y - 0.3_dp))
    call check(zero_error < 1.0e-14_dp, 'zero-input offset limit')
    print '(A,ES12.4)', 'analytic recorded gain/offset absolute error: ', analytic_error
    print '(A,ES12.4)', 'zero-input limit absolute error: ', zero_error
    print '(A,ES12.4)', 'recorded diagnostic FD relative error: ', err
    print '(A,ES12.4)', 'recorded diagnostic adjoint relative error: ', adj
    print '(A,ES12.4)', 'moving-point chain FD relative error: ', motion_error
    if (fails > 0) error stop 1
    print '(A)', 'recorded diagnostic and calibration latents: PASS'
contains
    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (ok) return
        write (error_unit, '(A)') 'FAIL: '//text
        fails = fails + 1
    end subroutine
    subroutine check_provider()
        type(circular_state_t) :: state
        type(kinetic_profiles_t) :: profiles
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: located(:)
        real(dp), allocatable :: z(:), zdot(:), zbar(:), f(:, :), df(:, :)
        real(dp), allocatable :: fbar(:, :), sbar(:), points(:, :), deltas(:, :, :)
        real(dp) :: p0(4), pdot(4), pbar(4), edot(4), ebar(4), h, yp(2), ym(2)
        real(dp) :: ybar(2), qbar(N_QUANTITIES, 2), lhs, rhs, error
        integer :: i
        call profiles%init([0.0_dp, 0.4_dp, 1.0_dp])
        call profiles%set_profile(P_TE, &
                                  [(200.0_dp*(1.0_dp - 0.7_dp*profiles%greville(i)), &
                  i=1, profiles%ncoef)], .true., 0.3_dp, 0.8_dp, 0.2_dp, 0.6_dp, 0.3_dp)
        allocate (z(profiles%n_latent), zdot(profiles%n_latent), &
                  zbar(profiles%n_latent), f(N_PROFILES, 2), df(N_PROFILES, 2), &
                  fbar(N_PROFILES, 2), sbar(2))
        z = 0.1_dp
        zdot = [(cos(real(i, dp)), i=1, size(z))]
        call state%get_params(p0)
        pdot = [0.2_dp, -0.1_dp, 0.05_dp, 0.3_dp]
        edot = [0.01_dp, -0.02_dp, 0.3_dp, -0.4_dp]
        call diag%sample_points(eta, points, deltas)
        call sample_state(state, points, smp, located)
        call profiles%evaluate(z, smp%q(Q_S, :), f, df)
        smp%q(Q_TE, :) = f(P_TE, :)
        smp%dq_dx(:, Q_TE, :) = smp%dq_dx(:, Q_S, :)* &
                                spread(df(P_TE, :), 1, 3)
        ybar = [0.6_dp, -0.2_dp]
        call diag%vjp(smp, eta, ybar, qbar, ebar)
        do i = 1, 4
            ebar(i) = ebar(i) + sum(qbar* &
                                    matmul_gradient_at(smp, deltas, i))
        end do
        fbar = 0.0_dp
        fbar(P_TE, :) = qbar(Q_TE, :)
        zbar = 0.0_dp
        call profiles%vjp(z, smp%q(Q_S, :), fbar, zbar, sbar)
        qbar(Q_S, :) = qbar(Q_S, :) + sbar
        pbar = 0.0_dp
        call sample_state_vjp(state, located, qbar, pbar)
        h = 1.0e-6_dp
        call state%set_params(p0 + h*pdot)
        call provider_prediction(state, profiles, z + h*zdot, eta + h*edot, yp)
        call state%set_params(p0 - h*pdot)
        call provider_prediction(state, profiles, z - h*zdot, eta - h*edot, ym)
        lhs = dot_product(ybar, (yp - ym)/(2*h))
        rhs = dot_product(pbar, pdot) + dot_product(zbar, zdot) + &
              dot_product(ebar, edot)
        error = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
        call check(error < 1.0e-7_dp, 'provider/profile/calibration VJP vs FD')
       print '(A,ES12.4)', 'provider/profile/calibration VJP relative FD error: ', error
    end subroutine

    function matmul_gradient_at(smp, deltas, k) result(g)
        type(plasma_sample_t), intent(in) :: smp
        real(dp), intent(in) :: deltas(:, :, :)
        integer, intent(in) :: k
        real(dp) :: g(N_QUANTITIES, 2)
        integer :: i
        do i = 1, 2
            g(:, i) = matmul(transpose(smp%dq_dx(:, :, i)), deltas(:, i, k))
        end do
    end function

    subroutine provider_prediction(state, profiles, z, eta_at, prediction)
        type(circular_state_t), intent(in) :: state
        type(kinetic_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: z(:), eta_at(:)
        real(dp), intent(out) :: prediction(:)
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: located(:)
        real(dp), allocatable :: points(:, :), deltas(:, :, :)
        real(dp) :: f(N_PROFILES, 2), df(N_PROFILES, 2)
        call diag%sample_points(eta_at, points, deltas)
        call sample_state(state, points, smp, located)
        call profiles%evaluate(z, smp%q(Q_S, :), f, df)
        smp%q(Q_TE, :) = f(P_TE, :)
        call diag%predict(smp, eta_at, prediction)
    end subroutine

    subroutine check_motion()
        type(plasma_sample_t) :: plus, minus
        real(dp) :: direction(4), qdot(N_QUANTITIES, 2), yd(2), yp(2), ym(2)
        real(dp) :: qb(N_QUANTITIES, 2), eb(4), yb(2), h, lhs, rhs, gradient
        real(dp), allocatable :: xp(:, :), xm(:, :), delta(:, :, :)
        integer :: i
        h = 1.0e-5_dp
        direction = [0.3_dp, -0.4_dp, 0.7_dp, -0.2_dp]
        sample%dq_dx(1, Q_TE, 1) = 4.0_dp
        sample%dq_dx(3, Q_B, 2) = -0.5_dp
        call diag%sample_points(eta + h*direction, xp, delta)
        call diag%sample_points(eta - h*direction, xm, delta)
        plus = sample
        minus = sample
        qdot = 0.0_dp
        do i = 1, 2
            plus%q(:, i) = sample%q(:, i) + &
                           matmul(transpose(sample%dq_dx(:, :, i)), xp(:, i) - x(:, i))
            minus%q(:, i) = sample%q(:, i) + &
                            matmul(transpose(sample%dq_dx(:, :, i)), xm(:, i) - x(:, i))
            qdot(:, i) = matmul(transpose(sample%dq_dx(:, :, i)), &
                                matmul(dx(:, i, :), direction))
        end do
        call diag%predict(plus, eta + h*direction, yp)
        call diag%predict(minus, eta - h*direction, ym)
        call diag%jvp(sample, eta, qdot, direction, yd)
        motion_error = maxval(abs((yp - ym)/(2*h) - yd))/maxval(abs(yd))
        call check(motion_error < 1.0e-8_dp, 'equilibrium/point-motion chain FD')
        yb = [0.7_dp, -0.2_dp]
        call diag%vjp(sample, eta, yb, qb, eb)
        do i = 1, 4
            gradient = sum(qb*matmul_gradient(i))
            eb(i) = eb(i) + gradient
        end do
        lhs = dot_product(yb, yd)
        rhs = dot_product(eb, direction)
        call check(abs(lhs - rhs) < 1.0e-12_dp, 'moving-point chain transpose')
    end subroutine
    function matmul_gradient(k) result(g)
        integer, intent(in) :: k
        real(dp) :: g(N_QUANTITIES, 2)
        integer :: i
        do i = 1, 2
            g(:, i) = matmul(transpose(sample%dq_dx(:, :, i)), dx(:, i, k))
        end do
    end function
end program
