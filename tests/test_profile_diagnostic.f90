program test_profile_diagnostic
    !! The point-probe reference diagnostic satisfies the contract: jvp agrees
    !! with central differences of predict, vjp is its adjoint, and the point
    !! motion of the alignment nuisance is reported through dx_deta.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_profile_diagnostic, only: plasma_sample_t, check_diagnostic_derivatives, &
        N_QUANTITIES, Q_NE, Q_TE, Q_TI
    use tiago_point_probe, only: point_probe_t
    implicit none
    type(point_probe_t) :: probe
    type(plasma_sample_t) :: sample
    real(dp), allocatable :: x(:, :), dx(:, :, :)
    real(dp) :: jvp_error, adjoint_error, eta(2)
    integer :: j, fails

    fails = 0
    probe%x0 = reshape([1.5_dp, 0.0_dp, 0.0_dp, 1.6_dp, 0.0_dp, 0.1_dp, 1.7_dp, 0.0_dp, 0.0_dp], [3, 3])
    probe%quantity = [Q_NE, Q_TE, Q_TI]
    probe%aligned = .false.
    eta = [0.01_dp, -0.02_dp]
    call probe%sample_points(eta, x, dx)
    if (abs(x(1, 2) - 1.61_dp) > 1.0e-14_dp .or. abs(x(3, 2) - 0.08_dp) > 1.0e-14_dp) then
        write(error_unit, *) 'FAIL: shifted points'
        fails = fails + 1
    end if
    if (any(dx(1, :, 1) /= 1.0_dp) .or. any(dx(3, :, 2) /= 1.0_dp) .or. any(dx(2, :, :) /= 0.0_dp)) then
        write(error_unit, *) 'FAIL: dx_deta'
        fails = fails + 1
    end if
    sample%n = 3
    sample%x = x
    sample%inside = [.true., .true., .false.]
    allocate(sample%q(N_QUANTITIES, 3), sample%dq_dx(3, N_QUANTITIES, 3))
    do j = 1, 3
        sample%q(:, j) = [0.2_dp * j, 2.5_dp, 3.0e19_dp / j, 2.0e3_dp / j, 1.5e3_dp / j, 1.5_dp]
    end do
    sample%q(3:5, 3) = 0.0_dp
    sample%dq_dx = 0.0_dp
    call check_diagnostic_derivatives(probe, sample, eta, 1.0e-6_dp, jvp_error, adjoint_error)
    if (jvp_error > 1.0e-8_dp .or. adjoint_error > 1.0e-12_dp) then
        write(error_unit, *) 'FAIL: derivatives', jvp_error, adjoint_error
        fails = fails + 1
    end if
    if (fails > 0) error stop 1
    print '(A)', 'profile diagnostic contract: point probe passes'
end program test_profile_diagnostic
