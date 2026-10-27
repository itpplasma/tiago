program benchmark_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_TE
    use tiago_point_probe, only: point_probe_t
    use tiago_recorded_diagnostic, only: calibrated_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_signal_source, only: signal_window_t
    use tiago_signal_reduction, only: linear_reduction_t
    implicit none
    integer, parameter :: n = 1000000, repeats = 30
    type(point_probe_t) :: probe
    type(calibrated_diagnostic_t) :: diag
    type(calibration_prior_t) :: prior
    type(signal_window_t) :: raw, reduced
    type(linear_reduction_t) :: reduction
    type(plasma_sample_t) :: sample
    real(dp), allocatable :: y(:), qdot(:, :), qbar(:, :), ybar(:)
    real(dp) :: eta(2), etadot(2), etabar(2), times(4), checksum
    real(dp) :: c(2), dc(2), expected(2), relative_error
    integer(int64) :: start, finish, rate
    integer :: status, i, mode
    allocate (probe%quantity(n), probe%x0(3, n))
    probe%quantity = Q_TE
    probe%x0 = 0.0_dp
    call prior%init([2.0_dp, 0.3_dp], &
                    reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2]), status)
    if (status /= 0) error stop 'prior'
    call diag%init(probe, prior, status)
    if (status /= 0) error stop 'diagnostic'
    sample%n = n
    allocate (sample%q(N_QUANTITIES, n), y(n), ybar(n), &
              qdot(N_QUANTITIES, n), qbar(N_QUANTITIES, n))
    sample%q = 3.0_dp
    qdot = 0.2_dp
    ybar = 0.7_dp
    eta = [0.1_dp, -0.2_dp]
    etadot = [0.3_dp, -0.1_dp]
    allocate(raw%time(n))
    do i = 1, n
        raw%time(i) = real(i, dp)
    end do
    allocate (raw%samples(n), raw%variance(n))
    raw%samples = 3.0_dp
    raw%variance = 0.1_dp
    raw%units = 'V'
    raw%provenance%source_id = 'synthetic-benchmark'
    raw%provenance%version = 'v1'
    raw%provenance%request_key = 'synthetic-benchmark'
    raw%provenance%record_id = 'synthetic-benchmark'
    raw%provenance%metadata = ''
    call raw%seal()
    reduction%id = 'synthetic-average-v1'
    reduction%units = 'V'
    allocate (reduction%weights(2, n))
    reduction%weights(1, :) = 1.0_dp/n
    reduction%weights(2, :) = 2.0_dp/n
    reduction%noise_covariance = reshape([0.01_dp, 0.0_dp, 0.0_dp, 0.02_dp], [2, 2])
    call system_clock(count_rate=rate)
    do mode = 1, 4
        call system_clock(start)
        do i = 0, repeats
            if (i == 1) call system_clock(start)
            select case (mode)
            case (1)
                call diag%predict(sample, eta, y)
            case (2)
                call diag%jvp(sample, eta, qdot, etadot, y)
            case (3)
                call diag%vjp(sample, eta, ybar, qbar, etabar)
            case (4)
                call reduction%apply(raw, reduced, status)
                if (status /= 0) error stop 'reduction'
            end select
        end do
        call system_clock(finish)
        times(mode) = real(finish - start, dp)/real(rate, dp)/repeats
    end do
    c = prior%physical(eta)
    dc = prior%jvp(etadot)
    expected = prior%vjp([2.1_dp*n, 0.7_dp*n])
    relative_error = maxval(abs(etabar - expected))/maxval(abs(expected))
    if (relative_error > 4*n*epsilon(1.0_dp)) error stop 'VJP reduction error'
    if (maxval(abs(qbar(Q_TE, :) - 0.7_dp*c(1))) > epsilon(1.0_dp)) &
        error stop 'VJP channel error'
    if (maxval(abs(y - (0.2_dp*c(1) + 3*dc(1) + dc(2)))) > 8*epsilon(1.0_dp)) &
        error stop 'JVP error'
    if (maxval(abs(reduced%samples - [3.0_dp, 6.0_dp])) > 4*n*epsilon(1.0_dp)) &
        error stop 'reduced value'
    checksum = sum(qbar) + sum(etabar) + sum(y)
   print '(A)', 'channels,predict_s,jvp_s,vjp_s,reduction_s,checksum,vjp_relative_error'
    print '(I0,6(",",ES24.16))', n, times, checksum, relative_error
end program
