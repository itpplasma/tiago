program benchmark_magnetics
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_max_threads
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, &
                                        Q_TE, Q_NE, Q_TI, Q_ZEFF
    use tiago_point_probe, only: point_probe_t
    use tiago_integrator_diagnostic, only: integrator_diagnostic_t
    use tiago_energy_consistency, only: energy_consistency_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    integer, parameter :: n = 1000000, reps = 30
    type(point_probe_t) :: probe
    type(integrator_diagnostic_t) :: diag
    type(energy_consistency_t) :: energy
    type(calibration_prior_t) :: prior
    type(plasma_sample_t) :: sample
    real(dp), allocatable :: y(:), dy(:), qb(:, :), qd(:, :), yb(:), time(:)
    real(dp) :: eta(3), ed(3), eb(3), ee(4), eed(4), eeb(4), ey(1), eyb(1)
    real(dp) :: cov(3, 3), elapsed, checksum, err, lhs, rhs
    integer(int64) :: first, last, rate
    integer :: status, i, mode, repeat
    allocate (probe%quantity(n), probe%x0(3, n), y(n), dy(n), yb(n), time(n))
    probe%quantity = Q_TE
    probe%x0 = 0.0_dp
    cov = 0.0_dp
    do i = 1, 3
        cov(i, i) = 0.01_dp
    end do
    do i = 1, n
        time(i) = real(i - 1, dp)/n
    end do
    call prior%init([1.2_dp, 0.1_dp, -0.03_dp], cov, status)
    if (status /= 0) error stop 'benchmark prior'
    call diag%init(probe, prior, time, 0.0_dp, status)
    if (status /= 0) error stop 'benchmark integrator'
    sample%n = n
    allocate (sample%q(N_QUANTITIES, n), qd(N_QUANTITIES, n), qb(N_QUANTITIES, n))
    sample%q = 0.0_dp
    sample%q(Q_TE, :) = 300.0_dp
    sample%q(Q_TI, :) = 200.0_dp
    sample%q(Q_NE, :) = 1.0e19_dp
    sample%q(Q_ZEFF, :) = 1.0_dp
    qd = 0.01_dp*sample%q
    yb = 0.7_dp
    eta = 0.1_dp
    ed = 0.2_dp
    energy%points = probe%x0
    allocate (energy%volume_weights(n))
    energy%volume_weights = 1.0_dp/n
    ee = [-0.001_dp, 2.0_dp, 10.0_dp, 0.0_dp]
    eed = [0.00001_dp, 0.1_dp, 0.01_dp, 0.0_dp]
    eyb = 0.7_dp
    call system_clock(count_rate=rate)
    print '(A,I0,A,I0)', 'threads=', omp_get_max_threads(), ',samples=', n
    print '(A)', 'operation,repeat,seconds_per_call,checksum'
    eeb = 0.0_dp
    do mode = 1, 6
        do repeat = 0, 3
            call system_clock(first)
            do i = 1, reps
                select case (mode)
                case (1)
                    call diag%predict(sample, eta, y)
                case (2)
                    call diag%jvp(sample, eta, qd, ed, dy)
                case (3)
                    call diag%vjp(sample, eta, yb, qb, eb)
                case (4)
                    call energy%predict(sample, ee, ey)
                case (5)
                    call energy%jvp(sample, ee, qd, eed, ey)
                case (6)
                    call energy%vjp(sample, ee, eyb, qb, eeb)
                end select
            end do
            call system_clock(last)
            elapsed = real(last - first, dp)/rate/reps
            select case (mode)
            case (1)
                checksum = sum(y)
            case (2)
                checksum = sum(dy)
            case (3)
                checksum = sum(qb) + sum(eb)
            case default
                checksum = sum(ey) + sum(qb) + sum(eeb)
            end select
            if (repeat > 0) print '(I0,A,I0,2(A,ES24.16))', &
                mode, ',', repeat, ',', elapsed, ',', checksum
        end do
    end do
    call diag%jvp(sample, eta, qd, ed, dy)
    call diag%vjp(sample, eta, yb, qb, eb)
    lhs = dot_product(yb, dy)
    rhs = sum(qb*qd) + dot_product(eb, ed)
    err = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
    if (err > 1.0e-9_dp) error stop 'benchmark integrator adjoint'
    print '(A,ES14.6)', 'integrator_adjoint_relative_error=', err
    call energy%jvp(sample, ee, qd, eed, ey)
    call energy%vjp(sample, ee, eyb, qb, eeb)
    lhs = dot_product(eyb, ey)
    rhs = sum(qb*qd) + dot_product(eeb, eed)
    err = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
    if (err > 1.0e-9_dp) error stop 'benchmark energy adjoint'
    print '(A,ES14.6)', 'energy_adjoint_relative_error=', err
end program
