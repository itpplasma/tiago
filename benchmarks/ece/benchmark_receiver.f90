program benchmark_ece_receiver
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_ece_receiver, only: ece_receiver_t, ece_response_channel_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    integer, parameter :: nc = 256, na = 4, nf = 4, nr = nc*na*nf, repeats = 1000
    type(ece_response_channel_t) :: channel(nc)
    type(ece_receiver_t) :: receiver
    type(calibration_prior_t) :: prior
    real(dp), allocatable :: covariance(:, :), mean(:), eta(:), temperature(:)
    real(dp) :: voltage(nc), variance(nc), elapsed, error
    integer(int64) :: start, finish, rate
    integer :: i, j, c, status
    allocate (covariance(3*nc, 3*nc), mean(3*nc), eta(3*nc), temperature(nr))
    covariance = 0.0_dp
    do i = 1, 3*nc
        covariance(i, i) = .01_dp
    end do
    do c = 1, nc
        mean(3*c - 2:3*c) = [.003_dp, .01_dp, 10.0_dp]
        channel(c)%rays = reshape([(i, i = (c - 1)*na*nf + 1, c*na*nf)], [na, nf])
        channel(c)%frequency = [(1.0e11_dp + 1.0e9_dp*j, j = 1, nf)]
        channel(c)%weights = [(1.0e9_dp, j = 1, nf)]
        channel(c)%response = [(1.0_dp, j = 1, nf)]
        channel(c)%antenna = [(1.0_dp, j = 1, na)]
        channel(c)%duration = .001_dp
    end do
    call prior%init(mean, covariance, status)
    if (status /= 0) error stop 'benchmark prior'
    call receiver%init(channel, nr, prior, status)
    if (status /= 0) error stop 'benchmark receiver'
    temperature = 500.0_dp
    eta = 0.0_dp
    call receiver%evaluate(temperature, eta, voltage, variance, status)
    if (status /= 0) error stop 'benchmark warmup'
    call system_clock(start, rate)
    do i = 1, repeats
        call receiver%evaluate(temperature, eta, voltage, variance, status)
        if (status /= 0) error stop 'benchmark evaluation'
    end do
    call system_clock(finish)
    elapsed = real(finish - start, dp)/rate
    error = maxval(abs(voltage - 1.51_dp))
    if (error > 1.0e-12_dp) error stop 'benchmark voltage accuracy'
    if (maxval(abs(variance/((.003_dp*510)**2/4.0e6_dp) - 1)) > 1.0e-12_dp) &
        error stop 'benchmark variance accuracy'
    print '(A)', 'channels,rays,repeats,seconds,channels_per_second,max_error_V'
    print '(I0,A,I0,A,I0,A,ES16.8,A,ES16.8,A,ES16.8)', &
        nc, ',', nr, ',', repeats, ',', elapsed, ',', nc*real(repeats, dp)/elapsed, &
        ',', error
end program
