program test_integrated_drift
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_signal_source, only: signal_window_t
    use tiago_integrated_drift, only: integrated_drift_t
    implicit none
    type(signal_window_t) :: input, output, plus, minus, op, om
    type(integrated_drift_t) :: reduction
    real(dp) :: covariance(6, 6), expected(2, 2), weights(2, 6), &
                direction(6), tangent(2), adjoint(6), ybar(2), err, h, fit(2, 2)
    integer :: i, j, status
    input%time = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    input%samples = 0.7_dp - 0.03_dp*input%time
    input%samples(3:4) = input%samples(3:4) + [1.2_dp, 2.1_dp]
    input%variance = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp]
    input%units = 'Wb'
    allocate (input%covariance(6, 6))
    ! Shared offset noise is correlated across baseline and active samples.
    input%covariance = 0.002_dp
    input%provenance%source_id = 'synthetic-integrator'
    input%provenance%version = '1'
    input%provenance%request_key = 'synthetic-drift'
    input%provenance%record_id = 'synthetic-drift-v1'
    input%provenance%metadata = 'seconds; zero physical flux at samples 1 and 6'
    call input%seal()
    call reduction%init(input, [1, 6], [3, 4], 0.0_dp, 'endpoint-baseline-v1', status)
    call check(status == 0, 'baseline setup')
    call reduction%apply(input, output, status)
    call check(status == 0, 'drift removal')
    err = maxval(abs(output%samples - [1.2_dp, 2.1_dp]))
    call check(err < 1.0e-14_dp, 'synthetic linear drift recovery')
    print '(A,ES14.6)', 'drift-removed signal absolute error [Wb]: ', err
    call check(all(output%time == [2.0_dp, 3.0_dp]), 'retain physical time tags')
    weights = 0.0_dp
    weights(1, :) = [-0.6_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -0.4_dp]
    weights(2, :) = [-0.4_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -0.6_dp]
    covariance = input%total_covariance()
    expected = matmul(weights, matmul(covariance, transpose(weights)))
    err = maxval(abs(output%total_covariance() - expected))
    call check(err < 1.0e-15_dp, 'covariance with baseline / signal cross terms')
    print '(A,ES14.6)', 'drift covariance max absolute error [Wb^2]: ', err
    fit = reduction%baseline_covariance()
    call check(abs(fit(1, 1) - 0.012_dp) < 1.0e-15_dp, 'offset fit variance')
    call check(abs(fit(2, 2) - 0.07_dp/25) < 1.0e-15_dp, 'drift fit variance')
    direction = [0.4_dp, -0.1_dp, 0.3_dp, 0.7_dp, -0.2_dp, -0.8_dp]
    h = 1.0e-5_dp
    plus = input
    minus = input
    plus%samples = input%samples + h*direction
    minus%samples = input%samples - h*direction
    call plus%seal()
    call minus%seal()
    call reduction%apply(plus, op, status)
    call check(status == 0, 'positive perturbed signal')
    call reduction%apply(minus, om, status)
    call check(status == 0, 'negative perturbed signal')
    call reduction%jvp(input, direction, tangent)
    err = maxval(abs(tangent - (op%samples - om%samples)/(2*h)))
    call check(err < 1.0e-10_dp, 'drift reduction FD')
    print '(A,ES14.6)', 'drift reduction FD absolute error: ', err
    ybar = [0.2_dp, -0.5_dp]
    call reduction%vjp(input, ybar, adjoint)
    err = abs(dot_product(ybar, tangent) - dot_product(adjoint, direction))
    call check(err < 1.0e-15_dp, 'drift reduction transpose')
    print '(A,ES14.6)', 'drift reduction transpose absolute error: ', err
    call check(index(output%provenance%metadata, input%provenance%checksum) > 0, &
               'reduction retains parent checksum')
    ! Overdetermined GLS must annihilate the same baseline and retain PSD covariance.
    call reduction%init(input, [1, 2, 5, 6], [3, 4], 2.5_dp, 'gls-baseline-v1', status)
    call check(status == 0, 'correlated overdetermined baseline')
    call reduction%apply(input, output, status)
    call check(status == 0, 'GLS drift removal')
 call check(maxval(abs(output%samples - [1.2_dp, 2.1_dp])) < 1.0e-14_dp, 'GLS recovery')
    ! Reconstruct the linear operator by basis probes for a separate propagation oracle.
    do j = 1, 6
        direction = 0.0_dp
        direction(j) = 1.0_dp
        call reduction%jvp(input, direction, weights(:, j))
    end do
    expected = 0.0_dp
    do i = 1, 6
        do j = 1, 6
            expected = expected + covariance(i, j)* &
                       spread(weights(:, i), 2, 2)*spread(weights(:, j), 1, 2)
        end do
    end do
    call check(maxval(abs(output%total_covariance() - expected)) < 1.0e-15_dp, &
               'GLS full correlated covariance')
    plus%time(1) = -1.0_dp
    call reduction%apply(plus, output, status)
    call check(status /= 0, 'reject changed time grid')
    call reduction%init(input, [1, 1], [3, 4], 0.0_dp, 'bad', status)
    call check(status /= 0, 'reject repeated baseline index')
    call reduction%init(input, [1, 6], [1, 4], 0.0_dp, 'bad', status)
    call check(status /= 0, 'reject baseline output overlap')
    plus = input
    plus%time = 0.0_dp
    call reduction%init(plus, [1, 6], [3, 4], 0.0_dp, 'bad', status)
    call check(status /= 0, 'reject rank deficient drift fit')
contains
    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (ok) return
        print '(A)', 'FAIL: '//text
        error stop 1
    end subroutine
end program
