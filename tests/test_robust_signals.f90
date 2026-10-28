program test_robust_signals
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_robust_signals, only: robust_rows, robust_normalization
    use tiago_calibration_campaign, only: campaign_factor, constrain_discrepancy
    implicit none
    real(dp)::r(7), original(7), slope(7), plus(7), minus(7), h, error
    real(dp)::l(4, 4), protected(4, 2), basis(4, 1), cov(3, 3), expected(3, 3)
    real(dp), allocatable::factor(:, :), constrained(:, :)
    integer::i, code
    original = [-100._dp, -1._dp, -1.e-9_dp, 0._dp, 1.e-9_dp, 1._dp, 100._dp]
    r = original
    call robust_rows(r, 4._dp, slope)
    error = maxval(abs(.5_dp*r**2 - 2.5_dp*log(1 + original**2/4)))
    print '(A,ES14.6)', 'Student log-density absolute error=', error
    if (error > 2.e-14_dp) error stop 'Student density reference'
    if (abs(slope(4) - sqrt(1.25_dp)) > 1.e-14_dp) error stop 'zero-residual derivative'
    h = 1.e-5_dp
    plus = original + h; minus = original - h
    call robust_rows(plus, 4._dp); call robust_rows(minus, 4._dp)
    error = maxval(abs((plus - minus)/(2*h) - slope))
    print '(A,ES14.6)', 'Student derivative absolute error=', error
    if (error > 1.e-9_dp) error stop 'Student derivative'
    if (abs(robust_normalization(4._dp) - log(.375_dp*sqrt(2*acos(-1._dp)))) &
        > 1.e-14_dp) error stop 'Student normalization'
    call campaign_factor([0._dp, 1._dp, 3._dp], .2_dp, .1_dp, factor, code)
    if (code /= 0) error stop 'campaign factor'
    cov = matmul(factor, transpose(factor))
    expected = reshape([.04_dp, .04_dp, .04_dp, .04_dp, .05_dp, .05_dp, &
                        .04_dp, .05_dp, .07_dp], [3, 3])
    error = maxval(abs(cov - expected))
    print '(A,ES14.6)', 'Brownian calibration covariance error=', error
    if (error > 2.e-17_dp) error stop 'campaign covariance'
    call campaign_factor([1._dp, 0._dp], .2_dp, .1_dp, factor, code)
    if (code == 0) error stop 'unsorted campaign accepted'
    l = 0
    do i = 1, 4
        l(i, i) = i
    end do
    l(3, 1) = .3_dp
    protected(:, 1) = [1._dp, 2._dp, 1._dp, 4._dp]
    protected(:, 2) = 2*protected(:, 1)
    basis(:, 1) = [3._dp, -1._dp, 2._dp, -3._dp]
    call constrain_discrepancy(l, protected, basis, constrained, code)
    if (code /= 0) error stop 'discrepancy projection'
    do i = 1, 4
        protected(i, :) = (protected(i, :) - matmul(l(i, :i - 1), protected(:i - 1, :)))/l(i, i)
        constrained(i, :) = (constrained(i, :) - matmul(l(i, :i - 1), constrained(:i - 1, :)))/l(i, i)
    end do
    error = maxval(abs(matmul(transpose(protected), constrained)))
    print '(A,ES14.6)', 'protected-space discrepancy overlap=', error
    if (error > 1.e-13_dp) error stop 'discrepancy absorbs protected signals'
    if (norm2(constrained) < .1_dp) error stop 'discrepancy vanished unnecessarily'
end program
