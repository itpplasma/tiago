program consumer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_recorded_affine, only: recorded_affine
    use tiago_student_residual, only: student_residual
    implicit none
    real(dp) :: signal, d_raw, d_gain, d_offset, value, derivative, expected

    call recorded_affine(2.0_dp, 3.0_dp, 1.0_dp, signal, d_raw, d_gain, d_offset)
    if (abs(signal - 7.0_dp) > 1.0e-12_dp) error stop 'affine signal'
    if (abs(d_raw - 3.0_dp) > 1.0e-12_dp) error stop 'affine raw derivative'
    if (abs(d_gain - 2.0_dp) > 1.0e-12_dp) error stop 'affine gain derivative'
    if (abs(d_offset - 1.0_dp) > 1.0e-12_dp) error stop 'affine offset derivative'

    call student_residual(2.0_dp, 4.0_dp, value, derivative)
    expected = sqrt(5.0_dp*log(2.0_dp))
    if (abs(value - expected) > 1.0e-12_dp) error stop 'student value'
    if (abs(derivative - 1.25_dp/expected) > 1.0e-12_dp) error stop 'student derivative'
end program consumer
