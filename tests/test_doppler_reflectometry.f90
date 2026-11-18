program test_doppler_reflectometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_doppler_reflectometry, only: doppler_reflectometry
    implicit none
    real(dp), parameter :: x(5) = [-20000.0_dp, 300.0_dp, 700.0_dp, 4.0_dp, -2.0_dp]
    real(dp) :: omega, row(5), trial(5), plus, minus, ignored(5), step, fd
    integer :: i

    call predict(x, omega, row)
    if (abs(omega - 7210000.0_dp) > 1.0e-8_dp) error stop 'Doppler shift and sign'
    do i = 1, size(x)
        step = 1.0e-5_dp*abs(x(i))
        trial = x
        trial(i) = x(i) + step
        call predict(trial, plus, ignored)
        trial(i) = x(i) - step
        call predict(trial, minus, ignored)
        fd = (plus - minus)/(2.0_dp*step)
        if (abs(fd - row(i)) > 1.0e-8_dp*abs(row(i))) error stop 'Doppler gradient'
    end do
    trial = x
    trial(3) = -x(3)
    call predict(trial, plus, ignored)
    if (abs(plus + omega) > 1.0e-8_dp) error stop 'wavevector reversal'
    trial = x
    trial(5) = 0.0_dp
    call predict(trial, plus, ignored)
    if (abs(plus - x(2)*x(3)) > 1.0e-8_dp) error stop 'phase velocity on blind channel'
contains
    subroutine predict(values, signal, gradient)
        real(dp), intent(in) :: values(5)
        real(dp), intent(out) :: signal, gradient(5)
        call doppler_reflectometry(values(1), values(2), values(3), values(4), values(5), &
                                    signal, gradient(1), gradient(2), gradient(3), gradient(4), gradient(5))
    end subroutine predict
end program test_doppler_reflectometry
