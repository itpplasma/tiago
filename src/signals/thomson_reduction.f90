module tiago_thomson_reduction
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_signal_reduction, only: linear_reduction_t
    use tiago_ts_noise, only: ts_noise
    implicit none
    private
    public :: pulse_integral, thomson_noise
contains
    subroutine pulse_integral(step, time, pulse, baseline, status, extra_variance)
        !! Trapezoidal pulse integral minus the time-averaged baseline.
        !! Gates name inclusive endpoint indices, with no shared intervals.
        type(linear_reduction_t), intent(out) :: step
        real(dp), intent(in) :: time(:)
        integer, intent(in) :: pulse(2), baseline(2)
        integer, intent(out) :: status
        real(dp), optional, intent(in) :: extra_variance
        real(dp) :: p(size(time)), b(size(time)), dt, noise
        integer :: i, n
        status = 1
        n = size(time)
        if (n < 2) return
        if (.not. all(ieee_is_finite(time))) return
        if (any(time(2:) <= time(:n - 1))) return
        if (any(pulse < 1) .or. any(pulse > n)) return
        if (any(baseline < 1) .or. any(baseline > n)) return
        if (pulse(2) <= pulse(1) .or. baseline(2) <= baseline(1)) return
        if (pulse(1) < baseline(2) .and. baseline(1) < pulse(2)) return
        noise = 0.0_dp
        if (present(extra_variance)) noise = extra_variance
        if (.not. ieee_is_finite(noise)) return
        if (noise < 0.0_dp) return
        p = 0.0_dp
        b = 0.0_dp
        do i = pulse(1), pulse(2) - 1
            dt = (time(i + 1) - time(i))/2
            p(i:i + 1) = p(i:i + 1) + dt
        end do
        do i = baseline(1), baseline(2) - 1
            dt = (time(i + 1) - time(i))/2
            b(i:i + 1) = b(i:i + 1) + dt
        end do
        allocate (step%weights(1, n), step%noise_covariance(1, 1))
        step%weights(1, :) = p - sum(p)*b/sum(b)
        step%noise_covariance(1, 1) = noise
        step%id = 'thomson-pulse-trapezoid-baseline-1'
        step%units = 'V s'
        status = 0
    end subroutine
    subroutine thomson_noise(photons, bg, gain, digitvar, excess, bgscale, variance, row)
        !! bg: expected background photoelectrons inside pulse gate.
        !! bgscale: pulse duration / independent baseline duration.
        !! digitvar: already integrated digitiser variance [V^2 s^2].
        real(dp), intent(in) :: photons, bg, gain, digitvar, excess, bgscale
        real(dp), intent(out) :: variance, row(6)
        if (.not. all(ieee_is_finite([photons, bg, gain, digitvar, excess, bgscale]))) &
            error stop 'Nonfinite Thomson noise input'
        if (min(photons, bg, digitvar, bgscale) < 0.0_dp .or. gain <= 0.0_dp .or. excess < 1.0_dp) &
            error stop 'Invalid Thomson noise input'
        call ts_noise(photons, bg, gain, digitvar, excess, bgscale, variance, &
                      row(1), row(2), row(3), row(4), row(5), row(6))
    end subroutine
end module
