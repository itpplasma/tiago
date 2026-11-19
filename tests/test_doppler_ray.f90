program test_doppler_ray
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ece_ray_fixture, only: polynomial_medium_t
    use tiago_ece_ray, only: ECE_O_MODE, ECE_X_MODE
    use tiago_doppler_ray, only: doppler_turning_t, find_doppler_turning, &
        doppler_ray_shift, doppler_ray_shift_jvp
    use fortnum_status, only: fortnum_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), light = 299792458.0_dp
    real(dp), parameter :: frequency = 9.0e10_dp
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: charge = 1.602176634e-19_dp, eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: critical = (2*pi*frequency)**2*eps0*electron_mass/charge**2
    real(dp), parameter :: rtol = 1.0e-11_dp, atol = 1.0e-13_dp
    real(dp), parameter :: radial_field = -2.0e4_dp, phase_velocity = 800.0_dp
    type(polynomial_medium_t) :: medium
    type(doppler_turning_t) :: turning
    type(fortnum_status_t) :: status
    real(dp) :: k0, sin_tilt, initial(6), omega, wavenumber, expected, bz, error
    integer :: i

    k0 = 2*pi*frequency/light

    ! Gusakov and Surkov (2004), slab O-mode: the poloidal index K/k0 = sin(theta)
    ! is conserved, the ray turns where X = cos^2(theta), and backscattering
    ! selects q = -2K, so the shift is q v with v = E_r/B + v_ph.
    medium%p = 0.0_dp
    medium%p(3) = 1.0_dp
    medium%p(5) = 0.6_dp*critical
    error = 0.0_dp
    do i = 1, 4
        sin_tilt = 0.1_dp*i
        initial = [0.0_dp, 0.0_dp, 0.0_dp, sqrt(1 - sin_tilt**2), sin_tilt, 0.0_dp]
        call find_doppler_turning(medium, frequency, ECE_O_MODE, initial, &
            3*initial(4)/0.6_dp, turning, status, rtol, atol)
        call check(status%code == 0, 'slab turning '//status%msg)
        error = max(error, abs(turning%state(1) - (1 - sin_tilt**2)/0.6_dp))
        error = max(error, abs(turning%state(2) - 2*sin_tilt*initial(4)/0.6_dp))
        call doppler_ray_shift(turning, radial_field, phase_velocity, omega, &
            wavenumber)
        bz = 1 + turning%state(1)/30
        call check(abs(wavenumber + 2*k0*sin_tilt) < 1.0e-8_dp*2*k0*sin_tilt, &
            'slab Bragg wavenumber')
        expected = -2*k0*sin_tilt*(radial_field/bz + phase_velocity)
        call check(abs(omega - expected) < 1.0e-8_dp*abs(expected), 'slab Doppler shift')
    end do
    call check(error < 1.0e-8_dp, 'slab turning point')
    print '(A,ES12.4)', 'slab turning-point error [m]: ', error

    call check_tangents(ECE_O_MODE, 4.0_dp)
    call check_tangents(ECE_X_MODE, 2.0_dp)

    call find_doppler_turning(medium, frequency, ECE_O_MODE, &
        [0.0_dp, 0.0_dp, 0.0_dp, sqrt(0.84_dp), 0.4_dp, 0.0_dp], 0.5_dp, turning, status)
    call check(status%code /= 0, 'short ray must not report a turning point')
contains
    subroutine check_tangents(mode, duration)
        integer, intent(in) :: mode
        real(dp), intent(in) :: duration
        type(polynomial_medium_t) :: base, trial
        type(doppler_turning_t) :: nominal, plus, minus
        real(dp) :: launch(6), direction(3), tangent(3), p_dot(8), y_dot(6)
        real(dp) :: omega_dot, k_dot, x_dot(3), worst, step
        real(dp) :: op, om, kp, km, fd_omega, fd_k, fd_x(3)
        integer :: j
        base%p = [0.3_dp, 0.2_dp, 2.0_dp, 0.0_dp, 0.5_dp*critical, &
            0.05_dp*critical, 0.02_dp*critical, 0.1_dp*critical]
        direction = [0.9_dp, 0.4_dp, 0.17_dp]
        direction = direction/norm2(direction)
        launch = [0.0_dp, 0.0_dp, 0.0_dp, direction]
        call find_doppler_turning(base, frequency, mode, launch, duration, nominal, &
            status, rtol, atol)
        call check(status%code == 0, 'oblique turning '//status%msg)
        worst = 0.0_dp
        ! Parameters 5..8 and the frequency keep the vacuum launch on shell.
        do j = 5, 9
            p_dot = 0.0_dp
            if (j <= 8) then
                p_dot(j) = 1.0_dp
                step = 1.0e-4_dp*base%p(j)
            else
                step = 1.0e-4_dp*frequency
            end if
            call doppler_ray_shift_jvp(base, nominal, radial_field, phase_velocity, &
                [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], p_dot, &
                merge(1.0_dp, 0.0_dp, j == 9), 0.0_dp, 0.0_dp, omega_dot, k_dot, &
                x_dot, status)
            call check(status%code == 0, 'Doppler tangent '//status%msg)
            trial = base
            trial%p = base%p + step*p_dot
            call shift_of(trial, frequency + merge(step, 0.0_dp, j == 9), launch, &
                mode, duration, plus, op, kp)
            trial%p = base%p - step*p_dot
            call shift_of(trial, frequency - merge(step, 0.0_dp, j == 9), launch, &
                mode, duration, minus, om, km)
            fd_omega = (op - om)/(2*step)
            fd_k = (kp - km)/(2*step)
            fd_x = (plus%state(:3) - minus%state(:3))/(2*step)
            worst = max(worst, abs(fd_omega - omega_dot)/abs(omega_dot))
            worst = max(worst, abs(fd_k - k_dot)/max(abs(k_dot), 1.0e-30_dp*k0))
            worst = max(worst, norm2(fd_x - x_dot)/norm2(x_dot))
        end do
        ! Launch direction, kept on the unit sphere.
        tangent = [-0.4_dp, 0.9_dp, 0.3_dp]
        tangent = tangent - dot_product(tangent, direction)*direction
        tangent = tangent/norm2(tangent)
        y_dot = [0.0_dp, 0.0_dp, 0.0_dp, tangent]
        p_dot = 0.0_dp
        call doppler_ray_shift_jvp(base, nominal, radial_field, phase_velocity, &
            y_dot, p_dot, 0.0_dp, 0.0_dp, 0.0_dp, omega_dot, k_dot, x_dot, status)
        call check(status%code == 0, 'direction tangent '//status%msg)
        step = 1.0e-6_dp
        call shift_of(base, frequency, [0.0_dp, 0.0_dp, 0.0_dp, &
            (direction + step*tangent)/norm2(direction + step*tangent)], mode, &
            duration, plus, op, kp)
        call shift_of(base, frequency, [0.0_dp, 0.0_dp, 0.0_dp, &
            (direction - step*tangent)/norm2(direction - step*tangent)], mode, &
            duration, minus, om, km)
        worst = max(worst, abs((op - om)/(2*step) - omega_dot)/abs(omega_dot))
        worst = max(worst, norm2((plus%state(:3) - minus%state(:3))/(2*step) - x_dot) &
            /norm2(x_dot))
        ! The shift is linear in the radial field and the phase velocity.
        call doppler_ray_shift(nominal, radial_field, phase_velocity, op, kp)
        call doppler_ray_shift(nominal, 0.0_dp, 0.0_dp, om, km)
        call doppler_ray_shift_jvp(base, nominal, radial_field, phase_velocity, &
            [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], p_dot, 0.0_dp, &
            radial_field, phase_velocity, omega_dot, k_dot, x_dot, status)
        call check(abs(om) == 0.0_dp, 'zero velocity gives zero shift')
        call check(abs(omega_dot - op) < 1.0e-12_dp*abs(op), 'linear velocity tangent')
        call check(worst < 1.0e-5_dp, 'Doppler tangent against finite differences')
        print '(A,I0,A,ES12.4)', 'mode ', mode, ' worst relative tangent error: ', worst
    end subroutine

    subroutine shift_of(trial, f, launch, mode, duration, turn, omega_value, k_value)
        type(polynomial_medium_t), intent(in) :: trial
        real(dp), intent(in) :: f, launch(6), duration
        integer, intent(in) :: mode
        type(doppler_turning_t), intent(out) :: turn
        real(dp), intent(out) :: omega_value, k_value
        call find_doppler_turning(trial, f, mode, launch, duration, turn, status, &
            rtol, atol)
        call check(status%code == 0, 'perturbed turning '//status%msg)
        call doppler_ray_shift(turn, radial_field, phase_velocity, omega_value, k_value)
    end subroutine

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(A)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program test_doppler_ray
