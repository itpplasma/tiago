module tiago_doppler_ray
    !! Doppler reflectometry along a cold ray: turning point, Bragg wavenumber
    !! and frequency shift, with tangents that include the moving turning point.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use tiago_ece_ray, only: cold_medium_t, cold_ray_t, trace_cold_ray, &
        cold_ray_jvp, ray_rhs
    use tiago_doppler_turning_shift, only: doppler_turning_shift
    implicit none
    private

    type, public :: doppler_turning_t
        logical :: valid = .false.
        !! Ray traced from the launch exactly to the turning point.
        type(cold_ray_t) :: ray
        real(dp) :: time = 0.0_dp
        real(dp) :: state(6) = 0.0_dp
        real(dp) :: field(3) = 0.0_dp
        real(dp) :: density_gradient(3) = 0.0_dp
        !! Signed Bragg wavenumber -2 k0 |N_perp| in 1/m.
        real(dp) :: wavenumber = 0.0_dp
    end type

    public :: find_doppler_turning, doppler_ray_shift, doppler_ray_shift_jvp
contains

    subroutine find_doppler_turning(medium, frequency, mode, initial, duration, &
            turning, status, rtol, atol)
        class(cold_medium_t), intent(in) :: medium
        real(dp), intent(in) :: frequency, initial(6), duration
        integer, intent(in) :: mode
        type(doppler_turning_t), intent(out) :: turning
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: rtol, atol
        type(cold_ray_t) :: ray
        real(dp) :: g, gdot, lo, hi, t, step, g_prev, g_k, ignored
        integer :: k, iteration, last
        logical :: bracketed

        call trace_cold_ray(medium, frequency, mode, initial, duration, ray, &
            status, rtol=rtol, atol=atol)
        if (status%code /= 0) return
        ! The index modulus is minimal where g = N.dN/dt changes sign upward.
        bracketed = .false.
        call turning_function(medium, frequency, mode, ray%trace%y(:, 1), g_prev, &
            ignored, status)
        if (status%code /= 0) return
        do k = 2, size(ray%trace%t)
            call turning_function(medium, frequency, mode, ray%trace%y(:, k), &
                g_k, ignored, status)
            if (status%code /= 0) return
            if (g_prev < 0.0_dp .and. g_k >= 0.0_dp) then
                lo = ray%trace%t(k - 1)
                hi = ray%trace%t(k)
                t = lo + (hi - lo)*g_prev/(g_prev - g_k)
                bracketed = .true.
                exit
            end if
            g_prev = g_k
        end do
        if (.not. bracketed) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'Doppler ray does not pass a turning point within the duration')
            return
        end if
        last = 0
        do iteration = 1, 60
            call trace_cold_ray(medium, frequency, mode, initial, t, ray, status, &
                rtol=rtol, atol=atol)
            if (status%code /= 0) return
            call turning_function(medium, frequency, mode, &
                ray%trace%y(:, size(ray%trace%t)), g, gdot, status)
            if (status%code /= 0) return
            if (g < 0.0_dp) lo = t
            if (g >= 0.0_dp) hi = t
            if (gdot > 0.0_dp) then
                step = -g/gdot
            else
                step = 0.5_dp*(lo + hi) - t
            end if
            if (t + step <= lo .or. t + step >= hi) step = 0.5_dp*(lo + hi) - t
            if (abs(step) <= 1.0e-12_dp*max(1.0_dp, t)) then
                last = iteration
                exit
            end if
            t = t + step
        end do
        if (last == 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'Doppler turning point did not converge')
            return
        end if
        turning%ray = ray
        turning%time = t
        turning%state = ray%trace%y(:, size(ray%trace%t))
        call sample_turning(medium, turning)
        turning%valid = .true.
    end subroutine

    subroutine doppler_ray_shift(turning, radial_field, phase_velocity, omega, &
            wavenumber)
        !! Received angular frequency shift in rad/s for an outward radial
        !! field in V/m and a fluctuation phase velocity in m/s along the
        !! perpendicular index direction.
        type(doppler_turning_t), intent(in) :: turning
        real(dp), intent(in) :: radial_field, phase_velocity
        real(dp), intent(out) :: omega, wavenumber
        real(dp) :: d_omega(12), d_wavenumber(12)
        call evaluate(turning, turning%ray%frequency, radial_field, phase_velocity, omega, &
            wavenumber, d_omega, d_wavenumber)
    end subroutine

    subroutine doppler_ray_shift_jvp(medium, turning, radial_field, &
            phase_velocity, initial_dot, parameter_dot, frequency_dot, &
            radial_field_dot, phase_velocity_dot, omega_dot, wavenumber_dot, &
            position_dot, status)
        class(cold_medium_t), intent(in) :: medium
        type(doppler_turning_t), intent(in) :: turning
        real(dp), intent(in) :: radial_field, phase_velocity
        real(dp), intent(in) :: initial_dot(6), parameter_dot(:), frequency_dot
        real(dp), intent(in) :: radial_field_dot, phase_velocity_dot
        real(dp), intent(out) :: omega_dot, wavenumber_dot, position_dot(3)
        type(fortnum_status_t), intent(out) :: status
        integer :: np
        real(dp) :: fixed_dot(6), state_dot(6), f(6), jac(6, 6), h, dgdy(6)
        real(dp) :: pj(6, medium%n_params() + 1), dgdp(medium%n_params() + 1)
        real(dp) :: gdot, time_dot, omega, wavenumber
        real(dp) :: q(4), dq(3, 4), d2q(3, 3, 4)
        real(dp) :: qp(4, medium%n_params()), dqp(3, 4, medium%n_params())
        real(dp) :: field_dot(3), gradient_dot(3), inputs_dot(12)
        real(dp) :: d_omega(12), d_wavenumber(12)
        logical :: valid
        np = medium%n_params()
        omega_dot = 0.0_dp
        wavenumber_dot = 0.0_dp
        position_dot = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid Doppler tangent input')
        if (.not. turning%valid) return
        if (size(parameter_dot) /= np) return
        call cold_ray_jvp(medium, turning%ray, initial_dot, parameter_dot, &
            frequency_dot, fixed_dot, status)
        if (status%code /= 0) return
        call ray_rhs(medium, np, turning%ray%frequency, turning%ray%mode, turning%state, f, &
            jac, pj, valid, h)
        if (.not. valid) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'Doppler turning state is singular')
            return
        end if
        ! Implicit turning time from g(y(t), p) = N.dN/dt = 0.
        dgdy = matmul(turning%state(4:), jac(4:, :))
        dgdy(4:) = dgdy(4:) + f(4:)
        dgdp = matmul(turning%state(4:), pj(4:, :))
        gdot = dot_product(dgdy, f)
        time_dot = -(dot_product(dgdy, fixed_dot) &
            + dot_product(dgdp(:np), parameter_dot) &
            + dgdp(np + 1)*frequency_dot)/gdot
        state_dot = fixed_dot + f*time_dot
        call medium%sample(turning%state(:3), q, dq, d2q, qp, dqp)
        field_dot = matmul(state_dot(:3), dq(:, :3)) + matmul(qp(:3, :), parameter_dot)
        gradient_dot = matmul(d2q(:, :, 4), state_dot(:3)) &
            + matmul(dqp(:, 4, :), parameter_dot)
        inputs_dot = [state_dot(4:), field_dot, gradient_dot, frequency_dot, &
            radial_field_dot, phase_velocity_dot]
        call evaluate(turning, turning%ray%frequency, radial_field, phase_velocity, &
            omega, wavenumber, d_omega, d_wavenumber)
        omega_dot = dot_product(d_omega, inputs_dot)
        wavenumber_dot = dot_product(d_wavenumber, inputs_dot)
        position_dot = state_dot(:3)
        call status_set(status, 0, '')
    end subroutine

    subroutine turning_function(medium, frequency, mode, y, g, gdot, status)
        class(cold_medium_t), intent(in) :: medium
        real(dp), intent(in) :: frequency, y(6)
        integer, intent(in) :: mode
        real(dp), intent(out) :: g, gdot
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: f(6), jac(6, 6), h, dgdy(6)
        logical :: valid
        call ray_rhs(medium, medium%n_params(), frequency, mode, y, f, jac, &
            valid=valid, h=h)
        g = 0.0_dp
        gdot = 0.0_dp
        if (.not. valid) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'Doppler ray state is singular')
            return
        end if
        g = dot_product(y(4:), f(4:))
        dgdy = matmul(y(4:), jac(4:, :))
        dgdy(4:) = dgdy(4:) + f(4:)
        gdot = dot_product(dgdy, f)
        call status_set(status, 0, '')
    end subroutine

    subroutine sample_turning(medium, turning)
        class(cold_medium_t), intent(in) :: medium
        type(doppler_turning_t), intent(inout) :: turning
        real(dp) :: q(4), dq(3, 4), d2q(3, 3, 4), omega, d_omega(12), d_k(12)
        call medium%sample(turning%state(:3), q, dq, d2q)
        turning%field = q(:3)
        turning%density_gradient = dq(:, 4)
        call evaluate(turning, turning%ray%frequency, 0.0_dp, 0.0_dp, omega, &
            turning%wavenumber, d_omega, d_k)
    end subroutine

    subroutine evaluate(turning, frequency, radial_field, phase_velocity, omega, &
            wavenumber, d_omega, d_wavenumber)
        type(doppler_turning_t), intent(in) :: turning
        real(dp), intent(in) :: frequency, radial_field, phase_velocity
        real(dp), intent(out) :: omega, wavenumber, d_omega(12), d_wavenumber(12)
        associate (n => turning%state(4:), b => turning%field, &
                   g => turning%density_gradient)
            call doppler_turning_shift(n(1), n(2), n(3), b(1), b(2), b(3), &
                g(1), g(2), g(3), frequency, radial_field, phase_velocity, omega, &
                wavenumber, d_omega, d_wavenumber)
        end associate
    end subroutine
end module tiago_doppler_ray
