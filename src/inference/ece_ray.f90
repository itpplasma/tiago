module tiago_ece_ray
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp, ode_integrate_parameter_vjp
    use tiago_ece_hamiltonian_o, only: ece_hamiltonian_o
    use tiago_ece_hamiltonian_x, only: ece_hamiltonian_x
    use tiago_ece_hamiltonian_unmagnetized, &
        only: ece_hamiltonian_unmagnetized
    implicit none
    private
    integer, parameter, public :: ECE_O_MODE = 1, ECE_X_MODE = 2

    type, abstract, public :: cold_medium_t
    contains
        procedure(sample_i), deferred :: sample
        procedure(count_i), deferred :: n_params
    end type
    abstract interface
        pure integer function count_i(self)
            import :: cold_medium_t
            class(cold_medium_t), intent(in) :: self
        end function
        subroutine sample_i(self, x, q, dq, d2q, qp, dqp)
            import :: cold_medium_t, dp
            class(cold_medium_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: q(4), dq(3, 4), d2q(3, 3, 4)
            real(dp), intent(out) :: qp(:, :), dqp(:, :, :)
        end subroutine
    end interface
    type, public :: cold_ray_t
        logical :: valid = .false.
        logical :: derivatives_valid = .false.
        type(ode_solution_t) :: trace
        real(dp) :: frequency = 0.0_dp
        integer :: mode = ECE_O_MODE
        real(dp) :: max_hamiltonian_error = 0.0_dp
        real(dp) :: max_index = 0.0_dp
    end type
    public :: cold_hamiltonian, trace_cold_ray, trace_cold_rays, cold_ray_jvp, &
        cold_ray_vjp
contains
    subroutine cold_hamiltonian(n, b, ne, frequency, mode, h, dh, d2h)
        real(dp), intent(in) :: n(3), b(3), ne, frequency
        integer, intent(in) :: mode
        real(dp), intent(out) :: h, dh(8), d2h(8, 8)
        if (norm2(b) == 0.0_dp) then
            call ece_hamiltonian_unmagnetized(n(1), n(2), n(3), &
                b(1), b(2), b(3), ne, frequency, h, dh, d2h)
        else if (mode == ECE_O_MODE) then
            call ece_hamiltonian_o(n(1), n(2), n(3), &
                b(1), b(2), b(3), ne, frequency, h, dh, d2h)
        else
            call ece_hamiltonian_x(n(1), n(2), n(3), &
                b(1), b(2), b(3), ne, frequency, h, dh, d2h)
        end if
    end subroutine

    subroutine ray_rhs(medium, np, frequency, mode, y, f, jac, pj, valid, h)
        class(cold_medium_t), intent(in) :: medium
        real(dp), intent(in) :: frequency, y(6)
        integer, intent(in) :: np, mode
        real(dp), intent(out) :: f(6), jac(6, 6)
        real(dp), intent(out) :: pj(6, np + 1), h
        logical, intent(out) :: valid
        real(dp) :: q(4), dq(3, 4), d2q(3, 3, 4), dh(8), d2h(8, 8)
        real(dp) :: qp(4, medium%n_params()), dqp(3, 4, medium%n_params())
        real(dp) :: a(8, 6), ap(8, medium%n_params() + 1)
        integer :: i, j
        call medium%sample(y(:3), q, dq, d2q, qp, dqp)
        valid = all(ieee_is_finite(q)) .and. all(ieee_is_finite(y))
        valid = valid .and. all(ieee_is_finite(dq))
        valid = valid .and. all(ieee_is_finite(d2q))
        valid = valid .and. all(ieee_is_finite(qp))
        valid = valid .and. all(ieee_is_finite(dqp))
        valid = valid .and. q(4) >= 0.0_dp
        f = 0.0_dp
        jac = 0.0_dp
        pj = 0.0_dp
        h = 0.0_dp
        if (.not. valid) return
        if (norm2(q(:3)) > 0.0_dp .and. norm2(y(4:)) == 0.0_dp) then
            valid = .false.
            return
        end if
        call cold_hamiltonian(y(4:), q(:3), q(4), frequency, mode, h, dh, d2h)
        valid = ieee_is_finite(h) .and. all(ieee_is_finite(dh))
        valid = valid .and. all(ieee_is_finite(d2h))
        if (.not. valid) return
        f(:3) = dh(:3)
        f(4:) = -matmul(dq, dh(4:7))
        a = 0.0_dp
        a(4:7, :3) = transpose(dq)
        do i = 1, 3
            a(i, i + 3) = 1.0_dp
        end do
        jac(:3, :) = matmul(d2h(:3, :), a)
        jac(4:, :) = -matmul(dq, matmul(d2h(4:7, :), a))
        do j = 1, 3
            do i = 1, 3
                jac(i + 3, j) = jac(i + 3, j) - dot_product(dh(4:7), d2q(i, j, :))
            end do
        end do
        ap = 0.0_dp
        ap(4:7, :np) = qp
        ap(8, np + 1) = 1.0_dp
        pj(:3, :) = matmul(d2h(:3, :), ap)
        pj(4:, :) = -matmul(dq, matmul(d2h(4:7, :), ap))
        do j = 1, np
            pj(4:, j) = pj(4:, j) - matmul(dqp(:, :, j), dh(4:7))
        end do
    end subroutine

    subroutine trace_cold_ray(medium, frequency, mode, initial, duration, &
            ray, status, rtol, atol, hmax)
        class(cold_medium_t), intent(in) :: medium
        real(dp), intent(in) :: frequency, initial(6), duration
        integer, intent(in) :: mode
        type(cold_ray_t), intent(out) :: ray
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: rtol, atol, hmax
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        real(dp) :: f(6), jac(6, 6), pj(6, medium%n_params() + 1), h
        real(dp) :: q(4), dq(3, 4), d2q(3, 3, 4)
        real(dp) :: qp(4, medium%n_params()), dqp(3, 4, medium%n_params())
        logical :: valid, failed
        integer :: k
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid cold-ray input')
        if (.not. all(ieee_is_finite(initial))) return
        if (.not. ieee_is_finite(frequency)) return
        if (.not. ieee_is_finite(duration)) return
        if (frequency <= 0.0_dp .or. duration <= 0.0_dp) return
        if (mode /= ECE_O_MODE .and. mode /= ECE_X_MODE) return
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-11_dp
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol
        if (present(hmax)) problem%hmax = hmax
        call ray_rhs(medium, medium%n_params(), frequency, mode, initial, &
            f, jac, pj, valid, h)
        if (.not. valid) return
        if (abs(h) > 1.0e-7_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray launch is off dispersion surface')
            return
        end if
        problem%rhs => rhs
        problem%t1 = duration
        problem%y0 = initial
        failed = .false.
        call ode_integrate(problem, workspace, ray%trace, status)
        ray%frequency = frequency
        ray%mode = mode
        if (failed) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray medium or dispersion is singular')
            return
        end if
        if (status%code /= 0) return
        ray%derivatives_valid = .true.
        do k = 1, size(ray%trace%t)
            call ray_rhs(medium, medium%n_params(), frequency, mode, &
                ray%trace%y(:, k), f, jac, pj, valid, h)
            if (.not. valid) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'cold-ray accepted state is singular')
                return
            end if
            call medium%sample(ray%trace%y(:3, k), q, dq, d2q, qp, dqp)
            if (norm2(q(:3)) == 0.0_dp) then
                if (any(qp(:3, :) /= 0.0_dp) .or. any(dq(:, :3) /= 0.0_dp)) then
                    ray%derivatives_valid = .false.
                end if
            end if
            ray%max_hamiltonian_error = max(ray%max_hamiltonian_error, abs(h))
            ray%max_index = max(ray%max_index, norm2(ray%trace%y(4:, k)))
        end do
        if (ray%max_index > 1.0_dp + max(1.0e-8_dp, 20*problem%rtol)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray fast branch requires |N| <= 1')
        end if
        if (ray%max_hamiltonian_error > max(1.0e-7_dp, 100*problem%rtol)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray dispersion residual exceeds tolerance')
        end if
        ray%trace%status = status
        ray%valid = status%code == 0
    contains
        subroutine rhs(t, y, dydt, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: a(6, 6), p(6, medium%n_params() + 1), hh
            logical :: ok
            call ray_rhs(medium, medium%n_params(), frequency, mode, y, &
                dydt, a, p, ok, hh)
            if (.not. ok) failed = .true.
        end subroutine
    end subroutine

    subroutine trace_cold_rays(medium, frequencies, modes, initial, &
            durations, &
            rays, statuses, status, rtol, atol, hmax)
        class(cold_medium_t), intent(in) :: medium
        real(dp), intent(in) :: frequencies(:), initial(:, :), durations(:)
        integer, intent(in) :: modes(:)
        type(cold_ray_t), allocatable, intent(out) :: rays(:)
        type(fortnum_status_t), allocatable, intent(out) :: statuses(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: rtol, atol, hmax
        integer :: i, n
        n = size(frequencies)
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'inconsistent cold-ray batch shapes')
        if (size(initial, 1) /= 6 .or. size(initial, 2) /= n) return
        if (size(modes) /= n .or. size(durations) /= n) return
        allocate (rays(n), statuses(n))
        !$omp parallel do default(shared) private(i) schedule(dynamic, 1)
        do i = 1, n
            call trace_cold_ray(medium, frequencies(i), modes(i), initial(:, i), &
                durations(i), rays(i), statuses(i), rtol, atol, hmax)
        end do
        !$omp end parallel do
        call status_set(status, 0, '')
        do i = 1, n
            if (statuses(i)%code == 0) cycle
            status = statuses(i)
            return
        end do
    end subroutine

    subroutine cold_ray_jvp(medium, ray, initial_dot, parameter_dot, &
            frequency_dot, final_dot, status)
        class(cold_medium_t), intent(in) :: medium
        type(cold_ray_t), intent(in) :: ray
        real(dp), intent(in) :: initial_dot(6), parameter_dot(:), frequency_dot
        real(dp), intent(out) :: final_dot(6)
        type(fortnum_status_t), intent(out) :: status
        type(ode_problem_t) :: problem
        real(dp), allocatable :: result(:)
        logical :: failed
        final_dot = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'invalid cold-ray tangent input')
        if (size(parameter_dot) /= medium%n_params()) return
        if (.not. ray%valid) return
        if (.not. ray%derivatives_valid) return
        if (ray%trace%status%code /= 0) return
        if (.not. allocated(ray%trace%y)) return
        problem%rhs => rhs
        failed = .false.
        call ode_integrate_jvp(problem, tangent, initial_dot, ray%trace, &
            result, status)
        if (failed) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray tangent encountered singular medium')
        end if
        if (status%code == 0) final_dot = result
    contains
        subroutine rhs(t, y, dydt, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: a(6, 6), p(6, medium%n_params() + 1), h
            logical :: ok
            call ray_rhs(medium, medium%n_params(), ray%frequency, ray%mode, &
                y, dydt, a, p, ok, h)
            if (.not. ok) failed = .true.
        end subroutine
        subroutine tangent(t, y, s, ds, ctx)
            real(dp), intent(in) :: t, y(:), s(:)
            real(dp), intent(out) :: ds(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: f(6), a(6, 6), p(6, medium%n_params() + 1), h
            logical :: ok
            call ray_rhs(medium, medium%n_params(), ray%frequency, ray%mode, &
                y, f, a, p, ok, h)
            ds = matmul(a, s) + matmul(p, [parameter_dot, frequency_dot])
            if (.not. ok) failed = .true.
        end subroutine
    end subroutine

    subroutine cold_ray_vjp(medium, ray, final_bar, initial_bar, &
            parameter_bar, frequency_bar, status)
        class(cold_medium_t), intent(in) :: medium
        type(cold_ray_t), intent(in) :: ray
        real(dp), intent(in) :: final_bar(6)
        real(dp), intent(out) :: initial_bar(6), parameter_bar(:), frequency_bar
        type(fortnum_status_t), intent(out) :: status
        type(ode_problem_t) :: problem
        real(dp), allocatable :: result(:), pbar(:)
        logical :: failed
        initial_bar = 0.0_dp
        parameter_bar = 0.0_dp
        frequency_bar = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'invalid cold-ray adjoint input')
        if (size(parameter_bar) /= medium%n_params()) return
        if (.not. ray%valid) return
        if (.not. ray%derivatives_valid) return
        if (ray%trace%status%code /= 0) return
        if (.not. allocated(ray%trace%y)) return
        problem%rhs => rhs
        failed = .false.
        call ode_integrate_parameter_vjp(problem, adjoint, parameters, &
            final_bar, &
            medium%n_params() + 1, ray%trace, result, pbar, status)
        if (failed) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'cold-ray adjoint encountered singular medium')
        end if
        if (status%code /= 0) return
        initial_bar = result
        parameter_bar = pbar(:medium%n_params())
        frequency_bar = pbar(medium%n_params() + 1)
    contains
        subroutine rhs(t, y, dydt, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: a(6, 6), p(6, medium%n_params() + 1), h
            logical :: ok
            call ray_rhs(medium, medium%n_params(), ray%frequency, ray%mode, &
                y, dydt, a, p, ok, h)
            if (.not. ok) failed = .true.
        end subroutine
        subroutine adjoint(t, y, s, ds, ctx)
            real(dp), intent(in) :: t, y(:), s(:)
            real(dp), intent(out) :: ds(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: f(6), a(6, 6), p(6, medium%n_params() + 1), h
            logical :: ok
            call ray_rhs(medium, medium%n_params(), ray%frequency, ray%mode, &
                y, f, a, p, ok, h)
            ds = matmul(transpose(a), s)
            if (.not. ok) failed = .true.
        end subroutine
        subroutine parameters(t, y, cotangent, gradient, ctx)
            real(dp), intent(in) :: t, y(:), cotangent(:)
            real(dp), intent(inout) :: gradient(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: f(6), a(6, 6), p(6, medium%n_params() + 1), h
            logical :: ok
            call ray_rhs(medium, medium%n_params(), ray%frequency, ray%mode, &
                y, f, a, p, ok, h)
            gradient = gradient + matmul(transpose(p), cotangent)
            if (.not. ok) failed = .true.
        end subroutine
    end subroutine
end module
