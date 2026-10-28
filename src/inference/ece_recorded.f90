module tiago_ece_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use tiago_ece_radiometer, only: ece_t
    use tiago_ece_opacity_coordinates, only: ece_opacity_coordinates
    use tiago_ece_ray, only: ECE_X_MODE
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_ray, only: cold_ray_t, trace_cold_ray, cold_ray_jvp, &
        cold_ray_vjp, &
        cold_hamiltonian
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_B, &
        Q_NE, Q_TE, Q_S, Q_TI, Q_ZEFF
    use tiago_real_space_state, only: located_t, locate_ok, located_vjp
    use tiago_kinetic_profiles, only: N_PROFILES, P_TE, P_TI, P_ZEFF
    use tiago_state_sampler, only: sample_state, cylinder_jacobian
    use tiago_sight_geometry, only: to_rpz
    use tiago_ece_field_norm, only: ece_field_norm
    use tiago_ece_launch_scale, only: ece_launch_scale
    use tiago_ece_parallel_fraction, only: ece_parallel_fraction
    use tiago_ece_segment_length, only: ece_segment_length
    implicit none
    private
    integer, parameter :: rows(3) = [Q_B, Q_NE, Q_TE]

    type :: ece_path_t
        type(cold_ray_t), allocatable :: segment(:)
        real(dp), allocatable :: launch_row(:, :)
    end type

    type, public :: recorded_ece_t
        private
        type(ece_t) :: diagnostic
        type(state_medium_t) :: medium
        type(ece_path_t), allocatable :: paths(:)
        type(plasma_sample_t) :: sample
        real(dp), allocatable :: eta(:), qx(:, :, :), qp(:, :, :), length_row(:, :)
        logical, public :: valid = .false.
        real(dp), allocatable, public :: voltage(:), variance(:)
    contains
        procedure :: jvp => recorded_jvp
        procedure :: vjp => recorded_vjp
        procedure :: sampled => recorded_sampled
    end type
    public :: record_ece
contains
    subroutine record_ece(diagnostic, medium, eta, tape, status)
        type(ece_t), intent(in) :: diagnostic
        type(state_medium_t), intent(in) :: medium
        real(dp), intent(in) :: eta(:)
        type(recorded_ece_t), intent(out) :: tape
        type(fortnum_status_t), intent(out) :: status
        type(fortnum_status_t) :: codes(size(diagnostic%channels))
        integer :: n, np, c, code
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'ECE: invalid composition input')
        if (.not. allocated(medium%state)) return
        if (.not. allocated(medium%z)) return
        if (size(eta) /= diagnostic%n_nuisance()) return
        if (.not. all(ieee_is_finite(eta))) return
        tape%diagnostic = diagnostic
        tape%medium = medium
        tape%eta = eta
        n = diagnostic%first(size(diagnostic%first)) - 1
        np = medium%n_params()
        allocate (tape%paths(size(diagnostic%channels)), tape%qx(3, 3, n), &
            tape%qp(3, np, n), tape%length_row(6, n))
        tape%sample%n = n
        allocate (tape%sample%x(3, n), tape%sample%inside(n), &
            tape%sample%q(N_QUANTITIES, n), tape%sample%dq_dx(3, &
            N_QUANTITIES, n))
        tape%sample%q = 0
        tape%sample%dq_dx = 0
        if (allocated(tape%diagnostic%points)) deallocate (tape%diagnostic%points, &
            tape%diagnostic%cell_length, &
            tape%diagnostic%parallel_cos2)
        allocate (tape%diagnostic%points(3, n), tape%diagnostic%cell_length(n), &
            tape%diagnostic%parallel_cos2(n))
        !$omp parallel do default(shared) private(c) schedule(dynamic, 1)
        do c = 1, size(diagnostic%channels)
            call build_path(tape, c, codes(c))
        end do
        !$omp end parallel do
        do c = 1, size(codes)
            if (codes(c)%code == 0) cycle
            status = codes(c)
            return
        end do
        tape%diagnostic%points = tape%sample%x
        allocate (tape%voltage(diagnostic%n_signals()), &
            tape%variance(diagnostic%n_signals()))
        call tape%diagnostic%record(tape%sample, eta, tape%voltage, &
            tape%variance, code)
        if (code == 3) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'ECE: emission ray exceeds perpendicular-opacity angle bound')
        else if (code /= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'ECE: unsupported opacity regime or invalid wall/receiver parameters')
        else
            tape%valid = .true.
            call status_set(status, 0, '')
        end if
    end subroutine

    subroutine build_path(tape, c, status)
        type(recorded_ece_t), intent(inout) :: tape
        integer, intent(in) :: c
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: y(6), previous(6), direction(3), mid(3), q(4), dq(3, 4), &
            d2q(3, 3, 4), qp(4, tape%medium%n_params()), &
            dqp(3, 4, tape%medium%n_params()), h, dh(8), d2h(8, 8), scale, &
            sr(1), &
            cos2, unused(6), duration, coordinates(4), coordinate_rows(4, 4)
        integer :: j, k, last
        associate (ch => tape%diagnostic%channels(c), path => tape%paths(c))
            direction = ch%direction/norm2(ch%direction)
            if (ch%mode == ECE_X_MODE .and. ch%harmonic == 1) then
                k = tape%diagnostic%first(c)
                call sample_cell(tape, k, ch%origin, q, status)
                if (status%code /= 0) return
                call ece_opacity_coordinates(tape%sample%q(Q_TE, k), q(4), &
                    tape%sample%q(Q_B, k), ch%frequency, coordinates, &
                    coordinate_rows)
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'ECE: first-harmonic X mode requires the tenuous density regime')
                if (coordinates(3) > .2_dp*coordinates(1)) return
            end if
            call tape%medium%sample(ch%origin, q, dq, d2q, qp, dqp)
            call cold_hamiltonian(direction, q(:3), q(4), ch%frequency, ch%mode, &
                h, dh, d2h)
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'ECE: nonpropagating launch')
            if (.not. ieee_is_finite(h)) return
            if (1 - 2*h <= 0) return
            call ece_launch_scale(h, scale, sr)
            previous = [ch%origin, scale*direction]
            allocate (path%launch_row(6, tape%medium%n_params()), &
                path%segment(ch%nscan))
            path%launch_row = 0
            do j = 1, 3
                path%launch_row(j + 3, :) = direction(j)*sr(1)*matmul(dh(4:7), qp)
            end do
            duration = ch%length/ch%nscan
            do j = 1, ch%nscan
                k = tape%diagnostic%first(c) + j - 1
                call trace_cold_ray(tape%medium, ch%frequency, ch%mode, previous, &
                    duration, path%segment(j), status, &
                    rtol=1.0e-10_dp, &
                    atol=1.0e-12_dp, hmax=duration)
                if (status%code /= 0) return
                if (.not. path%segment(j)%derivatives_valid) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        'ECE: ray derivatives are unavailable')
                    return
                end if
                last = size(path%segment(j)%trace%t)
                y = path%segment(j)%trace%y(:, last)
                mid = (previous(:3) + y(:3))/2
                call ece_segment_length(previous(1), previous(2), previous(3), &
                    y(1), y(2), y(3), tape%diagnostic%cell_length(k), &
                    tape%length_row(:, k))
                call sample_cell(tape, k, mid, q, status)
                if (status%code /= 0) return
                if (ch%mode == ECE_X_MODE .and. ch%harmonic == 1) then
                    call ece_opacity_coordinates(tape%sample%q(Q_TE, k), q(4), &
                        tape%sample%q(Q_B, k), ch%frequency, coordinates, &
                        coordinate_rows)
                    if (coordinates(3) > .2_dp*coordinates(1)) then
                        call status_set(status, FORTNUM_DOMAIN_ERROR, &
                       'ECE: first-harmonic X mode requires the tenuous density regime')
                        return
                    end if
                end if
                direction = (previous(4:) + y(4:))/2
                call ece_parallel_fraction(direction(1), direction(2), &
                    direction(3), &
                    q(1), q(2), q(3), cos2, unused)
                tape%diagnostic%parallel_cos2(k) = cos2
                direction = y(:3) - previous(:3)
                call ece_parallel_fraction(direction(1), direction(2), &
                    direction(3), &
                    q(1), q(2), q(3), cos2, unused)
                tape%diagnostic%parallel_cos2(k) = max(cos2, &
                    tape%diagnostic%parallel_cos2(k))
                previous = y
            end do
        end associate
    end subroutine

    subroutine sample_cell(tape, k, x, q, status)
        type(recorded_ece_t), intent(inout) :: tape
        integer, intent(in) :: k
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: q(4)
        type(fortnum_status_t), intent(out) :: status
        type(located_t) :: loc
        type(located_t), allocatable :: locations(:)
        type(plasma_sample_t) :: sampled
        real(dp) :: dq(3, 4), d2q(3, 3, 4), qp(4, tape%medium%n_params()), &
            dqp(3, 4, tape%medium%n_params()), b, br(3), &
            f(N_PROFILES, 1), df(N_PROFILES, 1), fb(N_PROFILES, 1), sb(1), &
            zb(tape%medium%profiles%n_latent), &
            pb(tape%medium%state%n_params()), m(3, 3)
        integer :: ns
        call tape%medium%sample(x, q, dq, d2q, qp, dqp)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'ECE: invalid medium cell')
        if (.not. all(ieee_is_finite(q))) return
        if (norm2(q(:3)) <= 0) return
        call sample_state(tape%medium%state, reshape(to_rpz(x), [3, 1]), sampled, &
            locations)
        loc = locations(1)
        if (loc%status /= locate_ok) return
        call ece_field_norm(q(1), q(2), q(3), b, br)
        tape%sample%x(:, k) = sampled%x(:, 1)
        tape%sample%q(:, k) = sampled%q(:, 1)
        tape%sample%dq_dx(:, :, k) = sampled%dq_dx(:, :, 1)
        tape%sample%inside(k) = loc%u(1) <= 1
        call tape%medium%profiles%evaluate(tape%medium%z, [loc%u(1)**2], f, df)
        tape%sample%q(Q_B, k) = b
        tape%sample%q(Q_NE, k) = q(4)
        tape%sample%q(Q_TE, k) = f(P_TE, 1)
        tape%sample%q(Q_TI, k) = f(P_TI, 1)
        tape%sample%q(Q_ZEFF, k) = f(P_ZEFF, 1)
        tape%qx(:, 1, k) = matmul(dq(:, :3), br)
        tape%qx(:, 2, k) = dq(:, 4)
        tape%qp(1, :, k) = matmul(br, qp(:3, :))
        tape%qp(2, :, k) = qp(4, :)
        fb = 0
        fb(P_TE, 1) = 1
        zb = 0
        call tape%medium%profiles%vjp(tape%medium%z, [loc%u(1)**2], fb, zb, sb)
        pb = 0
        call located_vjp(tape%medium%state, loc, [2*loc%u(1)*sb(1), 0.0_dp, &
            0.0_dp], 0.0_dp, pb)
        ns = size(pb)
        tape%qp(3, :ns, k) = pb
        tape%qp(3, ns + 1:, k) = zb
        tape%qx(:, 3, k) = 2*loc%u(1)*df(P_TE, 1)*loc%jinv(1, :)
        m = cylinder_jacobian(tape%sample%x(:, k))
        tape%sample%dq_dx(:, rows, k) = matmul(transpose(m), tape%qx(:, :, k))
        tape%sample%dq_dx(:, Q_TI, k) = df(P_TI, 1)*tape%sample%dq_dx(:, Q_S, k)
        tape%sample%dq_dx(:, Q_ZEFF, k) = df(P_ZEFF, 1)*tape%sample%dq_dx(:, Q_S, k)
        call status_set(status, 0, '')
    end subroutine

    subroutine recorded_sampled(self, diagnostic, sample)
        class(recorded_ece_t), intent(in) :: self
        type(ece_t), intent(out) :: diagnostic
        type(plasma_sample_t), intent(out) :: sample
        if (.not. self%valid) error stop 'ECE: invalid composition tape'
        diagnostic = self%diagnostic
        sample = self%sample
    end subroutine

    subroutine recorded_jvp(self, parameter_dot, eta_dot, voltage_dot, &
        variance_dot, status)
        class(recorded_ece_t), intent(in) :: self
        real(dp), intent(in) :: parameter_dot(:), eta_dot(:)
        real(dp), intent(out) :: voltage_dot(:), variance_dot(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: qdot(:, :), ld(:)
        type(fortnum_status_t), allocatable :: codes(:)
        integer :: c, code
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'ECE: invalid composition tangent')
        if (.not. self%valid) return
        if (size(parameter_dot) /= self%medium%n_params()) return
        if (size(eta_dot) /= size(self%eta)) return
        allocate (qdot(N_QUANTITIES, self%sample%n), ld(self%sample%n), &
            codes(size(self%paths)))
        qdot = 0
        !$omp parallel do default(shared) private(c) schedule(dynamic, 1)
        do c = 1, size(self%paths)
            call path_jvp(self, c, parameter_dot, qdot, ld, codes(c))
        end do
        !$omp end parallel do
        do c = 1, size(codes)
            if (codes(c)%code == 0) cycle
            status = codes(c)
            return
        end do
        call self%diagnostic%record_jvp(self%sample, self%eta, qdot, eta_dot, &
            ld, voltage_dot, variance_dot, code)
        if (code == 0) call status_set(status, 0, '')
    end subroutine

    subroutine path_jvp(self, c, pd, qdot, ld, status)
        class(recorded_ece_t), intent(in) :: self
        integer, intent(in) :: c
        real(dp), intent(in) :: pd(:)
        real(dp), intent(inout) :: qdot(:, :), ld(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: previous(6), final(6), xdot(3)
        integer :: j, k
        previous = matmul(self%paths(c)%launch_row, pd)
        do j = 1, size(self%paths(c)%segment)
            k = self%diagnostic%first(c) + j - 1
            call cold_ray_jvp(self%medium, self%paths(c)%segment(j), previous, pd, &
                0.0_dp, final, status)
            if (status%code /= 0) return
            xdot = (previous(:3) + final(:3))/2
            qdot(rows, k) = matmul(self%qp(:, :, k), pd) + &
                matmul(transpose(self%qx(:, :, k)), xdot)
            ld(k) = dot_product(self%length_row(:, k), [previous(:3), final(:3)])
            previous = final
        end do
    end subroutine

    subroutine recorded_vjp(self, voltage_bar, variance_bar, parameter_bar, &
        eta_bar, status)
        class(recorded_ece_t), intent(in) :: self
        real(dp), intent(in) :: voltage_bar(:), variance_bar(:)
        real(dp), intent(out) :: parameter_bar(:), eta_bar(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: qb(:, :), lb(:), partial(:, :)
        type(fortnum_status_t), allocatable :: codes(:)
        integer :: c, code
        parameter_bar = 0
        eta_bar = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'ECE: invalid composition adjoint')
        if (.not. self%valid) return
        if (size(parameter_bar) /= self%medium%n_params()) return
        if (size(eta_bar) /= size(self%eta)) return
        allocate (qb(N_QUANTITIES, self%sample%n), lb(self%sample%n), &
            partial(size(parameter_bar), size(self%paths)), &
            codes(size(self%paths)))
        call self%diagnostic%record_vjp(self%sample, self%eta, voltage_bar, &
            variance_bar, &
            qb, eta_bar, lb, code)
        if (code /= 0) return
        !$omp parallel do default(shared) private(c) schedule(dynamic, 1)
        do c = 1, size(self%paths)
            call path_vjp(self, c, qb, lb, partial(:, c), codes(c))
        end do
        !$omp end parallel do
        do c = 1, size(codes)
            if (codes(c)%code == 0) cycle
            status = codes(c)
            return
        end do
        parameter_bar = sum(partial, dim=2)
        call status_set(status, 0, '')
    end subroutine

    subroutine path_vjp(self, c, qb, lb, pb, status)
        class(recorded_ece_t), intent(in) :: self
        integer, intent(in) :: c
        real(dp), intent(in) :: qb(:, :), lb(:)
        real(dp), intent(out) :: pb(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: carry(6), initial(6), xb(3), direct(size(pb)), fb
        integer :: j, k
        carry = 0
        pb = 0
        do j = size(self%paths(c)%segment), 1, -1
            k = self%diagnostic%first(c) + j - 1
            xb = matmul(self%qx(:, :, k), qb(rows, k))/2
            carry(:3) = carry(:3) + xb + self%length_row(4:, k)*lb(k)
            call cold_ray_vjp(self%medium, self%paths(c)%segment(j), carry, &
                initial, &
                direct, fb, status)
            if (status%code /= 0) return
            pb = pb + direct + matmul(transpose(self%qp(:, :, k)), qb(rows, k))
            initial(:3) = initial(:3) + xb + self%length_row(:3, k)*lb(k)
            carry = initial
        end do
        pb = pb + matmul(transpose(self%paths(c)%launch_row), carry)
    end subroutine
end module
