module tiago_state_sampler
    !! Provider side of the profile-diagnostic contract: the flux label s and
    !! |B| of a real_space_state_t at the diagnostic sample points, their
    !! spatial derivatives, and their tangents and adjoints with respect to the
    !! provider parameters at fixed points.
    !!
    !! Points are (R, phi, Z). The state is inverted by tiago_real_space_state
    !! (Newton) and differentiated with the implicit-function theorem; s = rho^2.
    !! The kinetic rows Q_NE .. Q_ZEFF are left to the caller (kinetic profiles
    !! evaluated at s). A point is inside when it is located with rho <= 1;
    !! points outside the state's chart keep s and |B| of the last Newton
    !! iterate and zero derivatives.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok, &
        located_jvp, located_vjp
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, N_QUANTITIES
    implicit none
    private

    public :: sample_state, sample_state_jvp, sample_state_vjp

contains

    subroutine sample_state(state, x, sample, loc)
        class(real_space_state_t), intent(in) :: state
        real(dp), intent(in) :: x(:, :) !! (3, n): R, phi, Z
        type(plasma_sample_t), intent(out) :: sample
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp) :: m(3, 3), dudx(3, 3)
        integer :: k, n

        n = size(x, 2)
        sample%n = n
        sample%x = x
        allocate(sample%inside(n), sample%q(N_QUANTITIES, n), &
            sample%dq_dx(3, N_QUANTITIES, n), loc(n))
        sample%q = 0.0_dp
        sample%dq_dx = 0.0_dp
        !$omp parallel do private(m, dudx) schedule(dynamic, 4)
        do k = 1, n
            call state%locate(cartesian(x(:, k)), loc(k))
            sample%q(Q_S, k) = loc(k)%u(1)**2
            sample%q(Q_B, k) = loc(k)%b
            sample%inside(k) = loc(k)%status == locate_ok .and. loc(k)%u(1) <= 1.0_dp
            if (loc(k)%status /= locate_ok) cycle
            m = cylinder_jacobian(x(:, k))
            dudx = matmul(loc(k)%jinv, m)
            sample%dq_dx(:, Q_S, k) = 2.0_dp * loc(k)%u(1) * dudx(1, :)
            sample%dq_dx(:, Q_B, k) = matmul(loc(k)%dbdu, dudx)
        end do
        !$omp end parallel do
    end subroutine sample_state

    subroutine sample_state_jvp(state, loc, p_dot, q_dot)
        !! Rows Q_S and Q_B of q_dot for a provider tangent p_dot.
        class(real_space_state_t), intent(in) :: state
        type(located_t), intent(in) :: loc(:)
        real(dp), intent(in) :: p_dot(:)
        real(dp), intent(inout) :: q_dot(:, :)
        real(dp) :: du(3), dbx
        integer :: k

        !$omp parallel do private(du, dbx)
        do k = 1, size(loc)
            q_dot(Q_S, k) = 0.0_dp
            q_dot(Q_B, k) = 0.0_dp
            if (loc(k)%status /= locate_ok) cycle
            call located_jvp(state, loc(k), p_dot, du, dbx)
            q_dot(Q_S, k) = 2.0_dp * loc(k)%u(1) * du(1)
            q_dot(Q_B, k) = dbx
        end do
        !$omp end parallel do
    end subroutine sample_state_jvp

    subroutine sample_state_vjp(state, loc, q_bar, p_bar)
        !! p_bar += adjoint of sample_state_jvp for rows Q_S and Q_B of q_bar.
        class(real_space_state_t), intent(in) :: state
        type(located_t), intent(in) :: loc(:)
        real(dp), intent(in) :: q_bar(:, :)
        real(dp), intent(inout) :: p_bar(:)
        real(dp), allocatable :: acc(:)
        integer :: k

        allocate(acc(size(p_bar)))
        acc = 0.0_dp
        !$omp parallel do reduction(+:acc) schedule(dynamic, 4)
        do k = 1, size(loc)
            if (loc(k)%status /= locate_ok) cycle
            call located_vjp(state, loc(k), [2.0_dp * loc(k)%u(1) * q_bar(Q_S, k), &
                0.0_dp, 0.0_dp], q_bar(Q_B, k), acc)
        end do
        !$omp end parallel do
        p_bar = p_bar + acc
    end subroutine sample_state_vjp

    pure function cartesian(rpz) result(x)
        real(dp), intent(in) :: rpz(3)
        real(dp) :: x(3)
        x = [rpz(1) * cos(rpz(2)), rpz(1) * sin(rpz(2)), rpz(3)]
    end function cartesian

    pure function cylinder_jacobian(rpz) result(m)
        !! d (x, y, z) / d (R, phi, Z)
        real(dp), intent(in) :: rpz(3)
        real(dp) :: m(3, 3)
        m = 0.0_dp
        m(1, 1) = cos(rpz(2))
        m(2, 1) = sin(rpz(2))
        m(1, 2) = -rpz(1) * sin(rpz(2))
        m(2, 2) = rpz(1) * cos(rpz(2))
        m(3, 3) = 1.0_dp
    end function cylinder_jacobian
end module tiago_state_sampler
