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
    use tiago_real_space_state, only: real_space_state_t, located_t, vector_jets_t, &
        locate_ok, locate_singular, &
        located_jvp, located_vjp
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, N_QUANTITIES
    use tiago_kinetic_profiles, only: kinetic_profiles_t
    use tiago_medium_pullback, only: medium_pullback
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    public :: sample_medium, medium_sample_jvp, medium_sample_vjp
    public :: cylinder_jacobian
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

    subroutine sample_medium(state, profiles, z, x, q, dq, d2q, qp, dqp, status)
        class(real_space_state_t), intent(in) :: state
        type(kinetic_profiles_t), intent(in) :: profiles
        real(dp), intent(in) :: z(:), x(3)
        real(dp), intent(out) :: q(4), dq(3, 4), d2q(3, 3, 4)
        !! Without qp and dqp the parameter derivatives are skipped.
        real(dp), intent(out), optional :: qp(:, :), dqp(:, :, :)
        integer, intent(out) :: status
        type(located_t) :: loc
        type(vector_jets_t) :: jet
        real(dp) :: gu(3, 4), hu(3, 3, 4), ne, nr, nrr
        real(dp) :: nep(profiles%n_latent), nrp(profiles%n_latent)
        real(dp) :: xp(3), xup(3, 3), fp, gup(3), unused, unused3(3)
        real(dp) :: unused33(3, 3), zero3(3), zero33(3, 3)
        integer :: np, nz, i, j
        logical :: columns
        columns = present(qp) .and. present(dqp)
        q = ieee_value(0.0_dp, ieee_quiet_nan)
        dq = q(1)
        d2q = q(1)
        np = state%n_params()
        nz = profiles%n_latent
        status = locate_singular
        if (size(z) /= nz) return
        if (columns) then
            qp = q(1)
            dqp = q(1)
            if (size(qp, 1) /= 4 .or. size(qp, 2) /= np + nz) return
            if (size(dqp, 1) /= 3 .or. size(dqp, 2) /= 4 .or. size(dqp, &
                3) /= np + nz) return
        end if
        call state%locate(x, loc)
        status = loc%status
        if (status /= locate_ok) return
        call state%vector_jets(loc%u, jet, status, columns)
        if (status /= locate_ok) return
        call profiles%density_jet(z, loc%u(1), ne, nr, nrr, nep, nrp)
        q(:3) = jet%b
        q(4) = ne
        gu = 0
        hu = 0
        gu(:, :3) = transpose(jet%bu)
        do i = 1, 3
            hu(:, :, i) = jet%buu(i, :, :)
        end do
        gu(1, 4) = nr
        hu(1, 1, 4) = nrr
        zero3 = 0
        zero33 = 0
        do i = 1, 4
            call medium_pullback(loc%jinv, jet%xuu, gu(:, i), hu(:, :, i), &
                                 zero3, zero33, 0.0_dp, zero3, dq(:, i), d2q(:, :, i), &
                                 unused, unused3)
            if (.not. columns) cycle
            do j = 1, np + nz
                xp = 0
                xup = 0
                fp = 0
                gup = 0
                if (j <= np) then
                    xp = jet%xp(:, j)
                    xup = jet%xup(:, :, j)
                    if (i <= 3) then
                        fp = jet%bp(i, j)
                        gup = jet%bup(i, :, j)
                    end if
                else if (i == 4) then
                    fp = nep(j - np)
                    gup(1) = nrp(j - np)
                end if
                call medium_pullback(loc%jinv, jet%xuu, gu(:, i), hu(:, :, i), &
                                     xp, xup, fp, gup, unused3, unused33, qp(i, j), &
                                     dqp(:, i, j))
            end do
        end do
    end subroutine sample_medium

    pure subroutine medium_sample_jvp(dq, d2q, qp, dqp, xdot, pdot, qdot, dqdot)
        real(dp), intent(in) :: dq(3, 4), d2q(3, 3, 4), qp(:, :), dqp(:, :, :)
        real(dp), intent(in) :: xdot(3), pdot(:)
        real(dp), intent(out) :: qdot(4), dqdot(3, 4)
        integer :: i
        qdot = matmul(transpose(dq), xdot) + matmul(qp, pdot)
        do i = 1, 4
            dqdot(:, i) = matmul(d2q(:, :, i), xdot) + matmul(dqp(:, i, :), pdot)
        end do
    end subroutine

    pure subroutine medium_sample_vjp(dq, d2q, qp, dqp, qbar, dqbar, xbar, pbar)
        real(dp), intent(in) :: dq(3, 4), d2q(3, 3, 4), qp(:, :), dqp(:, :, :)
        real(dp), intent(in) :: qbar(4), dqbar(3, 4)
        real(dp), intent(out) :: xbar(3), pbar(:)
        integer :: i
        xbar = matmul(dq, qbar)
        pbar = matmul(transpose(qp), qbar)
        do i = 1, 4
            xbar = xbar + matmul(transpose(d2q(:, :, i)), dqbar(:, i))
            pbar = pbar + matmul(transpose(dqp(:, i, :)), dqbar(:, i))
        end do
    end subroutine

end module tiago_state_sampler
