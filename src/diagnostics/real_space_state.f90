module tiago_real_space_state
    !! Provider-neutral view of a plasma state in real space for profile
    !! diagnostics (Thomson, ECE, line integrals).
    !!
    !! A provider supplies a smooth map from labels u = (rho, theta, phi) to
    !! Cartesian x, the magnetic-field strength |B| on the same labels, and the
    !! derivatives of both with respect to its own parameters p at fixed u.
    !! rho is the flux label of the kinetic profiles (rho = 1 on the boundary).
    !! A VMEC nested state is one such map; a KIN6D real-space state can supply
    !! another (rho from its profile label, theta and phi any smooth chart).
    !!
    !! TIAGO inverts x = X(u; p) by Newton's method and differentiates the
    !! inverse with the implicit-function theorem,
    !!     du = J^{-1} (dx - dX/dp dp),       J = dX/du,
    !! so that a quantity q(u; p) at a fixed real-space point changes by
    !!     dq = dq/dp dp + grad_u q . du.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    integer, parameter, public :: locate_ok = 0, locate_outside = 1, &
        locate_no_convergence = 2, locate_singular = 3, locate_unsupported = 4

    type, abstract, public :: real_space_state_t
        !> Newton iterates beyond this rho count as outside the state's chart.
        real(dp) :: rho_chart = 1.5_dp
    contains
        procedure(n_params_i), deferred :: n_params
        procedure(get_params_i), deferred :: get_params
        procedure(set_params_i), deferred :: set_params
        procedure(eval_i), deferred :: eval
        procedure(param_jvp_i), deferred :: param_jvp
        procedure(param_vjp_i), deferred :: param_vjp
        procedure :: initial_guess => default_initial_guess
        procedure :: locate => state_locate
        procedure :: vector_jets => unsupported_vector_jets
    end type real_space_state_t

    type, public :: located_t
        !! A real-space point located in the state's labels.
        real(dp) :: x(3) = 0.0_dp
        real(dp) :: u(3) = 0.0_dp
        real(dp) :: jinv(3, 3) = 0.0_dp !! (dX/du)^{-1}
        real(dp) :: b = 0.0_dp, dbdu(3) = 0.0_dp
        integer :: status = locate_outside
    end type located_t

    type, public :: vector_jets_t
        real(dp) :: x(3), xu(3, 3), xuu(3, 3, 3)
        real(dp) :: b(3), bu(3, 3), buu(3, 3, 3)
        real(dp), allocatable :: xp(:, :), xup(:, :, :), bp(:, :), bup(:, :, :)
    end type

    abstract interface
        pure integer function n_params_i(self)
            import :: real_space_state_t
            class(real_space_state_t), intent(in) :: self
        end function n_params_i

        pure subroutine get_params_i(self, p)
            import :: real_space_state_t, dp
            class(real_space_state_t), intent(in) :: self
            real(dp), intent(out) :: p(:)
        end subroutine get_params_i

        subroutine set_params_i(self, p)
            import :: real_space_state_t, dp
            class(real_space_state_t), intent(inout) :: self
            real(dp), intent(in) :: p(:)
        end subroutine set_params_i

        pure subroutine eval_i(self, u, x, dxdu, b, dbdu)
            import :: real_space_state_t, dp
            class(real_space_state_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3), dxdu(3, 3), b, dbdu(3)
        end subroutine eval_i

        pure subroutine param_jvp_i(self, u, dp_, dx, db)
            !! dX/dp dp and d|B|/dp dp at fixed labels.
            import :: real_space_state_t, dp
            class(real_space_state_t), intent(in) :: self
            real(dp), intent(in) :: u(3), dp_(:)
            real(dp), intent(out) :: dx(3), db
        end subroutine param_jvp_i

        pure subroutine param_vjp_i(self, u, xbar, bbar, pbar)
            !! pbar += (dX/dp)^T xbar + (d|B|/dp)^T bbar at fixed labels.
            import :: real_space_state_t, dp
            class(real_space_state_t), intent(in) :: self
            real(dp), intent(in) :: u(3), xbar(3), bbar
            real(dp), intent(inout) :: pbar(:)
        end subroutine param_vjp_i
    end interface

    public :: located_jvp, located_vjp, invert3

contains

    subroutine unsupported_vector_jets(self, u, jet, status, parameters)
        !! parameters = .false. leaves the parameter columns (xp, xup, bp, bup)
        !! empty in every override.
        class(real_space_state_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        type(vector_jets_t), intent(out) :: jet
        integer, intent(out) :: status
        logical, intent(in), optional :: parameters
        status = locate_unsupported
        call self%eval(u, jet%x, jet%xu, jet%b(1), jet%bu(1, :))
    end subroutine

    pure subroutine default_initial_guess(self, x, u)
        !! Closest point of a coarse label grid at phi = atan2(y, x).
        class(real_space_state_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        real(dp) :: v(3), xv(3), jac(3, 3), b, db(3), best, d
        integer :: i, j

        best = huge(1.0_dp)
        u = [0.5_dp, 0.0_dp, atan2(x(2), x(1))]
        do i = 1, 4
            do j = 0, 15
                v = [0.25_dp * i, 2.0_dp * pi * j / 16.0_dp, atan2(x(2), x(1))]
                call self%eval(v, xv, jac, b, db)
                d = sum((xv - x)**2)
                if (d < best) then
                    best = d
                    u = v
                end if
            end do
        end do
    end subroutine default_initial_guess

    pure subroutine state_locate(self, x, loc)
        !! Damped Newton for X(u) = x. Status locate_outside when an iterate
        !! leaves the chart (rho > rho_chart); a converged point may still have
        !! rho > 1, which callers treat as outside the plasma.
        class(real_space_state_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        type(located_t), intent(out) :: loc
        real(dp) :: u(3), xu(3), jac(3, 3), jinv(3, 3), b, db(3), r(3), step(3)
        real(dp) :: res, trial(3), xt(3), jt(3, 3), bt, dbt(3), rt, scale, det
        integer :: iter, k

        loc%x = x
        call self%initial_guess(x, u)
        call self%eval(u, xu, jac, b, db)
        r = xu - x
        res = norm2(r)
        scale = max(1.0_dp, norm2(x))
        loc%status = locate_no_convergence
        do iter = 1, 80
            call invert3(jac, jinv, det)
            if (det == 0.0_dp) then
                loc%status = locate_singular
                return
            end if
            if (res <= 1.0e-13_dp * scale) then
                loc%status = locate_ok
                exit
            end if
            step = -matmul(jinv, r)
            do k = 1, 40
                trial = u + step
                if (trial(1) < 0.0_dp) trial(1:2) = [-trial(1), trial(2) + pi]
                call self%eval(trial, xt, jt, bt, dbt)
                rt = norm2(xt - x)
                if (rt < res) exit
                step = 0.5_dp * step
            end do
            if (rt >= res) then
                ! No decrease: accept a stagnated point at roundoff level.
                if (res <= 1.0e-10_dp * scale) loc%status = locate_ok
                exit
            end if
            u = trial
            xu = xt
            jac = jt
            b = bt
            db = dbt
            r = xu - x
            res = rt
            if (u(1) > self%rho_chart) then
                loc%status = locate_outside
                exit
            end if
        end do
        call invert3(jac, jinv, det)
        loc%u = u
        loc%jinv = jinv
        loc%b = b
        loc%dbdu = db
        if (loc%status == locate_ok .and. u(1) > self%rho_chart) loc%status = locate_outside
    end subroutine state_locate

    pure subroutine located_jvp(self, loc, dp_, du, dbx)
        !! Label and |B| changes at the fixed real-space point loc%x.
        class(real_space_state_t), intent(in) :: self
        type(located_t), intent(in) :: loc
        real(dp), intent(in) :: dp_(:)
        real(dp), intent(out) :: du(3), dbx
        real(dp) :: dxp(3), dbp

        call self%param_jvp(loc%u, dp_, dxp, dbp)
        du = -matmul(loc%jinv, dxp)
        dbx = dbp + dot_product(loc%dbdu, du)
    end subroutine located_jvp

    pure subroutine located_vjp(self, loc, ubar, bxbar, pbar)
        !! Adjoint of located_jvp: pbar += (du/dp)^T ubar + (d|B|_x/dp)^T bxbar.
        class(real_space_state_t), intent(in) :: self
        type(located_t), intent(in) :: loc
        real(dp), intent(in) :: ubar(3), bxbar
        real(dp), intent(inout) :: pbar(:)
        real(dp) :: utot(3)

        utot = ubar + bxbar * loc%dbdu
        call self%param_vjp(loc%u, -matmul(transpose(loc%jinv), utot), bxbar, pbar)
    end subroutine located_vjp

    pure subroutine invert3(a, ainv, det)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: ainv(3, 3), det

        ainv(1, 1) = a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)
        ainv(1, 2) = a(1, 3) * a(3, 2) - a(1, 2) * a(3, 3)
        ainv(1, 3) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
        ainv(2, 1) = a(2, 3) * a(3, 1) - a(2, 1) * a(3, 3)
        ainv(2, 2) = a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1)
        ainv(2, 3) = a(1, 3) * a(2, 1) - a(1, 1) * a(2, 3)
        ainv(3, 1) = a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1)
        ainv(3, 2) = a(1, 2) * a(3, 1) - a(1, 1) * a(3, 2)
        ainv(3, 3) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
        det = a(1, 1) * ainv(1, 1) + a(1, 2) * ainv(2, 1) + a(1, 3) * ainv(3, 1)
        if (det /= 0.0_dp) then
            ainv = ainv / det
        else
            ainv = 0.0_dp
        end if
    end subroutine invert3
end module tiago_real_space_state
