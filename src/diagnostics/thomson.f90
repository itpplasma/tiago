module tiago_thomson
    !! Thomson scattering: electron density and temperature of finite
    !! scattering volumes. Channel i averages over the laser segment of length
    !! L_i centred on its volume (Gauss-Legendre points, normalized weights w_k):
    !!     ne_i = sum_k w_k ne(rho_k),
    !!     Te_i = sum_k w_k ne(rho_k) Te(rho_k) / ne_i,
    !! since the scattered power, and therefore the spectral weight of each
    !! part of the volume, is proportional to ne. Points outside the plasma
    !! (rho > 1) contribute nothing; the boundary motion across such a point is
    !! not differentiated, so edge channels should lie inside the boundary.
    !! Outputs y = (ne_1 .. ne_n, Te_1 .. Te_n).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t
    use tiago_label_profile, only: radial_profile_t
    use tiago_los_geometry, only: sample_set_t, gauss_legendre01, locate_samples, &
        samples_jvp, samples_vjp
    implicit none
    private

    type, public :: thomson_channel_t
        character(len=:), allocatable :: label
        real(dp) :: center(3) = 0.0_dp
        real(dp) :: direction(3) = [0.0_dp, 0.0_dp, 1.0_dp] !! laser direction
        real(dp) :: length = 0.0_dp !! scattering length [m]
        integer :: npoints = 1
    end type thomson_channel_t

    public :: thomson_samples, thomson_forward, thomson_jvp, thomson_vjp

contains

    subroutine thomson_samples(state, channels, set)
        class(real_space_state_t), intent(in) :: state
        type(thomson_channel_t), intent(in) :: channels(:)
        type(sample_set_t), intent(out) :: set
        real(dp), allocatable :: xg(:), wg(:)
        real(dp) :: d(3)
        integer :: c, k, m, np

        set%n = 0
        do c = 1, size(channels)
            set%n = set%n + npts(channels(c))
        end do
        allocate(set%x(3, set%n), set%w(set%n), set%owner(set%n))
        m = 0
        do c = 1, size(channels)
            np = npts(channels(c))
            if (np == 1) then
                m = m + 1
                set%x(:, m) = channels(c)%center
                set%w(m) = 1.0_dp
                set%owner(m) = c
                cycle
            end if
            call gauss_legendre01(np, xg, wg)
            d = channels(c)%direction / norm2(channels(c)%direction)
            do k = 1, np
                m = m + 1
                set%x(:, m) = channels(c)%center + (xg(k) - 0.5_dp) * channels(c)%length * d
                set%w(m) = wg(k)
                set%owner(m) = c
            end do
        end do
        call locate_samples(state, set)
    contains
        pure integer function npts(ch)
            type(thomson_channel_t), intent(in) :: ch
            npts = merge(max(ch%npoints, 1), 1, ch%length > 0.0_dp)
        end function npts
    end subroutine thomson_samples

    subroutine sample_values(set, ne, te, n, dn, t, dt)
        type(sample_set_t), intent(in) :: set
        class(radial_profile_t), intent(in) :: ne, te
        real(dp), intent(out) :: n(:), dn(:), t(:), dt(:)
        integer :: k

        !$omp parallel do
        do k = 1, set%n
            if (set%inside(k)) then
                call ne%eval(set%loc(k)%u(1), n(k), dn(k))
                call te%eval(set%loc(k)%u(1), t(k), dt(k))
            else
                n(k) = 0.0_dp
                dn(k) = 0.0_dp
                t(k) = 0.0_dp
                dt(k) = 0.0_dp
            end if
        end do
        !$omp end parallel do
    end subroutine sample_values

    subroutine thomson_forward(set, nch, ne, te, y)
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nch
        class(radial_profile_t), intent(in) :: ne, te
        real(dp), intent(out) :: y(:)
        real(dp) :: n(set%n), dn(set%n), t(set%n), dt(set%n), nt(nch)
        integer :: k, c

        call sample_values(set, ne, te, n, dn, t, dt)
        y(:2 * nch) = 0.0_dp
        nt = 0.0_dp
        do k = 1, set%n
            c = set%owner(k)
            y(c) = y(c) + set%w(k) * n(k)
            nt(c) = nt(c) + set%w(k) * n(k) * t(k)
        end do
        where (y(:nch) /= 0.0_dp) y(nch + 1:2 * nch) = nt / y(:nch)
    end subroutine thomson_forward

    subroutine thomson_jvp(state, set, nch, ne, te, dp_state, dp_ne, dp_te, dy)
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nch
        class(radial_profile_t), intent(in) :: ne, te
        real(dp), intent(in) :: dp_state(:), dp_ne(:), dp_te(:)
        real(dp), intent(out) :: dy(:)
        real(dp) :: n(set%n), dn(set%n), t(set%n), dt(set%n), drho(set%n), dw(set%n)
        real(dp) :: y(2 * nch), dnt(nch), dnk, dtk
        integer :: k, c

        call sample_values(set, ne, te, n, dn, t, dt)
        call thomson_forward(set, nch, ne, te, y)
        call samples_jvp(state, set, dp_state, drho, dw)
        dy(:2 * nch) = 0.0_dp
        dnt = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            c = set%owner(k)
            call ne%param_jvp(set%loc(k)%u(1), dp_ne, dnk)
            call te%param_jvp(set%loc(k)%u(1), dp_te, dtk)
            dnk = dnk + dn(k) * drho(k)
            dtk = dtk + dt(k) * drho(k)
            dy(c) = dy(c) + set%w(k) * dnk
            dnt(c) = dnt(c) + set%w(k) * (dnk * t(k) + n(k) * dtk)
        end do
        do c = 1, nch
            if (y(c) /= 0.0_dp) dy(nch + c) = (dnt(c) - y(nch + c) * dy(c)) / y(c)
        end do
    end subroutine thomson_jvp

    subroutine thomson_vjp(state, set, nch, ne, te, ybar, pbar_state, pbar_ne, pbar_te)
        !! Accumulates the cotangents of the state and profile parameters.
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nch
        class(radial_profile_t), intent(in) :: ne, te
        real(dp), intent(in) :: ybar(:)
        real(dp), intent(inout) :: pbar_state(:), pbar_ne(:), pbar_te(:)
        real(dp) :: n(set%n), dn(set%n), t(set%n), dt(set%n), rhobar(set%n), wbar(set%n)
        real(dp) :: y(2 * nch), nbar, tbar
        integer :: k, c

        call sample_values(set, ne, te, n, dn, t, dt)
        call thomson_forward(set, nch, ne, te, y)
        rhobar = 0.0_dp
        wbar = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            c = set%owner(k)
            nbar = set%w(k) * ybar(c)
            tbar = 0.0_dp
            if (y(c) /= 0.0_dp) then
                nbar = nbar + set%w(k) * ybar(nch + c) * (t(k) - y(nch + c)) / y(c)
                tbar = set%w(k) * ybar(nch + c) * n(k) / y(c)
            end if
            call ne%param_vjp(set%loc(k)%u(1), nbar, pbar_ne)
            call te%param_vjp(set%loc(k)%u(1), tbar, pbar_te)
            rhobar(k) = nbar * dn(k) + tbar * dt(k)
        end do
        call samples_vjp(state, set, rhobar, wbar, pbar_state)
    end subroutine thomson_vjp
end module tiago_thomson
