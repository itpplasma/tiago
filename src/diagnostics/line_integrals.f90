module tiago_line_integrals
    !! Line-integrated profile diagnostics on chords (tiago_los_geometry):
    !!   interferometer   y_i = int ne dl   over the plasma part of chord i;
    !!   XICS-like        b_i = int eps dl  and  T_i = int eps T dl / b_i,
    !! with an emissivity profile eps(rho) (a latent or derived profile) and a
    !! temperature profile T (Te or Ti). XICS outputs are (b_1 .. b_n, T_1 .. T_n).
    !! The chord geometry, including the moving plasma boundary crossings, is
    !! differentiated exactly for the discrete Gauss rule (samples_jvp/_vjp).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t
    use tiago_label_profile, only: radial_profile_t
    use tiago_los_geometry, only: sample_set_t, samples_jvp, samples_vjp
    implicit none
    private

    public :: interferometer_forward, interferometer_jvp, interferometer_vjp
    public :: xics_forward, xics_jvp, xics_vjp

contains

    subroutine values(set, prof, f, df)
        type(sample_set_t), intent(in) :: set
        class(radial_profile_t), intent(in) :: prof
        real(dp), intent(out) :: f(:), df(:)
        integer :: k

        !$omp parallel do
        do k = 1, set%n
            f(k) = 0.0_dp
            df(k) = 0.0_dp
            if (set%inside(k)) call prof%eval(set%loc(k)%u(1), f(k), df(k))
        end do
        !$omp end parallel do
    end subroutine values

    subroutine interferometer_forward(set, nchord, ne, y)
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nchord
        class(radial_profile_t), intent(in) :: ne
        real(dp), intent(out) :: y(:)
        real(dp) :: n(set%n), dn(set%n)
        integer :: k

        call values(set, ne, n, dn)
        y(:nchord) = 0.0_dp
        do k = 1, set%n
            y(set%owner(k)) = y(set%owner(k)) + set%w(k) * n(k)
        end do
    end subroutine interferometer_forward

    subroutine interferometer_jvp(state, set, nchord, ne, dp_state, dp_ne, dy)
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nchord
        class(radial_profile_t), intent(in) :: ne
        real(dp), intent(in) :: dp_state(:), dp_ne(:)
        real(dp), intent(out) :: dy(:)
        real(dp) :: n(set%n), dn(set%n), drho(set%n), dw(set%n), dnk
        integer :: k

        call values(set, ne, n, dn)
        call samples_jvp(state, set, dp_state, drho, dw)
        dy(:nchord) = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            call ne%param_jvp(set%loc(k)%u(1), dp_ne, dnk)
            dy(set%owner(k)) = dy(set%owner(k)) + dw(k) * n(k) &
                + set%w(k) * (dnk + dn(k) * drho(k))
        end do
    end subroutine interferometer_jvp

    subroutine interferometer_vjp(state, set, ne, ybar, pbar_state, pbar_ne)
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        class(radial_profile_t), intent(in) :: ne
        real(dp), intent(in) :: ybar(:)
        real(dp), intent(inout) :: pbar_state(:), pbar_ne(:)
        real(dp) :: n(set%n), dn(set%n), rhobar(set%n), wbar(set%n), nbar
        integer :: k

        call values(set, ne, n, dn)
        rhobar = 0.0_dp
        wbar = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            nbar = set%w(k) * ybar(set%owner(k))
            call ne%param_vjp(set%loc(k)%u(1), nbar, pbar_ne)
            rhobar(k) = nbar * dn(k)
            wbar(k) = ybar(set%owner(k)) * n(k)
        end do
        call samples_vjp(state, set, rhobar, wbar, pbar_state)
    end subroutine interferometer_vjp

    subroutine xics_forward(set, nchord, emis, temp, y)
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nchord
        class(radial_profile_t), intent(in) :: emis, temp
        real(dp), intent(out) :: y(:)
        real(dp) :: e(set%n), de(set%n), t(set%n), dt(set%n), et(nchord)
        integer :: k, c

        call values(set, emis, e, de)
        call values(set, temp, t, dt)
        y(:2 * nchord) = 0.0_dp
        et = 0.0_dp
        do k = 1, set%n
            c = set%owner(k)
            y(c) = y(c) + set%w(k) * e(k)
            et(c) = et(c) + set%w(k) * e(k) * t(k)
        end do
        where (y(:nchord) /= 0.0_dp) y(nchord + 1:2 * nchord) = et / y(:nchord)
    end subroutine xics_forward

    subroutine xics_jvp(state, set, nchord, emis, temp, dp_state, dp_emis, dp_temp, dy)
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nchord
        class(radial_profile_t), intent(in) :: emis, temp
        real(dp), intent(in) :: dp_state(:), dp_emis(:), dp_temp(:)
        real(dp), intent(out) :: dy(:)
        real(dp) :: e(set%n), de(set%n), t(set%n), dt(set%n), drho(set%n), dw(set%n)
        real(dp) :: y(2 * nchord), det(nchord), dek, dtk
        integer :: k, c

        call values(set, emis, e, de)
        call values(set, temp, t, dt)
        call xics_forward(set, nchord, emis, temp, y)
        call samples_jvp(state, set, dp_state, drho, dw)
        dy(:2 * nchord) = 0.0_dp
        det = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            c = set%owner(k)
            call emis%param_jvp(set%loc(k)%u(1), dp_emis, dek)
            call temp%param_jvp(set%loc(k)%u(1), dp_temp, dtk)
            dek = dek + de(k) * drho(k)
            dtk = dtk + dt(k) * drho(k)
            dy(c) = dy(c) + dw(k) * e(k) + set%w(k) * dek
            det(c) = det(c) + dw(k) * e(k) * t(k) + set%w(k) * (dek * t(k) + e(k) * dtk)
        end do
        do c = 1, nchord
            if (y(c) /= 0.0_dp) dy(nchord + c) = (det(c) - y(nchord + c) * dy(c)) / y(c)
        end do
    end subroutine xics_jvp

    subroutine xics_vjp(state, set, nchord, emis, temp, ybar, pbar_state, pbar_emis, &
            pbar_temp)
        class(real_space_state_t), intent(in) :: state
        type(sample_set_t), intent(in) :: set
        integer, intent(in) :: nchord
        class(radial_profile_t), intent(in) :: emis, temp
        real(dp), intent(in) :: ybar(:)
        real(dp), intent(inout) :: pbar_state(:), pbar_emis(:), pbar_temp(:)
        real(dp) :: e(set%n), de(set%n), t(set%n), dt(set%n), rhobar(set%n), wbar(set%n)
        real(dp) :: y(2 * nchord), bbar(nchord), etbar(nchord), ebar, tbar
        integer :: k, c

        call values(set, emis, e, de)
        call values(set, temp, t, dt)
        call xics_forward(set, nchord, emis, temp, y)
        ! T = ET / B:  B gets ybar_b - ybar_T T / B, ET gets ybar_T / B.
        bbar = ybar(:nchord)
        etbar = 0.0_dp
        do c = 1, nchord
            if (y(c) /= 0.0_dp) then
                bbar(c) = bbar(c) - ybar(nchord + c) * y(nchord + c) / y(c)
                etbar(c) = ybar(nchord + c) / y(c)
            end if
        end do
        rhobar = 0.0_dp
        wbar = 0.0_dp
        do k = 1, set%n
            if (.not. set%inside(k)) cycle
            c = set%owner(k)
            ebar = set%w(k) * (bbar(c) + etbar(c) * t(k))
            tbar = set%w(k) * etbar(c) * e(k)
            call emis%param_vjp(set%loc(k)%u(1), ebar, pbar_emis)
            call temp%param_vjp(set%loc(k)%u(1), tbar, pbar_temp)
            rhobar(k) = ebar * de(k) + tbar * dt(k)
            wbar(k) = bbar(c) * e(k) + etbar(c) * e(k) * t(k)
        end do
        call samples_vjp(state, set, rhobar, wbar, pbar_state)
    end subroutine xics_vjp
end module tiago_line_integrals
