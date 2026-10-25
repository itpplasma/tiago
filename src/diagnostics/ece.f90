module tiago_ece
    !! Electron-cyclotron emission, cold-resonance model.
    !!
    !! Channel i at frequency f_i and harmonic n_i observes along its line of
    !! sight x(t) = origin + t d (from the antenna, t in [0, L]). Its cold
    !! resonance is the first point inside the plasma, counted from the
    !! antenna, where
    !!     n_i e |B(x(t))| / (2 pi m_e) = f_i,
    !! and the synthetic signal is the radiation temperature Te(rho_res) of an
    !! optically thick, unbroadened layer. The resonance moves with the state:
    !! with g(t, p) = |B|(x(t); p) - B_res(f) = 0,
    !!     dt = (dB_res - d|B|_x) / (d|B|/dt),
    !!     drho_res = drho_x + (drho/dt) dt.
    !!
    !! Extension point (documented, not implemented): a radiation-transport
    !! model (TRAVIS-like: relativistic, Doppler-broadened absorption and
    !! emission integrated along the ray, finite optical depth, refraction)
    !! would replace the single resonance by a weighted integral over the ray,
    !! reusing tiago_los_geometry samples and the same state/profile VJPs. See
    !! docs/profile-diagnostics.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_real_space_state, only: real_space_state_t, located_t, located_jvp, &
        located_vjp
    use tiago_label_profile, only: radial_profile_t
    use tiago_los_geometry, only: inside_plasma
    use tiago_cold_resonance, only: cold_resonance
    implicit none
    private

    type, public :: ece_channel_t
        character(len=:), allocatable :: label
        real(dp) :: frequency = 0.0_dp !! [Hz]
        integer :: harmonic = 2
        real(dp) :: origin(3) = 0.0_dp !! antenna
        real(dp) :: direction(3) = [-1.0_dp, 0.0_dp, 0.0_dp]
        real(dp) :: length = 1.0_dp
        integer :: nscan = 400
    end type ece_channel_t

    type, public :: ece_resonance_t
        logical :: found = .false.
        real(dp) :: t = 0.0_dp, b_res = 0.0_dp, db_res_df = 0.0_dp
        real(dp) :: dbdt = 0.0_dp, drho_dt = 0.0_dp
        type(located_t) :: loc
    end type ece_resonance_t

    public :: ece_resonances, ece_forward, ece_jvp, ece_vjp, ece_frequency_jacobian

contains

    subroutine ece_resonances(state, channels, res)
        class(real_space_state_t), intent(in) :: state
        type(ece_channel_t), intent(in) :: channels(:)
        type(ece_resonance_t), allocatable, intent(out) :: res(:)
        integer :: c

        allocate(res(size(channels)))
        !$omp parallel do schedule(dynamic, 1)
        do c = 1, size(channels)
            call find_resonance(state, channels(c), res(c))
        end do
        !$omp end parallel do
    end subroutine ece_resonances

    subroutine find_resonance(state, ch, r)
        class(real_space_state_t), intent(in) :: state
        type(ece_channel_t), intent(in) :: ch
        type(ece_resonance_t), intent(out) :: r
        type(located_t) :: prev, cur, mid
        real(dp) :: d(3), t, gp, gc, a, b, ga, tm
        integer :: j, it

        d = ch%direction / norm2(ch%direction)
        call cold_resonance(ch%frequency, real(ch%harmonic, dp), r%b_res, r%db_res_df)
        r%found = .false.
        call state%locate(ch%origin, prev)
        gp = prev%b - r%b_res
        do j = 1, ch%nscan
            t = ch%length * real(j, dp) / real(ch%nscan, dp)
            call state%locate(ch%origin + t * d, cur)
            gc = cur%b - r%b_res
            if (inside_plasma(prev) .and. inside_plasma(cur)) then
                if (sign(1.0_dp, gp) /= sign(1.0_dp, gc)) then
                    r%found = .true.
                    exit
                end if
            end if
            prev = cur
            gp = gc
        end do
        if (.not. r%found) return
        a = t - ch%length / real(ch%nscan, dp)
        b = t
        ga = gp
        do it = 1, 200
            tm = 0.5_dp * (a + b)
            if (abs(b - a) <= 4.0_dp * epsilon(1.0_dp) * max(ch%length, abs(tm))) exit
            call state%locate(ch%origin + tm * d, mid)
            if (sign(1.0_dp, mid%b - r%b_res) == sign(1.0_dp, ga)) then
                a = tm
                ga = mid%b - r%b_res
            else
                b = tm
            end if
        end do
        r%t = 0.5_dp * (a + b)
        call state%locate(ch%origin + r%t * d, r%loc)
        r%dbdt = dot_product(r%loc%dbdu, matmul(r%loc%jinv, d))
        r%drho_dt = dot_product(r%loc%jinv(1, :), d)
        r%found = inside_plasma(r%loc) .and. r%dbdt /= 0.0_dp
    end subroutine find_resonance

    subroutine ece_forward(res, te, y)
        !! y_i = Te(rho_res,i); channels without a resonance in the plasma give 0.
        type(ece_resonance_t), intent(in) :: res(:)
        class(radial_profile_t), intent(in) :: te
        real(dp), intent(out) :: y(:)
        real(dp) :: dt
        integer :: c

        do c = 1, size(res)
            y(c) = 0.0_dp
            if (res(c)%found) call te%eval(res(c)%loc%u(1), y(c), dt)
        end do
    end subroutine ece_forward

    subroutine ece_jvp(state, res, te, dp_state, dp_te, dy)
        class(real_space_state_t), intent(in) :: state
        type(ece_resonance_t), intent(in) :: res(:)
        class(radial_profile_t), intent(in) :: te
        real(dp), intent(in) :: dp_state(:), dp_te(:)
        real(dp), intent(out) :: dy(:)
        real(dp) :: du(3), dbx, dt, drho, f, dfdrho, dte
        integer :: c

        do c = 1, size(res)
            dy(c) = 0.0_dp
            if (.not. res(c)%found) cycle
            call located_jvp(state, res(c)%loc, dp_state, du, dbx)
            dt = -dbx / res(c)%dbdt
            drho = du(1) + res(c)%drho_dt * dt
            call te%eval(res(c)%loc%u(1), f, dfdrho)
            call te%param_jvp(res(c)%loc%u(1), dp_te, dte)
            dy(c) = dfdrho * drho + dte
        end do
    end subroutine ece_jvp

    subroutine ece_vjp(state, res, te, ybar, pbar_state, pbar_te)
        class(real_space_state_t), intent(in) :: state
        type(ece_resonance_t), intent(in) :: res(:)
        class(radial_profile_t), intent(in) :: te
        real(dp), intent(in) :: ybar(:)
        real(dp), intent(inout) :: pbar_state(:), pbar_te(:)
        real(dp) :: f, dfdrho, rhobar
        integer :: c

        do c = 1, size(res)
            if (.not. res(c)%found) cycle
            call te%eval(res(c)%loc%u(1), f, dfdrho)
            call te%param_vjp(res(c)%loc%u(1), ybar(c), pbar_te)
            rhobar = ybar(c) * dfdrho
            call located_vjp(state, res(c)%loc, [rhobar, 0.0_dp, 0.0_dp], &
                -rhobar * res(c)%drho_dt / res(c)%dbdt, pbar_state)
        end do
    end subroutine ece_vjp

    subroutine ece_frequency_jacobian(res, te, dydf)
        !! d y_i / d f_i (diagonal), e.g. for a frequency-calibration nuisance.
        type(ece_resonance_t), intent(in) :: res(:)
        class(radial_profile_t), intent(in) :: te
        real(dp), intent(out) :: dydf(:)
        real(dp) :: f, dfdrho
        integer :: c

        do c = 1, size(res)
            dydf(c) = 0.0_dp
            if (.not. res(c)%found) cycle
            call te%eval(res(c)%loc%u(1), f, dfdrho)
            dydf(c) = dfdrho * res(c)%drho_dt * res(c)%db_res_df / res(c)%dbdt
        end do
    end subroutine ece_frequency_jacobian
end module tiago_ece
