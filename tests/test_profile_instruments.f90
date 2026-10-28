program test_profile_instruments
    !! Thomson scattering, interferometer and XICS on the profile-diagnostic
    !! contract, sampled from a real-space state by tiago_state_sampler.
    !! Usage: test_profile_instruments            (analytic shifted-circle state)
    !!        test_profile_instruments <wout.nc>  (VMEC nested state)
    !! Oracles: closed-form integrals of profiles linear in s on concentric
    !! circles with |B| = B0 R0 / R; the contract's derivative check of every
    !! instrument; central finite differences of the complete chain
    !! state -> sample -> signals for the JVP; and the transposition identity
    !! <y_bar, J v> = <J^T y_bar, v> for the VJP.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok
    use tiago_analytic_state, only: circular_state_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_state_sampler, only: sample_state, sample_state_jvp, sample_state_vjp
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_S, &
        Q_NE, Q_TE, Q_TI, Q_ZEFF, N_QUANTITIES, check_diagnostic_derivatives
    use tiago_sight_geometry, only: chord_t
    use tiago_thomson_scattering, only: thomson_t, thomson_volume_t
    use tiago_line_integrals, only: interferometer_t, xics_t
    use tiago_dispersion_interferometer, only: dispersion_interferometer_t
    use tiago_calibration, only: calibration_prior_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: n0 = 5.0e19_dp, t0 = 3.0e3_dp, ti0 = 1.5e3_dp, zeff0 = 2.0_dp
    integer :: failures, i
    character(len=1024) :: wout
    class(real_space_state_t), allocatable :: state
    type(thomson_t) :: ts
    type(interferometer_t) :: interf
    type(xics_t) :: xics
    type(dispersion_interferometer_t) :: dispersion
    real(dp), allocatable :: p0(:), vstate(:)
    real(dp) :: tol_fd

    failures = 0
    if (command_argument_count() == 0) then
        allocate(circular_state_t :: state)
        call analytic_oracles()
        select type (state)
            type is (circular_state_t)
            state%delta = 0.15_dp
        end select
        tol_fd = 1.0e-6_dp
    else
        call get_command_argument(1, wout)
        call load_vmec(trim(wout))
        call locate_roundtrip()
        tol_fd = 1.0e-5_dp
    end if
    call geometry_from_state()
    allocate(p0(state%n_params()), vstate(state%n_params()))
    call state%get_params(p0)
    vstate = [(sin(1.7_dp * i + 0.3_dp), i = 1, size(p0))]
    vstate = vstate * max(abs(p0), 1.0e-2_dp * maxval(abs(p0)))
    call check_chain('Thomson', ts, 1.0e-6_dp)
    call check_chain('interferometer', interf, 1.0e-6_dp)
    call check_chain('XICS', xics, 1.0e-6_dp)
    call check_chain('dispersion I/Q', dispersion, 1.0e-4_dp)
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_profile_instruments passed'

contains

    subroutine load_vmec(path)
        character(len=*), intent(in) :: path
        type(vmec_state_t), allocatable :: v
        integer :: ierr
        character(len=:), allocatable :: message
        allocate(v)
        call v%load(path, ierr, message)
        if (ierr /= 0) then
            write(error_unit, '(A)') message
            error stop 1
        end if
        call move_alloc(v, state)
    end subroutine load_vmec

    subroutine sample(diag, smp, loc)
        !! Provider rows from the state; kinetic rows linear in s inside:
        !! ne = n0 (1 - s), Te = T0 (1 - s), Ti = Ti0 (1 - s), Zeff = 2.
        class(profile_diagnostic_t), intent(in) :: diag
        type(plasma_sample_t), intent(out) :: smp
        type(located_t), allocatable, intent(out) :: loc(:)
        real(dp), allocatable :: x(:, :), dx(:, :, :), eta(:)
        integer :: k

        allocate(eta(diag%n_nuisance()))
        eta = 0.0_dp
        call diag%sample_points(eta, x, dx)
        call sample_state(state, x, smp, loc)
        do k = 1, smp%n
            smp%q(Q_ZEFF, k) = 1.0_dp
            if (.not. smp%inside(k)) cycle
            smp%q(Q_NE, k) = n0 * (1.0_dp - smp%q(Q_S, k))
            smp%q(Q_TE, k) = t0 * (1.0_dp - smp%q(Q_S, k))
            smp%q(Q_TI, k) = ti0 * (1.0_dp - smp%q(Q_S, k))
            smp%q(Q_ZEFF, k) = zeff0
            smp%dq_dx(:, Q_NE, k) = -n0 * smp%dq_dx(:, Q_S, k)
            smp%dq_dx(:, Q_TE, k) = -t0 * smp%dq_dx(:, Q_S, k)
            smp%dq_dx(:, Q_TI, k) = -ti0 * smp%dq_dx(:, Q_S, k)
        end do
    end subroutine sample

    subroutine predict(diag, y)
        class(profile_diagnostic_t), intent(in) :: diag
        real(dp), allocatable, intent(out) :: y(:)
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: eta(:)
        allocate(eta(diag%n_nuisance()))
        eta = 0.0_dp
        call sample(diag, smp, loc)
        allocate(y(diag%n_signals()))
        call diag%predict(smp, eta, y)
    end subroutine predict

    subroutine analytic_oracles()
        !! Concentric circles (Delta = 0), R0 = 3, a = 1, B0 = 2.5.
        type(located_t) :: loc
        type(chord_t) :: chord(1)
        real(dp), allocatable :: y(:)
        real(dp) :: c, zl, zh, m1, m2, f, u0
        integer :: k

        do k = 1, 20
            f = 0.05_dp * k
            call state%locate([(3.0_dp + 0.7_dp * f * cos(2.0_dp * k)) * cos(0.3_dp * k), &
                (3.0_dp + 0.7_dp * f * cos(2.0_dp * k)) * sin(0.3_dp * k), &
                0.7_dp * f * sin(2.0_dp * k)], loc)
            call check(loc%status == locate_ok, 'analytic locate converges')
            call close('analytic locate rho', loc%u(1), 0.7_dp * f, 1.0e-12_dp)
        end do

        ! Thomson: vertical laser at R = 3.3 over z in [0.15, 0.25].
        call ts%init([thomson_volume_t('TS', [3.3_dp, 0.0_dp, 0.2_dp], &
            [0.0_dp, 0.0_dp, 1.0_dp], 0.1_dp, 4)])
        call predict(ts, y)
        c = 0.09_dp
        zl = 0.15_dp
        zh = 0.25_dp
        m1 = (anti1(zh, c) - anti1(zl, c)) / (zh - zl)
        m2 = (anti2(zh, c) - anti2(zl, c)) / (zh - zl)
        call close('Thomson ne volume average', y(1), n0 * m1, 1.0e-12_dp)
        call close('Thomson Te density-weighted average', y(2), t0 * m2 / m1, 1.0e-12_dp)

        ! Vertical chord at R = 3.4 between its boundary crossings.
        u0 = sqrt(1.0_dp - 0.16_dp)
        chord(1) = chord_t('VERT', [3.4_dp, 0.0_dp, -u0], [0.0_dp, 0.0_dp, 1.0_dp], &
            2.0_dp * u0, 20, 4)
        call interf%init(chord)
        call predict(interf, y)
        call close('interferometer int ne dl', y(1), n0 * 4.0_dp * u0**3 / 3.0_dp, 1.0e-12_dp)
        chord(1)%npanel = 400
        xics%e_line = 0.0_dp
        call xics%init(chord)
        call predict(xics, y)
        ! eps = n0^2 (Zeff - 1) / (Z (Z - 1) sqrt(T0)) (u0^2 - z^2)^(3/2)
        call close('XICS brightness', y(1), n0**2 * (zeff0 - 1.0_dp) / &
            (18.0_dp * 17.0_dp * sqrt(t0)) * 3.0_dp * pi / 8.0_dp * u0**4, 1.0e-5_dp)
        call close('XICS emissivity-weighted Ti', y(2), ti0 * 5.0_dp / 6.0_dp * u0**2, &
            1.0e-5_dp)
        xics%e_line = 3.1e3_dp
    end subroutine analytic_oracles

    pure real(dp) function anti1(z, c)
        real(dp), intent(in) :: z, c
        anti1 = (1.0_dp - c) * z - z**3 / 3.0_dp
    end function anti1

    pure real(dp) function anti2(z, c)
        real(dp), intent(in) :: z, c
        anti2 = (1.0_dp - c)**2 * z - 2.0_dp * (1.0_dp - c) * z**3 / 3.0_dp + z**5 / 5.0_dp
    end function anti2

    subroutine locate_roundtrip()
        type(located_t) :: loc
        real(dp) :: u(3), x(3), jac(3, 3), b, db(3)
        integer :: k

        do k = 1, 24
            u = [0.04_dp * k, 0.7_dp * k, 0.11_dp * k]
            call state%eval(u, x, jac, b, db)
            call state%locate(x, loc)
            call check(loc%status == locate_ok, 'VMEC locate converges')
            call close('VMEC locate roundtrip rho', loc%u(1), u(1), 1.0e-9_dp)
            call close('VMEC locate |B|', loc%b, b, 1.0e-9_dp)
        end do
    end subroutine locate_roundtrip

    subroutine geometry_from_state()
        !! Instruments placed relative to the state's own axis and boundary.
        type(thomson_volume_t) :: vol(4)
        type(chord_t) :: chords(3)
        type(calibration_prior_t) :: prior
        real(dp) :: cm(8), cov(8, 8), sd(8)
        integer :: status
        real(dp) :: xa(3), xb(3), xm(3), jac(3, 3), b, db(3), d(3), phi0, size0
        integer :: k

        phi0 = 0.2_dp
        call state%eval([0.0_dp, 0.0_dp, phi0], xa, jac, b, db)
        call state%eval([1.0_dp, 0.0_dp, phi0], xb, jac, b, db)
        size0 = norm2(xb - xa)
        do k = 1, 4
            call state%eval([0.2_dp * k, 0.3_dp + 0.4_dp * k, phi0], xm, jac, b, db)
            vol(k) = thomson_volume_t('TS', xm, [0.1_dp, 0.2_dp, 1.0_dp], 0.04_dp * size0, 3)
        end do
        call ts%init(vol, aligned=.false.)
        d = (xa - xb) / size0
        chords(1) = chord_t('C1', xa - 2.0_dp * size0 * [0.0_dp, 0.1_dp, 1.0_dp], &
            [0.0_dp, 0.1_dp, 1.0_dp], 4.0_dp * size0, 200, 4)
        chords(2) = chord_t('C2', xb + 0.5_dp * (xb - xa), d, 2.5_dp * size0, 200, 4)
        chords(3) = chord_t('C3', 0.5_dp * (xa + xb) + [0.0_dp, -2.0_dp, 0.05_dp] * size0, &
            [0.0_dp, 1.0_dp, 0.0_dp], 4.0_dp * size0, 200, 4)
        call interf%init(chords)
        call xics%init(chords)
        cm = [1.0e-5_dp, 1.0_dp, 1.2_dp, 0.1_dp, -0.2_dp, 0.3_dp, &
            1.0_dp, 0.2_dp]
        sd = [1.0e-8_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, &
            0.1_dp, 0.1_dp]
        cov = 0.0_dp
        do k = 1, 8
            cov(k, k) = sd(k)**2
        end do
        call prior%init(cm, cov, status)
        call check(status == 0, 'dispersion prior')
        call dispersion%init(chords, prior, [1.0e-8_dp, 2.0e-8_dp, 0.0_dp], &
            [0.2_dp, -0.3_dp, 0.5_dp], status)
        call check(status == 0, 'dispersion init')
    end subroutine geometry_from_state

    subroutine check_chain(name, diag, step)
        !! step: relative q step of the contract check. The ECE resonance
        !! fraction divides by sample differences of |B|, so its central-difference
        !! truncation error grows like (step |B| / dB)^2 and needs a smaller step;
        !! the linear line integrals need a larger one against roundoff.
        character(len=*), intent(in) :: name
        class(profile_diagnostic_t), intent(in) :: diag
        real(dp), intent(in) :: step
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: y(:), yp(:), ym(:), ydot(:), ybar(:), qdot(:, :), qbar(:, :)
        real(dp), allocatable :: eta(:), etadot(:), etabar(:), pbar(:)
        real(dp) :: h, jerr, aerr
        integer :: k

        allocate(eta(diag%n_nuisance()), etadot(diag%n_nuisance()), etabar(diag%n_nuisance()))
        eta = 0.0_dp
        etadot = 0.0_dp
        call sample(diag, smp, loc)
        call predict(diag, y)
        call check(any(abs(y) > 0.0_dp), name//' signals are non-zero')
        call check_diagnostic_derivatives(diag, smp, eta, step, jerr, aerr)
        print '(A,": contract JVP error ",ES10.3,", adjoint error ",ES10.3)', name, jerr, aerr
        call check(jerr <= 1.0e-6_dp .and. aerr <= 1.0e-12_dp, name//' contract derivatives')

        allocate(qdot(N_QUANTITIES, smp%n), qbar(N_QUANTITIES, smp%n), ydot(size(y)), &
            ybar(size(y)), pbar(size(p0)))
        qdot = 0.0_dp
        call sample_state_jvp(state, loc, vstate, qdot)
        do k = 1, smp%n
            if (.not. smp%inside(k)) cycle
            qdot(Q_NE, k) = -n0 * qdot(Q_S, k)
            qdot(Q_TE, k) = -t0 * qdot(Q_S, k)
            qdot(Q_TI, k) = -ti0 * qdot(Q_S, k)
        end do
        call diag%jvp(smp, eta, qdot, etadot, ydot)
        h = 1.0e-6_dp
        call state%set_params(p0 + h * vstate)
        call predict(diag, yp)
        call state%set_params(p0 - h * vstate)
        call predict(diag, ym)
        call state%set_params(p0)
        call compare(name//' chain JVP vs FD', ydot, (yp - ym) / (2.0_dp * h))

        ybar = [(cos(0.37_dp * k + 0.2_dp), k = 1, size(y))]
        call diag%vjp(smp, eta, ybar, qbar, etabar)
        do k = 1, smp%n
            if (.not. smp%inside(k)) cycle
            qbar(Q_S, k) = qbar(Q_S, k) - n0 * qbar(Q_NE, k) - t0 * qbar(Q_TE, k) &
                - ti0 * qbar(Q_TI, k)
        end do
        pbar = 0.0_dp
        call sample_state_vjp(state, loc, qbar, pbar)
        call close(name//' chain VJP transposition', dot_product(ybar, ydot), &
            dot_product(pbar, vstate), 1.0e-11_dp)
    end subroutine check_chain

    subroutine compare(name, analytic, fd)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: analytic(:), fd(:)
        real(dp) :: err
        err = maxval(abs(analytic - fd)) / max(maxval(abs(fd)), 1.0e-300_dp)
        print '(A,": max rel error ",ES10.3)', name, err
        if (.not. (err <= tol_fd)) then
            write(error_unit, '(A,ES12.4)') 'FAIL '//name//': ', err
            failures = failures + 1
        end if
    end subroutine compare

    subroutine close(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, tol
        if (.not. (abs(a - b) <= tol * max(1.0e-300_dp, abs(b)))) then
            write(error_unit, '(A,2ES24.16)') 'FAIL '//name//': ', a, b
            failures = failures + 1
        end if
    end subroutine close

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: name
        if (.not. ok) then
            write(error_unit, '(A)') 'FAIL '//name
            failures = failures + 1
        end if
    end subroutine check
end program test_profile_instruments
