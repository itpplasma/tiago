program test_soft_xray
    !! Soft X-ray pinhole cameras (tiago_soft_xray).
    !! Usage: test_soft_xray            (analytic shifted-circle state)
    !!        test_soft_xray <wout.nc>  (VMEC nested state)
    !! Oracles: the filter response against composite Simpson sums over photon
    !! energy and a step filter in closed form; the pinhole etendue by hand;
    !! chord integrals on concentric circles against closed forms (Gaussian
    !! emissivity profile, unfiltered bremsstrahlung) and against Simpson sums
    !! of the closed-form integrand (filtered bremsstrahlung); the contract's
    !! derivative check including the emissivity latents and per-camera log
    !! gains; central differences of the chain state -> sample -> signals and
    !! the transposition identity of its adjoint.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_real_space_state, only: real_space_state_t, located_t
    use tiago_analytic_state, only: circular_state_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_state_sampler, only: sample_state, sample_state_jvp, sample_state_vjp
    use tiago_profile_diagnostic, only: profile_diagnostic_t, plasma_sample_t, Q_S, &
        Q_NE, Q_TE, Q_ZEFF, N_QUANTITIES, check_diagnostic_derivatives
    use tiago_sight_geometry, only: chord_t
    use tiago_soft_xray, only: soft_xray_t, pinhole_camera_t, sxr_filter_t, camera_chords, &
        filter_response, SXR_BREMSSTRAHLUNG_MODEL
    use tiago_sxr_filter_segment, only: sxr_filter_segment
    use tiago_sxr_filter_segment_small, only: sxr_filter_segment_small
    use tiago_recorded_diagnostic, only: calibrated_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), four_pi = 4.0_dp * pi
    real(dp), parameter :: n0 = 5.0e19_dp, t0 = 3.0e3_dp, zeff0 = 2.0_dp
    real(dp), parameter :: brems = 1.69e-38_dp
    integer :: failures, i
    character(len=1024) :: wout
    class(real_space_state_t), allocatable :: state
    type(soft_xray_t) :: profile_sxr, brems_sxr
    type(calibrated_diagnostic_t) :: calibrated
    real(dp), allocatable :: p0(:), vstate(:), eta0(:)
    real(dp) :: tol_fd

    failures = 0
    call filter_oracles()
    call camera_oracle()
    if (command_argument_count() == 0) then
        allocate(circular_state_t :: state)
        call chord_oracles()
        select type (state)
        type is (circular_state_t)
            state%delta = 0.15_dp
        end select
        tol_fd = 1.0e-6_dp
    else
        call get_command_argument(1, wout)
        call load_vmec(trim(wout))
        tol_fd = 1.0e-5_dp
    end if
    call cameras_from_state()
    allocate(p0(state%n_params()), vstate(state%n_params()))
    call state%get_params(p0)
    vstate = [(sin(1.7_dp * i + 0.3_dp), i = 1, size(p0))]
    vstate = vstate * max(abs(p0), 1.0e-2_dp * maxval(abs(p0)))
    eta0 = [(0.4_dp * sin(1.3_dp * i), i = 1, profile_sxr%n_nuisance())]
    call check_chain('SXR emissivity profile', profile_sxr, eta0)
    call check_chain('SXR bremsstrahlung', brems_sxr, [real(dp) ::])
    call check_calibrated()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_soft_xray passed'

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

    subroutine filter_oracles()
        !! F(Te) = int t(E) exp(-E/Te) dE / Te for a beryllium-like table with
        !! one narrow interval (the series branch), against Simpson sums on
        !! every interval plus the closed-form tail; dF/dTe against central
        !! differences; the step filter exp(-E_c / Te); no filter F = 1.
        type(sxr_filter_t) :: filter, step, none
        real(dp) :: te, f, df, fp, fm, dfp, ref
        real(dp), parameter :: temps(6) = [30.0_dp, 200.0_dp, 1.0e3_dp, 2.0e3_dp, &
            5.0e3_dp, 2.0e4_dp]
        integer :: k

        filter%energy = [500.0_dp, 1500.0_dp, 3000.0_dp, 8000.0_dp, 8000.5_dp]
        filter%transmission = [0.0_dp, 0.2_dp, 0.7_dp, 0.95_dp, 0.96_dp]
        step%energy = [1.5e3_dp]
        step%transmission = [1.0_dp]
        none%energy = [0.0_dp]
        none%transmission = [1.0_dp]
        do k = 1, size(temps)
            te = temps(k)
            call filter_response(filter, te, f, df)
            ref = simpson_filter(filter, te)
            call close('filter response vs Simpson', f, ref, 1.0e-10_dp)
            call filter_response(filter, te * (1.0_dp + 1.0e-6_dp), fp, dfp)
            call filter_response(filter, te * (1.0_dp - 1.0e-6_dp), fm, dfp)
            call close('filter dF/dTe vs central differences', df, &
                (fp - fm) / (2.0e-6_dp * te), 1.0e-7_dp)
            call filter_response(step, te, f, df)
            call close('step filter', f, exp(-1.5e3_dp / te), 1.0e-14_dp)
            call close('step filter derivative', df, 1.5e3_dp / te**2 * exp(-1.5e3_dp / te), &
                1.0e-14_dp)
            call filter_response(none, te, f, df)
            call close('no filter', f, 1.0_dp, 0.0_dp)
        end do
        ! Series and closed-form kernels agree where they meet, x = 0.01; the
        ! closed form loses about 1e-16 / (x^2 / 2) = 2e-12 to cancellation.
        call sxr_filter_segment(1000.0_dp, 10.0_dp, 0.3_dp, 0.8_dp, 1000.0_dp, f, df)
        call sxr_filter_segment_small(1000.0_dp, 10.0_dp, 0.3_dp, 0.8_dp, 1000.0_dp, fp, dfp)
        call close('filter branches meet', f, fp, 1.0e-11_dp)
        call close('filter branch derivatives meet', df, dfp, 1.0e-9_dp)
    end subroutine filter_oracles

    real(dp) function simpson_filter(filter, te) result(f)
        type(sxr_filter_t), intent(in) :: filter
        real(dp), intent(in) :: te
        integer, parameter :: m = 20000
        real(dp) :: a, h, e, t, w
        integer :: j, k, n
        n = size(filter%energy)
        f = filter%transmission(n) * exp(-filter%energy(n) / te)
        do j = 1, n - 1
            a = filter%energy(j)
            h = (filter%energy(j + 1) - a) / m
            do k = 0, m
                e = a + k * h
                t = filter%transmission(j) + (filter%transmission(j + 1) &
                    - filter%transmission(j)) * k / real(m, dp)
                w = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == 0 .or. k == m)
                f = f + w * h / 3.0_dp * t * exp(-e / te) / te
            end do
        end do
    end function simpson_filter

    subroutine camera_oracle()
        !! Two channels behind one pinhole, normals along +-z.
        type(pinhole_camera_t) :: cam
        type(chord_t), allocatable :: chords(:)
        real(dp), allocatable :: g(:)
        real(dp) :: d(3), cosine

        cam%name = 'CAM'
        cam%aperture = [3.0_dp, 0.0_dp, 1.5_dp]
        cam%aperture_normal = [0.0_dp, 0.0_dp, 1.0_dp]
        cam%aperture_area = 2.0e-6_dp
        cam%detector = reshape([3.0_dp, 0.0_dp, 1.53_dp, 3.01_dp, 0.0_dp, 1.53_dp], [3, 2])
        cam%detector_normal = [0.0_dp, 0.0_dp, -2.0_dp]
        cam%detector_area = [4.0e-6_dp, 4.0e-6_dp]
        cam%length = 2.0_dp
        call camera_chords(cam, 10, chords, g)
        call check(chords(2)%label == 'CAM_2', 'camera chord label')
        call close('chord origin at the pinhole', norm2(chords(2)%origin - cam%aperture), &
            0.0_dp, 0.0_dp)
        d = cam%aperture - cam%detector(:, 2)
        call close('chord direction', dot_product(chords(2)%direction, d / norm2(d)), &
            1.0_dp, 1.0e-15_dp)
        call close('etendue on axis', g(1), 4.0e-6_dp * 2.0e-6_dp / 0.03_dp**2, 1.0e-14_dp)
        cosine = 0.03_dp / norm2(d)
        call close('oblique etendue', g(2), 4.0e-6_dp * 2.0e-6_dp * cosine**2 / norm2(d)**2, &
            1.0e-14_dp)
        cam%aperture_area = 0.0_dp
        call camera_chords(cam, 10, chords, g)
        call close('no etendue without areas', g(2), four_pi, 0.0_dp)
    end subroutine camera_oracle

    subroutine chord_oracles()
        !! Concentric circles (R0 = 3, a = 1); vertical chord at R = 3.4
        !! between its boundary crossings, s = 0.16 + z^2.
        type(soft_xray_t) :: sxr
        type(chord_t) :: chord(1)
        type(sxr_filter_t) :: step(1)
        real(dp), allocatable :: y(:), s(:)
        real(dp) :: u0, ref, z, h, w, te
        integer :: k, ierr

        u0 = sqrt(1.0_dp - 0.16_dp)
        chord(1) = chord_t('VERT', [3.4_dp, 0.0_dp, -u0], [0.0_dp, 0.0_dp, 1.0_dp], &
            2.0_dp * u0, 40, 4)
        ! eps = 2e4 exp(-3 s) W m^-3, exactly representable (log linear in s).
        call sxr%setup(chord, etendue=[1.0e-7_dp])
        call sxr%set_emissivity_basis([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp])
        s = [(0.05_dp * k, k = 0, 20)]
        call sxr%set_emissivity_table(s, 2.0e4_dp * exp(-3.0_dp * s), ierr)
        call check(ierr == 0, 'emissivity table fit')
        call predict(sxr, [real(dp) ::], y)
        ref = 2.0e4_dp * exp(-0.48_dp) * sqrt(pi / 3.0_dp) * erf(sqrt(3.0_dp) * u0)
        call close('Gaussian emissivity chord power', y(1), 1.0e-7_dp / four_pi * ref, 1.0e-12_dp)

        ! Bremsstrahlung, ne = n0 (1 - s), Te = T0 (1 - s), Zeff = 2:
        ! int (u0^2 - z^2)^(5/2) dz = 5 pi u0^6 / 16.
        chord(1)%npanel = 400
        call sxr%setup(chord)
        sxr%model = SXR_BREMSSTRAHLUNG_MODEL
        call predict(sxr, [real(dp) ::], y)
        ref = brems * n0**2 * zeff0 * sqrt(t0) * 5.0_dp * pi * u0**6 / 16.0_dp
        call close('unfiltered bremsstrahlung chord', y(1), ref, 1.0e-8_dp)

        ! Step filter at 1.5 keV: Simpson sums of the closed-form integrand.
        step(1)%energy = [1.5e3_dp]
        step(1)%transmission = [1.0_dp]
        call sxr%setup(chord, filters=step)
        sxr%model = SXR_BREMSSTRAHLUNG_MODEL
        call predict(sxr, [real(dp) ::], y)
        ref = 0.0_dp
        h = 2.0_dp * u0 / 20000
        do k = 0, 20000
            z = -u0 + k * h
            te = t0 * max(u0**2 - z**2, 0.0_dp)
            w = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == 0 .or. k == 20000)
            if (te > 0.0_dp) ref = ref + w * h / 3.0_dp * brems * (n0 * te / t0)**2 * zeff0 &
                * sqrt(te) * exp(-1.5e3_dp / te)
        end do
        call close('filtered bremsstrahlung chord vs Simpson', y(1), ref, 1.0e-8_dp)
    end subroutine chord_oracles

    subroutine cameras_from_state()
        !! Two cameras placed relative to the state's axis and boundary: one
        !! above the plasma looking down, one on the outboard side.
        type(pinhole_camera_t) :: cams(2)
        real(dp) :: xa(3), xb(3), jac(3, 3), b, db(3), size0, up(3)
        real(dp), allocatable :: mean(:)
        integer :: k

        call state%eval([0.0_dp, 0.0_dp, 0.2_dp], xa, jac, b, db)
        call state%eval([1.0_dp, 0.0_dp, 0.2_dp], xb, jac, b, db)
        size0 = norm2(xb - xa)
        up = [0.0_dp, 0.0_dp, 1.0_dp]
        cams(1)%name = 'TOP'
        cams(1)%aperture = xa + 1.6_dp * size0 * up
        allocate(cams(1)%detector(3, 5))
        do k = 1, 5
            cams(1)%detector(:, k) = cams(1)%aperture + 0.05_dp * size0 * &
                [0.15_dp * (k - 3), 0.0_dp, 1.0_dp]
        end do
        cams(1)%detector_area = [(1.0e-6_dp, k = 1, 5)]
        cams(1)%aperture_area = 1.0e-6_dp
        cams(1)%length = 3.2_dp * size0
        cams(2)%name = 'OUT'
        cams(2)%aperture = xb + 0.8_dp * (xb - xa)
        allocate(cams(2)%detector(3, 4))
        do k = 1, 4
            cams(2)%detector(:, k) = cams(2)%aperture + 0.05_dp * size0 * &
                ((xb - xa) / size0 + [0.0_dp, 0.0_dp, 0.2_dp * (k - 2.5_dp)])
        end do
        cams(2)%length = 2.6_dp * size0
        cams(2)%filter%energy = [800.0_dp, 2000.0_dp, 6000.0_dp]
        cams(2)%filter%transmission = [0.0_dp, 0.5_dp, 0.9_dp]
        call profile_sxr%set_cameras(cams, 100)
        call profile_sxr%set_emissivity_basis([0.0_dp, 0.3_dp, 0.6_dp, 1.0_dp])
        mean = 1.0e4_dp * (1.0_dp - 0.9_dp * profile_sxr%emissivity%greville)
        call profile_sxr%set_emissivity_prior(mean, 0.5_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call check(profile_sxr%n_nuisance() == 6, 'one emissivity latent per coefficient')
        call brems_sxr%set_cameras(cams, 100)
        brems_sxr%model = SXR_BREMSSTRAHLUNG_MODEL
        call check(all(brems_sxr%camera == [1, 1, 1, 1, 1, 2, 2, 2, 2]), 'camera index per chord')
    end subroutine cameras_from_state

    subroutine sample(diag, smp, loc)
        !! Provider rows from the state; ne = n0 (1 - s), Te = T0 (1 - s),
        !! Zeff = 2 inside.
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
            smp%q(Q_ZEFF, k) = zeff0
        end do
    end subroutine sample

    subroutine predict(diag, eta, y)
        class(profile_diagnostic_t), intent(in) :: diag
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: y(:)
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        call sample(diag, smp, loc)
        allocate(y(diag%n_signals()))
        call diag%predict(smp, eta, y)
    end subroutine predict

    subroutine check_chain(name, diag, eta)
        !! Contract derivatives (eta included), then the chain through the
        !! state parameters: JVP vs central differences, VJP transposition.
        character(len=*), intent(in) :: name
        class(profile_diagnostic_t), intent(in) :: diag
        real(dp), intent(in) :: eta(:)
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: y(:), yp(:), ym(:), ydot(:), ybar(:), qdot(:, :), qbar(:, :)
        real(dp), allocatable :: etadot(:), etabar(:), pbar(:)
        real(dp) :: h, jerr, aerr
        integer :: k

        allocate(etadot(size(eta)), etabar(size(eta)))
        etadot = 0.0_dp
        call sample(diag, smp, loc)
        call predict(diag, eta, y)
        call check(all(y > 0.0_dp), name//' signals are positive')
        call check_diagnostic_derivatives(diag, smp, eta, 1.0e-6_dp, jerr, aerr)
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
        end do
        call diag%jvp(smp, eta, qdot, etadot, ydot)
        h = 1.0e-6_dp
        call state%set_params(p0 + h * vstate)
        call predict(diag, eta, yp)
        call state%set_params(p0 - h * vstate)
        call predict(diag, eta, ym)
        call state%set_params(p0)
        call compare(name//' chain JVP vs FD', ydot, (yp - ym) / (2.0_dp * h))

        ybar = [(cos(0.37_dp * k + 0.2_dp), k = 1, size(y))]
        call diag%vjp(smp, eta, ybar, qbar, etabar)
        do k = 1, smp%n
            if (.not. smp%inside(k)) cycle
            qbar(Q_S, k) = qbar(Q_S, k) - n0 * qbar(Q_NE, k) - t0 * qbar(Q_TE, k)
        end do
        pbar = 0.0_dp
        call sample_state_vjp(state, loc, qbar, pbar)
        call close(name//' chain VJP transposition', dot_product(ybar, ydot), &
            dot_product(pbar, vstate), 1.0e-11_dp)
    end subroutine check_chain

    subroutine check_calibrated()
        !! One log gain per camera (calibrated_diagnostic_t groups): the signal
        !! is exp(g_camera) times the bare signal; contract derivatives.
        type(calibration_prior_t) :: prior
        type(plasma_sample_t) :: smp
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: raw(:), y(:), eta(:)
        real(dp) :: gains(2), jerr, aerr
        integer :: status

        call prior%init(log([2.0e5_dp, 3.0e5_dp]), reshape([0.04_dp, 0.0_dp, 0.0_dp, &
            0.01_dp], [2, 2]), status)
        call check(status == 0, 'per-camera log-gain prior')
        call calibrated%init(profile_sxr, prior, status, groups=profile_sxr%camera, &
            log_gains=.true.)
        call check(status == 0, 'calibrated soft X-ray setup')
        eta = [eta0, 0.5_dp, -1.0_dp]
        call predict(profile_sxr, eta0, raw)
        call predict(calibrated, eta, y)
        gains = [2.0e5_dp * exp(0.2_dp * 0.5_dp), 3.0e5_dp * exp(-0.1_dp)]
        call close('camera 1 gain', maxval(abs(y(:5) / raw(:5))), gains(1), 1.0e-12_dp)
        call close('camera 2 gain', minval(abs(y(6:) / raw(6:))), gains(2), 1.0e-12_dp)
        call sample(calibrated, smp, loc)
        call check_diagnostic_derivatives(calibrated, smp, eta, 1.0e-6_dp, jerr, aerr)
        print '(A,": contract JVP error ",ES10.3,", adjoint error ",ES10.3)', &
            'calibrated SXR', jerr, aerr
        call check(jerr <= 1.0e-6_dp .and. aerr <= 1.0e-12_dp, 'calibrated SXR derivatives')
    end subroutine check_calibrated

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
end program test_soft_xray
