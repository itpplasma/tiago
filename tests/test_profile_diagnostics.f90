program test_profile_diagnostics
    !! Thomson, ECE, interferometer and XICS forward models and derivatives.
    !! Usage: test_profile_diagnostics            (analytic shifted-circle state)
    !!        test_profile_diagnostics <wout.nc>  (VMEC nested state)
    !! Oracles: closed-form integrals of quadratic profiles on concentric
    !! circles with |B| = B0 R0 / R; central finite differences of the complete
    !! forward models (relocation and moving crossings included) for JVPs; and
    !! the transposition identity <ybar, J v> = <J^T ybar, v> for VJPs.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok
    use tiago_analytic_state, only: circular_state_t
    use tiago_vmec_state, only: vmec_state_t
    use tiago_label_profile, only: node_profile_t
    use tiago_los_geometry, only: chord_t, sample_set_t, chord_samples
    use tiago_thomson, only: thomson_channel_t, thomson_samples, thomson_forward, &
        thomson_jvp, thomson_vjp
    use tiago_ece, only: ece_channel_t, ece_resonance_t, ece_resonances, ece_forward, &
        ece_jvp, ece_vjp, ece_frequency_jacobian
    use tiago_line_integrals, only: interferometer_forward, interferometer_jvp, &
        interferometer_vjp, xics_forward, xics_jvp, xics_vjp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: n0 = 5.0_dp, t0 = 3.0_dp, e0 = 2.0_dp, ti0 = 1.5_dp
    real(dp), parameter :: me = 9.1093837015e-31_dp, qe = 1.602176634e-19_dp
    integer :: failures, nargs
    character(len=1024) :: wout
    type(node_profile_t) :: ne, te, emis, ti
    class(real_space_state_t), allocatable :: state
    type(thomson_channel_t), allocatable :: ts(:)
    type(ece_channel_t), allocatable :: ece(:)
    type(chord_t), allocatable :: chords(:)
    real(dp), allocatable :: p0(:), vstate(:)
    real(dp) :: tol_fd

    failures = 0
    call profile(ne, n0)
    call profile(te, t0)
    call profile(emis, e0)
    call profile(ti, ti0)
    nargs = command_argument_count()
    if (nargs == 0) then
        allocate(circular_state_t :: state)
        call analytic_oracles()
        select type (state)
            type is (circular_state_t)
            state%delta = 0.15_dp
        end select
        call geometry_from_state()
        tol_fd = 1.0e-6_dp
    else
        call get_command_argument(1, wout)
        call load_vmec(trim(wout))
        call locate_roundtrip()
        call geometry_from_state()
        tol_fd = 1.0e-5_dp
    end if
    allocate(p0(state%n_params()), vstate(state%n_params()))
    call state%get_params(p0)
    call random_direction(vstate)
    vstate = vstate * max(abs(p0), 1.0e-2_dp * maxval(abs(p0)))
    call check_thomson()
    call check_ece()
    call check_interferometer()
    call check_xics()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_profile_diagnostics passed'

contains

    subroutine profile(prof, scale)
        !! scale (1 - rho^2) on 11 uniform nodes; reproduced exactly.
        type(node_profile_t), intent(out) :: prof
        real(dp), intent(in) :: scale
        integer :: j
        prof%nodes = [(0.1_dp * j, j = 0, 10)]
        prof%values = scale * (1.0_dp - prof%nodes**2)
    end subroutine profile

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

    subroutine analytic_oracles()
        !! Concentric circles (Delta = 0), R0 = 3, a = 1, B0 = 2.5.
        type(sample_set_t) :: set
        type(ece_resonance_t), allocatable :: res(:)
        real(dp) :: y(4), c, z0, zl, zh, m1, m2, f, rres, u0
        type(located_t) :: loc
        integer :: k

        ! Locate: rho = sqrt((R - R0)^2 + Z^2) / a.
        do k = 1, 20
            f = 0.05_dp * k
            call state%locate([(3.0_dp + 0.7_dp * f * cos(2.0_dp * k)) * cos(0.3_dp * k), &
                (3.0_dp + 0.7_dp * f * cos(2.0_dp * k)) * sin(0.3_dp * k), &
                0.7_dp * f * sin(2.0_dp * k)], loc)
            call check(loc%status == locate_ok, 'analytic locate converges')
            call close('analytic locate rho', loc%u(1), 0.7_dp * f, 1.0e-12_dp)
        end do

        ! Thomson: vertical laser at R = 3.3 over z in [0.15, 0.25].
        allocate(ts(1))
        ts(1) = thomson_channel_t('TS', [3.3_dp, 0.0_dp, 0.2_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
            0.1_dp, 4)
        call thomson_samples(state, ts, set)
        call thomson_forward(set, 1, ne, te, y)
        c = 0.09_dp
        zl = 0.15_dp
        zh = 0.25_dp
        ! mean of q = 1 - c - z^2 and of q^2 over [zl, zh]
        m1 = (anti1(zh, c) - anti1(zl, c)) / (zh - zl)
        m2 = (anti2(zh, c) - anti2(zl, c)) / (zh - zl)
        call close('Thomson ne volume average', y(1), n0 * m1, 1.0e-12_dp)
        call close('Thomson Te density-weighted average', y(2), t0 * m2 / m1, 1.0e-12_dp)

        ! ECE: second harmonic resonant at R = 3.2 on a chord at z = 0.1.
        allocate(ece(1))
        rres = 3.2_dp
        f = 2.0_dp * qe * (2.5_dp * 3.0_dp / rres) / (2.0_dp * pi * me)
        ece(1) = ece_channel_t('ECE', f, 2, [4.5_dp, 0.0_dp, 0.1_dp], &
            [-1.0_dp, 0.0_dp, 0.0_dp], 3.0_dp, 300)
        call ece_resonances(state, ece, res)
        call check(res(1)%found, 'ECE resonance found')
        call close('ECE resonance position', 4.5_dp - res(1)%t, rres, 1.0e-12_dp)
        call ece_forward(res, te, y)
        call close('ECE Te at the resonance', y(1), t0 * (1.0_dp - 0.04_dp - 0.01_dp), &
            1.0e-12_dp)

        ! Vertical chord at R = 3.4 through the plasma.
        allocate(chords(1))
        chords(1) = chord_t('VERT', [3.4_dp, 0.0_dp, -2.0_dp], [0.0_dp, 0.0_dp, 1.0_dp], &
            4.0_dp, 200, 8)
        call chord_samples(state, chords, set)
        u0 = sqrt(1.0_dp - 0.16_dp)
        call interferometer_forward(set, 1, ne, y)
        call close('interferometer int ne dl', y(1), n0 * 4.0_dp * u0**3 / 3.0_dp, 1.0e-12_dp)
        call xics_forward(set, 1, emis, ti, y)
        call close('XICS brightness', y(1), e0 * 4.0_dp * u0**3 / 3.0_dp, 1.0e-12_dp)
        call close('XICS emissivity-weighted T', y(2), ti0 * 0.8_dp * u0**2, 1.0e-12_dp)
        deallocate(ts, ece, chords)
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
        !! Channels placed relative to the state's own axis and boundary at phi0.
        real(dp) :: xa(3), xb(3), xm(3), jac(3, 3), b, bm, db(3), d(3), phi0, f
        integer :: k

        phi0 = 0.2_dp
        call state%eval([0.0_dp, 0.0_dp, phi0], xa, jac, b, db)
        call state%eval([1.0_dp, 0.0_dp, phi0], xb, jac, b, db)
        allocate(ts(4), ece(3), chords(3))
        do k = 1, 4
            call state%eval([0.2_dp * k, 0.3_dp + 0.4_dp * k, phi0], xm, jac, b, db)
            ts(k) = thomson_channel_t('TS', xm, [0.1_dp, 0.2_dp, 1.0_dp], &
                0.04_dp * norm2(xb - xa), 3)
        end do
        d = (xa - xb) / norm2(xa - xb)
        do k = 1, 3
            call state%eval([0.25_dp * k, 0.05_dp, phi0], xm, jac, bm, db)
            f = 2.0_dp * qe * bm / (2.0_dp * pi * me)
            ece(k) = ece_channel_t('ECE', f, 2, xb - 0.6_dp * (xa - xb) &
                + [0.0_dp, 0.0_dp, 0.01_dp * k], d, 3.0_dp * norm2(xa - xb), 300)
        end do
        chords(1) = chord_t('C1', xa - 2.0_dp * norm2(xb - xa) * [0.0_dp, 0.1_dp, 1.0_dp], &
            [0.0_dp, 0.1_dp, 1.0_dp], 4.0_dp * norm2(xb - xa), 300, 12)
        chords(2) = chord_t('C2', xb + 0.5_dp * (xb - xa), d, 2.5_dp * norm2(xb - xa), 300, 12)
        chords(3) = chord_t('C3', 0.5_dp * (xa + xb) + [0.0_dp, -2.0_dp, 0.05_dp] &
            * norm2(xb - xa), [0.0_dp, 1.0_dp, 0.0_dp], 4.0_dp * norm2(xb - xa), 400, 12)
    end subroutine geometry_from_state

    subroutine check_thomson()
        integer :: n, nn, nt
        real(dp), allocatable :: y(:), yp(:), ym(:), dy(:), vn(:), vt(:), ybar(:)
        real(dp), allocatable :: ps(:), pn(:), pt(:)
        type(sample_set_t) :: set
        real(dp) :: h

        n = size(ts)
        nn = ne%n_params()
        nt = te%n_params()
        allocate(y(2 * n), yp(2 * n), ym(2 * n), dy(2 * n), vn(nn), vt(nt), ybar(2 * n))
        call random_direction(vn)
        call random_direction(vt)
        call thomson_samples(state, ts, set)
        call check(all(set%inside), 'Thomson samples inside the plasma')
        call thomson_jvp(state, set, n, ne, te, vstate, vn, vt, dy)
        h = 1.0e-6_dp
        call perturb(h, ne, n0, vn, te, t0, vt)
        call thomson_samples(state, ts, set)
        call thomson_forward(set, n, ne, te, yp)
        call perturb(-h, ne, n0, vn, te, t0, vt)
        call thomson_samples(state, ts, set)
        call thomson_forward(set, n, ne, te, ym)
        call perturb(0.0_dp, ne, n0, vn, te, t0, vt)
        call compare('Thomson JVP vs FD', dy, (yp - ym) / (2.0_dp * h))
        call thomson_samples(state, ts, set)
        call random_direction(ybar)
        allocate(ps(size(p0)), pn(nn), pt(nt))
        ps = 0.0_dp
        pn = 0.0_dp
        pt = 0.0_dp
        call thomson_vjp(state, set, n, ne, te, ybar, ps, pn, pt)
        call close('Thomson VJP transposition', dot_product(ybar, dy), &
            dot_product(ps, vstate) + dot_product(pn, vn) + dot_product(pt, vt), 1.0e-11_dp)
    end subroutine check_thomson

    subroutine check_ece()
        integer :: n, nt
        real(dp), allocatable :: y(:), yp(:), ym(:), dy(:), vt(:), ybar(:), ps(:), pt(:)
        real(dp), allocatable :: dydf(:)
        type(ece_resonance_t), allocatable :: res(:)
        type(ece_channel_t), allocatable :: shifted(:)
        real(dp) :: h

        n = size(ece)
        nt = te%n_params()
        allocate(y(n), yp(n), ym(n), dy(n), vt(nt), ybar(n), dydf(n))
        call random_direction(vt)
        call ece_resonances(state, ece, res)
        call check(all(res%found), 'ECE resonances found')
        call ece_jvp(state, res, te, vstate, vt, dy)
        call ece_frequency_jacobian(res, te, dydf)
        h = 1.0e-6_dp
        call perturb(h, te, t0, vt)
        call ece_resonances(state, ece, res)
        call ece_forward(res, te, yp)
        call perturb(-h, te, t0, vt)
        call ece_resonances(state, ece, res)
        call ece_forward(res, te, ym)
        call perturb(0.0_dp, te, t0, vt)
        call compare('ECE JVP vs FD', dy, (yp - ym) / (2.0_dp * h))
        shifted = ece
        shifted%frequency = ece%frequency * (1.0_dp + h)
        call ece_resonances(state, shifted, res)
        call ece_forward(res, te, yp)
        shifted%frequency = ece%frequency * (1.0_dp - h)
        call ece_resonances(state, shifted, res)
        call ece_forward(res, te, ym)
        call compare('ECE frequency derivative vs FD', dydf * ece%frequency, &
            (yp - ym) / (2.0_dp * h))
        call ece_resonances(state, ece, res)
        call random_direction(ybar)
        allocate(ps(size(p0)), pt(nt))
        ps = 0.0_dp
        pt = 0.0_dp
        call ece_vjp(state, res, te, ybar, ps, pt)
        call close('ECE VJP transposition', dot_product(ybar, dy), &
            dot_product(ps, vstate) + dot_product(pt, vt), 1.0e-11_dp)
    end subroutine check_ece

    subroutine check_interferometer()
        integer :: n, nn
        real(dp), allocatable :: yp(:), ym(:), dy(:), vn(:), ybar(:), ps(:), pn(:)
        type(sample_set_t) :: set
        real(dp) :: h

        n = size(chords)
        nn = ne%n_params()
        allocate(yp(n), ym(n), dy(n), vn(nn), ybar(n))
        call random_direction(vn)
        call chord_samples(state, chords, set)
        call check(size(set%seg_a) >= n, 'every chord crosses the plasma')
        call interferometer_jvp(state, set, n, ne, vstate, vn, dy)
        h = 1.0e-6_dp
        call perturb(h, ne, n0, vn)
        call chord_samples(state, chords, set)
        call interferometer_forward(set, n, ne, yp)
        call perturb(-h, ne, n0, vn)
        call chord_samples(state, chords, set)
        call interferometer_forward(set, n, ne, ym)
        call perturb(0.0_dp, ne, n0, vn)
        call compare('interferometer JVP vs FD', dy, (yp - ym) / (2.0_dp * h))
        call chord_samples(state, chords, set)
        call random_direction(ybar)
        allocate(ps(size(p0)), pn(nn))
        ps = 0.0_dp
        pn = 0.0_dp
        call interferometer_vjp(state, set, ne, ybar, ps, pn)
        call close('interferometer VJP transposition', dot_product(ybar, dy), &
            dot_product(ps, vstate) + dot_product(pn, vn), 1.0e-11_dp)
    end subroutine check_interferometer

    subroutine check_xics()
        integer :: n, ne_, nt
        real(dp), allocatable :: yp(:), ym(:), dy(:), ve(:), vt(:), ybar(:)
        real(dp), allocatable :: ps(:), pe(:), pt(:)
        type(sample_set_t) :: set
        real(dp) :: h

        n = size(chords)
        ne_ = emis%n_params()
        nt = ti%n_params()
        allocate(yp(2 * n), ym(2 * n), dy(2 * n), ve(ne_), vt(nt), ybar(2 * n))
        call random_direction(ve)
        call random_direction(vt)
        call chord_samples(state, chords, set)
        call xics_jvp(state, set, n, emis, ti, vstate, ve, vt, dy)
        h = 1.0e-6_dp
        call perturb(h, emis, e0, ve, ti, ti0, vt)
        call chord_samples(state, chords, set)
        call xics_forward(set, n, emis, ti, yp)
        call perturb(-h, emis, e0, ve, ti, ti0, vt)
        call chord_samples(state, chords, set)
        call xics_forward(set, n, emis, ti, ym)
        call perturb(0.0_dp, emis, e0, ve, ti, ti0, vt)
        call compare('XICS JVP vs FD', dy, (yp - ym) / (2.0_dp * h))
        call chord_samples(state, chords, set)
        call random_direction(ybar)
        allocate(ps(size(p0)), pe(ne_), pt(nt))
        ps = 0.0_dp
        pe = 0.0_dp
        pt = 0.0_dp
        call xics_vjp(state, set, n, emis, ti, ybar, ps, pe, pt)
        call close('XICS VJP transposition', dot_product(ybar, dy), &
            dot_product(ps, vstate) + dot_product(pe, ve) + dot_product(pt, vt), 1.0e-11_dp)
    end subroutine check_xics

    subroutine perturb(s, p1, s1, v1, p2, s2, v2)
        !! State and profiles at p0 + s v.
        real(dp), intent(in) :: s, s1, v1(:)
        type(node_profile_t), intent(inout) :: p1
        type(node_profile_t), intent(inout), optional :: p2
        real(dp), intent(in), optional :: s2, v2(:)
        call state%set_params(p0 + s * vstate)
        call p1%set_params(profile_values(s1) + s * v1)
        if (present(p2)) call p2%set_params(profile_values(s2) + s * v2)
    end subroutine perturb

    function profile_values(scale) result(v)
        real(dp), intent(in) :: scale
        real(dp) :: v(11)
        integer :: j
        v = scale * (1.0_dp - [(0.01_dp * j * j, j = 0, 10)])
    end function profile_values

    subroutine random_direction(v)
        real(dp), intent(out) :: v(:)
        integer :: k
        do k = 1, size(v)
            v(k) = sin(1.7_dp * k + 0.3_dp) + 0.5_dp * cos(3.1_dp * k * k)
        end do
    end subroutine random_direction

    subroutine compare(name, analytic, fd)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: analytic(:), fd(:)
        real(dp) :: err
        err = maxval(abs(analytic - fd)) / max(maxval(abs(fd)), 1.0e-30_dp)
        print '(A,": max rel error ",ES10.3)', name, err
        if (.not. (err <= tol_fd)) then
            write(error_unit, '(A,ES12.4)') 'FAIL '//name//': ', err
            failures = failures + 1
        end if
    end subroutine compare

    subroutine close(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b, tol
        if (.not. (abs(a - b) <= tol * max(1.0_dp, abs(b)))) then
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
end program test_profile_diagnostics
