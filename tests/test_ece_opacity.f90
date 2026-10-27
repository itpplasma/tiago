program test_ece_opacity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_ece_opacity
    use tiago_ece_opacity_transfer, only: thermal_ece_transfer, ece_opacity_transfer_t
    use tiago_ece_opacity_x1_flr, only: ece_opacity_x1_flr
    use fortnum_special_dnestrovskii, only: dnestrovskii
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), c = 299792458.0_dp, rest = 510998.95069_dp, &
                        e = 1.602176634e-19_dp, me = 9.1093837139e-31_dp, eps0 = 8.8541878188e-12_dp
    real(dp) :: maxfd = 0.0_dp, maxkirch = 0.0_dp
    integer :: mode, harmonic, regime
    do mode = ECE_O_MODE, ECE_X_MODE
        do harmonic = 1, 2
            regime = ECE_X1_FINITE
            call derivatives(mode, harmonic, regime)
            if (mode == ECE_X_MODE .and. harmonic == 1) &
                call derivatives(mode, harmonic, ECE_X1_TENUOUS)
        end do
    end do
    call slab(ECE_O_MODE, 1)
    call slab(ECE_O_MODE, 2)
    call slab(ECE_X_MODE, 1)
    call slab(ECE_X_MODE, 2)
    call limits()
    call coupled_transfer()
    print '(A,ES12.4,A,ES12.4)', 'maximum opacity FD relative error=', maxfd, &
        ' Kirchhoff relative error=', maxkirch
    call published_flr()
contains
    subroutine require(ok, why)
        logical, intent(in) :: ok
        character(*), intent(in) :: why
        if (.not. ok) then
            print *, why
            error stop 1
        end if
    end subroutine
    subroutine evaluate(p, mode, harmonic, regime, y, jac)
        real(dp), intent(in) :: p(4)
        integer, intent(in) :: mode, harmonic, regime
        real(dp), intent(out) :: y(3), jac(3, 4)
        integer :: status
        call ece_opacity(p(1), p(2), p(3), p(4), mode, harmonic, &
                         y(1), y(2), y(3), jac, status, regime)
        if (status /= 0) print *, 'rejected point', p, mode, harmonic, regime, status
        call require(status == 0, 'opacity rejected valid test point')
    end subroutine
    function parameters(beta, r, x, harmonic) result(p)
        real(dp), intent(in) :: beta, r, x
        integer, intent(in) :: harmonic
        real(dp) :: p(4)
        p(1) = beta*rest
        p(3) = 2.0_dp
        p(2) = r*eps0*p(3)**2/me
        p(4) = harmonic*e*p(3)/(2*pi*me*(1 + beta*x))
    end function
    subroutine derivatives(mode, harmonic, regime)
        integer, intent(in) :: mode, harmonic, regime
        real(dp) :: p(4), plus(4), minus(4), scales(4), y(3), yp(3), ym(3), jac(3, 4), &
                    unused(3, 4), fd(3), error, h, bb
        integer :: k, wing
        bb = 0.001_dp
        if (regime == ECE_X1_TENUOUS) then
            p = parameters(bb, 1.0e-5_dp, 3.1_dp, harmonic)
        else if (mode == ECE_X_MODE .and. harmonic == 1) then
            p = parameters(bb, 0.3_dp, 3.1_dp, harmonic)
        else
            p = parameters(bb, 0.05_dp, 3.1_dp, harmonic)
        end if
        do wing = 1, 2
            if (wing == 2) p(4) = harmonic*e*p(3)/(2*pi*me*(1 - bb*3.1_dp))
            scales = [p(1), p(2), p(3)*bb, p(4)*bb]
            call evaluate(p, mode, harmonic, regime, y, jac)
            do k = 1, 4
                h = 1.0e-5_dp*scales(k)
                plus = p
                minus = p
                plus(k) = plus(k) + h
                minus(k) = minus(k) - h
                call evaluate(plus, mode, harmonic, regime, yp, unused)
                call evaluate(minus, mode, harmonic, regime, ym, unused)
                fd = (yp - ym)/(2*h)
                error = maxval(abs(fd - jac(:, k))*scales(k)/max(abs(y), 1.0e-40_dp))
                maxfd = max(maxfd, error)
                if (error > 2.0e-7_dp) print *, 'FD', mode, harmonic, regime, k, error
                call require(error < 2.0e-7_dp, 'opacity finite differences')
            end do
            ! Per Hz and per single polarization: B_nu^vac = e Te f^2/c^2.
            if (wing == 1) then
                error = abs(y(2)/y(1)/(y(3)*e*p(1)*p(4)**2/c**2) - 1)
                maxkirch = max(maxkirch, error)
                call require(error < 2.0e-14_dp, 'Rayleigh-Jeans Kirchhoff law per Hz')
            else
                call require(y(1) == 0.0_dp .and. y(2) == 0.0_dp .and. y(3) > 0.0_dp, &
                             'blue-side opacity and propagation index')
            end if
            print '(A,3I3,A,ES12.4)', 'mode harmonic regime=', mode, harmonic, regime, &
                ' alpha [1/m]=', y(1)
        end do
    end subroutine
    subroutine published_flr()
        ! Bornatici (1983), Eq. (3.1.30), independent defining PV quadrature.
        ! Table III disagrees with this equation; see docs/ece-opacity.md.
        real(dp) :: p(4), v(1), row(1, 8), error, reference(3)
        complex(dp) :: f5, f7, d
        character(1024) :: path
        integer :: status, unit, ios, count
        call get_command_argument(1, path)
        open (newunit=unit, file=trim(path), status='old')
        error = 0
        count = 0
        do
            read (unit, *, iostat=ios) reference
            if (ios /= 0) exit
            p = parameters(.001_dp, reference(2), reference(1), 1)
            call dnestrovskii(5, -reference(1), f5, d, status)
            call require(status == 0, 'F5 oracle domain')
            call dnestrovskii(7, -reference(1), f7, d, status)
            call require(status == 0, 'F7 oracle domain')
            call ece_opacity_x1_flr(p(1), p(2), p(3), p(4), &
                                    real(f5, dp), aimag(f5), real(f7, dp), aimag(f7), v, row)
            error = max(error, abs(v(1) - reference(3)))
            count = count + 1
        end do
        close (unit)
        print '(A,I0,A,ES12.4)', 'Bornatici Eq. 3.1.30 PV samples=', count, &
            ' max absolute error=', error
        call require(count == 48, 'missing finite-Larmor-radius references')
        call require(error < 2.0e-13_dp, 'independent finite-Larmor-radius correction')
    end subroutine
    subroutine slab(mode, harmonic)
        integer, intent(in) :: mode, harmonic
        integer, parameter :: n = 65536
        real(dp), parameter :: r0 = 1.0e-10_dp, lb = 1.0_dp, dx = 50.0_dp/n
        real(dp) :: p(4), bb, bres, ne, tau, reference, error, y(3), jac(3, 4), x
        integer :: i, k, regime
        ! Bornatici (1983), Table IV, p.1192: tenuous perpendicular O1/X2.
        ! lambda_0=2*pi*c/omega_c; v_t^2=Te/m_e (Te in joules).
        ! tau_O1=(pi/2)*r*beta*omega_c*LB/c; tau_X2=2*tau_O1.
        regime = ECE_X1_FINITE
        if (mode == ECE_X_MODE .and. harmonic == 1) regime = ECE_X1_TENUOUS
        bres = 2.0_dp
        ne = r0*eps0*bres**2/me
        do k = 1, 3
            bb = 10.0_dp**(-5 - k)
            tau = 0
            do i = 1, n
                x = (i - .5_dp)*dx
                p = [bb*rest, ne, bres*(1 + bb*x), harmonic*e*bres/(2*pi*me)]
                call evaluate(p, mode, harmonic, regime, y, jac)
                tau = tau + y(1)*lb*bb*dx
            end do
            reference = pi/2*r0*bb*(e*bres/me)*lb/c
            if (mode == ECE_X_MODE .and. harmonic == 2) reference = 2*reference
            if (mode == ECE_O_MODE .and. harmonic == 2) reference = 2*bb*reference
            if (mode == ECE_X_MODE .and. harmonic == 1) reference = reference/bb
            error = abs(tau/reference - 1)
            print '(A,2I2,A,ES12.4,A,ES12.4)', 'slab mode/harmonic=', mode, harmonic, &
                ' Te [eV]=', bb*rest, ' asymptotic tau relative error=', error
            call require(error < 8*bb + 2.0e-9_dp, 'published slab optical depth')
        end do
    end subroutine
    subroutine limits()
        type(ece_opacity_transfer_t) :: trace
        real(dp) :: p(4), y(3), jac(3, 4), a, j, n2, thin_error, thick_error
        integer :: status, k, m, h, reg
        do m = ECE_O_MODE, ECE_X_MODE
            do h = 1, 2
                reg = ECE_X1_FINITE
                if (m == ECE_X_MODE .and. h == 1) reg = ECE_X1_TENUOUS
                p = parameters(.001_dp, 0.0_dp, 3.0_dp, h)
                call evaluate(p, m, h, reg, y, jac)
                call require(y(1) == 0.0_dp .and. y(2) == 0.0_dp .and. &
                             abs(y(3) - 1) < 2.0e-14_dp, 'all-mode vacuum')
                p = parameters(.001_dp, 1.0e-10_dp, 3.0_dp, h)
                p(1) = p(1)*1.0e-4_dp
                call evaluate(p, m, h, reg, y, jac)
                call require(y(1) < 1.0e-100_dp, 'all-mode cold wing')
            end do
        end do
        p = parameters(.001_dp, .3_dp, 3.0_dp, 1)
        p(1) = p(1)*1.0e-4_dp
        call evaluate(p, ECE_X_MODE, 1, ECE_X1_FINITE, y, jac)
        call require(y(1) < 1.0e-100_dp, 'finite-density X1 cold wing')
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        p(2) = 0.0_dp
        call evaluate(p, ECE_X_MODE, 2, ECE_X1_FINITE, y, jac)
        call require(y(1) == 0.0_dp .and. y(2) == 0.0_dp, 'vacuum emission and absorption')
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        p(4) = 1.01_dp*2*e*p(3)/(2*pi*me)
        call evaluate(p, ECE_X_MODE, 2, ECE_X1_FINITE, y, jac)
        call require(y(1) == 0.0_dp .and. y(2) == 0.0_dp, 'no blue-side perpendicular resonance')
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        do k = 1, 4
            p(1) = p(1)/10
            call evaluate(p, ECE_X_MODE, 2, ECE_X1_FINITE, y, jac)
        end do
        call require(y(1) < 1.0e-100_dp, 'cold off-resonant absorption vanishes')
        call ece_opacity(-1.0_dp, p(2), p(3), p(4), ECE_O_MODE, 1, a, j, n2, jac, status)
        call require(status /= 0, 'negative temperature accepted')
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        p(4) = 0.7_dp*2*e*p(3)/(2*pi*me)
        call ece_opacity(p(1), p(2), p(3), p(4), ECE_X_MODE, 2, a, j, n2, jac, status)
        call require(status == ECE_OPACITY_DOMAIN, 'large detuning accepted')

        p = parameters(.001_dp, .002_dp, 3.0_dp, 1)
        call ece_opacity(p(1), p(2), p(3), p(4), ECE_X_MODE, 1, a, j, n2, jac, status)
        call require(status == ECE_OPACITY_DOMAIN, 'X1 transition regime accepted')
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        call evaluate(p, ECE_X_MODE, 2, ECE_X1_FINITE, y, jac)
        call thermal_ece_transfer([p(1)], [p(2)], [p(3)], p(4), ECE_X_MODE, 2, &
                                  [1.0e-10_dp], 0.0_dp, 0.0_dp, trace, status)
        call require(status == 0, 'thin thermal opacity transport')
        thin_error = abs(trace%temperature/(p(1)*y(1)*1.0e-10_dp) - 1)
        call require(thin_error < 1.0e-8_dp, &
                     'optically thin physical emission')
        call thermal_ece_transfer([p(1)], [p(2)], [p(3)], p(4), ECE_X_MODE, 2, &
                                  [10.0_dp], 0.0_dp, 0.0_dp, trace, status)
        call require(status == 0, 'thick thermal opacity transport')
        thick_error = abs(trace%temperature/p(1) - 1)
        call require(thick_error < 2.0e-14_dp, &
                     'optically thick thermal equilibrium')
        print '(A,ES12.4,A,ES12.4)', 'thin-limit relative error=', thin_error, &
            ' thick-limit relative error=', thick_error
        print *, 'vacuum, blue side, cold wing, thin/thick and excluded domains: PASS'
    end subroutine
    subroutine coupled_transfer()
        integer, parameter :: n = 8
        type(ece_opacity_transfer_t) :: tape, plus, minus
        real(dp) :: p(4), te(n), ne(n), b(n), length(n), dt(n), dn(n), db(n), dl(n), &
                    tb(n), nb(n), bb(n), lb(n), fb, rb, wb, jvp, fd, dot, h, freq, rho, wall
        integer :: k, status
        p = parameters(.001_dp, .05_dp, 3.0_dp, 2)
        te = p(1)
        ne = p(2)
        b = p(3)
        do k = 1, n
            te(k) = te(k)*(1 + .02_dp*k)
        end do
        freq = p(4)
        length = .002_dp
        rho = .2_dp
        wall = 5.0_dp
        call thermal_ece_transfer(te, ne, b, freq, ECE_X_MODE, 2, length, rho, wall, tape, status)
        call require(status == 0 .and. tape%valid, 'coupled transfer failed')
        dt = .1_dp*te
        dn = .2_dp*ne
        db = .0001_dp*b
        dl = .1_dp*length
        call tape%jvp(dt, dn, db, freq*.0001_dp, dl, .01_dp, .3_dp, jvp)
        call tape%vjp(1.0_dp, tb, nb, bb, fb, lb, rb, wb)
        dot = sum(tb*dt + nb*dn + bb*db + lb*dl) + fb*freq*.0001_dp + rb*.01_dp + wb*.3_dp
        call require(abs(dot - jvp) < 1.0e-12_dp*max(1.0_dp, abs(jvp)), 'coupled transpose')
        h = 1.0e-5_dp
        call thermal_ece_transfer(te + h*dt, ne + h*dn, b + h*db, freq*(1 + h*.0001_dp), &
                         ECE_X_MODE, 2, length + h*dl, rho + h*.01_dp, wall + h*.3_dp, plus, status)
        call require(status == 0, 'positive transfer FD')
        call thermal_ece_transfer(te - h*dt, ne - h*dn, b - h*db, freq*(1 - h*.0001_dp), &
                        ECE_X_MODE, 2, length - h*dl, rho - h*.01_dp, wall - h*.3_dp, minus, status)
        call require(status == 0, 'negative transfer FD')
        fd = (plus%temperature - minus%temperature)/(2*h)
        print '(A,ES12.4,A,ES12.4)', 'coupled transfer FD relative error=', &
            abs(fd - jvp)/abs(jvp), ' transpose absolute error [eV]=', abs(dot - jvp)
        call require(abs(fd - jvp) < 2.0e-7_dp*abs(jvp), 'coupled transfer FD')
    end subroutine
end program
