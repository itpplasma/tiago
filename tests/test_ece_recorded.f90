program test_ece_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ece_recorded_fixture, only: configure
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_radiometer, only: ece_t
    use tiago_ece_recorded, only: recorded_ece_t, record_ece
    use tiago_ece_opacity, only: ece_opacity, ECE_X_MODE, ECE_X1_TENUOUS
    use tiago_ece_opacity_transfer, only: ece_opacity_transfer_t
    use tiago_profile_diagnostic, only: plasma_sample_t, &
        check_diagnostic_derivatives
    use tiago_sight_geometry, only: gauss_legendre01
    use tiago_ece_doppler_x1, only: ece_doppler_x1
    use tiago_ece_doppler_x2, only: ece_doppler_x2
    use tiago_kinetic_profiles, only: P_NE
    use fortnum_status, only: fortnum_status_t
    implicit none
    call thick()
    call derivatives()
    call domains()
    call vacuum()
    call refinement()
    call angle_reference()
contains
    subroutine require(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (.not. ok) then
            print *, text
            error stop 1
        end if
    end subroutine

    subroutine thick()
        type(state_medium_t) :: medium
        type(ece_t) :: diagnostic, prepared
        type(recorded_ece_t) :: tape
        type(ece_opacity_transfer_t) :: transport(1)
        type(plasma_sample_t) :: sample
        type(fortnum_status_t) :: status
        real(dp) :: temperature(1), error
        integer :: code
        call configure(medium, diagnostic, 1, 192, 100.0_dp, .false.)
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        if (status%code /= 0) print *, status%msg
        call require(status%code == 0, 'thick full-chain evaluation')
        call tape%sampled(prepared, sample)
        call prepared%brightness(sample, spread(0.0_dp, 1, 5), transport, &
            temperature, code)
        call require(code == 0, 'thick brightness')
        error = abs(temperature(1)/1000 - 1)
        print '(A,ES14.6,A,ES14.6)', 'thick optical depth=', &
            transport(1)%optical_depth, &
            ' T_rad/Te error=', error
        call require(transport(1)%optical_depth > 25, &
            'plasma must be optically thick')
        call require(error < 2.0e-11_dp, 'thick full-chain radiation temperature')
        call require(abs(tape%voltage(1) - 10.02_dp) < 2.0e-10_dp, &
            'thick recorded voltage')
    end subroutine

    subroutine derivatives()
        type(state_medium_t) :: medium, plus, minus
        type(ece_t) :: diagnostic, prepared
        type(recorded_ece_t) :: tape, tp, tm
        type(plasma_sample_t) :: sample
        type(fortnum_status_t) :: status
        real(dp), allocatable :: p(:), pd(:), pb(:), eta(:), ed(:), eb(:)
        real(dp) :: yd(2), vd(2), yb(2), vb(2), fd(2), fvd(2), h, lhs, rhs, &
            err, verr, adj, je, ae
        integer :: i, ns, np
        call configure(medium, diagnostic, 2, 48, 1.0_dp, .true.)
        ns = medium%state%n_params()
        np = medium%n_params()
        allocate (p(ns), pd(np), pb(np), eta(10), ed(10), eb(10))
        call medium%state%get_params(p)
        pd(:ns) = [(sin(.37_dp*i)*max(abs(p(i)), .001_dp), i=1, ns)]
        pd(ns + 1:) = [(cos(.43_dp*i), i=1, np - ns)]
        ! Toroidal-vector perturbations keep the ray perpendicular in the
        ! nominal state; tiny finite-difference perturbations remain in-domain.
        eta = .1_dp
        ed = [(cos(.17_dp*i), i=1, 10)]
        call record_ece(diagnostic, medium, eta, tape, status)
        if (status%code /= 0) print *, status%msg
        call require(status%code == 0, 'nominal derivative composition')
        call tape%sampled(prepared, sample)
        call check_diagnostic_derivatives(prepared, sample, eta, 1.0e-7_dp, je, ae)
        print '(A,2ES14.6)', 'W0 finite difference and transpose errors=', je, ae
        call require(je < 2.0e-6_dp .and. ae < 2.0e-12_dp, 'W0 sampled derivatives')
        call tape%jvp(pd, ed, yd, vd, status)
        call require(status%code == 0, 'composed tangent')
        h = 1.0e-7_dp
        plus = medium
        minus = medium
        call plus%state%set_params(p + h*pd(:ns))
        call minus%state%set_params(p - h*pd(:ns))
        plus%z = medium%z + h*pd(ns + 1:)
        minus%z = medium%z - h*pd(ns + 1:)
        call record_ece(diagnostic, plus, eta + h*ed, tp, status)
        if (status%code /= 0) print *, status%msg
        call require(status%code == 0, 'plus forward')
        call record_ece(diagnostic, minus, eta - h*ed, tm, status)
        if (status%code /= 0) print *, status%msg
        call require(status%code == 0, 'minus forward')
        fd = (tp%voltage - tm%voltage)/(2*h)
        fvd = (tp%variance - tm%variance)/(2*h)
        err = maxval(abs(fd - yd))/maxval(abs(fd))
        verr = maxval(abs(fvd - vd))/maxval(abs(fvd))
        yb = [.3_dp, -.7_dp]
        vb = [1.0e3_dp, -2.0e3_dp]
        call tape%vjp(yb, vb, pb, eb, status)
        call require(status%code == 0, 'composed adjoint')
        lhs = dot_product(yb, yd) + dot_product(vb, vd)
        rhs = dot_product(pb, pd) + dot_product(eb, ed)
        adj = abs(lhs - rhs)/max(abs(lhs), abs(rhs))
        print '(A,3ES14.6)', 'full-chain voltage/noise FD and transpose errors=', err, &
            verr, adj
        call require(err < 2.0e-6_dp .and. verr < 2.0e-6_dp, &
            'full-chain finite differences')
        call require(adj < 2.0e-11_dp, 'full-chain transpose identity')
    end subroutine

    subroutine domains()
        type(state_medium_t) :: medium
        type(ece_t) :: diagnostic
        type(recorded_ece_t) :: tape
        type(fortnum_status_t) :: status
        call configure(medium, diagnostic, 1, 32, 1.0_dp, .false.)
        diagnostic%channels(1)%direction = [-1.0_dp, .01_dp, 0.0_dp]
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        call require(status%code /= 0 .and. .not. tape%valid, &
            'finite obliquity must fail')
        call require(index(status%msg, 'angle bound') > 0, 'clear obliquity error')
        diagnostic%channels(1)%direction = [-1.0_dp, 0.0_dp, 0.0_dp]
        diagnostic%channels(1)%harmonic = 1
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        call require(status%code /= 0 .and. .not. tape%valid, &
            'excluded X1 density regime')
        call require(index(status%msg, 'tenuous density') > 0, &
            'clear X1 density error')
        call medium%profiles%set_profile(P_NE, &
            1.0e15_dp*exp(20*(.5_dp - medium%profiles%greville)), &
            .true., sigma_f=0.1_dp, l_core=0.4_dp, l_edge=0.2_dp, &
            s0=0.5_dp, sw=0.2_dp)
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        call require(status%code /= 0 .and. .not. tape%valid, &
            'excluded X1 density after valid launch')
        call require(index(status%msg, 'tenuous density') > 0, &
            'clear local X1 density error')
        print *, 'angle and excluded-domain rejection: PASS'
    end subroutine

    subroutine vacuum()
        type(state_medium_t) :: medium
        type(ece_t) :: diagnostic
        type(recorded_ece_t) :: tape
        type(fortnum_status_t) :: status
        call configure(medium, diagnostic, 1, 32, 1.0_dp, .false.)
        medium%profiles%p(1)%unit = 0
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        call require(status%code == 0, 'vacuum full-chain evaluation')
        print '(A,ES14.6)', 'vacuum recorded-voltage error [V]=', &
            abs(tape%voltage(1) - .12_dp)
        call require(abs(tape%voltage(1) - .12_dp) < 2.0e-14_dp, &
            'vacuum wall brightness')
    end subroutine

    subroutine refinement()
        type(state_medium_t) :: medium
        type(ece_t) :: diagnostic
        type(recorded_ece_t) :: tape
        type(fortnum_status_t) :: status
        real(dp) :: reference, errors(3)
        integer :: k
        call configure(medium, diagnostic, 1, 1024, 1.0_dp, .true.)
        call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
        call require(status%code == 0, 'refined reference')
        reference = tape%voltage(1)
        do k = 1, 3
            diagnostic%channels(1)%nscan = 16*2**k
            diagnostic%first(2) = 1 + diagnostic%channels(1)%nscan
            call record_ece(diagnostic, medium, spread(0.0_dp, 1, 5), tape, status)
            call require(status%code == 0, 'refinement forward')
            errors(k) = abs(tape%voltage(1) - reference)
        end do
        print '(A,3ES14.6)', '32/64/128-cell voltage errors [V]=', errors
        call require(errors(2) < .35_dp*errors(1), '32 to 64 cell convergence')
        call require(errors(3) < .35_dp*errors(2), '64 to 128 cell convergence')
        call require(errors(3) < 1.0e-4_dp, 'resolved source integration')
    end subroutine

    subroutine angle_reference()
        ! Complete the square in the weakly relativistic resonance
        ! u_perp^2+u_parallel^2-2*a*u_parallel+z=0. The Maxwellian shell
        ! integral is the tenuous-X reference, with u=v/sqrt(2*Te/m_e).
        ! Bornatici et al., Nucl. Fusion 23 (1983), Secs. 2.3 and 3.1.1.
        real(dp), allocatable :: x(:), w(:)
        real(dp) :: mu, z, a, v, row(3), oblique, perpendicular, error
        real(dp) :: alpha, emission, index2, jac(3, 4), alpha0, line0, shape_error, &
            frequency, b, te, ne, wc
        integer :: harmonic, j, k, code
        call gauss_legendre01(64, x, w)
        error = 0
        shape_error = 0
        b = 2
        te = 1.0e-5_dp*510998.95069_dp
        ne = 1.0e-12_dp*8.8541878188e-12_dp*b**2/9.1093837139e-31_dp
        wc = 1.602176634e-19_dp/9.1093837139e-31_dp*b
        a = 1.0e-4_dp/sqrt(2.0_dp)
        do harmonic = 1, 2
            frequency = harmonic*wc/(2*acos(-1.0_dp)*(1 + 1.0e-5_dp))
            call ece_opacity(te, ne, b, frequency, ECE_X_MODE, harmonic, &
                alpha0, emission, index2, jac, code, ECE_X1_TENUOUS)
            call require(code == 0, 'perpendicular shell normalization')
            line0 = 2.0_dp/3*exp(-1.0_dp)
            if (harmonic == 2) line0 = 8.0_dp/15*exp(-1.0_dp)
            do k = 1, 100
                z = -.1_dp*real(k, dp)
                oblique = 0
                perpendicular = 0
                do j = 1, 64
                    mu = 2*x(j) - 1
                    if (harmonic == 1) then
                        call ece_doppler_x1(a, z, mu, v, row)
                        oblique = oblique + w(j)*v
                        call ece_doppler_x1(0.0_dp, z, mu, v, row)
                    else
                        call ece_doppler_x2(a, z, mu, v, row)
                        oblique = oblique + w(j)*v
                        call ece_doppler_x2(0.0_dp, z, mu, v, row)
                    end if
                    perpendicular = perpendicular + w(j)*v
                end do
                error = max(error, abs(oblique/perpendicular - 1))
                frequency = harmonic*wc/(2*acos(-1.0_dp)*(1 - 1.0e-5_dp*z))
                call ece_opacity(te, ne, b, frequency, ECE_X_MODE, harmonic, &
                    alpha, emission, index2, jac, code, ECE_X1_TENUOUS)
                call require(code == 0, 'perpendicular shell opacity')
                shape_error = max(shape_error, &
                    abs((alpha/alpha0)/(perpendicular/line0) - 1))
            end do
        end do
        print '(A,ES14.6)', 'tenuous-X finite-angle core relative error=', error
        call require(error < 1.2e-7_dp, 'perpendicular core approximation bound')
        print '(A,ES14.6)', 'opacity/shell line-shape relative error=', shape_error
        call require(shape_error < 1.0e-8_dp, 'independent shell line shape')
    end subroutine
end program
