module island_oracle_field
    !! Exactly integrable 2/1 island in toroidal geometry. With x = R - R0,
    !! y = Z, r^2 = x^2 + y^2, theta = atan2(y, x),
    !!     B = B0 e_phi + grad(psi) x grad(phi),
    !!     psi/B0 = iota0 r^2/2 + iota2 r^4/(4 a^2) + eps r^2 cos(2 theta - phi),
    !! a field line obeys dx/dphi = -d(psi/B0)/dy, dy/dphi = d(psi/B0)/dx, so
    !! the rotational transform is iota0 + iota2 r^2/a^2 and
    !!     K = psi/B0 - r^2/4
    !! is conserved exactly along every line. The separatrix is the level set
    !! of K through the X point, which gives the width in closed form.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_field_lines, only: line_field_t
    use tiago_plasma_field_source, only: field_hint_t
    implicit none
    real(dp), parameter :: r0 = 3.0_dp, b0 = 1.0_dp, iota0 = 0.4_dp
    real(dp), parameter :: iota2 = 0.2_dp, a = 0.5_dp, eps = 0.0005_dp

    type, extends(line_field_t) :: oracle_t
        real(dp) :: axis_wobble = 0.0_dp, angular_warp = 0.0_dp
    contains
        procedure :: sample => oracle_sample
    end type oracle_t

    type, extends(line_field_t) :: hinted_oracle_t
        type(oracle_t) :: field
    contains
        procedure :: sample => require_hint
        procedure :: sample_with_hint => hinted_sample
    end type hinted_oracle_t
contains
    subroutine require_hint(self, x, b, status)
        class(hinted_oracle_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: status

        call self%field%sample(x, b, status)
        ! Make a tracer that bypasses the hinted entry point observable.
        status = 1
    end subroutine require_hint

    subroutine hinted_sample(self, x, b, status, hint)
        class(hinted_oracle_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: status
        type(field_hint_t), intent(inout) :: hint
        real(dp) :: dx, dy, conserved

        b = 0
        status = 1
        dx = x(1) - r0
        dy = x(3)
        conserved = invariant(dx*dx + dy*dy, &
                              cos(2*atan2(dy, dx) - x(2)))
        if (hint%volume == 0) then
            ! Each return starts at an integer turn; later stages need its hint.
            if (abs(sin(x(2)/2)) > 1.0e-10_dp) return
            hint%u = [conserved, 0.0_dp, 0.0_dp]
            hint%volume = 1
        else
            ! Different starts have different invariants, even on one thread.
            if (abs(hint%u(1) - conserved) > 1.0e-5_dp) return
        end if
        call self%field%sample(x, b, status)
    end subroutine hinted_sample

    subroutine oracle_sample(self, x, b, status)
        class(oracle_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: status
        real(dp) :: dx, dy, g1, px, py, radius, theta, alpha, radial, angular, metric
        integer :: k

        dx = x(1) - r0 - self%axis_wobble*cos(x(2))
        dy = x(3) - self%axis_wobble*sin(x(2))
        radius = norm2([dx, dy])
        theta = atan2(dy, dx)
        alpha = theta
        do k = 1, merge(12, 0, self%angular_warp /= 0)
            alpha = alpha - (alpha + self%angular_warp*cos(alpha) - theta) &
                            /(1 - self%angular_warp*sin(alpha))
        end do
        dx = radius*cos(alpha)
        dy = radius*sin(alpha)
        g1 = iota0/2 + iota2*(dx*dx + dy*dy)/(2*a*a)
        px = 2*g1*dx + eps*(2*dx*cos(x(2)) + 2*dy*sin(x(2)))
        py = 2*g1*dy + eps*(-2*dy*cos(x(2)) + 2*dx*sin(x(2)))
        metric = 1 - self%angular_warp*sin(alpha)
        radial = 0
        angular = 0
        if (radius > tiny(1.0_dp)) then
            radial = (-dx*py + dy*px)/radius
            angular = (dx*px + dy*py)/(radius*radius)
        end if
        b(2) = b0/metric
        b(1) = b(2)*(radial*cos(theta) - radius*sin(theta)*metric*angular &
                     - self%axis_wobble*sin(x(2)))/x(1)
        b(3) = b(2)*(radial*sin(theta) + radius*cos(theta)*metric*angular &
                     + self%axis_wobble*cos(x(2)))/x(1)
        status = 0
    end subroutine oracle_sample

    pure real(dp) function invariant(rho, c)
        !! K at r^2 = rho and cos(2 theta - phi) = c.
        real(dp), intent(in) :: rho, c
        invariant = iota0*rho/2 + iota2*rho*rho/(4*a*a) + eps*rho*c - rho/4
    end function invariant
end module island_oracle_field

program test_island_measure
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use island_oracle_field
    use tiago_field_lines, only: cartesian_line_field_t
    use tiago_coil_kernels, only: build_coil_set
    use tiago_island_measure, only: trace_options_t, island_point_t, &
        island_width_t, periodic_point, island_width, poincare_section
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(oracle_t) :: field
    type(trace_options_t) :: popt, wopt
    type(island_point_t) :: ax, op, xp
    type(island_width_t) :: wd
    real(dp) :: rho_o, rho_x, kx, qa, qb, disc, rho1, rho2, exact, dist, dir(2)
    real(dp) :: pts(2, 20, 1), kdrift
    integer :: failures, st(1), k

    failures = 0
    rho_o = (0.5_dp + 2*eps - iota0)*a*a/iota2
    rho_x = (0.5_dp - 2*eps - iota0)*a*a/iota2
    kx = invariant(rho_x, 1.0_dp)
    ! K(rho, cos u = -1) = K_X: iota2/(4a^2) rho^2 + (iota0/2 - 1/4 - eps) rho - K_X = 0
    qa = iota2/(4*a*a)
    qb = iota0/2 - 0.25_dp - eps
    disc = sqrt(qb*qb + 4*qa*kx)
    rho1 = (-qb - disc)/(2*qa)
    rho2 = (-qb + disc)/(2*qa)
    exact = sqrt(rho2) - sqrt(rho1)

    popt%rtol = 1.0e-11_dp
    popt%atol = 1.0e-13_dp
    popt%max_step = 2*pi/256
    popt%map_delta = 1.0e-4_dp
    popt%root_tol = 1.0e-8_dp
    call periodic_point(field, [r0 + 0.01_dp, 0.005_dp], 0.0_dp, 1, &
                        [r0 + 0.01_dp, 0.0_dp], 1, 0, popt, ax)
    call check('axis converges', ax%status == 0)
    call check('axis at (R0, 0)', norm2(ax%x - [r0, 0.0_dp]) < 1.0e-7_dp)
    call periodic_point(field, [r0, sqrt(rho_o) + 0.004_dp], 0.0_dp, 2, ax%x, 2, 1, &
                        popt, op)
    call periodic_point(field, [r0 + sqrt(rho_x) + 0.001_dp, 0.001_dp], 0.0_dp, 2, &
                        ax%x, 2, 1, popt, xp)
    call check('O point converges', op%status == 0)
    call check('X point converges', xp%status == 0)
    call check('O point position', norm2(op%x - [r0, sqrt(rho_o)]) < 1.0e-7_dp)
    call check('X point position', norm2(xp%x - [r0 + sqrt(rho_x), 0.0_dp]) < 1.0e-7_dp)
    call check('O residue in (0, 1)', op%residue > 0 .and. op%residue < 1)
    call check('X residue negative', xp%residue < 0)
    call check('area preserving', abs(op%determinant - 1) < 1.0e-6_dp)
    ! Helical phase 2 theta - phi of the O point is pi.
    call check('O phase', abs(cos(op%phase) + 1) < 1.0e-10_dp)
    call check('X phase', abs(cos(xp%phase) - 1) < 1.0e-10_dp)

    wopt%rtol = 1.0e-10_dp
    wopt%atol = 1.0e-12_dp
    wopt%max_step = 2*pi/128
    dist = norm2(op%x - ax%x)
    dir = (op%x - ax%x)/dist
    call island_width(field, [dist - 0.08_dp, dist + 0.08_dp], dist, ax%x, dir, &
                      0.0_dp, 256.0_dp, 2.0_dp, 2, wopt, wd, scans=16)
    call check('width converges', wd%status == 0)
    call check('width against the exact separatrix', &
               abs(wd%metres - exact) < 1.0e-3_dp*exact)
    call check('edges against the exact separatrix', &
               abs(wd%edge(1) - sqrt(rho1)) + abs(wd%edge(2) - sqrt(rho2)) &
               < 2.0e-3_dp*exact)

    ! A Poincare line keeps K: the tracer follows the exact invariant.
    call poincare_section(field, reshape([r0 + 0.3_dp, 0.0_dp], [2, 1]), 0.0_dp, &
                          1.0_dp, 20, pts, st, popt)
    kdrift = 0
    do k = 1, 20
        kdrift = max(kdrift, abs(k_at(pts(:, k, 1)) - k_at([r0 + 0.3_dp, 0.0_dp])))
    end do
    call check('Poincare section completes', st(1) == 0)
    call check('Poincare points conserve K', kdrift < 1.0e-9_dp)

    ! Half-turn sampling must advance phi in this phi-dependent field.
    call poincare_section(field, reshape([r0 + 0.3_dp, 0.0_dp], [2, 1]), 0.0_dp, &
                          0.5_dp, 20, pts, st, popt)
    kdrift = 0
    do k = 1, 20
        kdrift = max(kdrift, abs(k_at(pts(:, k, 1), k*pi) &
                                - k_at([r0 + 0.3_dp, 0.0_dp])))
    end do
    call check('half-turn section completes', st(1) == 0)
    call check('half-turn points conserve K', kdrift < 1.0e-9_dp)

    ! A periodic translation gives a helical axis with the same separatrix.
    field%axis_wobble = 0.25_dp
    call island_width(field, [sqrt(rho_o) - 0.08_dp, sqrt(rho_o) + 0.08_dp], &
                      sqrt(rho_o), [r0 + field%axis_wobble, 0.0_dp], &
                      [0.0_dp, 1.0_dp], 0.0_dp, 256.0_dp, 2.0_dp, 2, wopt, wd, &
                      scans=16)
    call check('helical-axis width converges', wd%status == 0)
    call check('helical-axis width against exact separatrix', &
               abs(wd%metres - exact) < 1.0e-3_dp*exact)

    ! An annular angular deformation preserves the chosen chord crossings.
    field%axis_wobble = 0
    field%angular_warp = -0.5_dp
    call island_width(field, [sqrt(rho_o) - 0.08_dp, sqrt(rho_o) + 0.08_dp], &
                      sqrt(rho_o), [r0, 0.0_dp], [0.0_dp, 1.0_dp], &
                      0.0_dp, 256.0_dp, 2.0_dp, 2, wopt, wd, scans=16)
    call check('deformed-angle width converges', wd%status == 0)
    call check('deformed-angle width against exact separatrix', &
               abs(wd%metres - exact) < 1.0e-3_dp*exact)

    call wire_check()
    call hint_check()

    print '(a,es12.4,a,es12.4,a,es10.2)', 'width ', wd%metres, ' exact ', exact, &
        ' duration error ', wd%duration_error
    if (failures > 0) then
        print '(i0,a)', failures, ' island check(s) failed'
        error stop 1
    end if
    print '(a)', 'island measurement checks passed'
contains
    subroutine check(name, ok)
        character(len=*), intent(in) :: name
        logical, intent(in) :: ok
        if (ok) return
        failures = failures + 1
        print '(2a)', 'FAIL: ', name
    end subroutine check

    real(dp) function k_at(x, phi)
        real(dp), intent(in) :: x(2)
        real(dp), intent(in), optional :: phi
        real(dp) :: dx, dy, angle
        dx = x(1) - r0
        dy = x(2)
        angle = 0
        if (present(phi)) angle = phi
        k_at = invariant(dx*dx + dy*dy, cos(2*atan2(dy, dx) - angle))
    end function k_at

    subroutine wire_check()
        !! A straight segment on the z axis from -L to L carrying I gives the
        !! purely toroidal B_phi = mu0 I L/(2 pi R sqrt(L^2 + R^2)) at Z = 0.
        type(cartesian_line_field_t) :: wire
        real(dp), parameter :: length = 50.0_dp, current = 1.0e5_dp
        real(dp) :: b(3), expected
        integer :: status
        wire%coils = build_coil_set([0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], &
                                    [-length, length], [current, 0.0_dp])
        wire%has_coils = .true.
        call wire%sample([0.5_dp, 0.7_dp, 0.0_dp], b, status)
        expected = 4*pi*1.0e-7_dp*current*length/(2*pi*0.5_dp*sqrt(length**2 + 0.25_dp))
        call check('wire sample status', status == 0)
        call check('wire B_phi', abs(b(2) - expected) < 1.0e-12_dp*expected)
        call check('wire B_R and B_Z vanish', abs(b(1)) + abs(b(3)) < 1.0e-12_dp*expected)
        call wire%sample([-0.1_dp, 0.0_dp, 0.0_dp], b, status)
        call check('negative R rejected', status /= 0)
    end subroutine wire_check

    subroutine hint_check()
        type(oracle_t) :: cold
        type(hinted_oracle_t) :: hinted
        type(trace_options_t) :: opt
        real(dp) :: starts(2, 3), ordinary(2, 2, 3), cached(2, 2, 3)
        integer :: ordinary_status(3), cached_status(3), schedule

        starts(:, 1) = [r0 + 0.2_dp, 0.03_dp]
        starts(:, 2) = [r0 + 0.27_dp, -0.02_dp]
        starts(:, 3) = [r0 + 0.32_dp, 0.01_dp]
        opt%rtol = 1.0e-10_dp
        opt%atol = 1.0e-12_dp
        do schedule = 1, 2
            opt%fixed_steps = merge(128, 0, schedule == 1)
            call poincare_section(cold, starts, 0.0_dp, 1.0_dp, 2, ordinary, &
                                  ordinary_status, opt)
            call poincare_section(hinted, starts, 0.0_dp, 1.0_dp, 2, cached, &
                                  cached_status, opt)
            call check('default sampling remains compatible', &
                       all(ordinary_status == 0))
            call check('tracer preserves separate trajectory hints', &
                       all(cached_status == 0))
            call check('hinted fixed/adaptive Poincare sections agree', &
                       maxval(abs(cached - ordinary)) < 1.0e-10_dp)
        end do
    end subroutine hint_check
end program test_island_measure
