module tiago_island_measure
    !! Poincare sections, periodic points and island width and phase on any
    !! line_field_t. Ported from MILCA's island diagnostics
    !! (itpplasma/milca src/mhd/island_diagnostics.f90 at f24dacdc), so that
    !! every level of the field hierarchy uses one measurement.
    !!
    !! A field line obeys dR/dphi = R B_R/B_phi, dZ/dphi = R B_Z/B_phi. The
    !! return map over `turns` toroidal turns is P; a periodic point solves
    !! P(x) = x. Greene's residue (2 - tr DP)/4 is in (0, 1) at an O point and
    !! negative at an X point. The helical phase of a point about the axis
    !! is m theta - n phi with theta = atan2(Z - Z_ax, R - R_ax).
    !!
    !! The width is the distance between the two libration edges along the
    !! straight chord from the axis through the O point: a line launched
    !! inside the island librates and one launched outside circulates. The
    !! helical phase is sampled where the line returns to the section, as m
    !! times the accumulated poloidal-angle change about the axis per return
    !! map (wrapped to (-pi, pi]); a line is classified as circulating once
    !! its poloidal angle spans a full turn within the trace. Sampling on the
    !! section keeps the angle
    !! meaningful for a helical magnetic axis, where an angle about one fixed
    !! point between sections is not (MILCA accumulates it continuously).
    !! Both edges are bisected. The result is a finite-time estimate, not an
    !! enclosure; the report carries its sensitivity to tolerance and trace
    !! duration.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_field_lines, only: line_field_t
    use tiago_plasma_field_source, only: field_hint_t
    use tiago_fieldline_slope, only: fieldline_slope
    use tiago_map_invariants, only: map_invariants
    use tiago_helical_phase, only: helical_phase
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    public :: trace_options_t, island_point_t, island_width_t
    public :: field_line_trace, map_stability, periodic_point
    public :: libration_span, libration_edge, island_width, poincare_section

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: trace_options_t
        real(dp) :: rtol = 1.0e-8_dp, atol = 1.0e-10_dp
        real(dp) :: max_step = 2*pi/128
        real(dp) :: map_delta = 1.0e-3_dp, root_tol = 1.0e-7_dp
        integer :: max_steps = 200000
        !! Nonzero selects a fixed Cash-Karp schedule of this many steps per
        !! turn, for explicit step ladders rather than tolerance control.
        integer :: fixed_steps = 0
    end type trace_options_t

    type :: island_point_t
        real(dp) :: x(2) = 0, closure = huge(1.0_dp)
        real(dp) :: residue = 0, determinant = 0, phase = 0
        real(dp) :: position_error = huge(1.0_dp), residue_error = huge(1.0_dp)
        real(dp) :: phase_error = huge(1.0_dp)
        !! 0: converged; 1: input or domain; 2: trace; 3: unresolved root.
        integer :: status = 1
    end type island_point_t

    type :: island_width_t
        !! edge: distances of the two libration edges from the chord origin
        !! [m]; metres: full width. trace_error: change under tightened
        !! tolerances; duration_error: change when the trace is doubled;
        !! scan_error: bisection interval left.
        real(dp) :: edge(2) = 0, metres = 0
        real(dp) :: trace_error = 0, duration_error = 0, scan_error = 0
        integer :: status = 1
    end type island_width_t

contains

    subroutine field_line_trace(field, start, phi0, turns, opt, finish, status)
        !! Follow a line from (R, Z) = start on phi0 over `turns` toroidal
        !! turns.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: start(2), phi0, turns
        type(trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: finish(2)
        integer, intent(out) :: status
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: work
        type(ode_solution_t) :: sol
        type(fortnum_status_t) :: stat
        type(field_hint_t) :: hint
        real(dp) :: k1(2), k2(2), k3(2), k4(2), k5(2), k6(2), tmp(2)
        real(dp) :: next(2), err(2), h, phi
        integer :: j, count, nfev, fault

        finish = start
        ! Each trajectory owns its seed, including concurrent return maps.
        hint = field_hint_t()
        status = 1
        fault = 0
        if (abs(turns) <= tiny(1.0_dp) .or. opt%max_step <= 0 .or. opt%rtol <= 0 &
            .or. opt%atol <= 0 .or. opt%fixed_steps < 0) return
        if (opt%fixed_steps > 0) then
            count = max(1, nint(abs(turns)*opt%fixed_steps))
            h = 2*pi*turns/count
            nfev = 0
            do j = 1, count
                phi = phi0 + h*(j - 1)
                call cash_karp_step(rhs, phi, finish, h, .false., k1, k2, &
                                    k3, k4, k5, k6, tmp, next, err, nfev)
                if (fault /= 0) exit
                finish = next
            end do
        else
            problem%rhs => rhs
            problem%t0 = phi0
            problem%t1 = phi0 + 2*pi*turns
            problem%y0 = start
            problem%rtol = opt%rtol
            problem%atol = opt%atol
            problem%hmax = opt%max_step
            problem%max_steps = opt%max_steps
            call ode_integrate(problem, work, sol, stat)
            if (stat%code /= 0) fault = 2
            if (fault == 0) finish = sol%y(:, size(sol%t))
        end if
        status = 0
        if (fault /= 0 .or. .not. all(ieee_is_finite(finish))) status = 2
    contains
        subroutine rhs(t, y, dy, ctx)
            real(dp), intent(in) :: t, y(:)
            real(dp), intent(out) :: dy(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: b(3), slope(2), jac(2, 4)
            integer :: st
            dy = 0
            if (fault /= 0) return
            call field%sample_with_hint([y(1), t, y(2)], b, st, hint)
            if (st /= 0 .or. abs(b(2)) <= tiny(1.0_dp)) then
                fault = 1
                return
            end if
            call fieldline_slope(y(1), b(1), b(2), b(3), slope, jac)
            dy = slope
        end subroutine rhs
    end subroutine field_line_trace

    subroutine libration_span(field, start, phi0, turns, map_turns, centre, m, &
                              opt, span, status)
        !! Range of the helical phase m (theta_k - theta_0) over the returns k
        !! to the section within `turns` toroidal turns, maps of `map_turns`
        !! turns each, theta about centre. Stops once theta spans a full turn.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: start(2), phi0, turns, map_turns, centre(2)
        integer, intent(in) :: m
        type(trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: span
        integer, intent(out) :: status
        real(dp) :: x(2), next(2), theta, old, angle, lo, hi
        integer :: k, maps

        span = 0
        status = 1
        if (map_turns <= 0 .or. m == 0 .or. &
            norm2(start - centre) <= tiny(1.0_dp)) return
        maps = max(1, nint(turns/map_turns))
        x = start
        old = atan2(x(2) - centre(2), x(1) - centre(1))
        angle = 0
        lo = 0
        hi = 0
        status = 0
        do k = 1, maps
            call field_line_trace(field, x, phi0 + 2*pi*(k - 1)*map_turns, &
                                  map_turns, opt, next, status)
            if (status /= 0) return
            x = next
            theta = atan2(x(2) - centre(2), x(1) - centre(1))
            angle = angle + atan2(sin(theta - old), cos(theta - old))
            old = theta
            lo = min(lo, m*angle)
            hi = max(hi, m*angle)
            if (hi - lo > 2*pi*abs(m)) exit
        end do
        span = hi - lo
    end subroutine libration_span

    subroutine map_stability(field, root, phi0, turns, opt, jac, closure, status)
        !! Central-difference Jacobian of the return map and its closure
        !! |P(root) - root| [m].
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: root(2), phi0, turns
        type(trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: jac(2, 2), closure
        integer, intent(out) :: status
        real(dp) :: finish(2), plus(2), minus(2), x(2)
        integer :: j, statuses(2)

        jac = 0
        closure = huge(1.0_dp)
        status = 1
        if (opt%map_delta <= 0) return
        call field_line_trace(field, root, phi0, turns, opt, finish, status)
        if (status /= 0) return
        closure = norm2(finish - root)
        statuses = 0
        !$omp parallel do default(shared) private(j, x, plus, minus) schedule(static)
        do j = 1, 2
            x = root
            x(j) = x(j) + opt%map_delta
            call field_line_trace(field, x, phi0, turns, opt, plus, statuses(j))
            if (statuses(j) /= 0) cycle
            x(j) = root(j) - opt%map_delta
            call field_line_trace(field, x, phi0, turns, opt, minus, statuses(j))
            if (statuses(j) /= 0) cycle
            jac(:, j) = (plus - minus)/(2*opt%map_delta)
        end do
        !$omp end parallel do
        status = maxval(statuses)
    end subroutine map_stability

    subroutine periodic_point(field, seed, phi0, turns, centre, m, n, opt, report)
        !! Newton (Powell hybrid) for P(x) = x from seed, solved twice: with
        !! opt and with tightened tolerances, half the step and half the
        !! difference step. The differences give the reported errors.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: seed(2), phi0, centre(2)
        integer, intent(in) :: turns, m, n
        type(trace_options_t), intent(in) :: opt
        type(island_point_t), intent(out) :: report
        type(trace_options_t) :: tight
        type(fortnum_status_t) :: stat
        real(dp) :: roots(2, 2), jac(2, 2), res(2), closure, residues(2), phases(2)
        real(dp) :: inv(2), pjac(1, 5), ijac(2, 4), phase(1)
        integer :: lev, st, fault

        report = island_point_t()
        tight = opt
        if (turns < 1 .or. opt%map_delta <= 0 .or. opt%root_tol <= 0) return
        roots(:, 1) = seed
        roots(:, 2) = seed
        do lev = 1, 2
            if (lev == 2) then
                tight%rtol = opt%rtol/8
                tight%atol = opt%atol/8
                tight%max_step = opt%max_step/2
                tight%map_delta = opt%map_delta/2
                tight%fixed_steps = 2*opt%fixed_steps
                roots(:, 2) = roots(:, 1)
            end if
            fault = 0
            call multiroot_hybrid(map_residual, 2, roots(:, lev), report%x, stat, &
                                  xtol=tight%root_tol/10, ftol=tight%root_tol, &
                                  max_iter=30)
            call map_residual(report%x, res, jac)
            if (fault /= 0 .or. norm2(res) > tight%root_tol) then
                report%status = 3
                return
            end if
            roots(:, lev) = report%x
            call map_stability(field, report%x, phi0, real(turns, dp), tight, jac, &
                               closure, st)
            if (st /= 0) then
                report%status = st
                return
            end if
            report%closure = closure
            call map_invariants(jac(1, 1), jac(1, 2), jac(2, 1), jac(2, 2), inv, ijac)
            residues(lev) = inv(1)
            report%determinant = inv(2)
            call helical_phase(report%x(1) - centre(1), report%x(2) - centre(2), &
                               real(m, dp), real(n, dp), phi0, phase, pjac)
            phases(lev) = phase(1)
        end do
        report%residue = residues(2)
        report%phase = atan2(sin(phases(2)), cos(phases(2)))
        report%position_error = norm2(roots(:, 2) - roots(:, 1)) + opt%root_tol
        report%residue_error = abs(residues(2) - residues(1))
        report%phase_error = abs(atan2(sin(phases(2) - phases(1)), &
                                       cos(phases(2) - phases(1))))
        if (norm2(report%x - centre) > report%position_error) then
            report%phase_error = report%phase_error + abs(m)*report%position_error &
                                 /(norm2(report%x - centre) - report%position_error)
        else
            report%phase_error = huge(1.0_dp)
        end if
        report%status = 0
        if (report%closure > opt%root_tol) report%status = 3
    contains
        subroutine map_residual(x, f, jacobian, ctx)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f(:), jacobian(:, :)
            class(*), intent(in), optional :: ctx
            real(dp) :: finish(2), defect
            integer :: j, ierr
            call field_line_trace(field, x, phi0, real(turns, dp), tight, finish, ierr)
            f = finish - x
            call map_stability(field, x, phi0, real(turns, dp), tight, jacobian, &
                               defect, st)
            if (ierr /= 0 .or. st /= 0) fault = 1
            do j = 1, 2
                jacobian(j, j) = jacobian(j, j) - 1
            end do
        end subroutine map_residual
    end subroutine periodic_point

    subroutine libration_edge(field, outside, inside, origin, direction, phi0, &
                              turns, map_turns, m, opt, edge, status, scans)
        !! Bisect between a circulating point at distance `outside` and a
        !! librating one at `inside` along origin + s direction (unit vector).
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: outside, inside, origin(2), direction(2)
        real(dp), intent(in) :: phi0, turns, map_turns
        integer, intent(in) :: m
        type(trace_options_t), intent(in) :: opt
        real(dp), intent(out) :: edge
        integer, intent(out) :: status
        integer, intent(in), optional :: scans
        real(dp) :: a, b, span
        integer :: j, count

        a = outside
        b = inside
        count = 12
        if (present(scans)) count = scans
        edge = a
        status = 1
        if (abs(norm2(direction) - 1) > 1.0e-12_dp .or. count < 1) return
        call libration_span(field, origin + a*direction, phi0, turns, map_turns, &
                            origin, m, opt, span, status)
        if (status /= 0) return
        status = 3
        if (span <= 2*pi*abs(m)) return
        call libration_span(field, origin + b*direction, phi0, turns, map_turns, &
                            origin, m, opt, span, status)
        if (status /= 0) return
        status = 3
        if (span > 2*pi*abs(m)) return
        do j = 1, count
            edge = (a + b)/2
            call libration_span(field, origin + edge*direction, phi0, turns, &
                                map_turns, origin, m, opt, span, status)
            if (status /= 0) return
            if (span > 2*pi*abs(m)) then
                a = edge
            else
                b = edge
            end if
        end do
        edge = (a + b)/2
    end subroutine libration_edge

    subroutine island_width(field, outside, inside, origin, direction, phi0, &
                            turns, map_turns, m, opt, report, scans)
        !! Both edges at three levels: opt, tightened tolerances, and
        !! tightened tolerances with twice the duration. The width is that of
        !! the tightened level; the others give trace and duration errors.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: outside(2), inside, origin(2), direction(2)
        real(dp), intent(in) :: phi0, turns, map_turns
        integer, intent(in) :: m
        type(trace_options_t), intent(in) :: opt
        type(island_width_t), intent(out) :: report
        integer, intent(in), optional :: scans
        type(trace_options_t) :: tights(3)
        real(dp) :: edges(2, 3), widths(3)
        integer :: lev, i, j, statuses(2, 3), count

        report = island_width_t()
        count = 12
        if (present(scans)) count = scans
        if (count < 1 .or. count > 24) return
        tights(1) = opt
        tights(2) = opt
        tights(2)%rtol = opt%rtol/8
        tights(2)%atol = opt%atol/8
        tights(2)%max_step = opt%max_step/2
        tights(2)%fixed_steps = 2*opt%fixed_steps
        tights(3) = tights(2)
        statuses = 0
        ! Six independent edge scans, longest first.
        !$omp parallel do default(shared) private(i, lev, j) schedule(dynamic, 1)
        do i = 1, 6
            lev = 3 - (i - 1)/2
            j = 1 + mod(i - 1, 2)
            call libration_edge(field, outside(j), inside, origin, direction, phi0, &
                                turns*merge(2, 1, lev == 3), map_turns, m, tights(lev), &
                                edges(j, lev), statuses(j, lev), scans=count)
        end do
        !$omp end parallel do
        do lev = 1, 3
            if (any(statuses(:, lev) /= 0)) then
                report%status = maxval(statuses(:, lev))
                return
            end if
            widths(lev) = abs(edges(2, lev) - edges(1, lev))
        end do
        report%edge = edges(:, 2)
        report%metres = widths(2)
        report%trace_error = abs(widths(2) - widths(1))
        report%duration_error = abs(widths(3) - widths(2))
        report%scan_error = sum(abs(outside - inside))/2**count
        report%status = 0
    end subroutine island_width

    subroutine poincare_section(field, starts, phi0, turns, crossings, points, &
                                statuses, opt)
        !! points(:, k, j): (R, Z) of start j after k maps of `turns` turns.
        !! A failed line keeps its last valid point and a nonzero status.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: starts(:, :), phi0, turns
        integer, intent(in) :: crossings
        real(dp), intent(out) :: points(:, :, :)
        integer, intent(out) :: statuses(:)
        type(trace_options_t), intent(in) :: opt
        real(dp) :: x(2), next(2)
        integer :: j, k

        statuses = 1
        points = 0
        if (size(points, 1) /= 2 .or. size(points, 2) /= crossings &
            .or. size(points, 3) /= size(starts, 2)) return
        !$omp parallel do default(shared) private(j, k, x, next) schedule(dynamic, 1)
        do j = 1, size(starts, 2)
            x = starts(:, j)
            statuses(j) = 0
            do k = 1, crossings
                call field_line_trace(field, x, phi0 + 2*pi*(k - 1)*turns, &
                                      turns, opt, next, statuses(j))
                if (statuses(j) /= 0) exit
                x = next
                points(:, k, j) = x
            end do
            if (statuses(j) /= 0) points(:, k:, j) = spread(x, 2, crossings - k + 1)
        end do
        !$omp end parallel do
    end subroutine poincare_section
end module tiago_island_measure
