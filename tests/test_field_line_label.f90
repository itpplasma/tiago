module test_helical_island_field
    !! Analytic field with a 2/1 island on concentric circles about (R0, 0).
    !! With r, theta polar about the axis and p = r^2 / 2, the lines follow
    !! the Hamiltonian H(theta, p, phi) = iota0 p + iota1 p^2 / 2
    !! + eps p cos(2 theta - phi):
    !!     dtheta/dphi = iota0 + iota1 p + eps cos(xi),
    !!     dp/dphi     = 2 eps p sin(xi),      xi = 2 theta - phi,
    !! so K = H - p / 2 is constant along each line; iota = 1/2 at
    !! p = (1/2 - iota0) / iota1. B_phi = B0 R0 / R and (B_R, B_Z) follow from
    !! dR/dphi = R B_R / B_phi = dx/dphi, dZ/dphi = R B_Z / B_phi = dy/dphi.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_field_line_label, only: line_field_t
    implicit none
    private

    type, extends(line_field_t), public :: helical_field_t
        real(dp) :: r0 = 3.0_dp, b0 = 2.0_dp, iota0 = 0.3_dp, iota1 = 0.8_dp
        real(dp) :: eps = 0.002_dp, a = 1.0_dp
    contains
        procedure :: eval => helical_eval
        procedure :: invariant
    end type helical_field_t

contains

    subroutine helical_eval(self, r, phi, z, b, inside)
        class(helical_field_t), intent(in) :: self
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: b(3)
        logical, intent(out) :: inside
        real(dp) :: x, y, rr, th, p, xi, dth, dpp, dr, dx, dy

        x = r - self%r0
        y = z
        rr = hypot(x, y)
        inside = rr < 1.2_dp * self%a
        b = 0.0_dp
        if (.not. inside) return
        th = atan2(y, x)
        p = 0.5_dp * rr**2
        xi = 2.0_dp * th - phi
        dth = self%iota0 + self%iota1 * p + self%eps * cos(xi)
        dpp = 2.0_dp * self%eps * p * sin(xi)
        dr = 0.0_dp
        if (rr > 0.0_dp) dr = dpp / rr
        dx = dr * cos(th) - rr * sin(th) * dth
        dy = dr * sin(th) + rr * cos(th) * dth
        b(2) = self%b0 * self%r0 / r
        b(1) = b(2) * dx / r
        b(3) = b(2) * dy / r
    end subroutine helical_eval

    pure real(dp) function invariant(self, p, xi)
        class(helical_field_t), intent(in) :: self
        real(dp), intent(in) :: p, xi
        invariant = (self%iota0 - 0.5_dp) * p + 0.5_dp * self%iota1 * p**2 &
            + self%eps * p * cos(xi)
    end function invariant

end module test_helical_island_field

program test_field_line_label
    !! Field-line labels (tiago_field_line_label) against independent oracles.
    !! Usage: test_field_line_label <workdir>
    !! 1. Without the perturbation every line stays on its circle: the label
    !!    equals s = r^2 / a^2 of the start point to the chart tolerance.
    !! 2. With the 2/1 island, the label of a circulating line equals the
    !!    invariant-measure average of s on its invariant curve, computed by
    !!    quadrature in the helical angle from K = const alone (no tracing).
    !! 3. Inside the island, the two lines started at mirror points of the
    !!    O point on the same invariant curve get one label, which lies closer
    !!    to the resonant s than their start points (flattening).
    !! 4. The same field tabulated on a grid (grid_field_t) reproduces the
    !!    analytic field to the tricubic interpolation error and its labels.
    !! 5. An mgrid file with two coil groups, read with currents (2, -0.5),
    !!    reproduces 2 B_1 - 0.5 B_2 at the grid nodes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, &
        nf90_close, nf90_clobber, nf90_double, nf90_int, nf90_char
    use tiago_analytic_state, only: circular_state_t
    use tiago_field_line_label, only: grid_field_t, line_labels_t, field_line_labels, &
        trace_returns
    use test_helical_island_field, only: helical_field_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), twopi = 2.0_dp * pi
    integer :: failures
    character(len=1024) :: workdir
    type(circular_state_t) :: chart
    type(helical_field_t) :: field

    failures = 0
    call get_command_argument(1, workdir)
    chart%r0 = field%r0
    chart%a = field%a
    chart%delta = 0.0_dp
    call nested_case()
    call circulating_case()
    call island_case()
    call grid_case()
    call mgrid_case()
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_field_line_label passed'

contains

    subroutine check(ok, what, value)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        real(dp), intent(in) :: value
        if (ok) then
            print '(A,ES10.3)', 'ok   '//what//' ', value
        else
            write(error_unit, '(A,ES10.3)') 'FAIL '//what//' ', value
            failures = failures + 1
        end if
    end subroutine check

    subroutine nested_case()
        type(helical_field_t) :: f
        type(line_labels_t) :: lab
        real(dp) :: r(4), z(4)

        f = field
        f%eps = 0.0_dp
        r = f%r0 + [0.2_dp, 0.45_dp, -0.6_dp, 0.0_dp]
        z = [0.0_dp, 0.1_dp, 0.0_dp, 0.75_dp]
        call field_line_labels(f, chart, r, z, 0.3_dp, 50, 360, lab)
        call check(.not. any(lab%lost), 'nested: no line lost', 0.0_dp)
        call check(maxval(abs(lab%label - lab%s_start)) < 1.0e-9_dp, &
            'nested: label = s of the start point', maxval(abs(lab%label - lab%s_start)))
        call check(maxval(lab%spread) < 1.0e-9_dp, 'nested: s constant on the line', &
            maxval(lab%spread))
    end subroutine nested_case

    real(dp) function oracle_label(k, branch_p) result(mean)
        !! Invariant-measure mean of s = 2 p / a^2 on the circulating curve
        !! K = k: p(xi) is the root of the quadratic nearest branch_p, and the
        !! measure is dxi / (dxi/dphi), dxi/dphi = 2 dtheta/dphi - 1.
        real(dp), intent(in) :: k, branch_p
        integer, parameter :: n = 4096
        real(dp) :: xi, qa, qb, disc, p1, p2, p, w, num, den
        integer :: j

        num = 0.0_dp
        den = 0.0_dp
        do j = 0, n - 1
            xi = twopi * real(j, dp) / real(n, dp)
            qa = 0.5_dp * field%iota1
            qb = field%iota0 - 0.5_dp + field%eps * cos(xi)
            disc = qb**2 + 4.0_dp * qa * k
            p1 = (-qb + sqrt(disc)) / (2.0_dp * qa)
            p2 = (-qb - sqrt(disc)) / (2.0_dp * qa)
            p = p1
            if (abs(p2 - branch_p) < abs(p1 - branch_p)) p = p2
            w = 1.0_dp / abs(2.0_dp * (field%iota0 + field%iota1 * p + field%eps * cos(xi)) &
                - 1.0_dp)
            num = num + w * 2.0_dp * p / field%a**2
            den = den + w
        end do
        mean = num / den
    end function oracle_label

    subroutine circulating_case()
        type(line_labels_t) :: lab
        real(dp), parameter :: phi0 = 0.0_dp
        real(dp) :: r(4), z(4), err, p0, k
        integer :: i

        ! Circles inside and outside the island (p_res = 0.25, half width ~0.05).
        r = field%r0 + [0.55_dp, 0.62_dp, 0.78_dp, 0.85_dp]
        z = 0.0_dp
        call field_line_labels(field, chart, r, z, phi0, 4000, 240, lab)
        err = 0.0_dp
        do i = 1, size(r)
            p0 = 0.5_dp * (r(i) - field%r0)**2
            k = field%invariant(p0, 2.0_dp * 0.0_dp - phi0)
            err = max(err, abs(lab%label(i) - oracle_label(k, p0)))
        end do
        call check(.not. any(lab%lost), 'circulating: no line lost', 0.0_dp)
        ! Finite-return errors are 9e-5 at 250 turns and 2e-6 at 4000 turns.
        call check(err < 1.0e-5_dp, 'circulating: label = invariant-curve mean of s', err)
        call check(maxval(abs(lab%label - lab%s_start)) > 1.0e-4_dp, &
            'circulating: perturbation moves the label', maxval(abs(lab%label - lab%s_start)))
    end subroutine circulating_case

    subroutine island_case()
        !! At phi = 0, the O point of eps > 0 lies where cos(2 theta) = -1,
        !! i.e. theta = +-pi/2; points at theta = pi/2 +- dtheta on the same
        !! circle have equal K (cos is even) and lie on one island tube.
        type(line_labels_t) :: lab
        real(dp) :: p_res, rr, r(2), z(2), th(2), s_res
        integer :: i

        p_res = (0.5_dp - field%iota0) / field%iota1
        s_res = 2.0_dp * p_res / field%a**2
        rr = sqrt(2.0_dp * (p_res + 0.02_dp))
        th = 0.5_dp * pi + [0.3_dp, -0.3_dp]
        do i = 1, 2
            r(i) = field%r0 + rr * cos(th(i))
            z(i) = rr * sin(th(i))
        end do
        call field_line_labels(field, chart, r, z, 0.0_dp, 1000, 240, lab)
        call check(abs(lab%label(1) - lab%label(2)) < 2.0e-4_dp, &
            'island: mirror points share the label', abs(lab%label(1) - lab%label(2)))
        call check(abs(lab%label(1) - s_res) < 0.5_dp * abs(lab%s_start(1) - s_res), &
            'island: label flattened towards the resonance', abs(lab%label(1) - s_res))
    end subroutine island_case

    subroutine fill_grid(nr, np, nz, rmin, rmax, zmin, zmax, br, bphi, bz)
        integer, intent(in) :: nr, np, nz
        real(dp), intent(in) :: rmin, rmax, zmin, zmax
        real(dp), intent(out) :: br(nr, np, nz), bphi(nr, np, nz), bz(nr, np, nz)
        real(dp) :: b(3), rr, ph, zz
        logical :: inside
        integer :: i, j, k

        do k = 1, nz
            do j = 1, np
                do i = 1, nr
                    rr = rmin + (rmax - rmin) * real(i - 1, dp) / real(nr - 1, dp)
                    ph = twopi * real(j - 1, dp) / real(np, dp)
                    zz = zmin + (zmax - zmin) * real(k - 1, dp) / real(nz - 1, dp)
                    call field%eval(rr, ph, zz, b, inside)
                    br(i, j, k) = b(1)
                    bphi(i, j, k) = b(2)
                    bz(i, j, k) = b(3)
                end do
            end do
        end do
    end subroutine fill_grid

    subroutine grid_case()
        integer, parameter :: nr = 61, np = 72, nz = 61
        real(dp), allocatable :: br(:, :, :), bphi(:, :, :), bz(:, :, :)
        type(grid_field_t) :: grid
        type(line_labels_t) :: exact, interp
        real(dp) :: b(3), bg(3), err, rr, ph, zz, r(3), z(3)
        logical :: inside, ins
        integer :: i

        allocate(br(nr, np, nz), bphi(nr, np, nz), bz(nr, np, nz))
        call fill_grid(nr, np, nz, 1.9_dp, 4.1_dp, -1.1_dp, 1.1_dp, br, bphi, bz)
        call grid%init(br, bphi, bz, 1.9_dp, 4.1_dp, -1.1_dp, 1.1_dp, 1)
        err = 0.0_dp
        do i = 1, 200
            rr = field%r0 + 0.9_dp * cos(0.37_dp * i) * sqrt(real(i, dp) / 200.0_dp)
            zz = 0.9_dp * sin(0.37_dp * i) * sqrt(real(i, dp) / 200.0_dp)
            ph = 0.113_dp * i
            call field%eval(rr, ph, zz, b, inside)
            call grid%eval(rr, ph, zz, bg, ins)
            err = max(err, maxval(abs(bg - b)) / abs(b(2)))
        end do
        call check(err < 1.0e-5_dp, 'grid: tricubic field against analytic', err)
        r = field%r0 + [0.6_dp, 0.8_dp, -0.4_dp]
        z = [0.0_dp, 0.0_dp, 0.2_dp]
        call field_line_labels(field, chart, r, z, 0.0_dp, 200, 240, exact)
        call field_line_labels(grid, chart, r, z, 0.0_dp, 200, 240, interp)
        err = maxval(abs(exact%label - interp%label))
        call check(err < 1.0e-4_dp, 'grid: labels against analytic field', err)
    end subroutine grid_case

    subroutine mgrid_case()
        integer, parameter :: nr = 21, np = 12, nz = 19, nfp = 5
        real(dp), parameter :: rmin = 2.0_dp, rmax = 4.0_dp, zmin = -1.0_dp, zmax = 1.0_dp
        real(dp) :: g(nr, nz, np, 3, 2), b(3), want(3), err, rr, ph, zz
        type(grid_field_t) :: grid
        character(len=:), allocatable :: path, message
        integer :: i, j, k, c, grp, ierr
        logical :: inside

        do grp = 1, 2
            do c = 1, 3
                do k = 1, np
                    do j = 1, nz
                        do i = 1, nr
                            g(i, j, k, c, grp) = real(grp, dp) * (1.0_dp + 0.1_dp * c) &
                                + 0.3_dp * i - 0.2_dp * j + 0.05_dp * k * c - 0.7_dp * grp * j * k
                        end do
                    end do
                end do
            end do
        end do
        path = trim(workdir)//'/mgrid_two_groups.nc'
        call write_mgrid(path, g, nfp, rmin, rmax, zmin, zmax)
        call grid%load_mgrid(path, [2.0_dp, -0.5_dp], ierr, message)
        call check(ierr == 0, 'mgrid: read two groups '//message, real(ierr, dp))
        if (ierr /= 0) return
        grid%margin = 0.0_dp
        err = 0.0_dp
        do k = 1, np
            do j = 2, nz - 1, 3
                do i = 2, nr - 1, 4
                    rr = rmin + (rmax - rmin) * real(i - 1, dp) / real(nr - 1, dp)
                    zz = zmin + (zmax - zmin) * real(j - 1, dp) / real(nz - 1, dp)
                    ph = twopi / real(nfp, dp) * real(k - 1, dp) / real(np, dp) &
                        + twopi / real(nfp, dp) * 3.0_dp
                    call grid%eval(rr, ph, zz, b, inside)
                    want = 2.0_dp * g(i, j, k, :, 1) - 0.5_dp * g(i, j, k, :, 2)
                    err = max(err, maxval(abs(b - want)))
                end do
            end do
        end do
        call check(err < 1.0e-10_dp, 'mgrid: nodes reproduce the weighted group sum', err)
    end subroutine mgrid_case

    subroutine write_mgrid(path, g, nfp, rmin, rmax, zmin, zmax)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: g(:, :, :, :, :), rmin, rmax, zmin, zmax
        integer, intent(in) :: nfp
        integer :: ncid, dr, dz, dph, dcur, dmode, vid(3, 2), v(7), grp, c
        character(len=6) :: name
        character(len=2), parameter :: comp(3) = ['br', 'bp', 'bz']

        call nc(nf90_create(path, nf90_clobber, ncid))
        call nc(nf90_def_dim(ncid, 'rad', size(g, 1), dr))
        call nc(nf90_def_dim(ncid, 'zee', size(g, 2), dz))
        call nc(nf90_def_dim(ncid, 'phi', size(g, 3), dph))
        call nc(nf90_def_dim(ncid, 'external_coils', 2, dcur))
        call nc(nf90_def_dim(ncid, 'dim_00001', 1, dmode))
        ! Fortran (rad, zee, phi) is netCDF (phi, zee, rad), the makegrid layout.
        do grp = 1, 2
            do c = 1, 3
                write(name, '(A,"_",I3.3)') comp(c), grp
                call nc(nf90_def_var(ncid, trim(name), nf90_double, [dr, dz, dph], vid(c, grp)))
            end do
        end do
        call nc(nf90_def_var(ncid, 'nfp', nf90_int, v(1)))
        call nc(nf90_def_var(ncid, 'nextcur', nf90_int, v(2)))
        call nc(nf90_def_var(ncid, 'rmin', nf90_double, v(3)))
        call nc(nf90_def_var(ncid, 'rmax', nf90_double, v(4)))
        call nc(nf90_def_var(ncid, 'zmin', nf90_double, v(5)))
        call nc(nf90_def_var(ncid, 'zmax', nf90_double, v(6)))
        call nc(nf90_def_var(ncid, 'mgrid_mode', nf90_char, [dmode], v(7)))
        call nc(nf90_enddef(ncid))
        do grp = 1, 2
            do c = 1, 3
                call nc(nf90_put_var(ncid, vid(c, grp), g(:, :, :, c, grp)))
            end do
        end do
        call nc(nf90_put_var(ncid, v(1), nfp))
        call nc(nf90_put_var(ncid, v(2), 2))
        call nc(nf90_put_var(ncid, v(3), rmin))
        call nc(nf90_put_var(ncid, v(4), rmax))
        call nc(nf90_put_var(ncid, v(5), zmin))
        call nc(nf90_put_var(ncid, v(6), zmax))
        call nc(nf90_put_var(ncid, v(7), 'S'))
        call nc(nf90_close(ncid))
    end subroutine write_mgrid

    subroutine nc(status)
        integer, intent(in) :: status
        if (status /= 0) then
            write(error_unit, '(A,I0)') 'netCDF error ', status
            error stop 2
        end if
    end subroutine nc

end program test_field_line_label
