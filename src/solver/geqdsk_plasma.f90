module tiago_geqdsk_plasma
    !! Plasma current and flux labels of an axisymmetric equilibrium from a
    !! G-EQDSK file, in the convention of docs/axisymmetric-free-boundary.md:
    !! (R, phi, Z) right-handed, currents positive along +phi, psi = R A_phi
    !! per radian. The file's current is taken along +phi. Its psi is converted
    !! as psi = sigma psi_file / (2 pi)^e, with e in {0, 1} chosen so that the
    !! integral of j_phi matches the file's current and sigma from the sign.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geqdsk_tools, only: geqdsk_t, geqdsk_read
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData2D, &
        construct_batch_splines_1d, construct_batch_splines_2d, &
        evaluate_batch_splines_1d, evaluate_batch_splines_2d_der
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: mu0 = 4.0e-7_dp * pi
    ! libneo stores G-EQDSK data in CGS.
    real(dp), parameter :: cm = 1.0e-2_dp, maxwell = 1.0e-8_dp, barye = 0.1_dp, &
        gauss = 1.0e-4_dp
    integer, parameter :: nflux = 2001

    type, public :: geqdsk_plasma_t
        integer :: nw = 0, nh = 0
        real(dp), allocatable :: r(:), z(:), psi(:, :) !! grid [m], psi [Wb/rad]
        real(dp), allocatable :: rb(:), zb(:) !! last closed surface [m]
        real(dp) :: psi_axis = 0, psi_b = 0 !! [Wb/rad]
        real(dp) :: r_axis = 0, z_axis = 0, current = 0 !! [m], file current [A]
        real(dp) :: sigma = 1 !! sign of the psi conversion
        integer :: exponent = 0 !! e of the 2 pi factor
        !> sign(F I_P) of the file, which no COCOS changes: B_phi = helicity |F| / R
        !> with the current along +phi.
        real(dp) :: helicity = 1
        real(dp), allocatable :: phi(:) !! toroidal flux on psi_N [Wb]
        real(dp), allocatable :: jphi(:, :) !! provider's j_phi on the grid [A/m^2]
        !> Integrate j_phi over grid cells instead of summing it at the nodes:
        !> 2 x 2 Gauss points in each cell inside the plasma and, in each cell the
        !> boundary cuts, one filament at the current centroid of n_sub x n_sub
        !> sub-cells weighted with their area fraction below psi_N = 1. The node
        !> sum reproduces a provider whose current lives on the nodes (FreeGS,
        !> FreeGSNKE); cells approximate a continuous current (analytic, CHEASE).
        logical :: cell_quadrature = .false.
        integer :: n_sub = 16
        type(BatchSplineData2D) :: psi_spline
        type(BatchSplineData1D) :: profile_spline !! p', F F', |q|, |F|, p on psi_N
    contains
        procedure :: label => plasma_label
        procedure :: rho_tor => plasma_rho_tor
        procedure :: inside => plasma_inside
        procedure :: filaments => plasma_filaments
        procedure :: toroidal_flux => plasma_toroidal_flux
        procedure :: modb => plasma_modb
        procedure :: field => plasma_field
        procedure :: pressure => plasma_pressure
        procedure :: midplane_jphi => plasma_midplane_jphi
    end type geqdsk_plasma_t

    public :: read_geqdsk_plasma

contains

    subroutine read_geqdsk_plasma(path, eq, ierr, message)
        character(*), intent(in) :: path
        type(geqdsk_plasma_t), intent(out) :: eq
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        type(geqdsk_t) :: g
        real(dp), allocatable :: prof(:, :), rf(:), zf(:), wf(:), x(:)
        real(dp) :: raw, ratio(0:1), dx
        integer :: i
        message = ''
        call geqdsk_read(g, path, ierr)
        if (ierr /= 0) then
            message = 'cannot read G-EQDSK file '//trim(path)
            return
        end if
        eq%nw = g%nw
        eq%nh = g%nh
        eq%r = g%r_eqd * cm
        eq%z = g%z_eqd * cm
        eq%psi = g%psirz * maxwell
        eq%psi_axis = g%simag * maxwell
        eq%psi_b = g%sibry * maxwell
        eq%r_axis = g%rmaxis * cm
        eq%z_axis = g%zmaxis * cm
        eq%current = g%current / (0.1_dp * 2.99792458e10_dp)
        eq%rb = g%rbbbs * cm
        eq%zb = g%zbbbs * cm
        allocate(prof(g%nw, 5))
        prof(:, 1) = g%pprime * barye / maxwell
        prof(:, 2) = g%ffprim * gauss
        prof(:, 3) = abs(g%qpsi)
        prof(:, 4) = abs(g%fpol) * gauss * cm
        prof(:, 5) = g%pres * barye
        call construct_batch_splines_1d(0.0_dp, 1.0_dp, prof, 3, .false., eq%profile_spline)
        call construct_batch_splines_2d([eq%r(1), eq%z(1)], [eq%r(eq%nw), eq%z(eq%nh)], &
            reshape(eq%psi, [eq%nw, eq%nh, 1]), [3, 3], [.false., .false.], eq%psi_spline)
        ! File units first, then the convention from the plasma current.
        call eq%filaments(rf, zf, wf)
        raw = sum(wf)
        if (raw == 0 .or. eq%current == 0) then
            ierr = 2
            message = 'no plasma current in '//trim(path)
            return
        end if
        ratio(0) = abs(raw / eq%current)
        ratio(1) = abs(2 * pi * raw / eq%current)
        eq%exponent = minloc(abs(log(ratio)), 1) - 1
        call read_jphi(path//'.jphi', eq)
        if (abs(log(ratio(eq%exponent))) > 0.1_dp .and. .not. allocated(eq%jphi)) then
            ierr = 3
            message = 'integral of j_phi does not match the current of '//trim(path)
            return
        end if
        eq%sigma = sign(1.0_dp, raw * eq%current)
        eq%helicity = sign(1.0_dp, g%fpol(g%nw) * eq%current)
        eq%psi = eq%sigma * eq%psi / (2 * pi)**eq%exponent
        eq%psi_axis = eq%sigma * eq%psi_axis / (2 * pi)**eq%exponent
        eq%psi_b = eq%sigma * eq%psi_b / (2 * pi)**eq%exponent
        prof(:, 1:2) = eq%sigma * (2 * pi)**eq%exponent * prof(:, 1:2)
        call construct_batch_splines_1d(0.0_dp, 1.0_dp, prof, 3, .false., eq%profile_spline)
        call construct_batch_splines_2d([eq%r(1), eq%z(1)], [eq%r(eq%nw), eq%z(eq%nh)], &
            reshape(eq%psi, [eq%nw, eq%nh, 1]), [3, 3], [.false., .false.], eq%psi_spline)
        ! Phi(psi_N) = 2 pi |psi_b - psi_axis| int_0^psi_N |q| dpsi_N, trapezoidal.
        allocate(eq%phi(nflux), x(nflux))
        dx = 1.0_dp / (nflux - 1)
        do i = 1, nflux
            x(i) = q_at(eq, (i - 1) * dx)
        end do
        eq%phi(1) = 0
        do i = 2, nflux
            eq%phi(i) = eq%phi(i - 1) + (x(i - 1) + x(i)) * dx / 2
        end do
        eq%phi = 2 * pi * abs(eq%psi_b - eq%psi_axis) * eq%phi
    end subroutine read_geqdsk_plasma

    subroutine read_jphi(path, eq)
        !! Optional sidecar of the provider's toroidal current density on the
        !! G-EQDSK grid, R fastest, one value per line [A/m^2]. With it the
        !! plasma current is the provider's own, independent of the boundary
        !! polygon and the profile tables.
        character(*), intent(in) :: path
        type(geqdsk_plasma_t), intent(inout) :: eq
        integer :: unit, ios
        logical :: present
        inquire(file=path, exist=present)
        if (.not. present) return
        allocate(eq%jphi(eq%nw, eq%nh))
        open(newunit=unit, file=path, status='old', action='read')
        read(unit, *, iostat=ios) eq%jphi
        close(unit)
        if (ios /= 0) deallocate(eq%jphi)
    end subroutine read_jphi

    real(dp) function q_at(eq, x)
        type(geqdsk_plasma_t), intent(in) :: eq
        real(dp), intent(in) :: x
        real(dp) :: y(5)
        call evaluate_batch_splines_1d(eq%profile_spline, x, y)
        q_at = y(3)
    end function q_at

    real(dp) function plasma_pressure(self, psin) result(p)
        !! Pressure [Pa] at psi_N from the file's table.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: psin
        real(dp) :: y(5)
        call evaluate_batch_splines_1d(self%profile_spline, min(max(psin, 0.0_dp), 1.0_dp), y)
        p = y(5)
    end function plasma_pressure

    subroutine plasma_midplane_jphi(self, psin, j, r)
        !! j_phi = R p' + F F' / (mu0 R) [A/m^2] at the outboard point r [m] of
        !! label psi_N on the horizontal line through the magnetic axis.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: psin
        real(dp), intent(out) :: j, r
        real(dp) :: lo, hi, x, a, b, y(5)
        integer :: it
        lo = self%r_axis
        hi = self%r(self%nw)
        do it = 1, 60
            r = (lo + hi) / 2
            call self%label(r, self%z_axis, x, a, b)
            if (x < psin) then
                lo = r
            else
                hi = r
            end if
        end do
        call evaluate_batch_splines_1d(self%profile_spline, psin, y)
        j = r * y(1) + y(2) / (mu0 * r)
    end subroutine plasma_midplane_jphi

    subroutine plasma_label(self, r, z, psin, dpsin_dr, dpsin_dz)
        !! psi_N = (psi - psi_axis) / (psi_b - psi_axis) and its gradient [1/m].
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: psin, dpsin_dr, dpsin_dz
        real(dp) :: y(1), dy(2, 1), scale
        call evaluate_batch_splines_2d_der(self%psi_spline, [r, z], y, dy)
        scale = 1 / (self%psi_b - self%psi_axis)
        psin = (y(1) - self%psi_axis) * scale
        dpsin_dr = dy(1, 1) * scale
        dpsin_dz = dy(2, 1) * scale
    end subroutine plasma_label

    subroutine plasma_toroidal_flux(self, psin, s, ds_dpsin)
        !! Normalized toroidal flux s = Phi(psi_N) / Phi(1): cubic Hermite
        !! interpolation of the tabulated Phi with the slopes
        !! dPhi/dpsi_N = 2 pi |psi_b - psi_axis| |q|, so that s and ds/dpsi_N
        !! are consistent. Beyond psi_N = 1 s continues linearly with the slope
        !! at the boundary, so that s > 1 outside the plasma.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: psin
        real(dp), intent(out) :: s, ds_dpsin
        real(dp) :: t, h, x0, m0, m1, c, h00, h10, h01, h11
        integer :: i
        c = 2 * pi * abs(self%psi_b - self%psi_axis)
        h = 1.0_dp / (nflux - 1)
        x0 = min(max(psin, 0.0_dp), 1.0_dp)
        i = min(int(x0 / h) + 1, nflux - 1)
        t = (x0 - (i - 1) * h) / h
        m0 = c * q_at(self, (i - 1) * h) * h
        m1 = c * q_at(self, i * h) * h
        h00 = (1 + 2 * t) * (1 - t)**2
        h10 = t * (1 - t)**2
        h01 = t**2 * (3 - 2 * t)
        h11 = t**2 * (t - 1)
        s = (h00 * self%phi(i) + h10 * m0 + h01 * self%phi(i + 1) + h11 * m1) / self%phi(nflux)
        ds_dpsin = (6 * t * (t - 1) * (self%phi(i) - self%phi(i + 1)) &
            + (1 - 4 * t + 3 * t**2) * m0 + (3 * t**2 - 2 * t) * m1) / (h * self%phi(nflux))
        if (psin > 1) s = 1 + (psin - 1) * ds_dpsin
        if (psin < 0) s = psin * ds_dpsin
    end subroutine plasma_toroidal_flux

    real(dp) function plasma_modb(self, r, z) result(b)
        !! |B| = sqrt(|grad psi|^2 + F^2) / R [T], with F of the boundary outside
        !! the plasma.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: psin, a, c, y(5), dpsi
        call self%label(r, z, psin, a, c)
        call evaluate_batch_splines_1d(self%profile_spline, min(max(psin, 0.0_dp), 1.0_dp), y)
        dpsi = self%psi_b - self%psi_axis
        b = sqrt((a**2 + c**2) * dpsi**2 + y(4)**2) / r
    end function plasma_modb

    subroutine plasma_field(self, r, z, b)
        !! B = (B_R, B_phi, B_Z) [T] with B_R = -dpsi/dZ / R, B_Z = dpsi/dR / R,
        !! B_phi = F / R, F of the boundary outside the plasma.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b(3)
        real(dp) :: psin, a, c, y(5), dpsi
        call self%label(r, z, psin, a, c)
        call evaluate_batch_splines_1d(self%profile_spline, min(max(psin, 0.0_dp), 1.0_dp), y)
        dpsi = self%psi_b - self%psi_axis
        b = [-c * dpsi, self%helicity * y(4), a * dpsi] / r
    end subroutine plasma_field

    real(dp) function plasma_rho_tor(self, psin) result(rho)
        !! sqrt(Phi(psi_N) / Phi(1)) of plasma_toroidal_flux, psi_N in [0, 1].
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: psin
        real(dp) :: s, ds
        call self%toroidal_flux(min(max(psin, 0.0_dp), 1.0_dp), s, ds)
        rho = sqrt(max(s, 0.0_dp))
    end function plasma_rho_tor

    logical function plasma_inside(self, r, z) result(in)
        !! Inside the last closed surface polygon and below psi_N = 1.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: psin, a, b
        integer :: i, j
        in = .false.
        j = size(self%rb)
        do i = 1, size(self%rb)
            if ((self%zb(i) > z) .neqv. (self%zb(j) > z)) then
                if (r < self%rb(i) + (z - self%zb(i)) * (self%rb(j) - self%rb(i)) &
                    / (self%zb(j) - self%zb(i))) in = .not. in
            end if
            j = i
        end do
        if (.not. in) return
        call self%label(r, z, psin, a, b)
        in = psin >= 0 .and. psin < 1
    end function plasma_inside

    subroutine plasma_filaments(self, rf, zf, wf)
        !! Grid nodes inside the plasma with the current j_phi dR dZ [A] of
        !! j_phi = R p'(psi) + F F'(psi) / (mu0 R), or of the provider's j_phi.
        class(geqdsk_plasma_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: rf(:), zf(:), wf(:)
        real(dp) :: area, psin, a, b, y(5), rbuf(self%nw * self%nh), &
            zbuf(self%nw * self%nh), wbuf(self%nw * self%nh)
        integer :: i, j, n
        area = (self%r(2) - self%r(1)) * (self%z(2) - self%z(1))
        n = 0
        if (self%cell_quadrature .and. .not. allocated(self%jphi)) then
            call cell_filaments(self, rf, zf, wf)
            return
        end if
        if (allocated(self%jphi)) then
            n = count(self%jphi /= 0)
            rf = pack(spread(self%r, 2, self%nh), self%jphi /= 0)
            zf = pack(spread(self%z, 1, self%nw), self%jphi /= 0)
            wf = pack(self%jphi * area, self%jphi /= 0)
            return
        end if
        do j = 1, self%nh
            do i = 1, self%nw
                if (.not. self%inside(self%r(i), self%z(j))) cycle
                call self%label(self%r(i), self%z(j), psin, a, b)
                call evaluate_batch_splines_1d(self%profile_spline, psin, y)
                n = n + 1
                rbuf(n) = self%r(i)
                zbuf(n) = self%z(j)
                wbuf(n) = (self%r(i) * y(1) + y(2) / (mu0 * self%r(i))) * area
            end do
        end do
        rf = rbuf(:n)
        zf = zbuf(:n)
        wf = wbuf(:n)
    end subroutine plasma_filaments

    subroutine cell_filaments(self, rf, zf, wf)
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), allocatable, intent(out) :: rf(:), zf(:), wf(:)
        real(dp), parameter :: g = 0.5_dp / sqrt(3.0_dp)
        real(dp), allocatable :: rbuf(:), zbuf(:), wbuf(:)
        real(dp) :: dr, dz, rc, zc, r, z, j, cur, mr, mz, psin, a, b, tol, f
        logical :: in(4)
        integer :: i, k, ia, ib, n, nmax
        tol = maxval(hypot(self%rb - cshift(self%rb, 1), self%zb - cshift(self%zb, 1)))
        dr = self%r(2) - self%r(1)
        dz = self%z(2) - self%z(1)
        nmax = 4 * (self%nw - 1) * (self%nh - 1)
        allocate(rbuf(nmax), zbuf(nmax), wbuf(nmax))
        n = 0
        do k = 1, self%nh - 1
            do i = 1, self%nw - 1
                in = [inside_fine(self, self%r(i), self%z(k), tol), &
                    inside_fine(self, self%r(i + 1), self%z(k), tol), &
                    inside_fine(self, self%r(i), self%z(k + 1), tol), &
                    inside_fine(self, self%r(i + 1), self%z(k + 1), tol)]
                rc = self%r(i) + dr / 2
                zc = self%z(k) + dz / 2
                if (all(in)) then
                    do ib = -1, 1, 2
                        do ia = -1, 1, 2
                            r = rc + ia * g * dr
                            z = zc + ib * g * dz
                            n = n + 1
                            rbuf(n) = r
                            zbuf(n) = z
                            wbuf(n) = current_density(self, r, z) * dr * dz / 4
                        end do
                    end do
                    cycle
                end if
                if (.not. any(in)) then
                    ! A boundary may still cross a cell with no corner inside.
                    call self%label(rc, zc, psin, a, b)
                    if (psin > 1 + 2 * (abs(a) * dr + abs(b) * dz)) cycle
                end if
                cur = 0
                mr = 0
                mz = 0
                do ib = 1, self%n_sub
                    do ia = 1, self%n_sub
                        r = self%r(i) + (ia - 0.5_dp) * dr / self%n_sub
                        z = self%z(k) + (ib - 0.5_dp) * dz / self%n_sub
                        if (.not. inside_region(self, r, z, tol + (dr + dz) / self%n_sub, &
                            (dr + dz) / self%n_sub)) cycle
                        call self%label(r, z, psin, a, b)
                        ! Area fraction of the sub-cell below psi_N = 1 for a
                        ! linear psi_N, so that the cut error is of second order.
                        f = min(max(0.5_dp + (1 - psin) / max(abs(a) * dr + abs(b) * dz, &
                            tiny(1.0_dp)) * self%n_sub, 0.0_dp), 1.0_dp)
                        if (f == 0 .or. psin < 0) cycle
                        j = f * current_density(self, r, z) * dr * dz / self%n_sub**2
                        cur = cur + j
                        mr = mr + j * r
                        mz = mz + j * z
                    end do
                end do
                if (cur == 0) cycle
                n = n + 1
                rbuf(n) = mr / cur
                zbuf(n) = mz / cur
                wbuf(n) = cur
            end do
        end do
        rf = rbuf(:n)
        zf = zbuf(:n)
        wf = wbuf(:n)
    end subroutine cell_filaments

    logical function inside_fine(self, r, z, tol) result(in)
        !! Inside the plasma: psi_N < 1 within the region of inside_region.
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z, tol
        real(dp) :: psin, a, b
        call self%label(r, z, psin, a, b)
        in = psin >= 0 .and. psin < 1
        if (in) in = inside_region(self, r, z, tol, 0.0_dp)
    end function inside_fine

    logical function inside_region(self, r, z, tol, skip) result(in)
        !! Inside the polygon of the last closed surface, or within tol of it and
        !! connected to the axis. psi_N then decides near the boundary, so the
        !! chords and the rounded X-point corner of a coarse polygon cut off no
        !! plasma.
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z, tol, skip
        real(dp) :: t, d, er, ez
        integer :: i, j
        in = polygon_inside(self, r, z)
        if (in) return
        j = size(self%rb)
        do i = 1, size(self%rb)
            er = self%rb(i) - self%rb(j)
            ez = self%zb(i) - self%zb(j)
            t = min(max(((r - self%rb(j)) * er + (z - self%zb(j)) * ez) / max(er**2 + ez**2, &
                tiny(1.0_dp)), 0.0_dp), 1.0_dp)
            d = hypot(r - self%rb(j) - t * er, z - self%zb(j) - t * ez)
            if (d < tol) then
                in = connected_to_axis(self, r, z, skip)
                return
            end if
            j = i
        end do
    end function inside_region

    logical function connected_to_axis(self, r, z, skip) result(in)
        !! The straight path to the magnetic axis stays below psi_N = 1. A point
        !! of the private flux region below an X-point crosses the separatrix
        !! on its way to the axis, possibly close to the X-point, where psi_N
        !! exceeds 1 only on a short piece. The largest sample is therefore
        !! refined by golden-section search. The first skip [m] of the path is
        !! not tested, so that points just outside the boundary keep their
        !! sub-cell area fraction.
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z, skip
        integer, parameter :: n = 64
        real(dp), parameter :: c = (sqrt(5.0_dp) - 1) / 2
        real(dp) :: f(n - 1), lo, hi, t1, t2, f1, f2, t0
        integer :: k, km, it
        t0 = min(skip / max(hypot(self%r_axis - r, self%z_axis - z), tiny(1.0_dp)), 0.5_dp)
        do k = 1, n - 1
            f(k) = path_label(t0 + (1 - t0) * real(k, dp) / n)
        end do
        in = .false.
        if (maxval(f) >= 1) return
        km = maxloc(f, 1)
        lo = t0 + (1 - t0) * real(km - 1, dp) / n
        hi = t0 + (1 - t0) * real(km + 1, dp) / n
        t1 = hi - c * (hi - lo)
        t2 = lo + c * (hi - lo)
        f1 = path_label(t1)
        f2 = path_label(t2)
        do it = 1, 40
            if (max(f1, f2) >= 1) return
            if (f1 > f2) then
                hi = t2
                t2 = t1
                f2 = f1
                t1 = hi - c * (hi - lo)
                f1 = path_label(t1)
            else
                lo = t1
                t1 = t2
                f1 = f2
                t2 = lo + c * (hi - lo)
                f2 = path_label(t2)
            end if
        end do
        in = .true.
    contains
        real(dp) function path_label(t)
            real(dp), intent(in) :: t
            real(dp) :: a, b
            call self%label(r + t * (self%r_axis - r), z + t * (self%z_axis - z), path_label, a, b)
        end function path_label
    end function connected_to_axis

    logical function polygon_inside(self, r, z) result(in)
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        integer :: i, j
        in = .false.
        j = size(self%rb)
        do i = 1, size(self%rb)
            if ((self%zb(i) > z) .neqv. (self%zb(j) > z)) then
                if (r < self%rb(i) + (z - self%zb(i)) * (self%rb(j) - self%rb(i)) &
                    / (self%zb(j) - self%zb(i))) in = .not. in
            end if
            j = i
        end do
    end function polygon_inside

    real(dp) function current_density(self, r, z) result(j)
        !! j_phi = R p'(psi) + F F'(psi) / (mu0 R) [A/m^2].
        type(geqdsk_plasma_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: psin, a, b, y(5)
        call self%label(r, z, psin, a, b)
        call evaluate_batch_splines_1d(self%profile_spline, min(max(psin, 0.0_dp), 1.0_dp), y)
        j = r * y(1) + y(2) / (mu0 * r)
    end function current_density

end module tiago_geqdsk_plasma
