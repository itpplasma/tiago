module tiago_point_table
    !! Equilibrium sampled on a cylindrical grid by any code (SPECTRE, HINT,
    !! KIN6D, VMEC), as a provider of the profile-diagnostic sample rows.
    !!
    !! File (netCDF-4, dimension order as written by netCDF4-python):
    !!   R(r) [m], Z(z) [m], phi(phi) [rad]: uniform grids; phi covers one
    !!     field period [0, 2 pi / nfp) without its end point; global
    !!     attribute nfp.
    !!   rho(phi, z, r): flux label, rho = 1 on the boundary, > 1 outside;
    !!   br, bphi, bz(phi, z, r) [T]: cylindrical field components;
    !!   status(phi, z, r): 0 inside the plasma, 1 outside the plasma (values
    !!     still defined), 2 no data;
    !!   optional island(phi, z, r) > 0: index of the island whose interior
    !!     contains the point, and island_rho(island): its plateau label.
    !! The island regions come from the equilibrium code's own field-line
    !! connectivity; TIAGO does not recompute them. Inside island k the
    !! profile label is island_rho(k), so kinetic profiles are flat there.
    !!
    !! Values between grid points are trilinear in (R, phi, Z), so a point
    !! on a grid node reproduces the tabulated values exactly.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int
    use tiago_netcdf_c, only: NC_NOERR, NC_NOWRITE, NC_GLOBAL, nc_open, nc_close, nc_inq_varid, &
        nc_get_var_double, nc_get_att_double, cstr, nc_shape
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_S, Q_B, N_QUANTITIES
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter, public :: TABLE_INSIDE = 0, TABLE_OUTSIDE = 1, TABLE_NO_DATA = 2

    type, public :: point_table_t
        integer :: nfp = 1
        real(dp), allocatable :: r(:), z(:), phi(:)
        real(dp), allocatable :: rho(:, :, :), babs(:, :, :) !! (r, z, phi), rho with plateaus
        integer, allocatable :: status(:, :, :)
    contains
        procedure :: load => table_load
        procedure :: sample => table_sample
    end type point_table_t

contains

    subroutine table_load(self, path, ierr, message)
        class(point_table_t), intent(inout) :: self
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer(c_int) :: ncid, varid
        real(dp), allocatable :: br(:, :, :), bphi(:, :, :), bz(:, :, :), st(:, :, :), isl(:, :, :), &
            plateau(:)
        real(dp) :: nfp(1)
        integer :: c, k

        ierr = 1
        if (nc_open(cstr(path), NC_NOWRITE, ncid) /= NC_NOERR) then
            message = 'cannot open point table '//path
            return
        end if
        message = path//': missing or malformed '
        if (nc_get_att_double(ncid, NC_GLOBAL, cstr('nfp'), nfp) /= NC_NOERR) then
            message = message//'attribute nfp'
            c = nc_close(ncid)
            return
        end if
        self%nfp = nint(nfp(1))
        if (.not. (axis('R', self%r) .and. axis('Z', self%z) .and. axis('phi', self%phi))) goto 9
        if (.not. (grid('rho', self%rho) .and. grid('br', br) .and. grid('bphi', bphi) .and. &
            grid('bz', bz) .and. grid('status', st))) goto 9
        self%babs = sqrt(br**2 + bphi**2 + bz**2)
        self%status = nint(st)
        if (nc_inq_varid(ncid, cstr('island'), varid) == NC_NOERR) then
            if (.not. (grid('island', isl) .and. axis('island_rho', plateau))) goto 9
            do k = 1, size(plateau)
                where (nint(isl) == k) self%rho = plateau(k)
            end do
            if (any(nint(isl) < 0 .or. nint(isl) > size(plateau))) then
                message = message//'island index'
                goto 9
            end if
        end if
        if (size(self%r) < 2 .or. size(self%z) < 2 .or. size(self%phi) < 1 .or. self%nfp < 1) then
            message = message//'grid size'
            goto 9
        end if
        if (any(self%status < TABLE_INSIDE .or. self%status > TABLE_NO_DATA)) then
            message = message//'status value'
            goto 9
        end if
        if (any(self%rho /= self%rho) .or. any(self%babs /= self%babs)) then
            message = message//'NaN where data is declared'
            if (any((self%rho /= self%rho .or. self%babs /= self%babs) .and. &
                self%status /= TABLE_NO_DATA)) goto 9
        end if
        ierr = 0
        message = ''
        9       c = nc_close(ncid)
    contains
        logical function axis(name, v)
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(out) :: v(:)
            integer, allocatable :: shape_(:)
            axis = .false.
            if (nc_inq_varid(ncid, cstr(name), varid) /= NC_NOERR) then
                message = message//name
                return
            end if
            if (nc_shape(ncid, varid, shape_) /= NC_NOERR) return
            if (size(shape_) /= 1) then
                message = message//name
                return
            end if
            allocate(v(shape_(1)))
            axis = nc_get_var_double(ncid, varid, v) == NC_NOERR
            if (.not. axis) message = message//name
        end function axis

        logical function grid(name, v)
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(out) :: v(:, :, :)
            integer, allocatable :: shape_(:)
            grid = .false.
            if (nc_inq_varid(ncid, cstr(name), varid) /= NC_NOERR) then
                message = message//name
                return
            end if
            if (nc_shape(ncid, varid, shape_) /= NC_NOERR) return
            if (size(shape_) /= 3) then
                message = message//name//' (needs phi, z, r)'
                return
            end if
            if (any(shape_ /= [size(self%phi), size(self%z), size(self%r)])) then
                message = message//name//' (needs phi, z, r)'
                return
            end if
            allocate(v(shape_(3), shape_(2), shape_(1)))
            grid = nc_get_var_double(ncid, varid, v) == NC_NOERR
            if (.not. grid) message = message//name
        end function grid
    end subroutine table_load

    subroutine table_sample(self, x, sample, valid)
        !! Rows Q_S = rho^2 and Q_B = |B| and their (R, phi, Z) gradients at
        !! points x(3, n) = (R, phi, Z). valid is false where a corner has no
        !! data or the point lies off the grid; inside needs every corner
        !! inside the plasma and rho <= 1.
        class(point_table_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        type(plasma_sample_t), intent(out) :: sample
        logical, allocatable, intent(out) :: valid(:)
        real(dp) :: rho, grho(3), b, gb(3)
        integer :: k, worst

        sample%n = size(x, 2)
        sample%x = x
        allocate(sample%inside(sample%n), sample%q(N_QUANTITIES, sample%n), &
            sample%dq_dx(3, N_QUANTITIES, sample%n), valid(sample%n))
        sample%q = 0.0_dp
        sample%dq_dx = 0.0_dp
        do k = 1, sample%n
            call trilinear(self, x(:, k), rho, grho, b, gb, worst)
            valid(k) = worst /= TABLE_NO_DATA
            sample%inside(k) = worst == TABLE_INSIDE .and. rho <= 1.0_dp
            if (.not. valid(k)) cycle
            sample%q(Q_S, k) = rho**2
            sample%q(Q_B, k) = b
            sample%dq_dx(:, Q_S, k) = 2.0_dp * rho * grho
            sample%dq_dx(:, Q_B, k) = gb
        end do
    end subroutine table_sample

    pure subroutine trilinear(self, x, rho, grho, b, gb, worst)
        type(point_table_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: rho, grho(3), b, gb(3)
        integer, intent(out) :: worst
        real(dp) :: t(3), h(3), period, p, w
        integer :: i(3), j(3), a, c, e, n(3), ip

        rho = 0.0_dp
        grho = 0.0_dp
        b = 0.0_dp
        gb = 0.0_dp
        worst = TABLE_NO_DATA
        n = [size(self%r), size(self%phi), size(self%z)]
        period = 2.0_dp * pi / self%nfp
        p = modulo(x(2), period)
        h = [self%r(2) - self%r(1), period / n(2), self%z(2) - self%z(1)]
        t = [(x(1) - self%r(1)) / h(1), (p - self%phi(1)) / h(2), (x(3) - self%z(1)) / h(3)]
        if (t(1) < 0.0_dp .or. t(1) > n(1) - 1 .or. t(3) < 0.0_dp .or. t(3) > n(3) - 1) return
        i = min(int(t), n - 2)
        i(2) = modulo(int(floor(t(2))), n(2))
        t = t - [real(i(1), dp), real(floor(t(2)), dp), real(i(3), dp)]
        worst = TABLE_INSIDE
        do a = 0, 1
            do c = 0, 1
                do e = 0, 1
                    j = [i(1) + a + 1, i(2) + c + 1, i(3) + e + 1]
                    ip = j(2)
                    if (ip > n(2)) ip = 1
                    worst = max(worst, self%status(j(1), j(3), ip))
                    w = weight(t(1), a) * weight(t(2), c) * weight(t(3), e)
                    rho = rho + w * self%rho(j(1), j(3), ip)
                    b = b + w * self%babs(j(1), j(3), ip)
                    grho = grho + self%rho(j(1), j(3), ip) * &
                        [dweight(a) * weight(t(2), c) * weight(t(3), e) / h(1), &
                        weight(t(1), a) * dweight(c) * weight(t(3), e) / h(2), &
                        weight(t(1), a) * weight(t(2), c) * dweight(e) / h(3)]
                    gb = gb + self%babs(j(1), j(3), ip) * &
                        [dweight(a) * weight(t(2), c) * weight(t(3), e) / h(1), &
                        weight(t(1), a) * dweight(c) * weight(t(3), e) / h(2), &
                        weight(t(1), a) * weight(t(2), c) * dweight(e) / h(3)]
                end do
            end do
        end do
        if (n(2) == 1) then
            grho(2) = 0.0_dp
            gb(2) = 0.0_dp
        end if
    contains
        pure real(dp) function weight(s, k)
            real(dp), intent(in) :: s
            integer, intent(in) :: k
            weight = merge(s, 1.0_dp - s, k == 1)
        end function weight

        pure real(dp) function dweight(k)
            integer, intent(in) :: k
            dweight = merge(1.0_dp, -1.0_dp, k == 1)
        end function dweight
    end subroutine trilinear
end module tiago_point_table
