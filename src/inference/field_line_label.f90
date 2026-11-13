module tiago_field_line_label
    !! Field-line profile labels in a magnetic field that need not have nested
    !! surfaces (vacuum islands, perturbed or non-nested equilibria).
    !!
    !! A field line is started at (R, Z) in the plane phi0 and followed with
    !! phi as the independent variable for a number of toroidal turns. Its
    !! label is the mean of the nested flux label s = rho^2 of a reference
    !! chart (any real_space_state_t) over its returns to the plane phi0. On a
    !! surface of the reference chart this is s itself; on a flux surface of a
    !! perturbed field it is the mean of s over the surface, and inside an
    !! island all points of one island tube share a label. Kinetic profiles
    !! that are constant on flux surfaces are then functions of this label.
    !!
    !! The field is any line_field_t. grid_field_t holds a field tabulated on a
    !! cylindrical grid, interpolated by tricubic splines periodic in phi
    !! (libneo batch splines), and reads makegrid/mgrid files with a vector of
    !! coil-group currents, so fields of VMEC++ free-boundary runs, coil sets
    !! and grids exported by other codes share one path.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_int
    use tiago_netcdf_c, only: NC_NOERR, NC_NOWRITE, nc_open, nc_close, nc_inq_varid, &
        nc_get_var_double, nc_shape, nc_read_text, cstr
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        evaluate_batch_splines_3d
    use tiago_real_space_state, only: real_space_state_t, located_t, locate_ok
    use tiago_field_line_rhs, only: field_line_rhs
    implicit none
    private

    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)

    type, abstract, public :: line_field_t
    contains
        procedure(field_eval_i), deferred :: eval
    end type line_field_t

    abstract interface
        subroutine field_eval_i(self, r, phi, z, b, inside)
            !! Cylindrical components b = (B_R, B_phi, B_Z) at (R, phi, Z);
            !! inside = .false. where the field is not defined.
            import :: line_field_t, dp
            class(line_field_t), intent(in) :: self
            real(dp), intent(in) :: r, phi, z
            real(dp), intent(out) :: b(3)
            logical, intent(out) :: inside
        end subroutine field_eval_i
    end interface

    type, extends(line_field_t), public :: grid_field_t
        integer :: nfp = 1
        real(dp) :: rmin = 0.0_dp, rmax = 0.0_dp, zmin = 0.0_dp, zmax = 0.0_dp
        !> Lines closer than this to the R or Z edge of the grid count as lost.
        real(dp) :: margin = 0.02_dp
        type(BatchSplineData3D) :: spline
    contains
        procedure :: init => grid_init
        procedure :: load_mgrid => grid_load_mgrid
        procedure :: eval => grid_eval
    end type grid_field_t

    type, public :: line_labels_t
        !! Per start point: s at the start, mean and standard deviation of s
        !! over the returns, the number of returns used, and whether the line
        !! left the field domain or the reference chart before the last turn.
        real(dp), allocatable :: s_start(:), label(:), spread(:)
        integer, allocatable :: returns(:)
        logical, allocatable :: lost(:)
    end type line_labels_t

    public :: trace_returns, field_line_labels, chart_s

contains

    subroutine grid_init(self, br, bphi, bz, rmin, rmax, zmin, zmax, nfp)
        !! Field components on (R, phi, Z) nodes: R and Z uniform including
        !! both ends, phi = (k - 1) 2 pi / (nfp nphi) without the end plane.
        class(grid_field_t), intent(inout) :: self
        real(dp), intent(in) :: br(:, :, :), bphi(:, :, :), bz(:, :, :)
        real(dp), intent(in) :: rmin, rmax, zmin, zmax
        integer, intent(in) :: nfp
        real(dp), allocatable :: y(:, :, :, :)
        integer :: nr, np, nz

        nr = size(br, 1)
        np = size(br, 2)
        nz = size(br, 3)
        self%nfp = nfp
        self%rmin = rmin
        self%rmax = rmax
        self%zmin = zmin
        self%zmax = zmax
        allocate(y(nr, np + 1, nz, 3))
        y(:, 1:np, :, 1) = br
        y(:, 1:np, :, 2) = bphi
        y(:, 1:np, :, 3) = bz
        y(:, np + 1, :, :) = y(:, 1, :, :)
        call construct_batch_splines_3d([rmin, 0.0_dp, zmin], &
            [rmax, twopi / real(nfp, dp), zmax], y, [3, 3, 3], &
            [.false., .true., .false.], self%spline)
    end subroutine grid_init

    subroutine grid_eval(self, r, phi, z, b, inside)
        class(grid_field_t), intent(in) :: self
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: b(3)
        logical, intent(out) :: inside

        b = 0.0_dp
        inside = r > self%rmin + self%margin .and. r < self%rmax - self%margin .and. &
            z > self%zmin + self%margin .and. z < self%zmax - self%margin
        if (.not. inside) return
        call evaluate_batch_splines_3d(self%spline, [r, phi, z], b)
    end subroutine grid_eval

    subroutine grid_load_mgrid(self, path, currents, ierr, message)
        !! makegrid file: br_NNN, bp_NNN, bz_NNN (phi, zee, rad) per coil
        !! group, one field period without its end plane. The field is
        !! sum_k currents(k) B_k; in raw mode ('R') each group is first
        !! divided by its raw_coil_cur. Groups with zero current are not read.
        !! Read through the netCDF C library, which reads any netCDF-4 layout.
        class(grid_field_t), intent(inout) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: currents(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer(c_int) :: ncid, varid
        integer :: nr, nz, np, nfp, nextcur, k, c
        integer, allocatable :: shape_(:)
        real(dp) :: rmin, rmax, zmin, zmax, scale, v(1)
        real(dp), allocatable :: raw(:), work(:, :, :), b(:, :, :, :)
        character(len=1) :: mode(1)
        character(len=16) :: tag
        character(len=2), parameter :: comp(3) = ['br', 'bp', 'bz']

        ierr = 0
        message = ''
        if (nc_open(cstr(trim(path)), NC_NOWRITE, ncid) /= NC_NOERR) then
            ierr = 1
            message = 'cannot open mgrid '//trim(path)
            return
        end if
        nfp = nint(scalar('nfp'))
        nextcur = nint(scalar('nextcur'))
        rmin = scalar('rmin')
        rmax = scalar('rmax')
        zmin = scalar('zmin')
        zmax = scalar('zmax')
        if (ierr == 0 .and. size(currents) /= nextcur) then
            ierr = 2
            write(tag, '(I0)') nextcur
            message = trim(path)//': expected '//trim(tag)//' coil-group currents'
        end if
        mode = 'S'
        if (ierr == 0) then
            allocate(raw(nextcur))
            raw = 1.0_dp
            if (nc_inq_varid(ncid, cstr('mgrid_mode'), varid) == NC_NOERR) then
                if (nc_read_text(ncid, varid, 1, 1, mode) /= NC_NOERR) mode = 'S'
            end if
            if (mode(1) == 'R' .or. mode(1) == 'r') then
                if (nc_inq_varid(ncid, cstr('raw_coil_cur'), varid) /= NC_NOERR) then
                    ierr = 3
                else if (nc_get_var_double(ncid, varid, raw) /= NC_NOERR) then
                    ierr = 3
                end if
                if (ierr /= 0) message = trim(path)//': missing variable raw_coil_cur'
            end if
        end if
        nr = 0
        if (ierr == 0) then
            do k = 1, nextcur
                if (currents(k) == 0.0_dp) cycle
                scale = currents(k)
                if (mode(1) == 'R' .or. mode(1) == 'r') scale = scale / raw(k)
                do c = 1, 3
                    write(tag, '(A,"_",I3.3)') comp(c), k
                    if (nc_inq_varid(ncid, cstr(trim(tag)), varid) /= NC_NOERR) then
                        ierr = 3
                    else if (nc_shape(ncid, varid, shape_) /= NC_NOERR) then
                        ierr = 3
                    else if (size(shape_) /= 3) then
                        ierr = 3
                    end if
                    if (ierr == 0 .and. nr == 0) then
                        np = shape_(1)
                        nz = shape_(2)
                        nr = shape_(3)
                        allocate(work(nr, nz, np), b(nr, np, nz, 3))
                        b = 0.0_dp
                    end if
                    if (ierr == 0) then
                        if (any(shape_ /= [np, nz, nr])) ierr = 3
                    end if
                    if (ierr == 0) then
                        if (nc_get_var_double(ncid, varid, work) /= NC_NOERR) ierr = 3
                    end if
                    if (ierr /= 0) then
                        message = trim(path)//': missing or malformed '//trim(tag)
                        exit
                    end if
                    b(:, :, :, c) = b(:, :, :, c) + scale * reshape(work, [nr, np, nz], &
                        order=[1, 3, 2])
                end do
                if (ierr /= 0) exit
            end do
            if (ierr == 0 .and. nr == 0) then
                ierr = 2
                message = trim(path)//': all coil-group currents are zero'
            end if
        end if
        if (nc_close(ncid) /= NC_NOERR .and. ierr == 0) then
            ierr = 1
            message = 'cannot close mgrid '//trim(path)
        end if
        if (ierr /= 0) return
        call self%init(b(:, :, :, 1), b(:, :, :, 2), b(:, :, :, 3), rmin, rmax, zmin, zmax, nfp)
    contains
        real(dp) function scalar(name)
            character(len=*), intent(in) :: name
            scalar = 0.0_dp
            if (ierr /= 0) return
            if (nc_inq_varid(ncid, cstr(name), varid) /= NC_NOERR) then
                ierr = 3
            else if (nc_get_var_double(ncid, varid, v) /= NC_NOERR) then
                ierr = 3
            end if
            if (ierr /= 0) then
                message = trim(path)//': missing variable '//name
            else
                scalar = v(1)
            end if
        end function scalar
    end subroutine grid_load_mgrid

    subroutine trace_returns(field, r0, z0, phi0, turns, steps, rz, alive)
        !! Classical fourth-order Runge-Kutta in phi from phi0 with steps per
        !! turn; rz(:, k) is (R, Z) after k full turns. alive is the number
        !! of turns completed before the line left the field domain.
        class(line_field_t), intent(in) :: field
        real(dp), intent(in) :: r0, z0, phi0
        integer, intent(in) :: turns, steps
        real(dp), intent(out) :: rz(2, turns)
        integer, intent(out) :: alive
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), h, phi
        logical :: ok
        integer :: turn, j

        rz = ieee_value(0.0_dp, ieee_quiet_nan)
        alive = 0
        y = [r0, z0]
        h = twopi / real(steps, dp)
        do turn = 1, turns
            do j = 0, steps - 1
                phi = phi0 + twopi * real(turn - 1, dp) + h * real(j, dp)
                call rhs(phi, y, k1, ok)
                if (ok) call rhs(phi + 0.5_dp * h, y + 0.5_dp * h * k1, k2, ok)
                if (ok) call rhs(phi + 0.5_dp * h, y + 0.5_dp * h * k2, k3, ok)
                if (ok) call rhs(phi + h, y + h * k3, k4, ok)
                if (.not. ok) return
                y = y + h * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
            end do
            rz(:, turn) = y
            alive = turn
        end do
    contains
        subroutine rhs(p, q, dq, inside)
            real(dp), intent(in) :: p, q(2)
            real(dp), intent(out) :: dq(2)
            logical, intent(out) :: inside
            real(dp) :: b(3)
            dq = 0.0_dp
            call field%eval(q(1), p, q(2), b, inside)
            if (.not. inside) return
            if (b(2) == 0.0_dp) then
                inside = .false.
                return
            end if
            call field_line_rhs(q(1), b(1), b(2), b(3), dq(1), dq(2))
        end subroutine rhs
    end subroutine trace_returns

    subroutine chart_s(chart, r, phi, z, s, ok)
        !! Reference flux label s = rho^2 of the point (R, phi, Z).
        class(real_space_state_t), intent(in) :: chart
        real(dp), intent(in) :: r, phi, z
        real(dp), intent(out) :: s
        logical, intent(out) :: ok
        type(located_t) :: loc

        call chart%locate([r * cos(phi), r * sin(phi), z], loc)
        ok = loc%status == locate_ok
        s = ieee_value(0.0_dp, ieee_quiet_nan)
        if (ok) s = loc%u(1)**2
    end subroutine chart_s

    subroutine field_line_labels(field, chart, r, z, phi0, turns, steps, out)
        !! Labels of the lines started at (r(i), z(i)) in the plane phi0. A
        !! line counts as lost if it leaves the field domain or a return
        !! cannot be located in the chart; its label is then the mean over
        !! the returns before that, NaN if there are none.
        class(line_field_t), intent(in) :: field
        class(real_space_state_t), intent(in) :: chart
        real(dp), intent(in) :: r(:), z(:), phi0
        integer, intent(in) :: turns, steps
        type(line_labels_t), intent(out) :: out
        real(dp) :: rz(2, turns), s(turns), mean
        integer :: i, k, alive, n
        logical :: ok

        n = size(r)
        allocate(out%s_start(n), out%label(n), out%spread(n), out%returns(n), out%lost(n))
        !$omp parallel do schedule(dynamic) default(shared) private(i, k, rz, s, alive, ok, mean)
        do i = 1, n
            call chart_s(chart, r(i), phi0, z(i), out%s_start(i), ok)
            call trace_returns(field, r(i), z(i), phi0, turns, steps, rz, alive)
            out%lost(i) = alive < turns
            do k = 1, alive
                call chart_s(chart, rz(1, k), phi0, rz(2, k), s(k), ok)
                if (.not. ok) then
                    alive = k - 1
                    out%lost(i) = .true.
                    exit
                end if
            end do
            out%returns(i) = alive
            if (alive == 0) then
                out%label(i) = ieee_value(0.0_dp, ieee_quiet_nan)
                out%spread(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                mean = sum(s(1:alive)) / real(alive, dp)
                out%label(i) = mean
                out%spread(i) = sqrt(sum((s(1:alive) - mean)**2) / real(alive, dp))
            end if
        end do
        !$omp end parallel do
    end subroutine field_line_labels

end module tiago_field_line_label
