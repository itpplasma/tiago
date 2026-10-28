module tiago_profile_geometry_io
    !! Readers for profile-diagnostic geometry and tabulated kinetic profiles.
    !!
    !! Geometry is one netCDF-4 file with an optional group per instrument
    !! (dimension order as written by netCDF4-python; Cartesian x, y, z in m):
    !!
    !!   interferometer, xics: start(chord, xyz), end(chord, xyz),
    !!                         label(chord, strlen)
    !!   thomson:              center(volume, xyz), direction(volume, xyz),
    !!                         length(volume) [m], label(volume, strlen)
    !!   ece:                  frequency(channel) [Hz], harmonic(channel),
    !!                         origin(channel, xyz), direction(channel, xyz),
    !!                         length(channel) [m], label(channel, strlen)
    !!   mse:                  origin(channel, xyz), direction(channel, xyz),
    !!                         beam_origin(channel, xyz),
    !!                         beam_direction(channel, xyz), label(channel, strlen)
    !!
    !! The profile table has a dimension s (normalized toroidal flux, 0 ... 1)
    !! and variables ne [m^-3], te [eV], ti [eV] and optionally zeff, each
    !! positive on s. A missing group or variable is reported, never filled in.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int
    use tiago_netcdf_c, only: NC_NOERR, NC_NOWRITE, nc_open, nc_close, nc_inq_ncid, &
        nc_inq_varid, nc_get_var_double, cstr, nc_shape, nc_read_text
    use tiago_sight_geometry, only: chord_t
    use tiago_thomson_scattering, only: thomson_volume_t
    implicit none
    private

    type, public :: ece_line_t
        !! One ECE channel: a straight line of sight from the antenna.
        character(len=:), allocatable :: label
        real(dp) :: frequency = 0.0_dp !! [Hz]
        integer :: harmonic = 2
        real(dp) :: origin(3) = 0.0_dp !! antenna position [m]
        real(dp) :: direction(3) = 0.0_dp !! into the plasma
        real(dp) :: length = 0.0_dp !! [m]
    end type ece_line_t

    type, public :: mse_line_t
        !! One MSE channel: sight line and heating beam, Cartesian [m].
        character(len=:), allocatable :: label
        real(dp) :: origin(3) = 0.0_dp, direction(3) = 0.0_dp
        real(dp) :: beam_origin(3) = 0.0_dp, beam_direction(3) = 0.0_dp
    end type mse_line_t

    type, public :: profile_table_t
        real(dp), allocatable :: s(:), ne(:), te(:), ti(:), zeff(:)
    end type profile_table_t

    public :: has_group, read_chords, read_thomson_volumes, read_ece_lines, &
        read_mse_lines, read_profile_table

contains

    logical function has_group(path, group)
        character(len=*), intent(in) :: path, group
        integer(c_int) :: ncid, gid
        has_group = .false.
        if (nc_open(cstr(path), NC_NOWRITE, ncid) /= NC_NOERR) return
        has_group = nc_inq_ncid(ncid, cstr(group), gid) == NC_NOERR
        if (nc_close(ncid) /= NC_NOERR) has_group = .false.
    end function has_group

    subroutine read_chords(path, group, npanel, chords, ierr, message)
        !! Chords from start to end; npanel panels of 4 Gauss points each.
        character(len=*), intent(in) :: path, group
        integer, intent(in) :: npanel
        type(chord_t), allocatable, intent(out) :: chords(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), allocatable :: a(:, :), b(:, :)
        character(len=:), allocatable :: labels(:)
        integer(c_int) :: ncid, gid
        integer :: c

        call open_group(path, group, ncid, gid, ierr, message)
        if (ierr /= 0) return
        call get_matrix(gid, 'start', a, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'end', b, ierr, message)
        if (ierr == 0) call get_labels(gid, size(a, 2), labels, ierr, message)
        if (ierr == 0) then
            if (size(b, 2) /= size(a, 2)) call fail('start and end differ in length')
        end if
        if (ierr == 0) then
            allocate(chords(size(a, 2)))
            do c = 1, size(a, 2)
                chords(c)%label = trim(labels(c))
                chords(c)%origin = a(:, c)
                chords(c)%direction = b(:, c) - a(:, c)
                chords(c)%length = norm2(b(:, c) - a(:, c))
                chords(c)%npanel = npanel
                chords(c)%ngauss = 4
                if (chords(c)%length <= 0.0_dp) call fail('chord '//trim(labels(c))//' has zero length')
            end do
        end if
        c = nc_close(ncid)
    contains
        subroutine fail(text)
            character(len=*), intent(in) :: text
            ierr = 1
            message = trim(path)//':'//group//': '//text
        end subroutine fail
    end subroutine read_chords

    subroutine read_thomson_volumes(path, npoints, volumes, ierr, message)
        character(len=*), intent(in) :: path
        integer, intent(in) :: npoints
        type(thomson_volume_t), allocatable, intent(out) :: volumes(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), allocatable :: center(:, :), direction(:, :), length(:)
        character(len=:), allocatable :: labels(:)
        integer(c_int) :: ncid, gid
        integer :: c

        call open_group(path, 'thomson', ncid, gid, ierr, message)
        if (ierr /= 0) return
        call get_matrix(gid, 'center', center, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'direction', direction, ierr, message)
        if (ierr == 0) call get_vector(gid, 'length', length, ierr, message)
        if (ierr == 0) call get_labels(gid, size(center, 2), labels, ierr, message)
        if (ierr == 0) then
            if (size(direction, 2) /= size(center, 2) .or. size(length) /= size(center, 2)) then
                ierr = 1
                message = trim(path)//':thomson: center, direction and length differ in length'
            end if
        end if
        if (ierr == 0) then
            allocate(volumes(size(center, 2)))
            do c = 1, size(center, 2)
                if (length(c) > 0.0_dp .and. norm2(direction(:, c)) <= 0.0_dp) then
                    ierr = 1
                    message = trim(path)//':thomson: volume '//trim(labels(c))//' has no direction'
                    exit
                end if
                volumes(c) = thomson_volume_t(trim(labels(c)), center(:, c), direction(:, c), &
                    length(c), npoints)
            end do
        end if
        c = nc_close(ncid)
    end subroutine read_thomson_volumes

    subroutine read_ece_lines(path, lines, ierr, message)
        character(len=*), intent(in) :: path
        type(ece_line_t), allocatable, intent(out) :: lines(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), allocatable :: origin(:, :), direction(:, :), length(:), frequency(:), &
            harmonic(:)
        character(len=:), allocatable :: labels(:)
        integer(c_int) :: ncid, gid
        integer :: c, n

        call open_group(path, 'ece', ncid, gid, ierr, message)
        if (ierr /= 0) return
        call get_vector(gid, 'frequency', frequency, ierr, message)
        if (ierr == 0) call get_vector(gid, 'harmonic', harmonic, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'origin', origin, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'direction', direction, ierr, message)
        if (ierr == 0) call get_vector(gid, 'length', length, ierr, message)
        if (ierr == 0) call get_labels(gid, size(frequency), labels, ierr, message)
        if (ierr == 0) then
            n = size(frequency)
            if (any([size(harmonic), size(origin, 2), size(direction, 2), size(length)] /= n)) then
                ierr = 1
                message = trim(path)//':ece: variables differ in length'
            end if
        end if
        if (ierr == 0) then
            allocate(lines(n))
            do c = 1, n
                lines(c)%label = trim(labels(c))
                lines(c)%frequency = frequency(c)
                lines(c)%harmonic = nint(harmonic(c))
                lines(c)%origin = origin(:, c)
                lines(c)%direction = direction(:, c)
                lines(c)%length = length(c)
                if (frequency(c) <= 0.0_dp .or. length(c) <= 0.0_dp .or. lines(c)%harmonic < 1 &
                    .or. norm2(direction(:, c)) <= 0.0_dp) then
                    ierr = 1
                    message = trim(path)//':ece: invalid channel '//trim(labels(c))
                    exit
                end if
            end do
        end if
        c = nc_close(ncid)
    end subroutine read_ece_lines

    subroutine read_mse_lines(path, lines, ierr, message)
        character(len=*), intent(in) :: path
        type(mse_line_t), allocatable, intent(out) :: lines(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), allocatable :: o(:, :), d(:, :), bo(:, :), bd(:, :)
        character(len=:), allocatable :: labels(:)
        integer(c_int) :: ncid, gid
        integer :: c, n

        call open_group(path, 'mse', ncid, gid, ierr, message)
        if (ierr /= 0) return
        call get_matrix(gid, 'origin', o, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'direction', d, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'beam_origin', bo, ierr, message)
        if (ierr == 0) call get_matrix(gid, 'beam_direction', bd, ierr, message)
        if (ierr == 0) call get_labels(gid, size(o, 2), labels, ierr, message)
        if (ierr == 0) then
            n = size(o, 2)
            if (any([size(d, 2), size(bo, 2), size(bd, 2)] /= n)) then
                ierr = 1
                message = trim(path)//':mse: variables differ in length'
            end if
        end if
        if (ierr == 0) then
            allocate(lines(n))
            do c = 1, n
                lines(c) = mse_line_t(trim(labels(c)), o(:, c), d(:, c), bo(:, c), bd(:, c))
            end do
        end if
        c = nc_close(ncid)
    end subroutine read_mse_lines

    subroutine read_profile_table(path, table, ierr, message)
        character(len=*), intent(in) :: path
        type(profile_table_t), intent(out) :: table
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer(c_int) :: ncid, varid
        integer :: c

        ierr = 1
        if (nc_open(cstr(path), NC_NOWRITE, ncid) /= NC_NOERR) then
            message = 'cannot open profile table '//trim(path)
            return
        end if
        call get_vector(ncid, 's', table%s, ierr, message)
        if (ierr == 0) call get_vector(ncid, 'ne', table%ne, ierr, message)
        if (ierr == 0) call get_vector(ncid, 'te', table%te, ierr, message)
        if (ierr == 0) call get_vector(ncid, 'ti', table%ti, ierr, message)
        if (ierr == 0 .and. nc_inq_varid(ncid, cstr('zeff'), varid) == NC_NOERR) &
            call get_vector(ncid, 'zeff', table%zeff, ierr, message)
        c = nc_close(ncid)
        if (ierr /= 0) then
            message = trim(path)//': '//message
            return
        end if
        if (any([size(table%ne), size(table%te), size(table%ti)] /= size(table%s))) then
            ierr = 1
            message = trim(path)//': profiles differ in length from s'
        else if (any(table%s < 0.0_dp) .or. any(table%s > 1.0_dp)) then
            ierr = 1
            message = trim(path)//': s outside [0, 1]'
        end if
        if (ierr == 0 .and. allocated(table%zeff)) then
            if (size(table%zeff) /= size(table%s)) then
                ierr = 1
                message = trim(path)//': zeff differs in length from s'
            end if
        end if
    end subroutine read_profile_table

    subroutine open_group(path, group, ncid, gid, ierr, message)
        character(len=*), intent(in) :: path, group
        integer(c_int), intent(out) :: ncid, gid
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer :: c
        ierr = 1
        if (nc_open(cstr(path), NC_NOWRITE, ncid) /= NC_NOERR) then
            message = 'cannot open geometry file '//trim(path)
            return
        end if
        if (nc_inq_ncid(ncid, cstr(group), gid) /= NC_NOERR) then
            message = trim(path)//': no group '//group
            c = nc_close(ncid)
            return
        end if
        ierr = 0
        message = ''
    end subroutine open_group

    subroutine var_shape(ncid, name, varid, shape_, ierr, message)
        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer(c_int), intent(out) :: varid
        integer, allocatable, intent(out) :: shape_(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        ierr = 1
        if (nc_inq_varid(ncid, cstr(name), varid) /= NC_NOERR) then
            message = 'missing variable '//name
            return
        end if
        if (nc_shape(ncid, varid, shape_) /= NC_NOERR) then
            message = 'cannot inquire variable '//name
            return
        end if
        ierr = 0
        message = ''
    end subroutine var_shape

    subroutine get_vector(ncid, name, v, ierr, message)
        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(out) :: v(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer, allocatable :: shape_(:)
        integer(c_int) :: varid
        call var_shape(ncid, name, varid, shape_, ierr, message)
        if (ierr /= 0) return
        ierr = 1
        if (size(shape_) /= 1) then
            message = name//' is not one-dimensional'
            return
        end if
        allocate(v(shape_(1)))
        if (nc_get_var_double(ncid, varid, v) /= NC_NOERR) then
            message = 'cannot read '//name
            return
        end if
        if (any(v /= v)) then
            message = name//' contains NaN'
            return
        end if
        ierr = 0
    end subroutine get_vector

    subroutine get_matrix(ncid, name, m, ierr, message)
        !! A (n, 3) variable of the file as m(3, n).
        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(out) :: m(:, :)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer, allocatable :: shape_(:)
        integer(c_int) :: varid
        call var_shape(ncid, name, varid, shape_, ierr, message)
        if (ierr /= 0) return
        ierr = 1
        if (size(shape_) /= 2) then
            message = name//' is not (n, 3)'
            return
        end if
        if (shape_(2) /= 3) then
            message = name//' is not (n, 3)'
            return
        end if
        allocate(m(3, shape_(1)))
        if (nc_get_var_double(ncid, varid, m) /= NC_NOERR) then
            message = 'cannot read '//name
            return
        end if
        if (any(m /= m)) then
            message = name//' contains NaN'
            return
        end if
        ierr = 0
    end subroutine get_matrix

    subroutine get_labels(ncid, n, labels, ierr, message)
        !! Character variable label(n, strlen); numbered labels if absent.
        integer(c_int), intent(in) :: ncid
        integer, intent(in) :: n
        character(len=:), allocatable, intent(out) :: labels(:)
        integer, intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer, allocatable :: shape_(:)
        integer(c_int) :: varid
        integer :: k
        ierr = 0
        message = ''
        if (nc_inq_varid(ncid, cstr('label'), varid) /= NC_NOERR) then
            allocate(character(len=12) :: labels(n))
            do k = 1, n
                write(labels(k), '(I0)') k
            end do
            return
        end if
        call var_shape(ncid, 'label', varid, shape_, ierr, message)
        if (ierr /= 0) return
        ierr = 1
        if (size(shape_) /= 2) then
            message = 'label is not (n, strlen) characters'
            return
        end if
        if (shape_(1) /= n) then
            message = 'label has the wrong length'
            return
        end if
        allocate(character(len=shape_(2)) :: labels(n))
        if (nc_read_text(ncid, varid, shape_(2), n, labels) /= NC_NOERR) then
            message = 'cannot read label'
            return
        end if
        ierr = 0
    end subroutine get_labels
end module tiago_profile_geometry_io
