program tiago_profile_points
    !! Ideal point profiles at supplied field-line labels; no spatial mapping.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, iostat_end, iostat_eor
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_kinetic_profiles, only: N_PROFILES, P_NE, P_TE, P_TI, P_ZEFF
    use tiago_profile_geometry_io, only: profile_table_t, read_profile_table
    use tiago_synthetic_forward, only: synth_state_t, set_profile_table
    use tiago_netcdf_c, only: NC_NOERR, NC_CLOBBER, NC_NETCDF4, NC_GLOBAL, NC_DOUBLE, &
        nc_create, nc_close, nc_def_dim, nc_def_var, nc_put_var_double, &
        nc_put_att_double, put_text_att, cstr
    use tiago_sha256, only: sha256_file
    use tiago_revision, only: tiago_source_revision
    implicit none

    character(len=4096) :: profiles = '', labels = '', output = 'points.nc', path
    integer :: profile_knots = 12
    real(dp) :: sol_decay = 0.05_dp
    namelist /profile_points/ profiles, labels, output, profile_knots, sol_decay
    type(synth_state_t) :: state
    type(profile_table_t) :: table
    real(dp), allocatable :: label(:), values(:, :), slopes(:, :)
    character(len=:), allocatable :: message
    character(len=64) :: profiles_hash, labels_hash
    integer :: unit, ios, ierr, i
    integer(c_int) :: ncid, dim

    if (command_argument_count() /= 1) call fail('usage: tiago_profile_points input.nml')
    call get_command_argument(1, path)
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('cannot open '//trim(path))
    read(unit, nml=profile_points, iostat=ios)
    close(unit)
    if (ios /= 0) call fail(trim(path)//': cannot read namelist profile_points')
    if (len_trim(profiles) == 0 .or. len_trim(labels) == 0) &
        call fail('profiles and labels are required')
    if (profile_knots < 1) call fail('profile_knots must be positive')
    if (.not. ieee_is_finite(sol_decay) .or. sol_decay <= 0.0_dp) &
        call fail('sol_decay must be finite and positive')
    call read_labels(trim(labels), label)
    call read_profile_table(trim(profiles), table, ierr, message)
    if (ierr /= 0) call fail(message)
    call state%profiles%init([(real(i, dp) / profile_knots, i = 0, profile_knots)])
    state%profiles%sol_decay = sol_decay
    call set_profile_table(state, table, ierr, message)
    if (ierr /= 0) call fail(trim(profiles)//': '//message)
    allocate(state%z(0), values(N_PROFILES, size(label)), slopes(N_PROFILES, size(label)))
    call state%profiles%evaluate(state%z, label, values, slopes)

    call reject_alias(trim(output), trim(profiles))
    call reject_alias(trim(output), trim(labels))
    call reject_alias(trim(output), trim(path))
    profiles_hash = sha256_file(trim(profiles))
    labels_hash = sha256_file(trim(labels))
    call check(nc_create(cstr(trim(output)), ior(NC_CLOBBER, NC_NETCDF4), ncid), 'create output')
    call check(nc_def_dim(ncid, cstr('point'), int(size(label), c_size_t), dim), 'define points')
    call check(put_text_att(ncid, NC_GLOBAL, 'tiago_revision', tiago_source_revision), 'revision')
    call check(put_text_att(ncid, NC_GLOBAL, 'profiles_path', trim(profiles)), 'profiles path')
    call check(put_text_att(ncid, NC_GLOBAL, 'profiles_sha256', profiles_hash), &
        'profiles hash')
    call check(put_text_att(ncid, NC_GLOBAL, 'labels_path', trim(labels)), 'labels path')
    call check(put_text_att(ncid, NC_GLOBAL, 'labels_sha256', labels_hash), &
        'labels hash')
    call check(nc_put_att_double(ncid, NC_GLOBAL, cstr('profile_knots'), NC_DOUBLE, &
        1_c_size_t, [real(profile_knots, dp)]), 'profile knots')
    call check(nc_put_att_double(ncid, NC_GLOBAL, cstr('sol_decay'), NC_DOUBLE, &
        1_c_size_t, [sol_decay]), 'scrape-off decay length')
    call check(put_text_att(ncid, NC_GLOBAL, 'observable', &
        'ideal local ne and Te; ECE Te assumes optically thick emission at the supplied point'), &
        'observable')
    call put('label', '1', label)
    call put('ne', 'm-3', values(P_NE, :))
    call put('te', 'eV', values(P_TE, :))
    call put('ti', 'eV', values(P_TI, :))
    call put('zeff', '1', values(P_ZEFF, :))
    call put('d_ne_dlabel', 'm-3', slopes(P_NE, :))
    call put('d_te_dlabel', 'eV', slopes(P_TE, :))
    call put('d_ti_dlabel', 'eV', slopes(P_TI, :))
    call put('d_zeff_dlabel', '1', slopes(P_ZEFF, :))
    call check(nc_close(ncid), 'close output')
    print '(A)', 'wrote '//trim(output)

contains

    subroutine fail(text)
        character(len=*), intent(in) :: text
        write(error_unit, '(A)') 'tiago_profile_points: '//text
        error stop 1
    end subroutine fail

    subroutine check(code, operation)
        integer(c_int), intent(in) :: code
        character(len=*), intent(in) :: operation
        if (code /= NC_NOERR) call fail(trim(output)//': cannot '//operation)
    end subroutine check

    subroutine put(name, units, data)
        character(len=*), intent(in) :: name, units
        real(dp), intent(in) :: data(:)
        integer(c_int) :: var
        call check(nc_def_var(ncid, cstr(name), NC_DOUBLE, 1_c_int, [dim], var), 'define '//name)
        call check(put_text_att(ncid, var, 'units', units), 'write units of '//name)
        call check(nc_put_var_double(ncid, var, data), 'write '//name)
    end subroutine put

    subroutine reject_alias(destination, source)
        character(len=*), intent(in) :: destination, source
        integer :: u, io
        logical :: connected

        ! A connected-file inquiry also detects links to the same input.
        open(newunit=u, file=source, status='old', access='stream', &
            form='unformatted', action='read', iostat=io)
        if (io /= 0) call fail('cannot open '//source)
        inquire(file=destination, opened=connected, iostat=io)
        close(u)
        if (io /= 0) call fail('cannot inspect '//destination)
        if (connected) call fail('output aliases input '//source)
    end subroutine reject_alias

    subroutine read_labels(filename, data)
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: data(:)
        real(dp) :: value
        integer :: u, io, n, k
        logical :: done

        open(newunit=u, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) call fail('cannot open '//filename)
        n = 0
        do
            call next_label(u, value, done)
            if (done) exit
            n = n + 1
        end do
        if (n == 0) call fail(filename//': no labels')
        rewind(u)
        allocate(data(n))
        do k = 1, n
            call next_label(u, data(k), done)
        end do
        close(u)
    end subroutine read_labels

    subroutine next_label(u, value, done)
        integer, intent(in) :: u
        real(dp), intent(out) :: value
        logical, intent(out) :: done
        character(len=4096) :: line
        integer :: io, width, first, last

        done = .false.
        do
            read(u, '(A)', advance='no', size=width, iostat=io) line
            if (io == iostat_end) then
                done = .true.
                return
            end if
            if (io /= iostat_eor) call fail(trim(labels)//': label row is too long or unreadable')
            first = 1
            do while (first <= width)
                if (line(first:first) /= ' ' .and. line(first:first) /= achar(9)) exit
                first = first + 1
            end do
            if (first > width) cycle
            last = width
            do while (last > first)
                if (line(last:last) /= ' ' .and. line(last:last) /= achar(9)) exit
                last = last - 1
            end do
            if (verify(line(first:last), '0123456789+-.eEdD') /= 0) &
                call fail(trim(labels)//': expected one numeric label per row')
            read(line(first:last), *, iostat=io) value
            if (io /= 0) call fail(trim(labels)//': malformed label')
            if (.not. ieee_is_finite(value)) call fail(trim(labels)//': non-finite label')
            if (value < 0.0_dp) call fail(trim(labels)//': negative label')
            return
        end do
    end subroutine next_label
end program tiago_profile_points
