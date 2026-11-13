program test_profile_points
    !! Independent exp(cubic) values and slopes, repeated island labels,
    !! and malformed inputs through the command-line reader.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tiago_profile_geometry_io, only: profile_table_t, read_profile_table
    use tiago_sha256, only: sha256_file
    use tiago_netcdf_c, only: NC_NOERR, NC_CLOBBER, NC_NETCDF4, NC_NOWRITE, NC_DOUBLE, &
        nc_create, nc_open, nc_close, nc_def_dim, nc_def_var, nc_put_var_double, &
        nc_inq_varid, nc_get_var_double, cstr
    implicit none

    real(dp), parameter :: points(9) = [0.0_dp, 0.11_dp, 0.25_dp, 0.25_dp, &
        0.25_dp, 0.63_dp, 0.9_dp, 1.0_dp, 1.04_dp]
    character(len=4096) :: executable, workdir
    character(len=:), allocatable :: prefix
    type(profile_table_t) :: table
    character(len=:), allocatable :: message
    integer :: failures, unit, i, bad, code, ierr
    character(len=64) :: saved_hash
    real(dp) :: ne(9), te(9), slope(9), values(9), expected

    failures = 0
    call get_command_argument(1, executable)
    call get_command_argument(2, workdir)
    prefix = trim(workdir)//'/point-profiles'
    call write_profiles(0)
    call write_labels()
    call run(.true.)
    call read_output('label', values)
    call check(all(values == points), 'point order and repeated labels')
    call read_output('ne', ne)
    call read_output('te', te)
    do i = 1, size(points)
        call close_to(ne(i), ne_of(points(i)), 2.0e-12_dp, 'density against exp(cubic)')
        call close_to(te(i), te_of(points(i)), 2.0e-12_dp, 'temperature against exp(cubic)')
    end do
    call check(ne(3) == ne(4) .and. ne(4) == ne(5), 'density plateau at one label')
    call check(te(3) == te(4) .and. te(4) == te(5), 'temperature plateau at one label')
    call read_output('d_ne_dlabel', slope)
    do i = 1, size(points)
        expected = ne_of(points(i)) * (-0.4_dp - 2.2_dp * points(i) + 0.9_dp * points(i)**2)
        if (points(i) > 1.0_dp) expected = -ne_of(points(i)) / 0.05_dp
        call close_to(slope(i), expected, 2.0e-11_dp, 'density label slope against analytic derivative')
    end do
    call read_output('d_te_dlabel', slope)
    do i = 1, size(points)
        expected = te_of(points(i)) * (-1.0_dp - 1.8_dp * points(i) + 0.6_dp * points(i)**2)
        if (points(i) > 1.0_dp) expected = -te_of(points(i)) / 0.05_dp
        call close_to(slope(i), expected, 2.0e-11_dp, 'temperature label slope against analytic derivative')
    end do
    call read_output('ti', values)
    call check(maxval(abs(values(:8) - 1500.0_dp)) < 1.0e-9_dp, 'constant ion temperature')
    call close_to(values(9), 1500.0_dp * exp(-0.04_dp / 0.05_dp), 2.0e-12_dp, 'ion SOL decay')
    call read_output('zeff', values)
    call check(maxval(abs(values - 2.0_dp)) < 1.0e-12_dp, 'constant effective charge')
    call read_output('d_ti_dlabel', values)
    call check(maxval(abs(values(:8))) < 1.0e-8_dp, 'zero ion-temperature slope')
    call close_to(values(9), -1500.0_dp * exp(-0.04_dp / 0.05_dp) / 0.05_dp, &
        2.0e-12_dp, 'ion SOL slope')
    call read_output('d_zeff_dlabel', values)
    call check(maxval(abs(values)) < 1.0e-11_dp, 'zero charge slope')
    do bad = 1, 4
        call write_profiles(bad)
        if (bad == 4) then
            call read_profile_table(prefix//'.nc', table, ierr, message)
            call check(ierr /= 0, 'profile reader rejects infinity')
        end if
        call run(.false.)
    end do
    call write_profiles(0)
    do bad = 1, 9
        open(newunit=unit, file=prefix//'.labels', status='replace')
        select case (bad)
        case (1)
            write(unit, '(A)') 'NaN'
        case (2)
            write(unit, '(A)') '-0.1'
        case (3)
            write(unit, '(A)') 'invalid'
        case (5)
            write(unit, '(A)') '0.25 0.75'
        case (6)
            write(unit, '(A)') '0.25 garbage'
        case (7)
            write(unit, '(A)') '/'
        case (8)
            write(unit, '(A)') ','
        case (9)
            write(unit, '(A)') '2*0.25'
        end select
        close(unit)
        call run(.false.)
    end do
    call write_labels()
    saved_hash = sha256_file(prefix//'.nc')
    call run(.false., prefix//'.nc')
    call check(sha256_file(prefix//'.nc') == saved_hash, 'profile collision preserves input')
    call execute_command_line('ln -sf "'//prefix//'.nc" "'//prefix//'-alias.nc"', exitstat=code)
    call check(code == 0, 'create symbolic input alias')
    call run(.false., prefix//'-alias.nc')
    call check(sha256_file(prefix//'.nc') == saved_hash, 'symbolic alias preserves input')
    call execute_command_line('rm "'//prefix//'-alias.nc"', exitstat=code)
    call execute_command_line('ln "'//prefix//'.nc" "'//prefix//'-alias.nc"', exitstat=code)
    call check(code == 0, 'create hard input alias')
    call run(.false., prefix//'-alias.nc')
    call check(sha256_file(prefix//'.nc') == saved_hash, 'hard alias preserves input')
    call execute_command_line('rm "'//prefix//'-alias.nc"', exitstat=code)
    saved_hash = sha256_file(prefix//'.labels')
    call run(.false., prefix//'.labels')
    call check(sha256_file(prefix//'.labels') == saved_hash, 'label collision preserves input')
    if (failures > 0) then
        write(error_unit, '(I0,A)') failures, ' check(s) failed'
        error stop 1
    end if
    print '(A)', 'test_profile_points passed'

contains

    pure real(dp) function ne_of(s)
        real(dp), intent(in) :: s
        real(dp) :: x
        x = min(s, 1.0_dp)
        ne_of = 5.0e19_dp * exp(-0.4_dp * x - 1.1_dp * x**2 + 0.3_dp * x**3 &
            - max(s - 1.0_dp, 0.0_dp) / 0.05_dp)
    end function ne_of

    pure real(dp) function te_of(s)
        real(dp), intent(in) :: s
        real(dp) :: x
        x = min(s, 1.0_dp)
        te_of = 3.0e3_dp * exp(-x - 0.9_dp * x**2 + 0.2_dp * x**3 &
            - max(s - 1.0_dp, 0.0_dp) / 0.05_dp)
    end function te_of

    subroutine check(ok, text)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: text
        if (ok) return
        failures = failures + 1
        write(error_unit, '(A)') 'FAIL: '//text
    end subroutine check

    subroutine close_to(a, b, tolerance, text)
        real(dp), intent(in) :: a, b, tolerance
        character(len=*), intent(in) :: text
        call check(abs(a - b) <= tolerance * abs(b), text)
    end subroutine close_to

    subroutine write_labels()
        integer :: u, k
        open(newunit=u, file=prefix//'.labels', status='replace')
        do k = 1, size(points)
            write(u, '(ES24.16E3)') points(k)
        end do
        close(u)
    end subroutine write_labels

    subroutine run(success, destination)
        logical, intent(in) :: success
        character(len=*), optional, intent(in) :: destination
        integer :: u, exit_code, command_code
        open(newunit=u, file=prefix//'.nml', status='replace')
        write(u, '(A)') '&profile_points'
        write(u, '(A)') 'profiles="'//prefix//'.nc", labels="'//prefix//'.labels"'
        if (present(destination)) then
            write(u, '(A)') 'output="'//destination//'", profile_knots=4'
        else
            write(u, '(A)') 'output="'//prefix//'-output.nc", profile_knots=4'
        end if
        write(u, '(A)') '/'
        close(u)
        call execute_command_line('"'//trim(executable)//'" "'//prefix//'.nml" > "'// &
            prefix//'.log" 2>&1', exitstat=exit_code, cmdstat=command_code)
        call check(command_code == 0, 'execute point-profile command')
        call check((exit_code == 0) .eqv. success, 'point-profile command exit status')
    end subroutine run

    subroutine write_profiles(bad)
        integer, intent(in) :: bad
        real(dp) :: s(33), n(33), t(33), ti(33), z(33)
        integer(c_int) :: id, dim, code
        integer :: k
        s = [(real(k, dp) / 32, k = 0, 32)]
        n = [(ne_of(s(k)), k = 1, size(s))]
        t = [(te_of(s(k)), k = 1, size(s))]
        ti = 1500.0_dp
        z = 2.0_dp
        select case (bad)
        case (1)
            t(7) = ieee_value(0.0_dp, ieee_quiet_nan)
        case (2)
            s(7) = ieee_value(0.0_dp, ieee_quiet_nan)
        case (3)
            z(7) = ieee_value(0.0_dp, ieee_quiet_nan)
        case (4)
            n(7) = ieee_value(0.0_dp, ieee_positive_inf)
        end select
        code = nc_create(cstr(prefix//'.nc'), ior(NC_CLOBBER, NC_NETCDF4), id)
        call check(code == NC_NOERR, 'create profile fixture')
        code = nc_def_dim(id, cstr('s'), int(size(s), c_size_t), dim)
        call put_fixture(id, dim, 's', s)
        call put_fixture(id, dim, 'ne', n)
        call put_fixture(id, dim, 'te', t)
        call put_fixture(id, dim, 'ti', ti)
        call put_fixture(id, dim, 'zeff', z)
        code = nc_close(id)
    end subroutine write_profiles

    subroutine put_fixture(id, dim, name, data)
        integer(c_int), intent(in) :: id, dim
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: data(:)
        integer(c_int) :: var, code
        code = nc_def_var(id, cstr(name), NC_DOUBLE, 1_c_int, [dim], var)
        call check(code == NC_NOERR, 'define fixture '//name)
        code = nc_put_var_double(id, var, data)
        call check(code == NC_NOERR, 'write fixture '//name)
    end subroutine put_fixture

    subroutine read_output(name, data)
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: data(:)
        integer(c_int) :: id, var, code
        code = nc_open(cstr(prefix//'-output.nc'), NC_NOWRITE, id)
        call check(code == NC_NOERR, 'open point-profile result')
        code = nc_inq_varid(id, cstr(name), var)
        call check(code == NC_NOERR, 'find '//name)
        code = nc_get_var_double(id, var, data)
        call check(code == NC_NOERR, 'read '//name)
        code = nc_close(id)
    end subroutine read_output
end program test_profile_points
