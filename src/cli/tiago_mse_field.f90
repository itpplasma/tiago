program tiago_mse_field
    !! Point MSE angles (docs/mse.md) from the interior field of an equilibrium:
    !! a VMEC wout, SPECTRE output, vacuum mgrid and/or filament coils.
    !! Writes label, gamma [rad], the channel midpoint [m] and Cartesian B [T].
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: read_mse_file, mse_geometry, mse_angles
    use tiago_vmec_state, only: vmec_state_t
    use tiago_mse_state, only: mse_state_angles
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_field_line_label, only: grid_field_t
    use tiago_grid_field_cartesian, only: grid_field_cartesian
#ifdef TIAGO_ENABLE_SPECTRE
    use tiago_spectre_source, only: spectre_source_t
#endif
    implicit none

    character(len=1024) :: mse_path = '', wout_path = '', spectre_path = '', &
        mgrid_path = '', coil_path = '', extcur_path = '', current_text = '', &
        output_path = 'tiago_mse_field.csv', jacobian_path = '', arg
    type(mse_channel_t), allocatable :: channels(:)
    real(dp), allocatable :: angles(:), field(:, :), points(:, :), field_jacobian(:, :)
    real(dp) :: coefficients(6)
    integer(i32) :: ierr
    character(len=:), allocatable :: message
    integer :: i, unit, ios, sources

    i = 1
    do while (i <= command_argument_count())
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ('--mse')
            call value(mse_path)
        case ('--wout')
            call value(wout_path)
        case ('--spectre')
            call value(spectre_path)
        case ('--mgrid')
            call value(mgrid_path)
        case ('--mgrid-currents')
            call value(current_text)
        case ('--coils')
            call value(coil_path)
        case ('--coil-extcur', '--extcur')
            call value(extcur_path)
        case ('--output')
            call value(output_path)
        case ('--field-jacobian')
            call value(jacobian_path)
        case ('-h', '--help')
            write (*, '(A)') 'usage: tiago_mse_field --mse CHANNELS '// &
                '[--wout WOUT.nc | --spectre SPECTRE.h5 | --mgrid MGRID.nc] '// &
                '[--mgrid-currents LIST] [--coils COILS [--coil-extcur EXTCUR]] '// &
                '[--output FILE.csv] [--field-jacobian JACOBIAN.csv]'
            stop
        case default
            call die('unknown argument '//trim(arg))
        end select
        i = i + 1
    end do
    if (len_trim(mse_path) == 0) call die('--mse is required')
    sources = count([len_trim(wout_path) > 0, len_trim(spectre_path) > 0, &
        len_trim(mgrid_path) > 0])
    if (sources > 1) call die('pass at most one of --wout, --spectre and --mgrid')
    if (sources == 0 .and. len_trim(coil_path) == 0) call die('a field source is required')
    if (len_trim(extcur_path) > 0 .and. len_trim(coil_path) == 0) &
        call die('--coil-extcur requires --coils')
    if ((len_trim(mgrid_path) > 0) .neqv. (len_trim(current_text) > 0)) &
        call die('--mgrid and --mgrid-currents are required together')

    call read_mse_file(trim(mse_path), channels, ierr, message)
    if (ierr /= 0) call die(message)
    allocate (points(3, size(channels)))
    do i = 1, size(channels)
        call mse_geometry(channels(i), points(:, i), coefficients, ierr, message)
        if (ierr /= 0) call die(trim(channels(i)%label)//': '//message)
    end do

    if (len_trim(wout_path) > 0) then
        call vmec_angles()
    else if (len_trim(spectre_path) > 0) then
        call spectre_angles()
    else if (len_trim(mgrid_path) > 0) then
        call grid_field()
    else
        allocate(field(3, size(channels)))
        field = 0.0_dp
    end if
    if (len_trim(coil_path) > 0) call add_coil_field()
    if (len_trim(jacobian_path) > 0) then
        call mse_angles(channels, field, angles, ierr, message, field_jacobian)
    else
        call mse_angles(channels, field, angles, ierr, message)
    end if
    if (ierr /= 0) call die(message)

    call reject_output_alias()
    open (newunit=unit, file=trim(output_path), action='write', status='replace', &
          iostat=ios)
    if (ios /= 0) call die('cannot write '//trim(output_path))
    write (unit, '(A)') 'label,gamma_rad,x_m,y_m,z_m,b_x_T,b_y_T,b_z_T'
    do i = 1, size(channels)
        write (unit, '(A,7(",",ES24.16E3))') trim(channels(i)%label), angles(i), &
            points(:, i), field(:, i)
    end do
    close (unit)
    if (len_trim(jacobian_path) > 0) then
        open (newunit=unit, file=trim(jacobian_path), action='write', status='replace', &
              iostat=ios)
        if (ios /= 0) call die('cannot write '//trim(jacobian_path))
        write (unit, '(A)') 'label,dgamma_dBx_rad_per_T,dgamma_dBy_rad_per_T,dgamma_dBz_rad_per_T'
        do i = 1, size(channels)
            write (unit, '(A,3(",",ES24.16E3))') trim(channels(i)%label), field_jacobian(:, i)
        end do
        close (unit)
    end if

contains

    subroutine vmec_angles()
        type(vmec_state_t) :: state
        integer :: status

        call state%load(trim(wout_path), status, message, vector=.true.)
        if (status /= 0) call die(message)
        call mse_state_angles(state, channels, angles, ierr, message, field=field)
        if (ierr /= 0) call die(message)
    end subroutine vmec_angles

    subroutine spectre_angles()
#ifdef TIAGO_ENABLE_SPECTRE
        type(spectre_source_t) :: source

        ! The casing grid is only used for exterior fields; keep it small.
        call source%init(trim(spectre_path), 8, 8)
        allocate (field(3, size(channels)))
        call source%sample_interior_bfield(points, field)
#else
        call die('--spectre requires a build with TIAGO_ENABLE_SPECTRE=ON')
#endif
    end subroutine spectre_angles

    subroutine add_coil_field()
        type(vacuum_solver_t) :: solver
        real(dp) :: added(3, size(channels))

        call solver%init(trim(coil_path), trim(extcur_path))
        call solver%sample_bfield(points, added)
        field = field + added
        call solver%finalize()
    end subroutine add_coil_field

    subroutine grid_field()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        type(grid_field_t) :: grid
        real(dp) :: currents(128), b(3), radius, phi
        integer :: n, j, code, first, last, comma
        character(len=:), allocatable :: token
        logical :: inside

        n = 0
        first = 1
        do
            comma = index(trim(current_text(first:)), ',')
            last = len_trim(current_text)
            if (comma > 0) last = first + comma - 2
            token = trim(adjustl(current_text(first:last)))
            if (len(token) == 0) call die('--mgrid-currents contains an empty entry')
            if (verify(token, '0123456789+-.eEdD') /= 0) &
                call die('--mgrid-currents requires comma-separated numbers')
            n = n + 1
            if (n > size(currents)) call die('--mgrid-currents has more than 128 entries')
            read(token, *, iostat=code) currents(n)
            if (code /= 0) call die('cannot read --mgrid-currents entry')
            if (.not. ieee_is_finite(currents(n))) call die('--mgrid-currents must be finite')
            if (comma == 0) exit
            first = last + 2
            if (first > len_trim(current_text)) call die('--mgrid-currents contains an empty entry')
        end do
        call grid%load_mgrid(trim(mgrid_path), currents(:n), code, message)
        if (code /= 0) call die(message)
        allocate(field(3, size(channels)))
        do j = 1, size(channels)
            radius = hypot(points(1, j), points(2, j))
            phi = atan2(points(2, j), points(1, j))
            call grid%eval(radius, phi, points(3, j), b, inside)
            if (.not. inside) call die(trim(channels(j)%label)//': midpoint outside mgrid')
            call grid_field_cartesian(phi, b(1), b(2), b(3), &
                field(1, j), field(2, j), field(3, j))
        end do
    end subroutine grid_field

    subroutine reject_output_alias()
        integer :: probe, code, j, k
        logical :: connected
        character(len=1024) :: inputs(6), outputs(2)

        inputs = [mse_path, wout_path, spectre_path, mgrid_path, coil_path, extcur_path]
        outputs = [output_path, jacobian_path]
        do j = 1, size(inputs)
            if (len_trim(inputs(j)) == 0) cycle
            open(newunit=probe, file=trim(inputs(j)), action='read', status='old', &
                 access='stream', form='unformatted', iostat=code)
            if (code /= 0) call die('cannot open '//trim(inputs(j)))
            do k = 1, size(outputs)
                if (len_trim(outputs(k)) == 0) cycle
                inquire(file=trim(outputs(k)), opened=connected)
                if (connected) then
                    close(probe)
                    call die('output aliases input '//trim(inputs(j)))
                end if
            end do
            close(probe)
        end do
        if (len_trim(jacobian_path) == 0) return
        if (trim(output_path) == trim(jacobian_path)) call die('output files must differ')
        ! Opening without replacement preserves any existing data for this check.
        open(newunit=probe, file=trim(output_path), action='readwrite', status='unknown', &
             access='stream', form='unformatted', iostat=code)
        if (code /= 0) call die('cannot open '//trim(output_path))
        inquire(file=trim(jacobian_path), opened=connected)
        close(probe)
        if (connected) call die('output files must differ')
    end subroutine reject_output_alias

    subroutine value(target)
        character(len=*), intent(out) :: target
        integer :: code

        i = i + 1
        if (i > command_argument_count()) call die(trim(arg)//' needs a value')
        call get_command_argument(i, target, status=code)
        if (code /= 0) call die(trim(arg)//' value is too long or unreadable')
    end subroutine value

    subroutine die(text)
        character(len=*), intent(in) :: text
        write (error_unit, '(A)') 'tiago_mse_field: '//text
        error stop 1
    end subroutine die
end program tiago_mse_field
