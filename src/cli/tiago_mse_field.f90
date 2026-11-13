program tiago_mse_field
    !! Point MSE angles (docs/mse.md) from the interior field of an equilibrium:
    !! a VMEC wout (stellarator-symmetric or not) or a SPECTRE output file.
    !! Writes label, gamma [rad], the channel midpoint [m] and Cartesian B [T].
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: read_mse_file, mse_geometry, mse_angles
    use tiago_vmec_state, only: vmec_state_t
    use tiago_mse_state, only: mse_state_angles
#ifdef TIAGO_ENABLE_SPECTRE
    use tiago_spectre_source, only: spectre_source_t
#endif
    implicit none

    character(len=1024) :: mse_path = '', wout_path = '', spectre_path = '', &
        output_path = 'tiago_mse_field.csv', arg
    type(mse_channel_t), allocatable :: channels(:)
    real(dp), allocatable :: angles(:), field(:, :), points(:, :)
    real(dp) :: coefficients(6)
    integer(i32) :: ierr
    character(len=:), allocatable :: message
    integer :: i, unit, ios

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
        case ('--output')
            call value(output_path)
        case ('-h', '--help')
            write (*, '(A)') 'usage: tiago_mse_field --mse CHANNELS '// &
                '(--wout WOUT.nc | --spectre SPECTRE.h5) [--output FILE.csv]'
            stop
        case default
            call die('unknown argument '//trim(arg))
        end select
        i = i + 1
    end do
    if (len_trim(mse_path) == 0) call die('--mse is required')
    if ((len_trim(wout_path) > 0) .eqv. (len_trim(spectre_path) > 0)) &
        call die('pass exactly one of --wout and --spectre')

    call read_mse_file(trim(mse_path), channels, ierr, message)
    if (ierr /= 0) call die(message)
    allocate (points(3, size(channels)))
    do i = 1, size(channels)
        call mse_geometry(channels(i), points(:, i), coefficients, ierr, message)
        if (ierr /= 0) call die(trim(channels(i)%label)//': '//message)
    end do

    if (len_trim(wout_path) > 0) then
        call vmec_angles()
    else
        call spectre_angles()
    end if

    open (newunit=unit, file=trim(output_path), action='write', status='replace', &
          iostat=ios)
    if (ios /= 0) call die('cannot write '//trim(output_path))
    write (unit, '(A)') 'label,gamma_rad,x_m,y_m,z_m,b_x_T,b_y_T,b_z_T'
    do i = 1, size(channels)
        write (unit, '(A,7(",",ES24.16E3))') trim(channels(i)%label), angles(i), &
            points(:, i), field(:, i)
    end do
    close (unit)

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
        call mse_angles(channels, field, angles, ierr, message)
        if (ierr /= 0) call die(message)
#else
        call die('--spectre requires a build with TIAGO_ENABLE_SPECTRE=ON')
#endif
    end subroutine spectre_angles

    subroutine value(target)
        character(len=*), intent(out) :: target
        i = i + 1
        if (i > command_argument_count()) call die(trim(arg)//' needs a value')
        call get_command_argument(i, target)
    end subroutine value

    subroutine die(text)
        character(len=*), intent(in) :: text
        write (error_unit, '(A)') 'tiago_mse_field: '//text
        error stop 1
    end subroutine die
end program tiago_mse_field
