program tiago_cli
    !! tiago diag lint <file> [--kind flux|segrog|sensors]
    !! Parses a DIAGNO diagnostic file, or a sensor-set manifest with all the
    !! geometry files it names, and checks it for problems.
    use, intrinsic :: iso_fortran_env, only: i32 => int32, error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file, lint_flux_loops
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file, &
        lint_segmented_rogowski
    use tiago_sensor_set, only: sensor_set_t, read_sensor_set, sensor_kind_name, &
        kind_combination
    implicit none

    character(len=512) :: arg, path
    character(len=:), allocatable :: kind
    integer :: i

    call get_command_argument(1, arg)
    if (trim(arg) == '--help' .or. trim(arg) == '-h') call usage(0)
    if (command_argument_count() < 3) call usage(1)
    if (trim(arg) /= 'diag') call usage(1)
    call get_command_argument(2, arg)
    if (trim(arg) /= 'lint') call usage(1)
    call get_command_argument(3, path)
    if (trim(path) == '--help' .or. trim(path) == '-h') call usage(0)

    kind = 'flux'
    i = 4
    do while (i <= command_argument_count())
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ('--kind')
            i = i + 1
            if (i > command_argument_count()) call die('--kind requires a value')
            call get_command_argument(i, arg)
            kind = trim(arg)
        case ('--help', '-h')
            call usage(0)
        case default
            call die('unknown option: '//trim(arg))
        end select
        i = i + 1
    end do

    call lint(trim(path), kind)

contains

    subroutine lint(path, kind)
        character(len=*), intent(in) :: path, kind
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        integer(i32) :: ierr
        character(len=:), allocatable :: message, report
        logical :: ok

        select case (kind)
        case ('flux', 'flux_loop')
            call read_flux_loop_file(path, loops, ierr, message)
            if (ierr /= 0_i32) call die(message)
            ok = lint_flux_loops(loops, report)
            call finish(ok, report, size(loops))
        case ('segrog', 'segmented_rogowski')
            call read_segmented_rogowski_file(path, segs, ierr, message)
            if (ierr /= 0_i32) call die(message)
            ok = lint_segmented_rogowski(segs, report)
            call finish(ok, report, size(segs))
        case ('sensors')
            call lint_sensor_set(path)
        case default
            call die('unknown diagnostic kind: '//kind)
        end select
    end subroutine lint

    subroutine lint_sensor_set(path)
        !! Counts per sensor kind, then the loop and Rogowski checks.
        character(len=*), intent(in) :: path
        type(sensor_set_t) :: set
        integer(i32) :: ierr
        character(len=:), allocatable :: message, report
        character(len=16) :: buffer
        logical :: ok
        integer :: k

        call read_sensor_set(path, set, ierr, message)
        if (ierr /= 0_i32) call die(message)
        do k = 1, kind_combination
            if (count(set%kinds == k) == 0) cycle
            write(buffer, '(I0)') count(set%kinds == k)
            write(*, '(A)') sensor_kind_name(k)//': '//trim(buffer)
        end do
        ok = .true.
        report = ''
        if (size(set%loops) > 0) ok = lint_flux_loops(set%loops, report)
        if (.not. ok) call finish(ok, report, 0)
        if (size(set%segs) > 0) ok = lint_segmented_rogowski(set%segs, report)
        call finish(ok, report, set%n_signals())
    end subroutine lint_sensor_set

    subroutine finish(ok, report, count)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: report
        integer, intent(in) :: count
        character(len=16) :: buffer

        if (.not. ok) then
            write(error_unit, '(A)') trim(adjustl(report))
            stop 2
        end if
        write(buffer, '(I0)') count
        write(*, '(A)') 'lint successful for '//trim(buffer)//' diagnostics'
    end subroutine finish

    subroutine usage(status)
        integer, intent(in) :: status
        write(*, '(A)') 'Usage: tiago_cli diag lint <file> [--kind flux|segrog|sensors]'
        if (status == 0) stop
        stop 1
    end subroutine usage

    subroutine die(message)
        character(len=*), intent(in) :: message
        write(error_unit, '(A)') trim(message)
        stop 1
    end subroutine die
end program tiago_cli
