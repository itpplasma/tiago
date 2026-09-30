program tiago_cli
    !! tiago diag lint <file> [--kind flux|segrog]
    !! Parses a DIAGNO diagnostic file and checks it for problems.
    use, intrinsic :: iso_fortran_env, only: i32 => int32, error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file, lint_flux_loops
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file, &
        lint_segmented_rogowski
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
        case default
            call die('unknown diagnostic kind: '//kind)
        end select
    end subroutine lint

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
        write(*, '(A)') 'Usage: tiago_cli diag lint <file> [--kind flux|segrog]'
        if (status == 0) stop
        stop 1
    end subroutine usage

    subroutine die(message)
        character(len=*), intent(in) :: message
        write(error_unit, '(A)') trim(message)
        stop 1
    end subroutine die
end program tiago_cli
