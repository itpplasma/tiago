program tiago_cli
    use, intrinsic :: iso_fortran_env, only: dp => real64, &
        i32 => int32, error_unit
    use tiago_diagnostic_types, only: flux_loop_t, &
        segmented_rogowski_t, metadata_entry_t
    use tiago_flux_loops, only: read_flux_loop_file, lint_flux_loops
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file, &
        lint_segmented_rogowski
    use tiago_metadata_registry, only: load_metadata_registry, &
        find_metadata_entry
    implicit none

    call handle_command()

contains
    subroutine handle_command()
        integer :: argc
        character(len=256) :: arg
        character(len=256) :: file_arg
        character(len=:), allocatable :: kind
        character(len=:), allocatable :: metadata_path
        character(len=:), allocatable :: family_name
        real(dp) :: area_override
        integer(i32) :: segment_override
        integer :: i
        logical :: have_file

        kind = 'flux'
        metadata_path = ''
        family_name = ''
        area_override = -1.0_dp
        segment_override = -1_i32
        have_file = .false.

        argc = command_argument_count()
        if (argc < 3) then
            call usage()
            stop 1
        end if

        call get_command_argument(1, arg)
        if (trim(arg) /= 'diag') then
            call usage()
            stop 1
        end if

        call get_command_argument(2, arg)
        if (trim(arg) /= 'lint') then
            call usage()
            stop 1
        end if

        call get_command_argument(3, file_arg)
        if (len_trim(file_arg) == 0) then
            call die('missing diagnostic file path')
        end if
        have_file = .true.

        i = 4
        do while (i <= argc)
            call get_command_argument(i, arg)
            select case (trim(arg))
            case ('--kind')
                i = i + 1
                if (i > argc) call die('--kind requires a value')
                call get_command_argument(i, arg)
                kind = trim(arg)
            case ('--metadata')
                i = i + 1
                if (i > argc) call die('--metadata requires a value')
                call get_command_argument(i, arg)
                metadata_path = trim(arg)
            case ('--family')
                i = i + 1
                if (i > argc) call die('--family requires a value')
                call get_command_argument(i, arg)
                family_name = trim(arg)
            case ('--area-m2')
                i = i + 1
                if (i > argc) call die('--area-m2 requires a numeric value')
                call get_command_argument(i, arg)
                read(arg, *) area_override
            case ('--segments')
                i = i + 1
                if (i > argc) call die('--segments requires an integer value')
                call get_command_argument(i, arg)
                read(arg, *) segment_override
            case ('--help', '-h')
                call usage()
                stop 0
            case default
                call die('unknown option: ' // trim(arg))
            end select
            i = i + 1
        end do

        if (.not. have_file) call die('no diagnostic file provided')

        call execute_lint(file_arg, kind, metadata_path, family_name, &
            area_override, segment_override)
    end subroutine handle_command

    subroutine execute_lint(path, kind, metadata_path, family_name, &
            area_override, segment_override)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: kind
        character(len=*), intent(in) :: metadata_path
        character(len=*), intent(in) :: family_name
        real(dp), intent(in) :: area_override
        integer(i32), intent(in) :: segment_override

        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(metadata_entry_t), allocatable :: registry(:)
        type(metadata_entry_t) :: entry
        logical :: found
        logical :: ok
        integer(i32) :: ierr
        character(len=:), allocatable :: message
        character(len=:), allocatable :: report
        real(dp) :: area_to_use
        integer(i32) :: segments_to_use
        character(len=:), allocatable :: registry_path

        ierr = 0_i32
        area_to_use = area_override
        segments_to_use = segment_override
        registry_path = trim(metadata_path)

        if (len_trim(registry_path) == 0 .and. len_trim(family_name) > 0) then
            registry_path = get_registry_from_env()
        end if

        if (len_trim(registry_path) > 0) then
            call load_metadata_registry(registry_path, registry, ierr, &
                message)
            if (ierr /= 0_i32) call die(message)
        elseif (len_trim(family_name) > 0) then
            call die('--family requires a registry (pass --metadata or set ' // &
                'TIAGO_DIAG_METADATA)')
        end if

        block
            character(len=:), allocatable :: selected_kind
            selected_kind = trim(kind)
            if (len_trim(family_name) > 0) then
                found = find_metadata_entry(registry, trim(family_name), entry)
                if (.not. found) then
                    call die('family '//trim(family_name)//' not found in ' // &
                        'registry')
                end if
                if (len_trim(entry%diag_type) > 0) then
                    selected_kind = trim(entry%diag_type)
                end if
                if (entry%effective_area > 0.0_dp) then
                    area_to_use = entry%effective_area
                end if
                if (entry%segments > 0_i32) segments_to_use = entry%segments
            end if

            select case (trim(selected_kind))
        case ('flux', 'flux_loop')
            call read_flux_loop_file(trim(path), loops, ierr, message)
            if (ierr /= 0_i32) call die(message)
            ok = lint_flux_loops(loops, report)
            call finish_lint(ok, report, size(loops))
        case ('segrog', 'segmented_rogowski')
            if (area_to_use <= 0.0_dp) then
                call die('segmented Rogowski lint requires --area-m2 or ' // &
                    'metadata entry')
            end if
            call read_segmented_rogowski_file(trim(path), segs, ierr, &
                message, area_to_use)
            if (ierr /= 0_i32) call die(message)
            if (segments_to_use > 0_i32) segs(1)%segments = segments_to_use
            call apply_segment_override(segs, segments_to_use)
            ok = lint_segmented_rogowski(segs, report)
            call finish_lint(ok, report, size(segs))
        case default
            call die('unknown diagnostic kind: '//trim(selected_kind))
        end select
        end block
    end subroutine execute_lint

    subroutine finish_lint(ok, report, count)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: report
        integer, intent(in) :: count

        if (ok) then
            write(*, '(A)') 'lint successful for ' // &
                trim(int_to_string(count)) // ' diagnostics'
        else
            write(error_unit, '(A)') trim(report)
            stop 2
        end if
    end subroutine finish_lint

    function get_registry_from_env() result(path)
        character(len=:), allocatable :: path
        character(len=512) :: buffer
        integer :: status
        integer :: length

        buffer = ''
        length = 0
        status = 1
        call get_environment_variable('TIAGO_DIAG_METADATA', buffer, length, status)
        if (status == 0 .and. length > 0) then
            path = buffer(1:length)
        else
            path = ''
        end if
    end function get_registry_from_env

    subroutine apply_segment_override(segments, value)
        type(segmented_rogowski_t), allocatable, intent(inout) :: segments(:)
        integer(i32), intent(in) :: value
        integer :: i

        if (value <= 0_i32) return
        do i = 1, size(segments)
            segments(i)%segments = value
        end do
    end subroutine apply_segment_override

    subroutine usage()
        write(*, '(A)') 'Usage: tiago diag lint <file> [--kind flux|segrog] ' // &
            '[--metadata registry.json]'
        write(*, '(A)') '                    [--family name] [--area-m2 value] ' // &
            '[--segments value]'
        write(*, '(A)') 'Environment: TIAGO_DIAG_METADATA can provide the ' // &
            'default registry path.'
    end subroutine usage

    subroutine die(message)
        character(len=*), intent(in) :: message
        write(error_unit, '(A)') trim(message)
        stop 1
    end subroutine die

    pure function int_to_string(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=32) :: buffer
        write(buffer, '(I0)') value
        text = trim(buffer)
    end function int_to_string
end program tiago_cli
