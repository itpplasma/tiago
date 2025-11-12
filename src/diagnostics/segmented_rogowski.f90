module tiago_segmented_rogowski
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file
    implicit none
    private

    public :: read_segmented_rogowski_file
    public :: lint_segmented_rogowski

contains
    subroutine read_segmented_rogowski_file(path, diagnostics, ierr, message, &
            default_area)
        character(len=*), intent(in) :: path
        type(segmented_rogowski_t), allocatable, intent(out) :: diagnostics(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), intent(in), optional :: default_area

        type(flux_loop_t), allocatable :: loops(:)
        integer(i32) :: i
        real(dp) :: area_value

        ierr = 0_i32
        message = ''
        if (allocated(diagnostics)) deallocate(diagnostics)

        call read_flux_loop_file(path, loops, ierr, message)
        if (ierr /= 0_i32) return

        allocate(diagnostics(size(loops)))
        area_value = 0.0_dp
        if (present(default_area)) area_value = default_area

        do i = 1, size(loops)
            diagnostics(i)%label = loops(i)%label
            diagnostics(i)%segments = max(1_i32, &
                size(loops(i)%points) - 1_i32)
            diagnostics(i)%effective_area = area_value
            diagnostics(i)%turn_scale = 1.0_dp
            allocate(diagnostics(i)%path(size(loops(i)%points)))
            diagnostics(i)%path = loops(i)%points
        end do
    end subroutine read_segmented_rogowski_file

    logical function lint_segmented_rogowski(diagnostics, report) result(ok)
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        character(len=:), allocatable, intent(out) :: report

        integer :: i

        ok = .true.
        report = ''

        if (.not. allocated(diagnostics)) then
            call append_line(report, 'segmented Rogowski array not allocated')
            ok = .false.
            return
        end if

        if (size(diagnostics) == 0) then
            call append_line(report, 'segmented Rogowski array is empty')
            ok = .false.
            return
        end if

        do i = 1, size(diagnostics)
            if (diagnostics(i)%segments <= 0) then
                call append_line(report, build_issue(i, &
                    'segment count must be positive'))
                ok = .false.
            end if
            if (.not. allocated(diagnostics(i)%path)) then
                call append_line(report, build_issue(i, 'path coordinates missing'))
                ok = .false.
            end if
            if (diagnostics(i)%effective_area <= 0.0_dp) then
                call append_line(report, build_issue(i, &
                    'effective area must be provided and positive'))
                ok = .false.
            end if
        end do
    end function lint_segmented_rogowski

    subroutine append_line(buffer, text)
        character(len=:), allocatable, intent(inout) :: buffer
        character(len=*), intent(in) :: text

        if (.not. allocated(buffer)) then
            buffer = text
        else
            buffer = buffer // new_line('a') // text
        end if
    end subroutine append_line

    pure function build_issue(index, text) result(line)
        integer, intent(in) :: index
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: line

        character(len=32) :: idx_buffer
        write(idx_buffer, '(I0)') index
        line = 'segrog ' // trim(idx_buffer) // ': ' // trim(text)
    end function build_issue
end module tiago_segmented_rogowski
