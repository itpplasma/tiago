module tiago_segmented_rogowski
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_diagnostic_types, only: loop_point_t, segmented_rogowski_t
    use tiago_flux_loops, only: parse_diag_header, check_no_trailing_data
    implicit none
    private

    public :: read_segmented_rogowski_file
    public :: lint_segmented_rogowski

contains
    subroutine read_segmented_rogowski_file(path, diagnostics, ierr, message, &
            default_area)
        !! DIAGNO segmented Rogowski file: header as for flux loops, then rows
        !! "x y z eff_area". Segment j is weighted with eff_area of point j. For
        !! rows without eff_area, default_area / (npts - 1) is used per segment;
        !! without default_area the areas stay zero (checked by the caller).
        character(len=*), intent(in) :: path
        type(segmented_rogowski_t), allocatable, intent(out) :: diagnostics(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        real(dp), intent(in), optional :: default_area

        integer :: unit, ios, i, j
        integer(i32) :: count, npts, flag1, flag2
        character(len=512) :: line
        character(len=:), allocatable :: label
        real(dp) :: row(4)
        logical :: has_area

        ierr = 0_i32
        message = ''
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1_i32
            message = 'unable to open segmented Rogowski file: '//trim(path)
            return
        end if
        read(unit, *, iostat=ios) count
        if (ios /= 0 .or. count <= 0) then
            ierr = 2_i32
            message = 'invalid diagnostic count in '//trim(path)
            close(unit)
            return
        end if

        allocate(diagnostics(count))
        do i = 1, count
            read(unit, '(A)', iostat=ios) line
            if (ios == 0) call parse_diag_header(line, npts, flag1, flag2, label, ios, message)
            if (ios /= 0 .or. npts < 2) then
                ierr = 4_i32
                if (len(message) == 0) message = 'invalid segmented Rogowski header: '//trim(line)
                close(unit)
                return
            end if
            diagnostics(i)%label = label
            allocate(diagnostics(i)%path(npts), diagnostics(i)%segment_area(npts - 1))
            diagnostics(i)%segment_area = 0.0_dp
            if (present(default_area)) diagnostics(i)%segment_area = default_area / real(npts - 1, dp)
            do j = 1, npts
                read(unit, '(A)', iostat=ios) line
                if (ios == 0) call read_row(line, row, has_area, ios)
                if (ios /= 0) then
                    ierr = 5_i32
                    message = 'insufficient coordinate rows for '//trim(label)
                    close(unit)
                    return
                end if
                diagnostics(i)%path(j) = loop_point_t(row(1), row(2), row(3))
                if (has_area .and. j < npts) diagnostics(i)%segment_area(j) = row(4)
            end do
        end do
        call check_no_trailing_data(unit, count, ierr, message)
        close(unit)
    end subroutine read_segmented_rogowski_file

    subroutine read_row(line, row, has_area, ios)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: row(4)
        logical, intent(out) :: has_area
        integer, intent(out) :: ios

        row = 0.0_dp
        read(line, *, iostat=ios) row
        has_area = (ios == 0)
        if (.not. has_area) read(line, *, iostat=ios) row(1:3)
    end subroutine read_row

    logical function lint_segmented_rogowski(diagnostics, report) result(ok)
        type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
        character(len=:), allocatable, intent(out) :: report

        integer :: i, j

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
            if (.not. allocated(diagnostics(i)%path)) then
                call append_line(report, build_issue(i, 'path coordinates missing'))
                ok = .false.
            end if
            if (.not. all(ieee_is_finite([diagnostics(i)%path%x, diagnostics(i)%path%y, &
                    diagnostics(i)%path%z]))) then
                call append_line(report, build_issue(i, 'coordinates contain NaN or Inf'))
                ok = .false.
            end if
            ! Files without an eff_area column get the area from --seg-area at run time.
            if (any(diagnostics(i)%segment_area < 0.0_dp)) then
                call append_line(report, build_issue(i, 'negative effective area'))
                ok = .false.
            end if
            if (any([(diagnostics(j)%label == diagnostics(i)%label, j = 1, i - 1)])) then
                call append_line(report, build_issue(i, 'duplicate label '// &
                    diagnostics(i)%label))
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
