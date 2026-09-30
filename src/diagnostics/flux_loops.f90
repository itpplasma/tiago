module tiago_flux_loops
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_diagnostic_types, only: flux_loop_t, loop_point_t
    implicit none
    private

    public :: read_flux_loop_file
    public :: lint_flux_loops
    public :: finalize_flux_signals

contains
    subroutine read_flux_loop_file(path, loops, ierr, message)
        character(len=*), intent(in) :: path
        type(flux_loop_t), allocatable, intent(out) :: loops(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: unit
        integer(i32) :: loop_count

        ierr = 0_i32
        message = ''
        if (allocated(loops)) deallocate(loops)

        call open_flux_unit(path, unit, ierr, message)
        if (ierr /= 0_i32) return

        call read_loop_count(unit, loop_count, ierr, message)
        if (ierr /= 0_i32) then
            call close_file(unit)
            return
        end if

        allocate(loops(loop_count))
        call populate_loops(unit, loops, ierr, message)
        call close_file(unit)
        if (ierr /= 0_i32) return
        call check_idia_references(loops, ierr, message)
    end subroutine read_flux_loop_file

    subroutine check_idia_references(loops, ierr, message)
        type(flux_loop_t), intent(in) :: loops(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message
        integer :: i

        ierr = 0_i32
        message = ''
        do i = 1, size(loops)
            if (-loops(i)%idia > size(loops) .or. -loops(i)%idia == i) then
                ierr = 6_i32
                message = 'loop '//trim(loops(i)%label)// &
                    ': idia refers to a non-existent loop or to itself'
                return
            end if
        end do
    end subroutine check_idia_references

    subroutine finalize_flux_signals(loops, fluxes)
        !! DIAGNO post-processing, applied in DIAGNO's order (diagno_flux.f90):
        !! for each loop in turn, subtract loop |idia| if idia < 0 (its value as
        !! already processed when |idia| < i), then apply the turn scale.
        !! idia = 1 would add the plasma's phiedge, which is zero in vacuum.
        type(flux_loop_t), intent(in) :: loops(:)
        real(dp), intent(inout) :: fluxes(:)
        integer :: i

        do i = 1, size(loops)
            if (loops(i)%idia < 0_i32) fluxes(i) = fluxes(i) - fluxes(-loops(i)%idia)
            fluxes(i) = fluxes(i) * loops(i)%turn_scale
        end do
    end subroutine finalize_flux_signals

    subroutine open_flux_unit(path, unit, ierr, message)
        character(len=*), intent(in) :: path
        integer, intent(out) :: unit
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: ios

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1_i32
            message = 'unable to open flux loop file: '//trim(path)
        else
            ierr = 0_i32
            message = ''
        end if
    end subroutine open_flux_unit

    subroutine read_loop_count(unit, loop_count, ierr, message)
        integer, intent(in) :: unit
        integer(i32), intent(out) :: loop_count
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: ios

        read(unit, *, iostat=ios) loop_count
        if (ios /= 0 .or. loop_count <= 0_i32) then
            ierr = 2_i32
            message = 'invalid loop count in diagnostic file'
        else
            ierr = 0_i32
            message = ''
        end if
    end subroutine read_loop_count

    subroutine populate_loops(unit, loops, ierr, message)
        integer, intent(in) :: unit
        type(flux_loop_t), intent(inout) :: loops(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: ios
        integer(i32) :: i
        character(len=256) :: line

        ierr = 0_i32
        message = ''

        do i = 1, size(loops)
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) then
                ierr = 3_i32
                message = 'missing loop header at index ' // trim(int_to_string(i))
                return
            end if
            call load_loop_header(line, loops(i), ios, message)
            if (ios /= 0) then
                ierr = 4_i32
                return
            end if
            call load_loop_points(unit, loops(i), ios, message)
            if (ios /= 0) then
                ierr = 5_i32
                return
            end if
        end do
    end subroutine populate_loops

    logical function lint_flux_loops(loops, report) result(ok)
        type(flux_loop_t), allocatable, intent(in) :: loops(:)
        character(len=:), allocatable, intent(out) :: report

        integer :: i
        logical :: local_ok
        character(len=:), allocatable :: local_report

        ok = .true.
        report = ''

        if (.not. allocated(loops)) then
            call append_line(report, 'no loops were provided to the linter')
            ok = .false.
            return
        end if

        if (size(loops) == 0) then
            call append_line(report, 'loop array is empty, nothing to validate')
            ok = .false.
            return
        end if

        do i = 1, size(loops)
            call lint_single_loop(loops(i), i, local_ok, local_report)
            if (.not. local_ok) then
                ok = .false.
                call append_line(report, trim(local_report))
            end if
        end do
    end function lint_flux_loops

    subroutine load_loop_header(line, loop, ios, message)
        character(len=*), intent(in) :: line
        type(flux_loop_t), intent(inout) :: loop
        integer, intent(out) :: ios
        character(len=:), allocatable, intent(out) :: message

        integer(i32) :: npts
        integer(i32) :: repeat_flag
        integer(i32) :: subtract_flag  ! DIAGNO idia column
        character(len=128) :: label_buffer

        message = ''
        read(line, *, iostat=ios) npts, repeat_flag, subtract_flag, label_buffer
        if (ios /= 0) then
            message = 'failed to parse loop header: '//trim(line)
            return
        end if

        if (npts <= 1_i32) then
            ios = -1
            message = 'loop must contain at least two points'
            return
        end if

        if (repeat_flag > 1_i32) then
            ios = -1
            message = 'unsupported iflflg > 1 (only 0 and 1 are defined): '//trim(line)
            return
        end if

        loop%idia = subtract_flag
        loop%one_period = (repeat_flag == 1_i32)
        loop%label = trim(label_buffer)
        if (len_trim(loop%label) == 0) then
            loop%label = 'loop_' // trim(adjustl(int_to_string(npts)))
        end if
        allocate(loop%points(npts))
    end subroutine load_loop_header

    subroutine load_loop_points(unit, loop, ios, message)
        integer, intent(in) :: unit
        type(flux_loop_t), intent(inout) :: loop
        integer, intent(out) :: ios
        character(len=:), allocatable, intent(out) :: message

        integer :: j
        real(dp) :: x, y, z

        message = ''
        ios = 0
        do j = 1, size(loop%points)
            read(unit, *, iostat=ios) x, y, z
            if (ios /= 0) then
                message = 'insufficient coordinate rows for loop '//trim(loop%label)
                return
            end if
            loop%points(j)%x = x
            loop%points(j)%y = y
            loop%points(j)%z = z
        end do
    end subroutine load_loop_points

    subroutine lint_single_loop(loop, index, ok, report)
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: index
        logical, intent(out) :: ok
        character(len=:), allocatable, intent(out) :: report

        integer :: npts
        logical :: has_issue

        ok = .true.
        report = ''
        has_issue = .false.
        npts = 0
        if (allocated(loop%points)) npts = size(loop%points)

        if (npts <= 1) then
            call append_line(report, format_issue(index, &
                'loop must contain at least two points'))
            ok = .false.
            has_issue = .true.
        end if

        if (.not. verify_point_coordinates(loop, index, report)) then
            ok = .false.
            has_issue = .true.
        end if

        if (.not. has_issue .and. len_trim(loop%label) == 0) then
            call append_line(report, format_issue(index, 'loop label is empty'))
            ok = .false.
        end if
    end subroutine lint_single_loop

    logical function verify_point_coordinates(loop, index, report) result(ok)
        type(flux_loop_t), intent(in) :: loop
        integer, intent(in) :: index
        character(len=:), allocatable, intent(inout) :: report

        integer :: j
        ok = .true.
        do j = 1, size(loop%points)
            if (.not. coordinates_are_finite(loop%points(j))) then
                ok = .false.
                call append_line(report, format_issue(index, 'point ' // &
                    trim(int_to_string(j)) // ' contains NaN or Inf'))
            end if
        end do
    end function verify_point_coordinates

    logical function coordinates_are_finite(point) result(ok)
        type(loop_point_t), intent(in) :: point

        ok = ieee_is_finite(point%x) .and. ieee_is_finite(point%y) .and. &
            ieee_is_finite(point%z)
    end function coordinates_are_finite



    subroutine append_line(buffer, line)
        character(len=:), allocatable, intent(inout) :: buffer
        character(len=*), intent(in) :: line

        if (.not. allocated(buffer)) then
            buffer = line
        else
            buffer = buffer // new_line('a') // line
        end if
    end subroutine append_line

    pure function format_issue(index, text) result(line)
        integer, intent(in) :: index
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: line

        line = 'loop ' // trim(int_to_string(index)) // ': ' // trim(text)
    end function format_issue

    pure function int_to_string(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text

        character(len=32) :: buffer
        write(buffer, '(I0)') value
        text = trim(buffer)
    end function int_to_string


    subroutine close_file(unit)
        integer, intent(in) :: unit
        integer :: ios
        close(unit, iostat=ios)
    end subroutine close_file
end module tiago_flux_loops
