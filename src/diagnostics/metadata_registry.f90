module tiago_metadata_registry
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32
    use tiago_diagnostic_types, only: metadata_entry_t
    implicit none
    private

    public :: load_metadata_registry
    public :: find_metadata_entry

contains
    subroutine load_metadata_registry(path, entries, ierr, message)
        character(len=*), intent(in) :: path
        type(metadata_entry_t), allocatable, intent(out) :: entries(:)
        integer(i32), intent(out) :: ierr
        character(len=:), allocatable, intent(out) :: message

        integer :: unit
        integer :: ios
        character(len=256) :: line
        character(len=256) :: trimmed
        integer :: active_index

        ierr = 0_i32
        message = ''
        active_index = 0
        if (allocated(entries)) deallocate(entries)
        allocate(entries(0))

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 10_i32
            message = 'unable to open registry file: '//trim(path)
            return
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            trimmed = adjustl(line)
            if (len_trim(trimmed) == 0) cycle

            if (index(trimmed, '"name"') > 0) then
                call push_entry(entries)
                active_index = size(entries)
                call assign_name(entries(active_index), trimmed)
            else if (active_index == 0) then
                cycle
            else if (index(trimmed, '"type"') > 0) then
                entries(active_index)%diag_type = extract_string_value(trimmed)
            else if (index(trimmed, '"effective_area_m2"') > 0) then
                entries(active_index)%effective_area = extract_real_value(trimmed)
            else if (index(trimmed, '"segments"') > 0) then
                entries(active_index)%segments = int( &
                    extract_real_value(trimmed), kind=i32)
            else if (index(trimmed, '"reference_orientation"') > 0) then
                entries(active_index)%reference_orientation = &
                    extract_orientation(trimmed)
            end if
        end do

        close(unit)

        call prune_incomplete(entries)
    end subroutine load_metadata_registry

    logical function find_metadata_entry(entries, name, entry) result(found)
        type(metadata_entry_t), allocatable, intent(in) :: entries(:)
        character(len=*), intent(in) :: name
        type(metadata_entry_t), intent(out) :: entry

        integer :: i

        found = .false.
        entry = metadata_entry_t()
        if (.not. allocated(entries)) return
        do i = 1, size(entries)
            if (trim(entries(i)%name) == trim(name)) then
                entry = entries(i)
                found = .true.
                return
            end if
        end do
    end function find_metadata_entry

    pure function extract_string_value(line) result(value)
        character(len=*), intent(in) :: line
        character(len=:), allocatable :: value
        integer :: first
        integer :: second
        integer :: colon

        colon = index(line, ':')
        first = index(line(colon + 1:), '"') + colon
        second = index(line(first + 1:), '"') + first
        if (first <= colon .or. second <= first) then
            value = ''
        else
            value = line(first + 1:second - 1)
        end if
    end function extract_string_value

    pure function extract_real_value(line) result(val)
        character(len=*), intent(in) :: line
        real(dp) :: val
        character(len=64) :: buffer
        integer :: colon
        integer :: stop

        colon = index(line, ':')
        stop = len_trim(line)
        if (index(line, ',') > 0) stop = index(line, ',') - 1
        buffer = line(colon + 1:stop)
        read(buffer, *) val
    end function extract_real_value

    pure function extract_orientation(line) result(vec)
        character(len=*), intent(in) :: line
        real(dp) :: vec(3)
        character(len=:), allocatable :: raw
        integer :: start_idx
        integer :: end_idx

        vec = 0.0_dp
        start_idx = index(line, '[')
        end_idx = index(line, ']')
        if (start_idx <= 0 .or. end_idx <= start_idx) return
        raw = line(start_idx + 1:end_idx - 1)
        read(raw, *) vec
    end function extract_orientation

    subroutine push_entry(entries)
        type(metadata_entry_t), allocatable, intent(inout) :: entries(:)
        type(metadata_entry_t), allocatable :: tmp(:)

        if (.not. allocated(entries)) then
            allocate(entries(1))
            call reset_entry(entries(1))
        else
            allocate(tmp(size(entries) + 1))
            tmp(1:size(entries)) = entries
            call move_alloc(tmp, entries)
            call reset_entry(entries(size(entries)))
        end if
    end subroutine push_entry

    subroutine reset_entry(entry)
        type(metadata_entry_t), intent(inout) :: entry
        entry%name = ''
        entry%diag_type = ''
        entry%effective_area = 0.0_dp
        entry%segments = 0_i32
        entry%reference_orientation = 0.0_dp
    end subroutine reset_entry

    subroutine assign_name(entry, line)
        type(metadata_entry_t), intent(inout) :: entry
        character(len=*), intent(in) :: line
        entry%name = extract_string_value(line)
    end subroutine assign_name

    subroutine prune_incomplete(entries)
        type(metadata_entry_t), allocatable, intent(inout) :: entries(:)
        type(metadata_entry_t), allocatable :: tmp(:)
        integer :: i
        integer :: kept

        if (.not. allocated(entries)) return
        kept = 0
        if (size(entries) == 0) return
        allocate(tmp(size(entries)))
        do i = 1, size(entries)
            if (len_trim(entries(i)%name) == 0) cycle
            if (len_trim(entries(i)%diag_type) == 0) cycle
            kept = kept + 1
            tmp(kept) = entries(i)
        end do
        if (kept == size(entries)) then
            deallocate(tmp)
            return
        end if
        if (kept == 0) then
            deallocate(entries)
            allocate(entries(0))
            deallocate(tmp)
            return
        end if
        block
            type(metadata_entry_t), allocatable :: compact(:)
            allocate(compact(kept))
            compact = tmp(1:kept)
            call move_alloc(compact, entries)
        end block
        deallocate(tmp)
    end subroutine prune_incomplete
end module tiago_metadata_registry
