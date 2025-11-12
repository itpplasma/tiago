module tiago_coil_loader
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, error_unit
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_init, coils_deinit
    implicit none
    private

    real(dp), parameter :: tiny_current = 1.0e-12_dp

    public :: load_coils_into_field

contains

    subroutine load_coils_into_field(field, coil_path, extcur_path)
        type(biotsavart_field_t), intent(inout) :: field
        character(len=*), intent(in) :: coil_path
        character(len=*), intent(in), optional :: extcur_path

        logical :: is_stellopt
        character(len=:), allocatable :: trimmed_extcur
        real(dp), allocatable :: x(:), y(:), z(:), current(:)

        call determine_format(trim(coil_path), is_stellopt)

        if (.not. is_stellopt) then
            call field%biotsavart_field_init(trim(coil_path))
            return
        end if

        if (present(extcur_path)) then
            if (len_trim(extcur_path) > 0) trimmed_extcur = trim(extcur_path)
        end if

        if (allocated(trimmed_extcur)) then
            call read_stellopt_coils(trim(coil_path), x, y, z, current, trimmed_extcur)
        else
            call read_stellopt_coils(trim(coil_path), x, y, z, current)
        end if

        if (allocated(field%coils%x)) then
            call coils_deinit(field%coils)
        end if
        call coils_init(x, y, z, current, field%coils)
        deallocate(x, y, z, current)
    end subroutine load_coils_into_field

    subroutine determine_format(path, is_stellopt)
        character(len=*), intent(in) :: path
        logical, intent(out) :: is_stellopt

        character(len=512) :: line
        integer :: unit, ios
        character(len=:), allocatable :: trimmed
        integer(int32) :: dummy_int

        is_stellopt = .false.
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A)') 'failed to open coil file: '//trim(path)
            stop 1
        end if
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            trimmed = adjustl(line)
            if (len_trim(trimmed) == 0) cycle
            if (trimmed(1:1) == '!') cycle
            call to_lower_inplace(trimmed)
            if (index(trimmed, 'periods') == 1) then
                is_stellopt = .true.
            else
                read(trimmed, *, iostat=ios) dummy_int
                if (ios /= 0) then
                    write(error_unit, '(A)') 'unrecognised coil header: '//trim(trimmed)
                    stop 1
                end if
                is_stellopt = .false.
            end if
            exit
        end do
        close(unit)
    end subroutine determine_format

    subroutine read_stellopt_coils(path, x, y, z, current, extcur_path)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable, intent(out) :: y(:)
        real(dp), allocatable, intent(out) :: z(:)
        real(dp), allocatable, intent(out) :: current(:)
        character(len=:), allocatable, intent(in), optional :: extcur_path

        integer :: unit, ios, n_points, capacity, coil_start
        real(dp), allocatable :: tmp_x(:), tmp_y(:), tmp_z(:), tmp_current(:)
        integer, allocatable :: group_ids(:)
        character(len=512) :: line
        character(len=:), allocatable :: trimmed
        real(dp) :: px, py, pz, icurrent
        integer :: group_id, max_group
        logical :: has_group
        real(dp), allocatable :: ref_current(:)
        real(dp), allocatable :: target_extcur(:)
        integer :: i

        capacity = 0
        n_points = 0
        coil_start = 0
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A)') 'failed to open STELLOPT coil: '//trim(path)
            stop 1
        end if

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            trimmed = adjustl(line)
            if (len_trim(trimmed) == 0) cycle
            call to_lower_inplace(trimmed)
            if (index(trimmed, 'periods') == 1) cycle
            if (index(trimmed, 'begin') == 1) cycle
            if (index(trimmed, 'mirror') == 1) cycle
            if (index(trimmed, 'end') == 1) exit
            call parse_coil_line(line, px, py, pz, icurrent, has_group, group_id)
            call append_point(px, py, pz, icurrent, tmp_x, tmp_y, tmp_z, tmp_current, &
                group_ids, n_points, capacity)
            if (coil_start == 0) coil_start = n_points
            if (has_group) then
                call assign_group_range(group_ids, coil_start, n_points, group_id)
                coil_start = 0
            end if
        end do
        close(unit)

        if (n_points == 0) then
            write(error_unit, '(A)') 'STELLOPT coil file had no points: '//trim(path)
            stop 1
        end if

        if (coil_start /= 0) then
            write(error_unit, '(A)') 'unterminated coil in: '//trim(path)
            stop 1
        end if

        call trim_storage(n_points, tmp_x, tmp_y, tmp_z, tmp_current, group_ids)
        max_group = maxval(group_ids)
        if (max_group <= 0) then
            write(error_unit, '(A)') 'could not detect coil groups in: '//trim(path)
            stop 1
        end if

        allocate(ref_current(max_group))
        ref_current = 0.0_dp
        do i = 1, n_points
            if (group_ids(i) < 1) cycle
            if (ref_current(group_ids(i)) == 0.0_dp) then
                ref_current(group_ids(i)) = abs(tmp_current(i))
            end if
        end do

        if (present(extcur_path)) then
            if (allocated(extcur_path)) then
                call load_extcur_values(extcur_path, max_group, target_extcur)
            else
                allocate(target_extcur(max_group))
                target_extcur = 0.0_dp
            end if
        else
            allocate(target_extcur(max_group))
            target_extcur = 0.0_dp
        end if

        do i = 1, max_group
            if (target_extcur(i) == 0.0_dp) then
                target_extcur(i) = ref_current(i)
            end if
        end do

        do i = 1, n_points
            if (group_ids(i) < 1) cycle
            if (ref_current(group_ids(i)) /= 0.0_dp) then
                tmp_current(i) = tmp_current(i) * target_extcur(group_ids(i)) / &
                    ref_current(group_ids(i))
            end if
        end do

        call move_alloc(tmp_x, x)
        call move_alloc(tmp_y, y)
        call move_alloc(tmp_z, z)
        call move_alloc(tmp_current, current)
        deallocate(group_ids, ref_current, target_extcur)
    end subroutine read_stellopt_coils

    subroutine parse_coil_line(line, x, y, z, current, has_group, group_id)
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: x, y, z, current
        logical, intent(out) :: has_group
        integer, intent(out) :: group_id

        real(dp) :: tmp_x, tmp_y, tmp_z, tmp_current
        integer :: ios, tmp_group

        read(line, *, iostat=ios) x, y, z, current
        if (ios /= 0) then
            write(error_unit, '(A)') 'invalid coil line: '//trim(line)
            stop 1
        end if

        read(line, *, iostat=ios) tmp_x, tmp_y, tmp_z, tmp_current, tmp_group
        if (ios == 0) then
            has_group = .true.
            group_id = tmp_group
        else
            has_group = .false.
            group_id = -1
        end if
    end subroutine parse_coil_line

    subroutine append_point(px, py, pz, pcurr, x, y, z, current, groups, n_points, capacity)
        real(dp), intent(in) :: px, py, pz, pcurr
        real(dp), allocatable, intent(inout) :: x(:), y(:), z(:), current(:)
        integer, allocatable, intent(inout) :: groups(:)
        integer, intent(inout) :: n_points
        integer, intent(inout) :: capacity

        call ensure_capacity(n_points + 1, x, y, z, current, groups, capacity)
        n_points = n_points + 1
        x(n_points) = px
        y(n_points) = py
        z(n_points) = pz
        current(n_points) = pcurr
        groups(n_points) = 0
    end subroutine append_point

    subroutine ensure_capacity(required, x, y, z, current, groups, capacity)
        integer, intent(in) :: required
        real(dp), allocatable, intent(inout) :: x(:), y(:), z(:), current(:)
        integer, allocatable, intent(inout) :: groups(:)
        integer, intent(inout) :: capacity

        integer :: new_cap

        if (required <= capacity) return
        new_cap = max(required, merge(1024, capacity * 2, capacity > 0))
        call grow_array(x, new_cap)
        call grow_array(y, new_cap)
        call grow_array(z, new_cap)
        call grow_array(current, new_cap)
        call grow_int_array(groups, new_cap)
        capacity = new_cap
    end subroutine ensure_capacity

    subroutine grow_array(array, new_cap)
        real(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: new_cap
        real(dp), allocatable :: tmp(:)
        integer :: old_size

        old_size = 0
        if (allocated(array)) then
            old_size = size(array)
            call move_alloc(array, tmp)
        end if
        allocate(array(new_cap))
        if (old_size > 0) array(1:old_size) = tmp(1:old_size)
        if (new_cap > old_size) array(old_size+1:new_cap) = 0.0_dp
        if (allocated(tmp)) deallocate(tmp)
    end subroutine grow_array

    subroutine grow_int_array(array, new_cap)
        integer, allocatable, intent(inout) :: array(:)
        integer, intent(in) :: new_cap
        integer, allocatable :: tmp(:)
        integer :: old_size

        old_size = 0
        if (allocated(array)) then
            old_size = size(array)
            call move_alloc(array, tmp)
        end if
        allocate(array(new_cap))
        if (old_size > 0) array(1:old_size) = tmp(1:old_size)
        if (new_cap > old_size) array(old_size+1:new_cap) = 0
        if (allocated(tmp)) deallocate(tmp)
    end subroutine grow_int_array

    subroutine assign_group_range(groups, first_idx, last_idx, group_id)
        integer, intent(inout) :: groups(:)
        integer, intent(in) :: first_idx, last_idx, group_id

        if (first_idx < 1 .or. last_idx > size(groups)) then
            write(error_unit, '(A)') 'invalid coil indexing while assigning groups'
            stop 1
        end if
        groups(first_idx:last_idx) = group_id
    end subroutine assign_group_range

    subroutine trim_storage(n_points, x, y, z, current, groups)
        integer, intent(in) :: n_points
        real(dp), allocatable, intent(inout) :: x(:), y(:), z(:), current(:)
        integer, allocatable, intent(inout) :: groups(:)

        if (size(x) > n_points) call shrink_real_array(x, n_points)
        if (size(y) > n_points) call shrink_real_array(y, n_points)
        if (size(z) > n_points) call shrink_real_array(z, n_points)
        if (size(current) > n_points) call shrink_real_array(current, n_points)
        if (size(groups) > n_points) call shrink_int_array(groups, n_points)
    end subroutine trim_storage

    subroutine shrink_real_array(array, new_size)
        real(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: new_size
        real(dp), allocatable :: tmp(:)

        call move_alloc(array, tmp)
        allocate(array(new_size))
        array = tmp(1:new_size)
        if (allocated(tmp)) deallocate(tmp)
    end subroutine shrink_real_array

    subroutine shrink_int_array(array, new_size)
        integer, allocatable, intent(inout) :: array(:)
        integer, intent(in) :: new_size
        integer, allocatable :: tmp(:)

        call move_alloc(array, tmp)
        allocate(array(new_size))
        array = tmp(1:new_size)
        if (allocated(tmp)) deallocate(tmp)
    end subroutine shrink_int_array

    subroutine load_extcur_values(path, max_group, values)
        character(len=*), intent(in) :: path
        integer, intent(in) :: max_group
        real(dp), allocatable, intent(out) :: values(:)

        character(len=512) :: line
        character(len=:), allocatable :: lowered
        integer :: unit, ios, idx, start_pos, end_pos
        logical :: found_keyword
        real(dp) :: value

        allocate(values(max_group))
        values = 0.0_dp
        if (len_trim(path) == 0) return
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A)') 'failed to open EXTCUR file: '//trim(path)
            stop 1
        end if
        found_keyword = .false.
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            lowered = adjustl(line)
            call to_lower_inplace(lowered)
            if (index(lowered, 'extcur') > 0) then
                found_keyword = .true.
                start_pos = index(lowered, '(')
                end_pos = index(lowered, ')')
                if (start_pos > 0 .and. end_pos > start_pos) then
                    read(lowered(start_pos+1:end_pos-1), *, iostat=ios) idx
                else
                    idx = -1
                end if
                start_pos = index(line, '=')
                if (start_pos > 0) then
                    read(line(start_pos+1:), *, iostat=ios) value
                else
                    ios = -1
                end if
                if (ios == 0 .and. idx >= 1 .and. idx <= max_group) then
                    values(idx) = value
                end if
            end if
        end do
        rewind(unit)
        if (.not. found_keyword) then
            do idx = 1, max_group
                read(unit, *, iostat=ios) value
                if (ios /= 0) exit
                values(idx) = value
            end do
        end if
        close(unit)
    end subroutine load_extcur_values

    subroutine to_lower_inplace(text)
        character(len=:), allocatable, intent(inout) :: text
        integer :: i, ia
        do i = 1, len(text)
            ia = iachar(text(i:i))
            if (ia >= iachar('A') .and. ia <= iachar('Z')) then
                text(i:i) = achar(ia + 32)
            end if
        end do
    end subroutine to_lower_inplace

end module tiago_coil_loader
