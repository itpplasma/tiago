module tiago_coil_loader
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32, error_unit
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_biotsavart, only: coils_init, coils_deinit
    implicit none
    private

    real(dp), parameter :: tiny_current = 1.0e-12_dp

    public :: load_coils_into_field

contains

    subroutine load_coils_into_field(field, coil_path, extcur_path, groups, unit_current)
        !! groups/unit_current (optional): coil group of every point and its
        !! current per unit EXTCUR of that group, for response matrices. Files
        !! without groups (libneo format) form a single group with EXTCUR = 1.
        type(biotsavart_field_t), intent(inout) :: field
        character(len=*), intent(in) :: coil_path
        character(len=*), intent(in), optional :: extcur_path
        integer, allocatable, intent(out), optional :: groups(:)
        real(dp), allocatable, intent(out), optional :: unit_current(:)

        logical :: is_stellopt
        character(len=:), allocatable :: trimmed_extcur
        real(dp), allocatable :: x(:), y(:), z(:), current(:), unit(:)
        integer, allocatable :: group_ids(:)

        call determine_format(trim(coil_path), is_stellopt)

        if (.not. is_stellopt) then
            call field%biotsavart_field_init(trim(coil_path))
            if (present(groups)) then
                allocate(groups(size(field%coils%current)))
                groups = 1
            end if
            if (present(unit_current)) unit_current = field%coils%current
            return
        end if

        if (present(extcur_path)) then
            if (len_trim(extcur_path) > 0) trimmed_extcur = trim(extcur_path)
        end if

        if (allocated(trimmed_extcur)) then
            call read_stellopt_coils(trim(coil_path), x, y, z, current, group_ids, unit, &
                trimmed_extcur)
        else
            call read_stellopt_coils(trim(coil_path), x, y, z, current, group_ids, unit)
        end if

        if (allocated(field%coils%x)) then
            call coils_deinit(field%coils)
        end if
        call coils_init(x, y, z, current, field%coils)
        if (present(groups)) call move_alloc(group_ids, groups)
        if (present(unit_current)) call move_alloc(unit, unit_current)
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

    subroutine read_stellopt_coils(path, x, y, z, current, groups, unit_current, extcur_path)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable, intent(out) :: y(:)
        real(dp), allocatable, intent(out) :: z(:)
        real(dp), allocatable, intent(out) :: current(:)
        integer, allocatable, intent(out) :: groups(:)
        real(dp), allocatable, intent(out) :: unit_current(:)
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
        logical, allocatable :: given(:)
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
                ref_current(group_ids(i)) = tmp_current(i)
            end if
        end do

        allocate(target_extcur(max_group), given(max_group))
        target_extcur = 0.0_dp
        given = .false.
        if (present(extcur_path)) then
            if (allocated(extcur_path)) then
                call load_extcur_values(extcur_path, target_extcur, given)
            end if
        end if

        ! Current per unit EXTCUR of the point's group (response matrices).
        allocate(unit_current(n_points))
        unit_current = 0.0_dp
        do i = 1, n_points
            if (group_ids(i) < 1) cycle
            if (ref_current(group_ids(i)) /= 0.0_dp) then
                unit_current(i) = tmp_current(i) / ref_current(group_ids(i))
            end if
        end do

        ! EXTCUR(g) replaces the file current of group g, keeping relative
        ! currents within the group; groups without EXTCUR keep the file currents.
        do i = 1, n_points
            if (group_ids(i) < 1) cycle
            if (.not. given(group_ids(i))) cycle
            if (ref_current(group_ids(i)) /= 0.0_dp) then
                tmp_current(i) = tmp_current(i) * target_extcur(group_ids(i)) / &
                    ref_current(group_ids(i))
            end if
        end do

        call move_alloc(tmp_x, x)
        call move_alloc(tmp_y, y)
        call move_alloc(tmp_z, z)
        call move_alloc(tmp_current, current)
        call move_alloc(group_ids, groups)
        deallocate(ref_current, target_extcur, given)
    end subroutine read_stellopt_coils

    subroutine parse_coil_line(line, x, y, z, current, has_group, group_id)
        !! "x y z I [group name]": only the closing line of a coil has a group.
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: x, y, z, current
        logical, intent(out) :: has_group
        integer, intent(out) :: group_id

        real(dp) :: values(4)
        integer :: ios

        read(line, *, iostat=ios) x, y, z, current
        if (ios /= 0) then
            write(error_unit, '(A)') 'invalid coil line: '//trim(line)
            stop 1
        end if

        has_group = .false.
        group_id = -1
        if (count_tokens(line) < 5) return    ! cheap test avoids a second read per line
        read(line, *, iostat=ios) values, group_id
        has_group = ios == 0
        if (.not. has_group) group_id = -1
    end subroutine parse_coil_line

    pure integer function count_tokens(line)
        character(len=*), intent(in) :: line
        integer :: i
        logical :: in_token

        count_tokens = 0
        in_token = .false.
        do i = 1, len_trim(line)
            if (line(i:i) == ' ' .or. line(i:i) == achar(9) .or. line(i:i) == ',') then
                in_token = .false.
            else if (.not. in_token) then
                in_token = .true.
                count_tokens = count_tokens + 1
            end if
        end do
    end function count_tokens

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
        new_cap = max(required, merge(capacity * 2, 1024, capacity > 0))
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

    subroutine load_extcur_values(path, values, given)
        !! EXTCUR from a VMEC &INDATA file (EXTCUR(i) = v, EXTCUR = a, b, ...,
        !! EXTCUR(i) = a b, repeat counts n*v, D exponents, ! comments) or,
        !! without any EXTCUR keyword, a plain list of numbers.
        character(len=*), intent(in) :: path
        real(dp), intent(inout) :: values(:)
        logical, intent(inout) :: given(:)

        character(len=:), allocatable :: text
        integer :: pos, start_idx, found

        text = read_lowercase_without_comments(path)
        found = 0
        pos = 1
        do
            pos = find_keyword(text, 'extcur', pos)
            if (pos == 0) exit
            found = found + 1
            pos = pos + len('extcur')
            start_idx = 1
            call skip_blanks(text, pos)
            if (pos <= len(text)) then
                if (text(pos:pos) == '(') then
                    call read_index(text, pos, start_idx, path)
                    call skip_blanks(text, pos)
                end if
            end if
            if (pos > len(text)) call extcur_error(path, 'missing "=" after EXTCUR')
            if (text(pos:pos) /= '=') call extcur_error(path, 'missing "=" after EXTCUR')
            pos = pos + 1
            call read_value_list(text, pos, start_idx, values, given, path)
        end do

        if (found == 0) call read_plain_list(text, values, given, path)
    end subroutine load_extcur_values

    function read_lowercase_without_comments(path) result(text)
        character(len=*), intent(in) :: path
        character(len=:), allocatable :: text
        character(len=1024) :: line
        character(len=:), allocatable :: lowered
        integer :: unit, ios, bang

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A)') 'failed to open EXTCUR file: '//trim(path)
            stop 1
        end if
        text = ''
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            lowered = trim(line)
            bang = index(lowered, '!')
            if (bang > 0) lowered = lowered(:bang - 1)
            call to_lower_inplace(lowered)
            text = text // ' ' // lowered
        end do
        close(unit)
        text = text // ' '
    end function read_lowercase_without_comments

    integer function find_keyword(text, key, from) result(pos)
        !! Next occurrence of key as a whole identifier (so LEXTCUR does not match).
        character(len=*), intent(in) :: text, key
        integer, intent(in) :: from
        integer :: k, after

        pos = from
        do
            k = index(text(pos:), key)
            if (k == 0) then
                pos = 0
                return
            end if
            pos = pos + k - 1
            after = pos + len(key)
            if (.not. is_ident_char(text, pos - 1) .and. .not. is_ident_char(text, after)) return
            pos = pos + 1
        end do
    end function find_keyword

    logical function is_ident_char(text, i)
        character(len=*), intent(in) :: text
        integer, intent(in) :: i

        is_ident_char = .false.
        if (i < 1 .or. i > len(text)) return
        is_ident_char = verify(text(i:i), 'abcdefghijklmnopqrstuvwxyz0123456789_') == 0
    end function is_ident_char

    subroutine skip_blanks(text, pos)
        character(len=*), intent(in) :: text
        integer, intent(inout) :: pos

        do while (pos <= len(text))
            if (text(pos:pos) /= ' ' .and. text(pos:pos) /= achar(9)) exit
            pos = pos + 1
        end do
    end subroutine skip_blanks

    subroutine read_index(text, pos, idx, path)
        character(len=*), intent(in) :: text, path
        integer, intent(inout) :: pos
        integer, intent(out) :: idx
        integer :: close_pos, ios

        close_pos = index(text(pos:), ')')
        if (close_pos == 0) call extcur_error(path, 'unterminated EXTCUR(')
        read(text(pos + 1:pos + close_pos - 2), *, iostat=ios) idx
        if (ios /= 0 .or. idx < 1) call extcur_error(path, 'invalid EXTCUR index')
        pos = pos + close_pos
    end subroutine read_index

    subroutine read_value_list(text, pos, start_idx, values, given, path)
        !! Values after "=" up to the next identifier, "/" or "&".
        character(len=*), intent(in) :: text, path
        integer, intent(inout) :: pos
        integer, intent(in) :: start_idx
        real(dp), intent(inout) :: values(:)
        logical, intent(inout) :: given(:)
        integer :: idx, tok_end, star, repeat, ios, k
        real(dp) :: v
        character(len=:), allocatable :: token

        idx = start_idx
        do
            do while (pos <= len(text))
                if (index(' ,'//achar(9), text(pos:pos)) == 0) exit
                pos = pos + 1
            end do
            if (pos > len(text)) exit
            if (index('/&', text(pos:pos)) > 0) exit
            if (verify(text(pos:pos), 'abcdefghijklmnopqrstuvwxyz_') == 0) exit
            tok_end = scan(text(pos:), ' ,/&'//achar(9)) + pos - 2
            token = text(pos:tok_end)
            pos = tok_end + 1
            repeat = 1
            star = index(token, '*')
            if (star > 0) then
                read(token(:star - 1), *, iostat=ios) repeat
                if (ios /= 0) call extcur_error(path, 'invalid repeat count: '//token)
                token = token(star + 1:)
            end if
            read(token, *, iostat=ios) v
            if (ios /= 0) call extcur_error(path, 'invalid EXTCUR value: '//token)
            do k = 1, repeat
                call store(idx, v, values, given, path)
                idx = idx + 1
            end do
        end do
        if (idx == start_idx) call extcur_error(path, 'EXTCUR without values')
    end subroutine read_value_list

    subroutine read_plain_list(text, values, given, path)
        character(len=*), intent(in) :: text, path
        real(dp), intent(inout) :: values(:)
        logical, intent(inout) :: given(:)
        integer :: pos

        pos = 1
        call read_value_list(text, pos, 1, values, given, path)
    end subroutine read_plain_list

    subroutine store(idx, v, values, given, path)
        integer, intent(in) :: idx
        real(dp), intent(in) :: v
        real(dp), intent(inout) :: values(:)
        logical, intent(inout) :: given(:)
        character(len=*), intent(in) :: path

        if (idx > size(values)) then
            write(error_unit, '(A,I0,A,I0,A)') 'WARNING: EXTCUR(', idx, ') ignored; coil file has ', &
                size(values), ' groups ('//trim(path)//')'
            return
        end if
        values(idx) = v
        given(idx) = .true.
    end subroutine store

    subroutine extcur_error(path, message)
        character(len=*), intent(in) :: path, message
        write(error_unit, '(A)') trim(path)//': '//message
        stop 1
    end subroutine extcur_error

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
