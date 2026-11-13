program tiago_field_line_labels
    !! Field-line profile labels of start points in one toroidal plane, for a
    !! scan of one coil-group current (docs/field-line-labels.md).
    !!
    !! Usage: tiago_field_line_labels input.nml
    !!
    !! &field_line_labels
    !!   mgrid = 'mgrid.nc'          ! makegrid file
    !!   currents = 1.0, 1.0, ...    ! one current per coil group
    !!   scan_group = 0              ! group whose current is replaced by each
    !!   scan_currents = 0.0, 1000.0 !   scan value; 0: a single case
    !!   wout = 'wout.nc'            ! reference chart (VMEC), s = rho^2
    !!   radii = 'radii.txt'         ! start radii [m], one per line, at Z = 0
    !!   phi_deg = 0.0               ! plane of the starts and returns
    !!   turns = 200, steps = 720    ! toroidal turns, RK4 steps per turn
    !!   output = 'labels.json'
    !! /
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_field_line_label, only: grid_field_t, line_labels_t, compute_labels => field_line_labels
    use tiago_vmec_state, only: vmec_state_t
    use tiago_sha256, only: sha256_file
    use tiago_revision, only: tiago_source_revision
    implicit none

    integer, parameter :: max_groups = 128
    real(dp), parameter :: unset = -huge(1.0_dp), pi = acos(-1.0_dp)
    character(len=1024) :: mgrid = '', wout = '', radii = '', output = 'labels.json', nml
    real(dp) :: currents(max_groups) = unset, scan_currents(max_groups) = unset
    real(dp) :: phi_deg = 0.0_dp
    integer :: scan_group = 0, turns = 200, steps = 720
    namelist /field_line_labels/ mgrid, currents, scan_group, scan_currents, wout, radii, &
        phi_deg, turns, steps, output

    type(vmec_state_t) :: chart
    type(grid_field_t), allocatable :: field
    type(line_labels_t), allocatable :: cases(:)
    real(dp), allocatable :: r(:), z(:), group(:), scan(:)
    character(len=:), allocatable :: message
    integer :: unit, ios, ierr, ngroup, nscan, k, n
    integer(8) :: t0, t1, rate

    call system_clock(t0, rate)
    if (command_argument_count() /= 1) then
        write(error_unit, '(A)') 'usage: tiago_field_line_labels input.nml'
        error stop 2
    end if
    call get_command_argument(1, nml)
    open(newunit=unit, file=trim(nml), status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('cannot open '//trim(nml))
    read(unit, nml=field_line_labels, iostat=ios)
    close(unit)
    if (ios /= 0) call fail('cannot read namelist field_line_labels from '//trim(nml))
    ngroup = count(currents > unset)
    group = currents(1:ngroup)
    nscan = count(scan_currents > unset)
    if (scan_group == 0) then
        scan = [0.0_dp]
    else
        if (scan_group < 1 .or. scan_group > ngroup .or. nscan == 0) &
            call fail('scan_group must name a group with scan_currents')
        scan = scan_currents(1:nscan)
    end if
    if (turns < 1 .or. steps < 4) call fail('turns >= 1 and steps >= 4 required')

    call chart%load(trim(wout), ierr, message)
    if (ierr /= 0) call fail(message)
    call read_radii(trim(radii), r)
    allocate(z(size(r)))
    z = 0.0_dp

    allocate(cases(size(scan)))
    do k = 1, size(scan)
        if (scan_group > 0) group(scan_group) = scan(k)
        allocate(field)
        call field%load_mgrid(trim(mgrid), group, ierr, message)
        if (ierr /= 0) call fail(message)
        call compute_labels(field, chart, r, z, phi_deg * pi / 180.0_dp, turns, steps, &
            cases(k))
        deallocate(field)
        call system_clock(t1)
        write(output_unit, '(A,ES12.4,A,F0.1,A)') 'case ', scan(k), ' A done after ', &
            real(t1 - t0, dp) / real(rate, dp), ' s'
        flush(output_unit)
    end do
    call system_clock(t1)
    call write_json(trim(output), real(t1 - t0, dp) / real(rate, dp))

contains

    subroutine fail(text)
        character(len=*), intent(in) :: text
        write(error_unit, '(A)') 'tiago_field_line_labels: '//text
        error stop 1
    end subroutine fail

    subroutine read_radii(path, values)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:)
        real(dp) :: v
        integer :: u, io, m

        open(newunit=u, file=path, status='old', action='read', iostat=io)
        if (io /= 0) call fail('cannot open radii file '//path)
        m = 0
        do
            read(u, *, iostat=io) v
            if (io /= 0) exit
            m = m + 1
        end do
        rewind(u)
        allocate(values(m))
        do m = 1, size(values)
            read(u, *) values(m)
        end do
        close(u)
        if (size(values) == 0) call fail('no radii in '//path)
    end subroutine read_radii

    function num(x) result(text)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: text
        character(len=32) :: buf
        if (ieee_is_finite(x)) then
            write(buf, '(ES24.16E3)') x
            text = trim(adjustl(buf))
        else
            text = 'null'
        end if
    end function num

    function vector(x) result(text)
        real(dp), intent(in) :: x(:)
        character(len=:), allocatable :: text
        integer :: i
        text = '['
        do i = 1, size(x)
            text = text//num(x(i))
            if (i < size(x)) text = text//', '
        end do
        text = text//']'
    end function vector

    function flags(x) result(text)
        logical, intent(in) :: x(:)
        character(len=:), allocatable :: text
        integer :: i
        text = '['
        do i = 1, size(x)
            if (x(i)) then
                text = text//'true'
            else
                text = text//'false'
            end if
            if (i < size(x)) text = text//', '
        end do
        text = text//']'
    end function flags

    function ints(x) result(text)
        integer, intent(in) :: x(:)
        character(len=:), allocatable :: text
        character(len=16) :: buf
        integer :: i
        text = '['
        do i = 1, size(x)
            write(buf, '(I0)') x(i)
            text = text//trim(buf)
            if (i < size(x)) text = text//', '
        end do
        text = text//']'
    end function ints

    subroutine write_json(path, seconds)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: seconds
        integer :: u, k
        character(len=16) :: buf

        open(newunit=u, file=path, status='replace', action='write')
        write(u, '(A)') '{'
        write(u, '(A)') '"tiago_revision": "'//trim(tiago_source_revision)//'",'
        write(u, '(A)') '"phi_section_rad": '//num(phi_deg * pi / 180.0_dp)//','
        write(buf, '(I0)') turns
        write(u, '(A)') '"turns": '//trim(buf)//','
        write(buf, '(I0)') steps
        write(u, '(A)') '"steps_per_turn": '//trim(buf)//','
        write(buf, '(I0)') scan_group
        write(u, '(A)') '"scan_group": '//trim(buf)//','
        write(u, '(A)') '"currents": '//vector(currents(1:ngroup))//','
        write(u, '(A)') '"R_m": '//vector(r)//','
        write(u, '(A)') '"Z_m": '//vector(z)//','
        write(u, '(A)') '"s_own": '//vector(cases(1)%s_start)//','
        write(u, '(A)') '"cases": ['
        do k = 1, size(scan)
            write(u, '(A)') '{"current_A": '//num(scan(k))//', "label": '// &
                vector(cases(k)%label)//', "label_spread": '//vector(cases(k)%spread)// &
                ', "returns": '//ints(cases(k)%returns)//', "lost": '// &
                flags(cases(k)%lost)//'}'
            if (k < size(scan)) write(u, '(A)') ','
        end do
        write(u, '(A)') '],'
        write(u, '(A)') '"inputs_sha256": {"'//trim(mgrid)//'": "'//sha256_file(trim(mgrid))// &
            '", "'//trim(wout)//'": "'//sha256_file(trim(wout))//'", "'//trim(radii)// &
            '": "'//sha256_file(trim(radii))//'"},'
        write(u, '(A)') '"wall_seconds": '//num(seconds)
        write(u, '(A)') '}'
        close(u)
    end subroutine write_json

end program tiago_field_line_labels
