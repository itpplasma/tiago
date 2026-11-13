program tiago_island
    !! Magnetic axis, island O and X points, width and phase, and an optional
    !! Poincare section, for a field of coils, a SPECTRE state and/or a nested
    !! VMEC state (inside its plasma).
    !! Usage: tiago_island input.nml   (namelist &island, see docs/islands.md)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_field_lines, only: cartesian_line_field_t, read_scaled_coils
    use tiago_vmec_state, only: vmec_state_t
    use tiago_island_measure, only: trace_options_t, island_point_t, &
        island_width_t, periodic_point, island_width, poincare_section
#ifdef TIAGO_ENABLE_SPECTRE
    use tiago_spectre_source, only: spectre_source_t
#endif
    implicit none
    integer, parameter :: max_files = 16, path_len = 1024
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(len=path_len) :: coil_files(max_files) = '', spectre_file = ''
    character(len=path_len) :: coil_extcur(max_files) = '', vmec_wout = ''
    character(len=path_len) :: output = 'island', input
    real(dp) :: coil_scales(max_files) = 1.0_dp
    real(dp) :: phi0 = 0.0_dp, axis_seed(2) = 0, o_seed(2) = 0, x_seed(2) = 0
    real(dp) :: half_span = 0.05_dp, width_turns = 32.0_dp, r_min = 0.0_dp
    real(dp) :: point_rtol = 1.0e-11_dp, point_atol = 1.0e-13_dp
    real(dp) :: point_max_step = 2*pi/512, map_delta = 1.0e-4_dp
    real(dp) :: root_tol = 2.0e-6_dp
    real(dp) :: width_rtol = 1.0e-10_dp, width_atol = 1.0e-12_dp
    real(dp) :: width_max_step = 2*pi/256
    integer :: map_turns = 2, m = 2, n = 1, scans = 14, width_fixed_steps = 0
    integer :: poincare_lines = 0, poincare_crossings = 0, max_steps = 2000000
    logical :: measure_x = .true., measure_width = .true., axis_given = .false.
    namelist /island/ coil_files, coil_scales, coil_extcur, spectre_file, vmec_wout, output, phi0, &
        axis_seed, o_seed, x_seed, half_span, width_turns, r_min, point_rtol, &
        point_atol, point_max_step, map_delta, root_tol, width_rtol, width_atol, &
        width_max_step, map_turns, m, n, scans, width_fixed_steps, &
        poincare_lines, poincare_crossings, max_steps, measure_x, measure_width, &
        axis_given
    type(cartesian_line_field_t) :: field
    type(trace_options_t) :: popt, wopt
    type(island_point_t) :: ax, op, xp
    type(island_width_t) :: wd
    real(dp) :: dist, dir(2), tm(5), type_margin
    integer :: unit, ios, nfiles, f, poincare_status = 0

    if (command_argument_count() /= 1) call die('usage: tiago_island input.nml')
    call get_command_argument(1, input)
    open (newunit=unit, file=trim(input), status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open '//trim(input))
    read (unit, nml=island, iostat=ios)
    if (ios /= 0) call die('cannot read namelist &island from '//trim(input))
    close (unit)

    nfiles = count(len_trim(coil_files) > 0)
    if (nfiles > 0) then
        call read_scaled_coils(coil_files(:nfiles), coil_scales(:nfiles), field%coils, &
                               coil_extcur(:nfiles))
        field%has_coils = .true.
    end if
    if (len_trim(spectre_file) > 0) call attach_spectre()
    if (len_trim(vmec_wout) > 0) call attach_vmec()
    if (nfiles == 0 .and. len_trim(spectre_file) == 0 .and. len_trim(vmec_wout) == 0) &
        call die('no field given')
    field%r_min = r_min

    popt%rtol = point_rtol
    popt%atol = point_atol
    popt%max_step = point_max_step
    popt%map_delta = map_delta
    popt%root_tol = root_tol
    popt%max_steps = max_steps
    wopt%rtol = width_rtol
    wopt%atol = width_atol
    wopt%max_step = width_max_step
    wopt%fixed_steps = width_fixed_steps
    wopt%max_steps = max_steps

    open (newunit=f, file=trim(output)//'.json', status='replace')
    write (f, '(a)') '{'
    write (f, '(a,i0,a,i0,a,i0,a,es24.15e3,a)') '"m": ', m, ', "n": ', n, &
        ', "map_turns": ', map_turns, ', "phi0": ', phi0, ','
    tm(1) = wall()
    if (axis_given) then
        ! A nested state's chart is singular on its axis, so its axis is given.
        ax = island_point_t(x=axis_seed, closure=0, position_error=0, status=0)
    else
        call periodic_point(field, axis_seed, phi0, 1, axis_seed, 1, 0, popt, ax)
    end if
    call put_point('axis', ax)
    tm(2) = wall()
    if (ax%status == 0) then
        call periodic_point(field, o_seed, phi0, map_turns, ax%x, m, n, popt, op)
        if (op%status == 0) then
            dist = norm2(op%x - ax%x)
            if (dist <= tiny(1.0_dp)) then
                op%status = 4
            else
                dir = (op%x - ax%x)/dist
            end if
            if (measure_width .or. measure_x) then
                type_margin = residue_margin(op)
                if (op%residue <= type_margin .or. &
                    op%residue >= 1 - type_margin) op%status = 4
            end if
        end if
        if (measure_width .or. measure_x) then
            call put_point('O', op)
        else
            call put_point('periodic_point', op)
        end if
        if (measure_x) then
            call periodic_point(field, x_seed, phi0, map_turns, ax%x, m, n, popt, xp)
            if (xp%status == 0) then
                if (xp%residue >= -residue_margin(xp)) xp%status = 4
            end if
            call put_point('X', xp)
        end if
        tm(3) = wall()
        if (op%status == 0 .and. measure_width) then
            call island_width(field, [dist - half_span, dist + half_span], dist, &
                              ax%x, dir, phi0, width_turns, real(map_turns, dp), m, wopt, wd, &
                              scans=scans)
            write (f, '(a,i0,6(a,es24.15e3),a)') '"width": {"status": ', wd%status, &
                ', "metres": ', wd%metres, ', "trace_error": ', wd%trace_error, &
                ', "duration_error": ', wd%duration_error, ', "scan_error": ', &
                wd%scan_error, ', "edge": [', wd%edge(1), ',', wd%edge(2), ']},'
        end if
        tm(4) = wall()
        if (op%status == 0 .and. poincare_lines > 0 .and. poincare_crossings > 0) &
            call write_poincare()
        tm(5) = wall()
        write (f, '(4(a,f12.3),a)') '"seconds_axis": ', tm(2) - tm(1), &
            ', "seconds_points": ', tm(3) - tm(2), ', "seconds_width": ', &
            tm(4) - tm(3), ', "seconds_poincare": ', tm(5) - tm(4), ','
    end if
    write (f, '(a)') '"end": true}'
    close (f)
    if (ax%status /= 0) call die('magnetic axis not found')
    if (op%status /= 0) call die('periodic point not found or wrong type')
    if (measure_x) then
        if (xp%status /= 0) call die('X point not found or wrong type')
    end if
    if (measure_width) then
        if (wd%status /= 0) call die('island width not resolved')
    end if
    if (poincare_status /= 0) call die('Poincare trace left field domain')
contains
    subroutine attach_spectre()
#ifdef TIAGO_ENABLE_SPECTRE
        type(spectre_source_t), allocatable :: source
        allocate (source)
        call source%init(trim(spectre_file), 8, 16)
        call move_alloc(source, field%plasma)
#else
        call die('spectre_file requires a build with TIAGO_ENABLE_SPECTRE=ON')
#endif
    end subroutine attach_spectre

    subroutine attach_vmec()
        type(vmec_state_t), allocatable :: state
        character(len=:), allocatable :: message
        integer :: ierr
        allocate (state)
        call state%load(trim(vmec_wout), ierr, message, vector=.true.)
        if (ierr /= 0) call die(message)
        call move_alloc(state, field%state)
    end subroutine attach_vmec

    subroutine write_poincare()
        !! Lines started on the axis-O chord, recorded every map_turns turns.
        real(dp), allocatable :: starts(:, :), pts(:, :, :)
        integer, allocatable :: st(:)
        integer :: j, k, u
        allocate (starts(2, poincare_lines), st(poincare_lines))
        allocate (pts(2, poincare_crossings, poincare_lines))
        do j = 1, poincare_lines
            starts(:, j) = ax%x + dir*(dist - half_span + 2*half_span*(j - 1) &
                                       /max(poincare_lines - 1, 1))
        end do
        call poincare_section(field, starts, phi0, real(map_turns, dp), &
                              poincare_crossings, pts, st, popt)
        open (newunit=u, file=trim(output)//'-poincare.csv', status='replace')
        write (u, '(a)') 'line,crossing,R_m,Z_m,status'
        do j = 1, poincare_lines
            do k = 1, poincare_crossings
                write (u, '(i0,a,i0,2(a,es23.15),a,i0)') j, ',', k, ',', pts(1, k, j), &
                    ',', pts(2, k, j), ',', st(j)
            end do
        end do
        close (u)
        poincare_status = maxval(st)
    end subroutine write_poincare

    subroutine put_point(name, p)
        character(*), intent(in) :: name
        type(island_point_t), intent(in) :: p
        write (f, '(3a,i0,a)') '"', name, '": {"status": ', p%status, &
            ', "R": '//trim(number(p%x(1)))//', "Z": '//trim(number(p%x(2)))// &
            ', "closure": '//trim(number(p%closure))// &
            ', "residue": '//trim(number(p%residue))// &
            ', "residue_error": '//trim(number(p%residue_error))// &
            ', "phase": '//trim(number(p%phase))// &
            ', "phase_error": '//trim(number(p%phase_error))// &
            ', "position_error": '//trim(number(p%position_error))//'},'
    end subroutine put_point

    real(dp) function residue_margin(point)
        type(island_point_t), intent(in) :: point

        ! Central differences lose coordinate roundoff divided by map_delta.
        residue_margin = max(2*point%residue_error, &
                             16*epsilon(1.0_dp)*max(1.0_dp, norm2(point%x))/map_delta)
    end function residue_margin

    function number(value) result(text)
        real(dp), intent(in) :: value
        character(len=32) :: text

        text = 'null'
        if (.not. ieee_is_finite(value)) return
        if (abs(value) >= huge(value)/2) return
        write (text, '(es24.15e3)') value
    end function number

    real(dp) function wall()
        integer(8) :: c, r
        call system_clock(c, r)
        wall = real(c, dp)/r
    end function wall

    subroutine die(message)
        character(*), intent(in) :: message
        write (error_unit, '(a)') trim(message)
        error stop 1
    end subroutine die
end program tiago_island
