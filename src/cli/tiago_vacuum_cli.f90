program tiago_vacuum_cli
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, &
        error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    implicit none

    character(len=512) :: coil_path
    character(len=512) :: flux_path
    character(len=512) :: segrog_path
    character(len=512) :: output_dir
    character(len=512) :: flux_out_path
    character(len=512) :: segrog_out_path
    character(len=512) :: flux_turn_path
    character(len=512) :: segrog_turn_path
    integer(i32) :: samples_per_segment
    integer(i32) :: nfp_value
    real(dp) :: seg_area
    logical :: have_flux
    logical :: have_seg
    integer :: argc
    character(len=512) :: coil_extcur_path

    argc = command_argument_count()
    if (argc < 3) call usage_and_stop()

    call get_command_argument(1, coil_path)
    call get_command_argument(2, flux_path)
    call get_command_argument(3, segrog_path)

    output_dir = '.'
    flux_out_path = 'tiago_flux.csv'
    segrog_out_path = 'tiago_segrog.csv'
    flux_turn_path = ''
    segrog_turn_path = ''
    samples_per_segment = 6_i32
    nfp_value = 1_i32
    seg_area = -1.0_dp
    coil_extcur_path = ''
    have_flux = len_trim(flux_path) > 0
    have_seg = len_trim(segrog_path) > 0

    call parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        flux_turn_path, segrog_turn_path, samples_per_segment, seg_area, &
        nfp_value, coil_extcur_path)
    call ensure_paths(output_dir, flux_out_path, segrog_out_path)
    call run_solver(trim(coil_path), trim(flux_path), trim(segrog_path), &
        trim(output_dir), trim(flux_out_path), trim(segrog_out_path), &
        trim(flux_turn_path), trim(segrog_turn_path), samples_per_segment, &
        seg_area, have_flux, have_seg, nfp_value, trim(coil_extcur_path))
contains

subroutine parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        flux_turn_path, segrog_turn_path, samples_per_segment, seg_area, &
        nfp_value, coil_extcur_path)
    use, intrinsic :: iso_fortran_env, only: i32 => int32, dp => real64
    implicit none
    integer, intent(in) :: argc
    character(len=*), intent(inout) :: output_dir
    character(len=*), intent(inout) :: flux_out_path
    character(len=*), intent(inout) :: segrog_out_path
    character(len=*), intent(inout) :: flux_turn_path
    character(len=*), intent(inout) :: segrog_turn_path
    integer(i32), intent(inout) :: samples_per_segment
    real(dp), intent(inout) :: seg_area
    integer(i32), intent(inout) :: nfp_value
    character(len=*), intent(inout) :: coil_extcur_path

    integer :: i
    character(len=512) :: arg

    i = 4
    do while (i <= argc)
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ('--output-dir')
            i = i + 1
            call ensure_arg(argc, i, '--output-dir')
            call get_command_argument(i, output_dir)
        case ('--flux-out')
            i = i + 1
            call ensure_arg(argc, i, '--flux-out')
            call get_command_argument(i, flux_out_path)
        case ('--segrog-out')
            i = i + 1
            call ensure_arg(argc, i, '--segrog-out')
            call get_command_argument(i, segrog_out_path)
        case ('--flux-turns')
            i = i + 1
            call ensure_arg(argc, i, '--flux-turns')
            call get_command_argument(i, flux_turn_path)
        case ('--segrog-turns')
            i = i + 1
            call ensure_arg(argc, i, '--segrog-turns')
            call get_command_argument(i, segrog_turn_path)
        case ('--samples')
            i = i + 1
            call ensure_arg(argc, i, '--samples')
            call get_command_argument(i, arg)
            read(arg, *) samples_per_segment
        case ('--nfp')
            i = i + 1
            call ensure_arg(argc, i, '--nfp')
            call get_command_argument(i, arg)
            read(arg, *) nfp_value
        case ('--seg-area')
            i = i + 1
            call ensure_arg(argc, i, '--seg-area')
            call get_command_argument(i, arg)
            read(arg, *) seg_area
        case ('--coil-extcur')
            i = i + 1
            call ensure_arg(argc, i, '--coil-extcur')
            call get_command_argument(i, coil_extcur_path)
        case ('--help', '-h')
            call usage_and_stop()
        case default
            call die('unknown option: '//trim(arg))
        end select
        i = i + 1
    end do
end subroutine parse_options

subroutine ensure_arg(argc, position, flag)
    use, intrinsic :: iso_fortran_env, only: error_unit
    integer, intent(in) :: argc
    integer, intent(in) :: position
    character(len=*), intent(in) :: flag
    if (position > argc) then
        write(error_unit, '(A)') trim(flag)//' requires an argument'
        stop 1
    end if
end subroutine ensure_arg

subroutine ensure_paths(output_dir, flux_out_path, segrog_out_path)
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(inout) :: flux_out_path
    character(len=*), intent(inout) :: segrog_out_path

    if (len_trim(output_dir) > 0) then
        if (index(flux_out_path, '/') == 0) then
            flux_out_path = trim(output_dir)//'/'//trim(flux_out_path)
        end if
        if (index(segrog_out_path, '/') == 0) then
            segrog_out_path = trim(output_dir)//'/'//trim(segrog_out_path)
        end if
    end if
end subroutine ensure_paths

subroutine run_solver(coil_path, flux_path, segrog_path, output_dir, &
        flux_out_path, segrog_out_path, flux_turn_path, segrog_turn_path, &
        samples_per_segment, seg_area, have_flux, have_seg, nfp_value, &
        coil_extcur_path)
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, &
        error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    character(len=*), intent(in) :: coil_path
    character(len=*), intent(in) :: flux_path
    character(len=*), intent(in) :: segrog_path
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(in) :: flux_out_path
    character(len=*), intent(in) :: segrog_out_path
    character(len=*), intent(in) :: flux_turn_path
    character(len=*), intent(in) :: segrog_turn_path
    integer(i32), intent(in) :: samples_per_segment
    real(dp), intent(in) :: seg_area
    logical, intent(in) :: have_flux
    logical, intent(in) :: have_seg
    integer(i32), intent(in) :: nfp_value
    character(len=*), intent(in) :: coil_extcur_path

    type(vacuum_solver_t) :: solver
    type(quadrature_rule_t) :: rule
    type(flux_loop_t), allocatable :: loops(:)
    type(segmented_rogowski_t), allocatable :: segs(:)
    real(dp), allocatable :: fluxes(:)
    real(dp), allocatable :: voltages(:)
    integer(i32) :: ierr
    character(len=:), allocatable :: message

    call ensure_directory(output_dir)

    rule%samples_per_segment = samples_per_segment
    if (len_trim(coil_extcur_path) > 0) then
        call solver%init(coil_path, coil_extcur_path)
    else
        call solver%init(coil_path)
    end if
    call solver%set_nfp(nfp_value)

    if (have_flux) then
        call read_flux_loop_file(flux_path, loops, ierr, message)
        if (ierr /= 0_i32) call die('flux parse failed: '//trim(message))
        call apply_flux_turns(trim(flux_turn_path), loops)
    end if

    if (have_seg) then
        if (seg_area <= 0.0_dp) then
            call die('--seg-area is required for segmented Rogowski diagnostics')
        end if
        call read_segmented_rogowski_file(segrog_path, segs, ierr, message, &
            seg_area)
        if (ierr /= 0_i32) call die('segrog parse failed: '//trim(message))
        call apply_segrog_turns(trim(segrog_turn_path), segs)
    end if

    if (have_flux .and. have_seg) then
        call solver%flux_and_segrog(loops, fluxes, segs, voltages, rule)
        call scale_flux_turns(loops, fluxes)
        call scale_segrog_turns(segs, voltages)
        call write_result(flux_out_path, loops, fluxes)
        call write_segrog(segrog_out_path, segs, voltages)
    else if (have_flux) then
        call solver%flux_loops(loops, fluxes, rule)
        call scale_flux_turns(loops, fluxes)
        call write_result(flux_out_path, loops, fluxes)
    else if (have_seg) then
        call solver%segrog(segs, voltages, rule)
        call scale_segrog_turns(segs, voltages)
        call write_segrog(segrog_out_path, segs, voltages)
    end if

    call solver%finalize()
end subroutine run_solver

subroutine ensure_directory(path)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: path
    character(len=1024) :: cmd
    integer :: status
    if (len_trim(path) == 0) return
    write(cmd, '(A)') 'mkdir -p '//trim(path)
    call execute_command_line(trim(cmd), exitstat=status)
    if (status /= 0) then
        write(error_unit, '(A)') 'failed to create directory: '//trim(path)
        stop 1
    end if
end subroutine ensure_directory

subroutine write_result(path, loops, values)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    character(len=*), intent(in) :: path
    type(flux_loop_t), allocatable, intent(in) :: loops(:)
    real(dp), allocatable, intent(in) :: values(:)

    integer :: unit
    integer :: i

    open(newunit=unit, file=path, action='write', status='replace')
    write(unit, '(A)') 'label,value'
    do i = 1, size(values)
        write(unit, '(A,",",ES24.16)') trim(loops(i)%label), values(i)
    end do
    close(unit)
end subroutine write_result

subroutine write_segrog(path, diagnostics, values)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    character(len=*), intent(in) :: path
    type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
    real(dp), allocatable, intent(in) :: values(:)

    integer :: unit
    integer :: i

    open(newunit=unit, file=path, action='write', status='replace')
    write(unit, '(A)') 'label,value'
    do i = 1, size(values)
        write(unit, '(A,",",ES24.16)') trim(diagnostics(i)%label), values(i)
    end do
    close(unit)
end subroutine write_segrog

subroutine usage_and_stop()
    use, intrinsic :: iso_fortran_env, only: error_unit
    write(error_unit, '(A)') 'Usage: tiago_vacuum_cli <coil> <flux> <segrog>'
    write(error_unit, '(A)') '       [--output-dir dir] [--flux-out file]'
    write(error_unit, '(A)') '       [--segrog-out file] [--samples N]'
    write(error_unit, '(A)') '       [--seg-area value] [--nfp value]'
    write(error_unit, '(A)') '       [--coil-extcur vmec_input_or_list]'
    write(error_unit, '(A)') '       [--flux-turns file] [--segrog-turns file]'
    stop 1
end subroutine usage_and_stop

subroutine die(message)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: message
    write(error_unit, '(A)') trim(message)
    stop 1
end subroutine die

subroutine apply_flux_turns(path, loops)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    character(len=*), intent(in) :: path
    type(flux_loop_t), allocatable, intent(inout) :: loops(:)

    if (len_trim(path) == 0) return
    if (.not. allocated(loops)) return
    call read_turn_file_flux(path, loops)
end subroutine apply_flux_turns

subroutine apply_segrog_turns(path, diagnostics)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    character(len=*), intent(in) :: path
    type(segmented_rogowski_t), allocatable, intent(inout) :: diagnostics(:)

    if (len_trim(path) == 0) return
    if (.not. allocated(diagnostics)) return
    call read_turn_file_seg(path, diagnostics)
end subroutine apply_segrog_turns

subroutine read_turn_file_flux(path, loops)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    character(len=*), intent(in) :: path
    type(flux_loop_t), allocatable, intent(inout) :: loops(:)

    integer :: unit
    integer :: ios
    character(len=256) :: line
    character(len=128) :: label
    real(dp) :: value
    character(len=:), allocatable :: trimmed

    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A)') 'unable to open turn file: '//trim(path)
        stop 1
    end if

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        trimmed = adjustl(line)
        if (len_trim(trimmed) == 0) cycle
        if (trimmed(1:1) == '#') cycle
        read(trimmed, *, iostat=ios) label, value
        if (ios /= 0) cycle
        call assign_flux_turn(trim(label), value, loops)
    end do
    close(unit)
end subroutine read_turn_file_flux

subroutine assign_flux_turn(label, value, loops)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value
    type(flux_loop_t), allocatable, intent(inout) :: loops(:)

    integer :: i

    do i = 1, size(loops)
        if (trim(loops(i)%label) == trim(label)) then
            loops(i)%turn_scale = value
            return
        end if
    end do
end subroutine assign_flux_turn

subroutine read_turn_file_seg(path, diagnostics)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    character(len=*), intent(in) :: path
    type(segmented_rogowski_t), allocatable, intent(inout) :: diagnostics(:)

    integer :: unit
    integer :: ios
    character(len=256) :: line
    character(len=128) :: label
    real(dp) :: value
    character(len=:), allocatable :: trimmed

    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A)') 'unable to open turn file: '//trim(path)
        stop 1
    end if

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        trimmed = adjustl(line)
        if (len_trim(trimmed) == 0) cycle
        if (trimmed(1:1) == '#') cycle
        read(trimmed, *, iostat=ios) label, value
        if (ios /= 0) cycle
        call assign_seg_turn(trim(label), value, diagnostics)
    end do
    close(unit)
end subroutine read_turn_file_seg

subroutine assign_seg_turn(label, value, diagnostics)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value
    type(segmented_rogowski_t), allocatable, intent(inout) :: diagnostics(:)

    integer :: i

    do i = 1, size(diagnostics)
        if (trim(diagnostics(i)%label) == trim(label)) then
            diagnostics(i)%turn_scale = value
            return
        end if
    end do
end subroutine assign_seg_turn

subroutine scale_flux_turns(loops, fluxes)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    type(flux_loop_t), allocatable, intent(in) :: loops(:)
    real(dp), allocatable, intent(inout) :: fluxes(:)
    integer :: i

    if (.not. allocated(loops)) return
    if (.not. allocated(fluxes)) return
    do i = 1, min(size(loops), size(fluxes))
        fluxes(i) = fluxes(i) * loops(i)%turn_scale
    end do
end subroutine scale_flux_turns

subroutine scale_segrog_turns(diagnostics, values)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    type(segmented_rogowski_t), allocatable, intent(in) :: diagnostics(:)
    real(dp), allocatable, intent(inout) :: values(:)
    integer :: i

    if (.not. allocated(diagnostics)) return
    if (.not. allocated(values)) return
    do i = 1, min(size(diagnostics), size(values))
        values(i) = values(i) * diagnostics(i)%turn_scale
    end do
end subroutine scale_segrog_turns
end program tiago_vacuum_cli
