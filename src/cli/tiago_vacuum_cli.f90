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
    integer(i32) :: samples_per_segment
    real(dp) :: seg_area
    logical :: have_flux
    logical :: have_seg
    integer :: argc

    argc = command_argument_count()
    if (argc < 3) call usage_and_stop()

    call get_command_argument(1, coil_path)
    call get_command_argument(2, flux_path)
    call get_command_argument(3, segrog_path)

    output_dir = '.'
    flux_out_path = 'tiago_flux.csv'
    segrog_out_path = 'tiago_segrog.csv'
    samples_per_segment = 6_i32
    seg_area = -1.0_dp
    have_flux = len_trim(flux_path) > 0
    have_seg = len_trim(segrog_path) > 0

    call parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        samples_per_segment, seg_area)
    call ensure_paths(output_dir, flux_out_path, segrog_out_path)
    call run_solver(trim(coil_path), trim(flux_path), trim(segrog_path), &
        trim(output_dir), trim(flux_out_path), trim(segrog_out_path), &
        samples_per_segment, seg_area, have_flux, have_seg)
contains

subroutine parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        samples_per_segment, seg_area)
    use, intrinsic :: iso_fortran_env, only: i32 => int32, dp => real64
    implicit none
    integer, intent(in) :: argc
    character(len=*), intent(inout) :: output_dir
    character(len=*), intent(inout) :: flux_out_path
    character(len=*), intent(inout) :: segrog_out_path
    integer(i32), intent(inout) :: samples_per_segment
    real(dp), intent(inout) :: seg_area

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
        case ('--samples')
            i = i + 1
            call ensure_arg(argc, i, '--samples')
            call get_command_argument(i, arg)
            read(arg, *) samples_per_segment
        case ('--seg-area')
            i = i + 1
            call ensure_arg(argc, i, '--seg-area')
            call get_command_argument(i, arg)
            read(arg, *) seg_area
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
        flux_out_path, segrog_out_path, samples_per_segment, seg_area, &
        have_flux, have_seg)
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
    integer(i32), intent(in) :: samples_per_segment
    real(dp), intent(in) :: seg_area
    logical, intent(in) :: have_flux
    logical, intent(in) :: have_seg

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
    call solver%init(coil_path)

    if (have_flux) then
        call read_flux_loop_file(flux_path, loops, ierr, message)
        if (ierr /= 0_i32) call die('flux parse failed: '//trim(message))
        call solver%flux_loops(loops, fluxes, rule)
        call write_result(flux_out_path, loops, fluxes)
    end if

    if (have_seg) then
        if (seg_area <= 0.0_dp) then
            call die('--seg-area is required for segmented Rogowski diagnostics')
        end if
        call read_segmented_rogowski_file(segrog_path, segs, ierr, message, &
            seg_area)
        if (ierr /= 0_i32) call die('segrog parse failed: '//trim(message))
        call solver%segrog(segs, voltages, rule)
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
    write(error_unit, '(A)') '       [--seg-area value]'
    stop 1
end subroutine usage_and_stop

subroutine die(message)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: message
    write(error_unit, '(A)') trim(message)
    stop 1
end subroutine die
end program tiago_vacuum_cli
