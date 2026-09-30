program tiago_vacuum_cli
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, &
        error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file, finalize_flux_signals
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    use tiago_build_config, only: simsopt_sample_path, simsopt_sample_url, &
        download_script_path
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
    character(len=512) :: plasma_wout
    integer(i32) :: plasma_nphi
    integer(i32) :: plasma_ntheta
    integer :: argc
    character(len=512) :: coil_extcur_path
    logical :: use_plasma_sample

    argc = command_argument_count()
    call check_help(argc)
    if (argc < 3) call usage_and_stop(1)

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
    seg_area = 0.0_dp   ! 0: not given (areas from the file)
    coil_extcur_path = ''
    have_flux = len_trim(flux_path) > 0
    have_seg = len_trim(segrog_path) > 0
    plasma_wout = ''
    plasma_nphi = 64_i32
    plasma_ntheta = 64_i32
    use_plasma_sample = .false.

    call parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        flux_turn_path, segrog_turn_path, samples_per_segment, seg_area, &
        nfp_value, coil_extcur_path, plasma_wout, plasma_nphi, plasma_ntheta, &
        use_plasma_sample)
    call validate_options(samples_per_segment, nfp_value, seg_area, plasma_nphi, plasma_ntheta)
    call prepare_plasma_support(plasma_wout, use_plasma_sample)
    call ensure_paths(output_dir, flux_out_path, segrog_out_path)
    call run_solver(trim(coil_path), trim(flux_path), trim(segrog_path), &
        trim(output_dir), trim(flux_out_path), trim(segrog_out_path), &
        trim(flux_turn_path), trim(segrog_turn_path), samples_per_segment, &
        seg_area, have_flux, have_seg, nfp_value, trim(coil_extcur_path), &
        trim(plasma_wout), plasma_nphi, plasma_ntheta)
contains

subroutine parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        flux_turn_path, segrog_turn_path, samples_per_segment, seg_area, &
        nfp_value, coil_extcur_path, plasma_wout, plasma_nphi, plasma_ntheta, &
        use_plasma_sample)
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
    character(len=*), intent(inout) :: plasma_wout
    integer(i32), intent(inout) :: plasma_nphi
    integer(i32), intent(inout) :: plasma_ntheta
    logical, intent(inout) :: use_plasma_sample

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
            samples_per_segment = parse_int(arg, '--samples')
        case ('--nfp')
            i = i + 1
            call ensure_arg(argc, i, '--nfp')
            call get_command_argument(i, arg)
            nfp_value = parse_int(arg, '--nfp')
        case ('--seg-area')
            i = i + 1
            call ensure_arg(argc, i, '--seg-area')
            call get_command_argument(i, arg)
            seg_area = parse_real(arg, '--seg-area')
        case ('--coil-extcur')
            i = i + 1
            call ensure_arg(argc, i, '--coil-extcur')
            call get_command_argument(i, coil_extcur_path)
        case ('--plasma-wout')
            i = i + 1
            call ensure_arg(argc, i, '--plasma-wout')
            call get_command_argument(i, plasma_wout)
        case ('--plasma-nphi')
            i = i + 1
            call ensure_arg(argc, i, '--plasma-nphi')
            call get_command_argument(i, arg)
            plasma_nphi = parse_int(arg, '--plasma-nphi')
        case ('--plasma-ntheta')
            i = i + 1
            call ensure_arg(argc, i, '--plasma-ntheta')
            call get_command_argument(i, arg)
            plasma_ntheta = parse_int(arg, '--plasma-ntheta')
        case ('--plasma-sample')
            use_plasma_sample = .true.
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

subroutine prepare_plasma_support(plasma_wout, use_plasma_sample)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(inout) :: plasma_wout
    logical, intent(in) :: use_plasma_sample

    if (.not. use_plasma_sample) then
        call ensure_sample_file(plasma_wout)
        return
    end if


    if (len_trim(plasma_wout) > 0) then
        write(error_unit, '(A)') 'WARNING: --plasma-sample ignored because '// &
            '--plasma-wout was provided'
        call ensure_sample_file(plasma_wout)
        return
    end if

    if (len(simsopt_sample_path) == 0) then
        call die('--plasma-sample: no sample path configured in this build')
    end if
    plasma_wout = simsopt_sample_path
    call ensure_sample_file(plasma_wout)
end subroutine prepare_plasma_support

subroutine ensure_sample_file(plasma_wout)
    character(len=*), intent(in) :: plasma_wout
    logical :: exists

    if (.not. needs_sample_download(plasma_wout)) return

    call ensure_parent_directory(simsopt_sample_path)
    inquire (file=trim(simsopt_sample_path), exist=exists)
    if (exists) return
    call download_sample_vmec()
end subroutine ensure_sample_file

logical function needs_sample_download(plasma_wout)
    character(len=*), intent(in) :: plasma_wout
    character(len=:), allocatable :: trimmed

    trimmed = trim(plasma_wout)
    needs_sample_download = len_trim(trimmed) > 0 .and. &
        trimmed == trim(simsopt_sample_path)
end function needs_sample_download

subroutine ensure_parent_directory(path)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: trimmed
    character(len=:), allocatable :: directory
    integer :: slash

    trimmed = trim(path)
    slash = scan(trimmed, '/', back=.true.)
    if (slash <= 1) return
    directory = trimmed(:slash - 1)
    call ensure_directory(directory)
end subroutine ensure_parent_directory

subroutine download_sample_vmec()
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=4096) :: command
    integer :: status

    write(error_unit, '(A)') 'INFO: downloading Simsopt VMEC sample to '// &
        trim(simsopt_sample_path)
    write(command, '(A)') 'cmake "-DURL='//trim(simsopt_sample_url)//'" '// &
        '"-DDEST='//trim(simsopt_sample_path)//'" -P "'// &
        trim(download_script_path)//'"'
    call execute_command_line(trim(command), exitstat=status)
    if (status /= 0) then
        call die('failed to download Simsopt VMEC sample; ensure cmake is available')
    end if
end subroutine download_sample_vmec

subroutine run_solver(coil_path, flux_path, segrog_path, output_dir, &
        flux_out_path, segrog_out_path, flux_turn_path, segrog_turn_path, &
        samples_per_segment, seg_area, have_flux, have_seg, nfp_value, &
        coil_extcur_path, plasma_wout, plasma_nphi, plasma_ntheta)
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, &
        error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t
    use tiago_flux_loops, only: read_flux_loop_file, finalize_flux_signals
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
    character(len=*), intent(in) :: plasma_wout
    integer(i32), intent(in) :: plasma_nphi
    integer(i32), intent(in) :: plasma_ntheta

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

    if (len_trim(plasma_wout) > 0) then
        call solver%enable_plasma_from_vmec(plasma_wout, plasma_nphi, plasma_ntheta)
    end if

    if (have_flux) then
        call read_flux_loop_file(flux_path, loops, ierr, message)
        if (ierr /= 0_i32) call die('flux parse failed: '//trim(message))
        call apply_flux_turns(trim(flux_turn_path), loops)
    end if

    if (have_seg) then
        if (seg_area > 0.0_dp) then
            call read_segmented_rogowski_file(segrog_path, segs, ierr, message, seg_area)
        else
            call read_segmented_rogowski_file(segrog_path, segs, ierr, message)
        end if
        if (ierr /= 0_i32) call die('segrog parse failed: '//trim(message))
        if (any_missing_area(segs)) then
            call die('segmented Rogowski file has no eff_area column; pass --seg-area')
        end if
        call apply_segrog_turns(trim(segrog_turn_path), segs)
    end if

    if (have_flux .and. have_seg) then
        call solver%flux_and_segrog(loops, fluxes, segs, voltages, rule)
        call finalize_flux_signals(loops, fluxes, solver%plasma%diamagnetic_flux)
        call scale_segrog_turns(segs, voltages)
        call write_result(flux_out_path, loops, fluxes)
        call write_segrog(segrog_out_path, segs, voltages)
    else if (have_flux) then
        call solver%flux_loops(loops, fluxes, rule)
        call finalize_flux_signals(loops, fluxes, solver%plasma%diamagnetic_flux)
        call write_result(flux_out_path, loops, fluxes)
    else if (have_seg) then
        call solver%segrog(segs, voltages, rule)
        call scale_segrog_turns(segs, voltages)
        call write_segrog(segrog_out_path, segs, voltages)
    end if

    call solver%finalize()

    if (have_flux) call check_finite('flux', loops_labels(loops), fluxes)
    if (have_seg) call check_finite('segrog', segrog_labels(segs), voltages)
end subroutine run_solver

logical function any_missing_area(segs)
    type(segmented_rogowski_t), intent(in) :: segs(:)
    integer :: i

    any_missing_area = .false.
    do i = 1, size(segs)
        if (all(segs(i)%segment_area == 0.0_dp)) any_missing_area = .true.
    end do
end function any_missing_area

subroutine check_finite(kind, labels, values)
    !! Outputs are written first; a non-finite value (e.g. a sensor point exactly on a
    !! coil filament) then turns into a non-zero exit status.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*), intent(in) :: kind
    character(len=128), intent(in) :: labels(:)
    real(dp), intent(in) :: values(:)
    integer :: i
    logical :: bad

    bad = .false.
    do i = 1, size(values)
        if (.not. ieee_is_finite(values(i))) then
            write(error_unit, '(A)') 'non-finite '//kind//' signal: '//trim(labels(i))
            bad = .true.
        end if
    end do
    if (bad) stop 2
end subroutine check_finite

function loops_labels(loops) result(labels)
    type(flux_loop_t), intent(in) :: loops(:)
    character(len=128) :: labels(size(loops))
    integer :: i
    do i = 1, size(loops)
        labels(i) = loops(i)%label
    end do
end function loops_labels

function segrog_labels(segs) result(labels)
    type(segmented_rogowski_t), intent(in) :: segs(:)
    character(len=128) :: labels(size(segs))
    integer :: i
    do i = 1, size(segs)
        labels(i) = segs(i)%label
    end do
end function segrog_labels

subroutine ensure_directory(path)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: path
    character(len=1024) :: cmd
    integer :: status
    if (len_trim(path) == 0) return
    write(cmd, '(A)') 'mkdir -p "'//trim(path)//'"'
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

subroutine check_help(argc)
    integer, intent(in) :: argc
    character(len=64) :: arg
    integer :: i

    do i = 1, argc
        call get_command_argument(i, arg)
        if (trim(arg) == '--help' .or. trim(arg) == '-h') call usage_and_stop(0)
    end do
end subroutine check_help

integer(i32) function parse_int(text, flag) result(value)
    character(len=*), intent(in) :: text, flag
    integer :: ios
    read(text, *, iostat=ios) value
    if (ios /= 0 .or. verify(trim(adjustl(text)), '+-0123456789') /= 0) then
        call die(flag//' expects an integer, got: '//trim(text))
    end if
end function parse_int

real(dp) function parse_real(text, flag) result(value)
    character(len=*), intent(in) :: text, flag
    integer :: ios
    read(text, *, iostat=ios) value
    if (ios /= 0) call die(flag//' expects a number, got: '//trim(text))
end function parse_real

subroutine validate_options(samples_per_segment, nfp_value, seg_area, plasma_nphi, plasma_ntheta)
    integer(i32), intent(in) :: samples_per_segment, nfp_value, plasma_nphi, plasma_ntheta
    real(dp), intent(in) :: seg_area

    if (samples_per_segment < 1) call die('--samples must be at least 1')
    if (nfp_value < 1) call die('--nfp must be at least 1')
    if (seg_area < 0.0_dp) call die('--seg-area must be positive')
    if (plasma_nphi < 4 .or. plasma_ntheta < 4) then
        call die('--plasma-nphi and --plasma-ntheta must be at least 4')
    end if
end subroutine validate_options

subroutine usage_and_stop(status)
    use, intrinsic :: iso_fortran_env, only: error_unit
    integer, intent(in) :: status
    write(error_unit, '(A)') 'Usage: tiago_vacuum_cli <coil> <flux> <segrog>   (pass "" to omit one)'
    write(error_unit, '(A)') '       [--output-dir dir] [--flux-out file]'
    write(error_unit, '(A)') '       [--segrog-out file] [--samples N]'
    write(error_unit, '(A)') '       [--seg-area value] [--nfp value]'
    write(error_unit, '(A)') '       [--coil-extcur vmec_input_or_list]'
    write(error_unit, '(A)') '       [--flux-turns file] [--segrog-turns file]'
    write(error_unit, '(A)') '       [--plasma-wout file] [--plasma-sample]'
    write(error_unit, '(A)') '       [--plasma-nphi N_per_period] [--plasma-ntheta N] (default 64 64)'
    if (status == 0) stop
    stop 1
end subroutine usage_and_stop

subroutine die(message)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: message
    write(error_unit, '(A)') trim(message)
    stop 1
end subroutine die

subroutine apply_flux_turns(path, loops)
    character(len=*), intent(in) :: path
    type(flux_loop_t), allocatable, intent(inout) :: loops(:)
    real(dp), allocatable :: scales(:)
    integer :: i

    if (len_trim(path) == 0 .or. .not. allocated(loops)) return
    call read_turn_file(path, loops_labels(loops), scales)
    do i = 1, size(loops)
        loops(i)%turn_scale = scales(i)
    end do
end subroutine apply_flux_turns

subroutine apply_segrog_turns(path, diagnostics)
    character(len=*), intent(in) :: path
    type(segmented_rogowski_t), allocatable, intent(inout) :: diagnostics(:)
    real(dp), allocatable :: scales(:)
    integer :: i

    if (len_trim(path) == 0 .or. .not. allocated(diagnostics)) return
    call read_turn_file(path, segrog_labels(diagnostics), scales)
    do i = 1, size(diagnostics)
        diagnostics(i)%turn_scale = scales(i)
    end do
end subroutine apply_segrog_turns

subroutine read_turn_file(path, labels, scales)
    !! Lines "label scale" (the scale is the last token, so labels may contain
    !! blanks; "label, scale" works too). '#' starts a comment line. Labels
    !! without an entry keep scale 1; entries matching no diagnostic are errors.
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: labels(:)
    real(dp), allocatable, intent(out) :: scales(:)

    integer :: unit, ios, cut, i
    character(len=512) :: line
    character(len=:), allocatable :: text, label
    real(dp) :: value
    logical :: found

    allocate(scales(size(labels)))
    scales = 1.0_dp
    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) call die('unable to open turn file: '//trim(path))
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        text = trim(adjustl(line))
        if (len(text) == 0) cycle
        if (text(1:1) == '#') cycle
        cut = scan(text, ' ,'//achar(9), back=.true.)
        if (cut == 0) call die(trim(path)//': expected "label scale": '//text)
        read(text(cut + 1:), *, iostat=ios) value
        if (ios /= 0) call die(trim(path)//': invalid scale: '//text)
        label = trim(text(:cut - 1))
        if (len(label) > 0) then
            if (label(len(label):len(label)) == ',') label = trim(label(:len(label) - 1))
        end if
        found = .false.
        do i = 1, size(labels)
            if (trim(labels(i)) == label) then
                scales(i) = value
                found = .true.
            end if
        end do
        if (.not. found) call die(trim(path)//': no diagnostic labelled '//label)
    end do
    close(unit)
end subroutine read_turn_file

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
