program tiago_vacuum_cli
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, &
        error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_flux_loops, only: read_flux_loop_file, finalize_flux_signals
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_bprobes, only: read_bprobe_file
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
    logical :: covariant_plasma = .false., conservative_plasma = .false.
    ! Magnetic probes and response matrices (set in parse_options by host association)
    character(len=512) :: bprobe_path = '', bprobe_out_path = 'tiago_bprobes.csv'
    character(len=512) :: bprobe_turn_path = '', response_out_path = ''
    character(len=512) :: plasma_response_out_path = '', plasma_shape_out_path = ''
    logical :: rphiz = .false.
    logical :: use_gauss = .false.   ! Gauss-Legendre instead of midpoint samples

    argc = command_argument_count()
    call check_help(argc)
    if (argc < 1) call usage_and_stop(1)

    coil_path = ''
    flux_path = ''
    segrog_path = ''
    output_dir = '.'
    flux_out_path = 'tiago_flux.csv'
    segrog_out_path = 'tiago_segrog.csv'
    flux_turn_path = ''
    segrog_turn_path = ''
    samples_per_segment = 6_i32
    nfp_value = 1_i32
    seg_area = 0.0_dp   ! 0: not given (areas from the file)
    coil_extcur_path = ''
    plasma_wout = ''
    plasma_nphi = 64_i32
    plasma_ntheta = 64_i32
    use_plasma_sample = .false.

    call parse_options(argc, output_dir, flux_out_path, segrog_out_path, &
        flux_turn_path, segrog_turn_path, samples_per_segment, seg_area, &
        nfp_value, coil_extcur_path, plasma_wout, plasma_nphi, plasma_ntheta, &
        use_plasma_sample)
    have_flux = len_trim(flux_path) > 0
    have_seg = len_trim(segrog_path) > 0
    if (.not. (have_flux .or. have_seg .or. len_trim(bprobe_path) > 0)) then
        call die('nothing to evaluate: pass --flux, --segrog and/or --bprobes')
    end if
    if (len_trim(coil_path) == 0 .and. len_trim(plasma_wout) == 0 .and. &
        .not. use_plasma_sample) then
        call die('no field source: pass --coils and/or --plasma-wout')
    end if
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

    i = 1
    do while (i <= argc)
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ('--coils')
            i = i + 1
            call ensure_arg(argc, i, '--coils')
            call get_command_argument(i, coil_path)
        case ('--flux')
            i = i + 1
            call ensure_arg(argc, i, '--flux')
            call get_command_argument(i, flux_path)
        case ('--segrog')
            i = i + 1
            call ensure_arg(argc, i, '--segrog')
            call get_command_argument(i, segrog_path)
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
        case ('--plasma-conservative')
            covariant_plasma = .true.
            conservative_plasma = .true.
        case ('--plasma-covariant')
            covariant_plasma = .true.
            conservative_plasma = .false.
        case ('--plasma-contravariant')
            covariant_plasma = .false.
            conservative_plasma = .false.
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
        case ('--bprobes')
            i = i + 1
            call ensure_arg(argc, i, '--bprobes')
            call get_command_argument(i, bprobe_path)
        case ('--bprobe-out')
            i = i + 1
            call ensure_arg(argc, i, '--bprobe-out')
            call get_command_argument(i, bprobe_out_path)
        case ('--bprobe-turns')
            i = i + 1
            call ensure_arg(argc, i, '--bprobe-turns')
            call get_command_argument(i, bprobe_turn_path)
        case ('--rphiz')
            rphiz = .true.
        case ('--gauss')
            use_gauss = .true.
        case ('--plasma-response-out')
            i = i + 1
            call ensure_arg(argc, i, '--plasma-response-out')
            call get_command_argument(i, plasma_response_out_path)
        case ('--plasma-shape-response-out')
            i = i + 1
            call ensure_arg(argc, i, '--plasma-shape-response-out')
            call get_command_argument(i, plasma_shape_out_path)
        case ('--response-out')
            i = i + 1
            call ensure_arg(argc, i, '--response-out')
            call get_command_argument(i, response_out_path)
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

    call in_output_dir(output_dir, flux_out_path)
    call in_output_dir(output_dir, segrog_out_path)
    call in_output_dir(output_dir, bprobe_out_path)
    if (len_trim(response_out_path) > 0) call in_output_dir(output_dir, response_out_path)
    if (len_trim(plasma_response_out_path) > 0) then
        call in_output_dir(output_dir, plasma_response_out_path)
    end if
    if (len_trim(plasma_shape_out_path) > 0) call in_output_dir(output_dir, plasma_shape_out_path)
end subroutine ensure_paths

subroutine in_output_dir(output_dir, path)
    !! Bare file names go into the output directory.
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(inout) :: path
    if (len_trim(output_dir) > 0 .and. index(path, '/') == 0) then
        path = trim(output_dir)//'/'//trim(path)
    end if
end subroutine in_output_dir

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
    use tiago_bprobes, only: read_bprobe_file
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
    type(bprobe_t), allocatable :: probes(:)
    real(dp), allocatable :: probe_values(:)
    integer(i32) :: ierr
    character(len=:), allocatable :: message

    call ensure_directory(output_dir)

    rule%samples_per_segment = samples_per_segment
    rule%gauss = use_gauss
    if (len_trim(coil_extcur_path) > 0) then
        call solver%init(coil_path, coil_extcur_path)
    else
        call solver%init(coil_path)
    end if
    call solver%set_nfp(nfp_value)

    if (len_trim(plasma_wout) > 0) then
        call solver%enable_plasma_from_vmec(plasma_wout, plasma_nphi, &
            plasma_ntheta, covariant_plasma, conservative_plasma)
        call write_plasma_provenance(solver, trim(output_dir)//'/plasma_model.txt')
        if (solver%plasma%covariant) then
            write(error_unit, '(A,ES12.4)') &
                'Plasma Fourier curl norm [T m]: ', solver%plasma%curl_norm
            write(error_unit, '(A,ES12.4)') &
                'Plasma toroidal current ripple bound [A]: ', &
                solver%plasma%current_ripple
            if (solver%plasma%conservative) then
                write(error_unit, '(A,ES12.4)') &
                    'Plasma coefficient projection norm [T m]: ', &
                    solver%plasma%projection_norm
            end if
        end if
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

    if (len_trim(bprobe_path) > 0) then
        call read_bprobe_file(trim(bprobe_path), probes, ierr, message, rphiz)
        if (ierr /= 0_i32) call die('B-probe parse failed: '//trim(message))
        call apply_bprobe_turns(trim(bprobe_turn_path), probes)
        call solver%bprobes(probes, probe_values)
        probe_values = probe_values * probes%turn_scale
        call write_values(bprobe_out_path, probe_labels(probes), probe_values)
    end if

    if (len_trim(response_out_path) > 0) then
        call write_response(solver, rule, loops, segs, probes, trim(response_out_path))
    end if

    if (len_trim(plasma_response_out_path) + len_trim(plasma_shape_out_path) > 0) then
        if (len_trim(plasma_wout) == 0) then
            call die('--plasma-response-out/--plasma-shape-response-out need --plasma-wout')
        end if
        call write_plasma_response(solver, rule, loops, segs, probes, &
            trim(plasma_response_out_path), trim(plasma_shape_out_path))
    end if

    call solver%finalize()

    if (have_flux) call check_finite('flux', loops_labels(loops), fluxes)
    if (have_seg) call check_finite('segrog', segrog_labels(segs), voltages)
    if (allocated(probes)) call check_finite('B-probe', probe_labels(probes), probe_values)
end subroutine run_solver

subroutine write_response(solver, rule, loops, segs, probes, path)
    !! kind,label,group,value: signal per unit EXTCUR of each coil group, with
    !! the same turns/idia post-processing as the signals (no plasma part).
    type(vacuum_solver_t), intent(in) :: solver
    type(quadrature_rule_t), intent(in) :: rule
    type(flux_loop_t), allocatable, intent(in) :: loops(:)
    type(segmented_rogowski_t), allocatable, intent(in) :: segs(:)
    type(bprobe_t), allocatable, intent(in) :: probes(:)
    character(len=*), intent(in) :: path
    real(dp), allocatable :: flux_resp(:, :), seg_resp(:, :), probe_resp(:, :)
    integer :: unit, g, i

    call solver%response(loops, segs, probes, flux_resp, seg_resp, probe_resp, rule)
    open(newunit=unit, file=path, action='write', status='replace')
    write(unit, '(A)') 'kind,label,group,value'
    do g = 1, solver%n_groups()
        if (allocated(loops)) then
            call finalize_flux_signals(loops, flux_resp(:, g))
            do i = 1, size(loops)
                write(unit, '(A,",",I0,",",ES24.16)') 'flux,'//loops(i)%label, g, flux_resp(i, g)
            end do
        end if
        if (allocated(segs)) then
            do i = 1, size(segs)
                write(unit, '(A,",",I0,",",ES24.16)') 'segrog,'//segs(i)%label, g, &
                    seg_resp(i, g) * segs(i)%turn_scale
            end do
        end if
        if (allocated(probes)) then
            do i = 1, size(probes)
                write(unit, '(A,",",I0,",",ES24.16)') 'bprobe,'//probes(i)%label, g, &
                    probe_resp(i, g) * probes(i)%turn_scale
            end do
        end if
    end do
    close(unit)
end subroutine write_response

subroutine write_plasma_provenance(solver, path)
    type(vacuum_solver_t), intent(in) :: solver
    character(len=*), intent(in) :: path
    integer :: unit, status

    open(newunit=unit, file=path, status='replace', iostat=status)
    if (status /= 0) call die('cannot write plasma model provenance: '//path)
    write(unit, '(A,L1)') 'covariant = ', solver%plasma%covariant
    write(unit, '(A,L1)') 'conservative_projection = ', solver%plasma%conservative
    write(unit, '(A,ES24.16)') 'Fourier_curl_norm_T_m = ', solver%plasma%curl_norm
    write(unit, '(A,ES24.16)') 'toroidal_current_ripple_bound_A = ', &
        solver%plasma%current_ripple
    write(unit, '(A,ES24.16)') 'coefficient_projection_norm_T_m = ', &
        solver%plasma%projection_norm
    close(unit)
end subroutine write_plasma_provenance

subroutine write_plasma_response(solver, rule, loops, segs, probes, path, shape_path)
    !! kind,label,coefficient,m,n,value: derivative of each signal's plasma part
    !! with respect to the VMEC boundary (s = 1) field coefficients (path) and,
    !! if shape_path is given, the boundary geometry coefficients rmnc, zmns, ...
    !! (at fixed field coefficients). After turns and idia < 0 differences; the
    !! idia = 1 phiedge term does not depend on them.
    type(vacuum_solver_t), intent(in) :: solver
    type(quadrature_rule_t), intent(in) :: rule
    type(flux_loop_t), allocatable, intent(in) :: loops(:)
    type(segmented_rogowski_t), allocatable, intent(in) :: segs(:)
    type(bprobe_t), allocatable, intent(in) :: probes(:)
    character(len=*), intent(in) :: path, shape_path
    real(dp), allocatable :: flux_resp(:, :), seg_resp(:, :), probe_resp(:, :)
    real(dp), allocatable :: flux_shape(:, :), seg_shape(:, :), probe_shape(:, :)
    character(len=:), allocatable :: coefficient
    integer :: c, m, n, unit

    if (len_trim(shape_path) > 0) then
        call solver%plasma_response(loops, segs, probes, flux_resp, seg_resp, probe_resp, rule, &
            flux_shape, seg_shape, probe_shape)
    else
        call solver%plasma_response(loops, segs, probes, flux_resp, seg_resp, probe_resp, rule)
    end if
    if (len_trim(path) > 0) then
        open(newunit=unit, file=path, action='write', status='replace')
        write(unit, '(A)') 'kind,label,coefficient,m,n,value'
        do c = 1, solver%plasma%n_mode_columns()
            call solver%plasma%mode_column_name(c, coefficient, m, n)
            call write_response_column(unit, coefficient, m, n, loops, segs, probes, &
                flux_resp, seg_resp, probe_resp, c)
        end do
        close(unit)
    end if
    if (len_trim(shape_path) > 0) then
        open(newunit=unit, file=shape_path, action='write', status='replace')
        write(unit, '(A)') 'kind,label,coefficient,m,n,value'
        do c = 1, solver%plasma%n_shape_columns()
            call solver%plasma%shape_column_name(c, coefficient, m, n)
            call write_response_column(unit, coefficient, m, n, loops, segs, probes, &
                flux_shape, seg_shape, probe_shape, c)
        end do
        close(unit)
    end if
end subroutine write_plasma_response

subroutine write_response_column(unit, coefficient, m, n, loops, segs, probes, &
        flux_resp, seg_resp, probe_resp, c)
    integer, intent(in) :: unit, m, n, c
    character(len=*), intent(in) :: coefficient
    type(flux_loop_t), allocatable, intent(in) :: loops(:)
    type(segmented_rogowski_t), allocatable, intent(in) :: segs(:)
    type(bprobe_t), allocatable, intent(in) :: probes(:)
    real(dp), intent(inout) :: flux_resp(:, :)
    real(dp), intent(in) :: seg_resp(:, :), probe_resp(:, :)
    integer :: i

    if (allocated(loops)) then
        call finalize_flux_signals(loops, flux_resp(:, c))
        do i = 1, size(loops)
            write(unit, '(A,2(",",I0),",",ES24.16)') 'flux,'//loops(i)%label//','// &
                coefficient, m, n, flux_resp(i, c)
        end do
    end if
    if (allocated(segs)) then
        do i = 1, size(segs)
            write(unit, '(A,2(",",I0),",",ES24.16)') 'segrog,'//segs(i)%label//','// &
                coefficient, m, n, seg_resp(i, c) * segs(i)%turn_scale
        end do
    end if
    if (allocated(probes)) then
        do i = 1, size(probes)
            write(unit, '(A,2(",",I0),",",ES24.16)') 'bprobe,'//probes(i)%label//','// &
                coefficient, m, n, probe_resp(i, c) * probes(i)%turn_scale
        end do
    end if
end subroutine write_response_column

subroutine apply_bprobe_turns(path, probes)
    character(len=*), intent(in) :: path
    type(bprobe_t), intent(inout) :: probes(:)
    real(dp), allocatable :: scales(:)

    if (len_trim(path) == 0) return
    call read_turn_file(path, probe_labels(probes), scales)
    probes%turn_scale = scales
end subroutine apply_bprobe_turns

function probe_labels(probes) result(labels)
    type(bprobe_t), intent(in) :: probes(:)
    character(len=128) :: labels(size(probes))
    integer :: i
    do i = 1, size(probes)
        labels(i) = probes(i)%label
    end do
end function probe_labels

subroutine write_values(path, labels, values)
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: labels(:)
    real(dp), intent(in) :: values(:)
    integer :: unit, i

    open(newunit=unit, file=path, action='write', status='replace')
    write(unit, '(A)') 'label,value'
    do i = 1, size(values)
        write(unit, '(A,",",ES24.16)') trim(labels(i)), values(i)
    end do
    close(unit)
end subroutine write_values

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
    write(error_unit, '(A)') 'Usage: tiago_vacuum_cli [--coils file] [--flux file] [--segrog file]'
    write(error_unit, '(A)') '       (at least one field source: --coils and/or --plasma-wout)'
    write(error_unit, '(A)') '       [--output-dir dir] [--flux-out file]'
    write(error_unit, '(A)') '       [--segrog-out file] [--samples N] [--gauss]'
    write(error_unit, '(A)') '       [--seg-area value] [--nfp value]'
    write(error_unit, '(A)') '       [--coil-extcur vmec_input_or_list]'
    write(error_unit, '(A)') '       [--flux-turns file] [--segrog-turns file]'
    write(error_unit, '(A)') '       [--plasma-wout file] [--plasma-sample]'
    write(error_unit, '(A)') '       [--plasma-covariant | --plasma-conservative | '// &
        '--plasma-contravariant]'
    write(error_unit, '(A)') '       [--plasma-nphi N_per_period] [--plasma-ntheta N] (default 64 64)'
    write(error_unit, '(A)') '       [--bprobes file [--rphiz] [--bprobe-out file] [--bprobe-turns file]]'
    write(error_unit, '(A)') '       [--response-out file]  (signals per unit EXTCUR per coil group)'
    write(error_unit, '(A)') '       [--plasma-response-out file]  (d signal / d VMEC boundary B coefficients)'
    write(error_unit, '(A)') '       [--plasma-shape-response-out file]  (d signal / d boundary rmnc, zmns)'
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
