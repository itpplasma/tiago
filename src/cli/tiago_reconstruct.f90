program tiago_reconstruct
    !! Equilibrium reconstruction from magnetic diagnostics with VMEC++.
    !! Usage: tiago_reconstruct input.nml
    !! See README ("Equilibrium reconstruction") for the namelist.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_flux_loops, only: read_flux_loop_file, finalize_flux_signals
    use tiago_segmented_rogowski, only: read_segmented_rogowski_file
    use tiago_bprobes, only: read_bprobe_file
    use tiago_reconstruction, only: rec, n_rows, model_signals, residual_at, jacobian_at, &
        fit, covariance, synthesize, signal_value, KIND_FLUX, KIND_SEGROG, KIND_BPROBE
    use tiago_reconstruction_plots, only: equilibrium_profiles_t, read_profiles, plot_signals, &
        plot_residuals, plot_convergence, plot_profiles, plot_boundary, plot_jacobian_check
    implicit none

    integer, parameter :: maxp = 64
    ! namelist
    character(len=512) :: vmec_input = '', output_dir = 'reconstruction'
    character(len=512) :: coils = '', coil_extcur = ''
    character(len=512) :: flux = '', segrog = '', bprobes = '', measurements = ''
    character(len=512) :: consistency_points = ''
    character(len=64) :: parameters(maxp) = ''
    real(dp) :: start_values(maxp), truth_values(maxp), parameter_scale(maxp)
    real(dp) :: prior_values(maxp), prior_sigma(maxp)
    real(dp) :: consistency_sigma = 1.0e-3_dp, seg_area = 0.0_dp
    real(dp) :: sigma_relative = 0.01_dp, sigma_flux = 0.0_dp, sigma_segrog = 0.0_dp
    real(dp) :: sigma_bprobe = 0.0_dp, fd_step = 1.0e-4_dp, vmec_ftol = 1.0e-14_dp
    integer :: samples = 4, plasma_nphi = 32, plasma_ntheta = 32, max_iterations = 20, seed = 1
    integer :: vmec_niter = 20000
    logical :: gauss = .true., synthesize_measurements = .false., check_jacobian = .false., &
        add_noise = .true.
    logical :: rphiz = .false.
    namelist /reconstruction/ vmec_input, output_dir, coils, coil_extcur, &
        flux, segrog, bprobes, measurements, consistency_points, parameters, start_values, &
        truth_values, parameter_scale, prior_values, prior_sigma, consistency_sigma, seg_area, &
        sigma_relative, sigma_flux, sigma_segrog, sigma_bprobe, fd_step, vmec_ftol, vmec_niter, samples, plasma_nphi, &
        plasma_ntheta, max_iterations, seed, gauss, synthesize_measurements, add_noise, check_jacobian, rphiz

    character(len=512) :: input_file
    real(dp), allocatable :: x(:), x_start(:), x_truth(:), defaults(:), cov(:, :), jac(:, :)
    real(dp), allocatable :: r_start(:), r_fit(:)
    real(dp), allocatable :: f0(:), s0(:), p0(:), f1(:), s1(:), p1(:)
    character(len=:), allocatable :: message
    integer :: unit, ios, k, steps, info, n
    real(dp) :: t_start, t_end
    logical :: have_truth

    start_values = huge(1.0_dp)
    truth_values = huge(1.0_dp)
    parameter_scale = 0.0_dp
    prior_values = huge(1.0_dp)
    prior_sigma = 0.0_dp
    if (command_argument_count() /= 1) then
        write(error_unit, '(A)') 'Usage: tiago_reconstruct input.nml'
        stop 1
    end if
    call get_command_argument(1, input_file)
    open(newunit=unit, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) call die('cannot open '//trim(input_file))
    read(unit, nml=reconstruction, iostat=ios)
    if (ios /= 0) call die('error in namelist &reconstruction of '//trim(input_file))
    close(unit)
    if (len_trim(vmec_input) == 0) call die('vmec_input is required')
    call execute_command_line('mkdir -p "'//trim(output_dir)//'"')
    call cpu_time(t_start)

    call setup()
    n = rec%n
    ! parameter values: VMEC input / coil file, overridden by start_values
    allocate(defaults(n), x_start(n), x_truth(n))
    call parameter_defaults(defaults)
    do k = 1, n
        x_start(k) = merge(defaults(k), start_values(k), start_values(k) == huge(1.0_dp))
        x_truth(k) = truth_values(k)
    end do
    have_truth = all(x_truth /= huge(1.0_dp))
    do k = 1, n
        rec%scale(k) = parameter_scale(k)
        if (rec%scale(k) <= 0.0_dp) rec%scale(k) = max(abs(x_start(k)), abs(defaults(k)), 1.0e-3_dp)
        rec%prior(k) = merge(x_start(k), prior_values(k), prior_values(k) == huge(1.0_dp))
        rec%prior_sigma(k) = prior_sigma(k)
    end do
    rec%prior_param = pack([(k, k = 1, n)], rec%prior_sigma > 0.0_dp)
    rec%n_prior = size(rec%prior_param)

    if (synthesize_measurements) then
        if (.not. have_truth) call die('synthesize_measurements needs truth_values for every parameter')
        call synthesize(x_truth, sigma_relative, sigma_flux, sigma_segrog, sigma_bprobe, seed, &
            rec%meas_value, rec%meas_sigma, add_noise)
        call write_measurements(trim(output_dir)//'/measurements.csv')
        call rec%eq%write_wout(trim(output_dir)//'/wout_truth.nc')
    end if

    x = x_start
    allocate(r_start(n_rows()), r_fit(n_rows()))
    call residual_at(x, r_start)
    call model_signals(x, f0, s0, p0)
    call rec%eq%write_wout(trim(output_dir)//'/wout_initial.nc')

    if (check_jacobian) call jacobian_check(x)

    call fit(x, max_iterations, steps, message)
    call residual_at(x, r_fit)
    call model_signals(x, f1, s1, p1)
    call rec%eq%write_wout(trim(output_dir)//'/wout_fit.nc')
    call covariance(x, cov, jac, info)
    call cpu_time(t_end)

    call write_results()
    call make_plots()
    print '(A)', 'reconstruction written to '//trim(output_dir)

contains

    subroutine setup()
        use tiago_vacuum_forward, only: vacuum_solver_t
        integer(i32) :: ierr
        character(len=:), allocatable :: msg
        type(bprobe_t), allocatable :: measured(:), fb(:)
        integer :: k, neq, npar
        character(len=64) :: name
        character(len=64), allocatable :: eq_names(:)
        character(len=:), allocatable :: workdir

        npar = count(len_trim(parameters) > 0)
        if (npar == 0) call die('no parameters given')
        rec%n = npar
        allocate(rec%names(npar), rec%is_coil(npar), rec%group(npar), rec%eq_index(npar))
        allocate(rec%x0(npar), rec%scale(npar), rec%prior(npar), rec%prior_sigma(npar))
        allocate(rec%history(0))
        neq = 0
        allocate(eq_names(0))
        do k = 1, npar
            name = adjustl(parameters(k))
            rec%names(k) = name
            rec%is_coil(k) = index(upper(name), 'EXTCUR(') == 1
            rec%group(k) = 0
            rec%eq_index(k) = 0
            if (rec%is_coil(k)) then
                read(name(8:index(name, ')') - 1), *, iostat=ios) rec%group(k)
                if (ios /= 0 .or. rec%group(k) < 1) call die('bad parameter '//trim(name))
            else
                neq = neq + 1
                rec%eq_index(k) = neq
                eq_names = [character(len=64) :: eq_names, name]
            end if
        end do
        workdir = trim(output_dir)//'/work'
        call execute_command_line('mkdir -p "'//workdir//'"')
        call rec%eq%init(trim(vmec_input), workdir, eq_names)
        call rec%eq%vmec%set_input('ftol', vmec_ftol)
        call rec%eq%vmec%set_input('niter', real(vmec_niter, dp))
        rec%rule%samples_per_segment = samples
        rec%rule%gauss = gauss
        rec%nphi = plasma_nphi
        rec%ntheta = plasma_ntheta

        if (len_trim(flux) > 0) then
            call read_flux_loop_file(trim(flux), rec%loops, ierr, msg)
            if (ierr /= 0) call die('flux: '//msg)
        end if
        if (len_trim(segrog) > 0) then
            if (seg_area > 0.0_dp) then
                call read_segmented_rogowski_file(trim(segrog), rec%segs, ierr, msg, seg_area)
            else
                call read_segmented_rogowski_file(trim(segrog), rec%segs, ierr, msg)
            end if
            if (ierr /= 0) call die('segrog: '//msg)
        end if
        allocate(measured(0))
        if (len_trim(bprobes) > 0) then
            call read_bprobe_file(trim(bprobes), measured, ierr, msg, rphiz)
            if (ierr /= 0) call die('bprobes: '//msg)
        end if
        call read_consistency_points(fb)
        rec%n_diag_probes = size(measured)
        rec%n_fb = size(fb)
        rec%fb_sigma = consistency_sigma
        rec%probes = [measured, fb]

        ! coils: response matrices per unit EXTCUR, computed once
        if (len_trim(coils) > 0) then
            if (len_trim(coil_extcur) > 0) then
                call rec%coils%init(trim(coils), trim(coil_extcur))
            else
                call rec%coils%init(trim(coils))
            end if
            call rec%coils%response(rec%loops, rec%segs, rec%probes, rec%flux_m, rec%seg_m, &
                rec%probe_m, rec%rule)
            rec%extcur = rec%coils%extcur
        else
            call rec%coils%init('')
            allocate(rec%extcur(0))
            allocate(rec%flux_m(0, 0), rec%seg_m(0, 0), rec%probe_m(0, 0))
            if (allocated(rec%loops)) then
                deallocate(rec%flux_m)
                allocate(rec%flux_m(size(rec%loops), 0))
            end if
            if (allocated(rec%segs)) then
                deallocate(rec%seg_m)
                allocate(rec%seg_m(size(rec%segs), 0))
            end if
            deallocate(rec%probe_m)
            allocate(rec%probe_m(size(rec%probes), 0))
        end if
        do k = 1, npar
            if (rec%is_coil(k) .and. rec%group(k) > size(rec%extcur)) then
                call die('parameter '//trim(rec%names(k))//': no such coil group')
            end if
        end do
        call rec%plasma%init('')

        if (synthesize_measurements) then
            call all_signals_measured()
        else
            if (len_trim(measurements) == 0) call die('measurements file required')
            call read_measurements(trim(measurements))
        end if
    end subroutine setup

    subroutine parameter_defaults(values)
        real(dp), intent(out) :: values(:)
        real(dp), allocatable :: eq(:)
        integer :: k

        eq = rec%eq%values()
        do k = 1, rec%n
            if (rec%is_coil(k)) then
                values(k) = rec%extcur(rec%group(k))
            else
                values(k) = eq(rec%eq_index(k))
            end if
        end do
        rec%x0 = values
    end subroutine parameter_defaults

    subroutine read_consistency_points(fb)
        !! Rows "x y z" [m]; each point gives three virtual probes (x, y, z).
        type(bprobe_t), allocatable, intent(out) :: fb(:)
        real(dp) :: p(3)
        character(len=512) :: line
        character(len=32) :: label
        integer :: unit, ios, n, axis

        allocate(fb(0))
        if (len_trim(consistency_points) == 0) return
        open(newunit=unit, file=trim(consistency_points), status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open '//trim(consistency_points))
        n = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat=ios) p
            if (ios /= 0) call die('bad consistency point: '//trim(line))
            n = n + 1
            do axis = 1, 3
                write(label, '(A,I4.4,A)') 'FB_', n, '_'//'XYZ'(axis:axis)
                fb = [fb, bprobe_t(label=trim(label), position=p, normal=unit_vector(axis), &
                    eff_area=1.0_dp, turn_scale=1.0_dp)]
            end do
        end do
        close(unit)
    end subroutine read_consistency_points

    function unit_vector(axis) result(e)
        integer, intent(in) :: axis
        real(dp) :: e(3)
        e = 0.0_dp
        e(axis) = 1.0_dp
    end function unit_vector

    subroutine all_signals_measured()
        integer :: i, m

        m = 0
        if (allocated(rec%loops)) m = m + size(rec%loops)
        if (allocated(rec%segs)) m = m + size(rec%segs)
        m = m + rec%n_diag_probes
        allocate(rec%meas_kind(m), rec%meas_index(m), rec%meas_value(m), rec%meas_sigma(m), &
            rec%meas_label(m))
        m = 0
        if (allocated(rec%loops)) then
            do i = 1, size(rec%loops)
                m = m + 1
                call set_meas(m, KIND_FLUX, i, rec%loops(i)%label)
            end do
        end if
        if (allocated(rec%segs)) then
            do i = 1, size(rec%segs)
                m = m + 1
                call set_meas(m, KIND_SEGROG, i, rec%segs(i)%label)
            end do
        end if
        do i = 1, rec%n_diag_probes
            m = m + 1
            call set_meas(m, KIND_BPROBE, i, rec%probes(i)%label)
        end do
        rec%n_meas = m
        rec%meas_value = 0.0_dp
        rec%meas_sigma = 1.0_dp
    end subroutine all_signals_measured

    subroutine set_meas(m, kind, index, label)
        integer, intent(in) :: m, kind, index
        character(len=*), intent(in) :: label
        rec%meas_kind(m) = kind
        rec%meas_index(m) = index
        rec%meas_label(m) = label
    end subroutine set_meas

    subroutine read_measurements(path)
        !! CSV rows kind,label,value,sigma with kind flux, segrog or bprobe.
        character(len=*), intent(in) :: path
        character(len=512) :: line
        character(len=16) :: kind
        character(len=64) :: label
        real(dp) :: value, sigma
        integer :: unit, ios, c1, c2, c3, idx, kd

        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) call die('cannot open '//path)
        allocate(rec%meas_kind(0), rec%meas_index(0), rec%meas_value(0), rec%meas_sigma(0), &
            rec%meas_label(0))
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#' .or. line(1:5) == 'kind,') cycle
            c1 = index(line, ',')
            c2 = c1 + index(line(c1 + 1:), ',')
            c3 = c2 + index(line(c2 + 1:), ',')
            kind = line(:c1 - 1)
            label = line(c1 + 1:c2 - 1)
            read(line(c2 + 1:c3 - 1), *, iostat=ios) value
            if (ios == 0) read(line(c3 + 1:), *, iostat=ios) sigma
            if (ios /= 0 .or. sigma <= 0.0_dp) call die('bad measurement: '//trim(line))
            call locate(trim(kind), trim(label), kd, idx)
            rec%meas_kind = [rec%meas_kind, kd]
            rec%meas_index = [rec%meas_index, idx]
            rec%meas_value = [rec%meas_value, value]
            rec%meas_sigma = [rec%meas_sigma, sigma]
            rec%meas_label = [character(len=64) :: rec%meas_label, label]
        end do
        close(unit)
        rec%n_meas = size(rec%meas_kind)
        if (rec%n_meas == 0) call die('no measurements in '//path)
    end subroutine read_measurements

    subroutine locate(kind, label, kd, idx)
        character(len=*), intent(in) :: kind, label
        integer, intent(out) :: kd, idx
        integer :: i

        idx = 0
        select case (kind)
        case ('flux')
            kd = KIND_FLUX
            if (allocated(rec%loops)) then
                do i = 1, size(rec%loops)
                    if (rec%loops(i)%label == label) idx = i
                end do
            end if
        case ('segrog')
            kd = KIND_SEGROG
            if (allocated(rec%segs)) then
                do i = 1, size(rec%segs)
                    if (rec%segs(i)%label == label) idx = i
                end do
            end if
        case ('bprobe')
            kd = KIND_BPROBE
            do i = 1, rec%n_diag_probes
                if (rec%probes(i)%label == label) idx = i
            end do
        case default
            call die('unknown measurement kind '//kind)
        end select
        if (idx == 0) call die('measurement of unknown '//kind//' '//label)
    end subroutine locate

    subroutine write_measurements(path)
        character(len=*), intent(in) :: path
        integer :: unit, i
        character(len=8), parameter :: kinds(3) = [character(len=8) :: 'flux', 'segrog', 'bprobe']

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') 'kind,label,value,sigma'
        do i = 1, rec%n_meas
            write(unit, '(A,",",A,2(",",ES24.16))') trim(kinds(rec%meas_kind(i))), &
                trim(rec%meas_label(i)), rec%meas_value(i), rec%meas_sigma(i)
        end do
        close(unit)
    end subroutine write_measurements

    subroutine jacobian_check(x)
        !! Adjoint Jacobian vs central differences of re-solved equilibria.
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: j_ad(:, :), j_fd(:, :), rp(:), rm(:), xp(:)
        real(dp) :: h, err
        integer :: k, unit

        call jacobian_at(x, j_ad)
        allocate(j_fd(n_rows(), rec%n), rp(n_rows()), rm(n_rows()))
        open(newunit=unit, file=trim(output_dir)//'/jacobian_check.csv', status='replace')
        write(unit, '(A)') 'parameter,step,max_abs_error_over_max_abs'
        do k = 1, rec%n
            h = fd_step * rec%scale(k)
            xp = x
            xp(k) = x(k) + h
            call residual_at(xp, rp)
            xp(k) = x(k) - h
            call residual_at(xp, rm)
            j_fd(:, k) = (rp - rm) / (2.0_dp * h)
            err = maxval(abs(j_ad(:, k) - j_fd(:, k))) / max(maxval(abs(j_fd(:, k))), tiny(1.0_dp))
            write(unit, '(A,",",ES12.4,",",ES12.4)') trim(rec%names(k)), h, err
            print '(A,A20,A,ES10.3)', 'Jacobian check ', trim(rec%names(k)), &
                ': max |adjoint - FD| / max |FD| = ', err
        end do
        close(unit)
        call plot_jacobian_check(trim(output_dir)//'/jacobian_check.png', &
            reshape(j_ad, [size(j_ad)]), reshape(j_fd, [size(j_fd)]))
    end subroutine jacobian_check

    subroutine write_results()
        integer :: unit, k, i, dof
        real(dp) :: chi2
        character(len=8), parameter :: kinds(3) = [character(len=8) :: 'flux', 'segrog', 'bprobe']

        chi2 = sum(r_fit**2)
        dof = max(n_rows() - rec%n, 1)
        open(newunit=unit, file=trim(output_dir)//'/parameters.csv', status='replace')
        write(unit, '(A)') 'name,start,fitted,sigma,truth,prior,prior_sigma'
        do k = 1, rec%n
            write(unit, '(A,6(",",ES24.16))') '"'//trim(rec%names(k))//'"', x_start(k), x(k), &
                sqrt(max(diag(k), 0.0_dp)), merge(x_truth(k), 0.0_dp, have_truth), rec%prior(k), &
                rec%prior_sigma(k)
        end do
        close(unit)
        open(newunit=unit, file=trim(output_dir)//'/signals.csv', status='replace')
        write(unit, '(A)') 'kind,label,measured,sigma,initial,fitted,normalized_residual'
        do i = 1, rec%n_meas
            write(unit, '(A,",",A,5(",",ES24.16))') trim(kinds(rec%meas_kind(i))), &
                trim(rec%meas_label(i)), rec%meas_value(i), rec%meas_sigma(i), &
                signal_value(rec%meas_kind(i), rec%meas_index(i), f0, s0, p0), &
                signal_value(rec%meas_kind(i), rec%meas_index(i), f1, s1, p1), r_fit(i)
        end do
        close(unit)
        open(newunit=unit, file=trim(output_dir)//'/covariance.csv', status='replace')
        do k = 1, rec%n
            write(unit, '(*(ES24.16,:,","))') cov(k, :)
        end do
        close(unit)
        open(newunit=unit, file=trim(output_dir)//'/history.csv', status='replace')
        write(unit, '(A)') 'evaluation,chi2'
        do i = 1, size(rec%history)
            write(unit, '(I0,",",ES24.16)') i, rec%history(i)
        end do
        close(unit)

        open(newunit=unit, file=trim(output_dir)//'/summary.txt', status='replace')
        write(unit, '(A)') 'Tiago equilibrium reconstruction'
        write(unit, '(A,I0,A,I0,A,I0,A,I0)') 'residuals: ', n_rows(), ' (measurements ', &
            rec%n_meas, ', consistency ', rec%n_fb, ', priors ', rec%n_prior
        write(unit, '(A,ES12.4,A,ES12.4,A,I0,A,F8.3)') 'chi^2 start ', sum(r_start**2), &
            '  fitted ', chi2, '  dof ', dof, '  chi^2/dof ', chi2 / dof
        write(unit, '(A,I0,A,I0,A,I0,A,I0)') 'LM steps ', steps, ', residual evaluations ', &
            rec%evaluations, ', VMEC++ solves ', rec%eq%solves, &
            ', adjoint right-hand sides ', rec%eq%adjoint_solves
        write(unit, '(A,F9.1,A,F9.1,A,F9.1,A)') 'time: VMEC++ solves ', rec%eq%t_solve, &
            ' s, Jacobians ', rec%eq%t_adjoint, ' s, total ', t_end - t_start, ' s'
        if (len_trim(message) > 0) write(unit, '(A)') 'optimizer: '//message
        if (info /= 0) write(unit, '(A)') 'WARNING: J^T J not positive definite; no covariance'
        write(unit, '(A)') ''
        write(unit, '(A12,5A16)') 'parameter', 'start', 'fitted', 'sigma', 'truth', &
            '(fit-truth)/sig'
        do k = 1, rec%n
            if (have_truth) then
                write(unit, '(A12,4ES16.7,F16.2)') trim(rec%names(k)), x_start(k), x(k), &
                    sqrt(max(diag(k), 0.0_dp)), x_truth(k), &
                    (x(k) - x_truth(k)) / max(sqrt(max(diag(k), 0.0_dp)), tiny(1.0_dp))
            else
                write(unit, '(A12,3ES16.7)') trim(rec%names(k)), x_start(k), x(k), &
                    sqrt(max(diag(k), 0.0_dp))
            end if
        end do
        close(unit)
        call execute_command_line('cat "'//trim(output_dir)//'/summary.txt"')
    end subroutine write_results

    real(dp) function diag(k)
        integer, intent(in) :: k
        diag = 0.0_dp
        if (info == 0) diag = cov(k, k)
    end function diag

    subroutine make_plots()
        type(equilibrium_profiles_t) :: p_initial, p_fit, p_truth
        real(dp), allocatable :: meas(:), sig(:), ini(:), fitv(:), band(:)
        integer, allocatable :: idx(:)
        integer :: kind
        character(len=8), parameter :: kinds(3) = [character(len=8) :: 'flux', 'segrog', 'bprobe']
        character(len=24), parameter :: units(3) = [character(len=24) :: 'flux [Wb]', &
            'signal [T m^3]', 'eff_area B.n [T m^2]']
        integer :: i

        do kind = 1, 3
            meas = pack(rec%meas_value, rec%meas_kind == kind)
            sig = pack(rec%meas_sigma, rec%meas_kind == kind)
            if (size(meas) == 0) cycle
            idx = pack(rec%meas_index(:rec%n_meas), rec%meas_kind(:rec%n_meas) == kind)
            ini = [(signal_value(kind, idx(i), f0, s0, p0), i = 1, size(idx))]
            fitv = [(signal_value(kind, idx(i), f1, s1, p1), i = 1, size(idx))]
            call plot_signals(trim(output_dir)//'/signals_'//trim(kinds(kind))//'.png', &
                'Measured vs. modelled: '//trim(kinds(kind)), trim(units(kind)), meas, sig, ini, fitv)
        end do
        call plot_residuals(trim(output_dir)//'/residuals.png', &
            r_start(:rec%n_meas + rec%n_fb), r_fit(:rec%n_meas + rec%n_fb))
        call plot_convergence(trim(output_dir)//'/convergence.png', rec%history)
        call read_profiles(trim(output_dir)//'/wout_initial.nc', p_initial)
        call read_profiles(trim(output_dir)//'/wout_fit.nc', p_fit)
        if (synthesize_measurements) call read_profiles(trim(output_dir)//'/wout_truth.nc', p_truth)
        if (pressure_band(p_fit%s, band)) then
            call plot_profiles(trim(output_dir)//'/profiles.png', p_initial, p_fit, p_truth, band)
        else
            call plot_profiles(trim(output_dir)//'/profiles.png', p_initial, p_fit, p_truth)
        end if
        call plot_boundary(trim(output_dir)//'/boundary.png', p_initial, p_fit, p_truth)
    end subroutine make_plots

    logical function pressure_band(s, sigma) result(ok)
        !! 1-sigma band of p(s) = pres_scale sum_i AM(i) s^i from the covariance
        !! of the AM(i) and PRES_SCALE fit parameters (power-series pressure).
        real(dp), intent(in) :: s(:)
        real(dp), allocatable, intent(out) :: sigma(:)
        character(len=64) :: names(12)
        real(dp), allocatable :: coeff(:), g(:, :)
        real(dp) :: pres_scale
        integer :: i, k, j

        ok = .false.
        if (info /= 0) return
        if (.not. any([(index(upper(rec%names(k)), 'AM(') == 1 .or. &
            trim(upper(rec%names(k))) == 'PRES_SCALE', k = 1, rec%n)])) return
        names(1) = 'PRES_SCALE'
        do i = 0, 10
            write(names(i + 2), '(A,I0,A)') 'AM(', i, ')'
        end do
        coeff = [rec%eq%value_of('pres_scale'), (rec%eq%value_of('am', i + 1), i = 0, 10)]
        do k = 1, rec%n
            do j = 1, size(names)
                if (trim(upper(rec%names(k))) == trim(names(j))) coeff(j) = x(k)
            end do
        end do
        pres_scale = coeff(1)
        allocate(g(size(s), rec%n))
        g = 0.0_dp
        do k = 1, rec%n
            if (trim(upper(rec%names(k))) == 'PRES_SCALE') then
                g(:, k) = [(sum([(coeff(j + 2) * s(i)**j, j = 0, 10)]), i = 1, size(s))]
            else if (index(upper(rec%names(k)), 'AM(') == 1) then
                read(rec%names(k)(4:index(rec%names(k), ')') - 1), *) j
                g(:, k) = pres_scale * s**j
            end if
        end do
        allocate(sigma(size(s)))
        do i = 1, size(s)
            sigma(i) = sqrt(max(dot_product(g(i, :), matmul(cov, g(i, :))), 0.0_dp))
        end do
        ok = .true.
    end function pressure_band



    pure function upper(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i
        out = text
        do i = 1, len(text)
            if (text(i:i) >= 'a' .and. text(i:i) <= 'z') out(i:i) = achar(iachar(text(i:i)) - 32)
        end do
    end function upper

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(error_unit, '(A)') trim(msg)
        stop 1
    end subroutine die
end program tiago_reconstruct
