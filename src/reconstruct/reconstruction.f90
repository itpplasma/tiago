module tiago_reconstruction
    !! Equilibrium reconstruction: find parameters x (VMEC profile and boundary
    !! inputs, PHIEDGE, coil-group currents) such that Tiago's signals match the
    !! measurements, by Levenberg-Marquardt on the whitened residual
    !!
    !!     r = [ (S(x) - m) / sigma ;  (B_coil + B_sheet)(p_j) / sigma_fb ;
    !!           (x - x_prior) / sigma_prior ].
    !!
    !! S = coil part (linear in EXTCUR, response matrices computed once) plus the
    !! plasma part of a fixed-boundary VMEC++ equilibrium. The second block is
    !! the free-boundary condition: inside the plasma the boundary sheet current
    !! produces minus the external field, so B_coil + B_sheet = 0 at interior
    !! points p_j exactly when the equilibrium is consistent with the coils.
    !!
    !! Jacobian: coil columns from the response matrices; equilibrium columns
    !! J = C (dy/dx) with C = dS/dy from Tiago (boundary field and shape
    !! Jacobians) and dy/dx from VMEC++'s implicit adjoint (tiago_equilibrium); the
    !! idia = 1 phiedge term is added analytically.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i32 => int32, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, FORTNUM_CONVERGENCE_ERROR
    use fortopt_least_squares, only: least_squares_t, levenberg_marquardt_t, &
        least_squares_options_t, least_squares_result_t
    use tiago_diagnostic_types, only: flux_loop_t, segmented_rogowski_t, bprobe_t
    use tiago_flux_loops, only: finalize_flux_signals
    use tiago_vacuum_forward, only: vacuum_solver_t, quadrature_rule_t
    use tiago_equilibrium, only: equilibrium_t
    use tiago_plasma_support, only: vmec_boundary_t
    implicit none
    private

    integer, parameter, public :: KIND_FLUX = 1, KIND_SEGROG = 2, KIND_BPROBE = 3

    type, public :: reconstruction_t
        ! parameters
        integer :: n = 0
        character(len=64), allocatable :: names(:)
        logical, allocatable :: is_coil(:)          !! EXTCUR(g) parameter
        integer, allocatable :: group(:)            !! coil group of a coil parameter
        integer, allocatable :: eq_index(:)         !! position among the equilibrium parameters
        real(dp), allocatable :: x0(:), scale(:), prior(:), prior_sigma(:)
        ! diagnostics and measurements
        type(flux_loop_t), allocatable :: loops(:)
        type(segmented_rogowski_t), allocatable :: segs(:)
        type(bprobe_t), allocatable :: probes(:)    !! measured probes, then consistency probes
        integer :: n_diag_probes = 0
        integer :: n_meas = 0, n_fb = 0, n_prior = 0
        integer, allocatable :: meas_kind(:), meas_index(:)
        real(dp), allocatable :: meas_value(:), meas_sigma(:)
        character(len=64), allocatable :: meas_label(:)
        real(dp) :: fb_sigma = 1.0_dp
        integer, allocatable :: prior_param(:)
        ! models
        type(vacuum_solver_t) :: coils                   !! coils (response matrices)
        type(vacuum_solver_t) :: plasma                  !! plasma-only solver
        type(quadrature_rule_t) :: rule
        integer :: nphi = 32, ntheta = 32
        real(dp), allocatable :: extcur(:)               !! all coil groups
        real(dp), allocatable :: flux_m(:, :), seg_m(:, :), probe_m(:, :)
        type(equilibrium_t) :: eq
        ! bookkeeping
        integer :: evaluations = 0
        real(dp), allocatable :: history(:)              !! chi^2 of every evaluation
        real(dp), allocatable :: eq_x(:)                 !! its equilibrium parameters
    end type reconstruction_t

    type(reconstruction_t), public, target, save :: rec

    ! cache of the last residual and Jacobian (fortopt callbacks carry no context)
    real(dp), allocatable, save :: cache_x(:), cache_r(:), cache_jx(:), cache_j(:, :)

    public :: n_rows, parameters_of, model_signals, residual_at, jacobian_at
    public :: fit, covariance, synthesize, signal_value

contains

    integer function n_rows()
        n_rows = rec%n_meas + rec%n_fb + rec%n_prior
    end function n_rows

    function parameters_of(z) result(x)
        !! fortopt works on scaled parameters z = x / scale.
        real(dp), intent(in) :: z(:)
        real(dp) :: x(size(z))
        x = z * rec%scale
    end function parameters_of

    function eq_values(x) result(values)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: values(:)
        integer :: k

        allocate(values(size(rec%eq%names)))
        do k = 1, rec%n
            if (.not. rec%is_coil(k)) values(rec%eq_index(k)) = x(k)
        end do
    end function eq_values

    function extcur_of(x) result(extcur)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: extcur(:)
        integer :: k

        extcur = rec%extcur
        do k = 1, rec%n
            if (rec%is_coil(k)) extcur(rec%group(k)) = x(k)
        end do
    end function extcur_of

    subroutine model_signals(x, flux, seg, probe, ok)
        !! All signals at x: flux loops (DIAGNO post-processing), Rogowskis and
        !! probes (measured probes, then consistency probes), after turns.
        !! ok = .false. if VMEC++ found no equilibrium (e.g. invalid parameters).
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: flux(:), seg(:), probe(:)
        logical, intent(out), optional :: ok
        real(dp), allocatable :: extcur(:), pf(:), ps(:), pp(:)
        logical :: solved

        extcur = extcur_of(x)
        call use_equilibrium(x, solved)
        if (present(ok)) ok = solved
        if (.not. solved) then
            if (.not. present(ok)) error stop 'VMEC++ found no equilibrium'
            return
        end if
        allocate(flux(0), seg(0), probe(0))
        ! (no matmul with a zero inner extent: gfortran's inlined matmul sizes
        ! the sum from it and overruns the result)
        if (allocated(rec%loops)) then
            call rec%plasma%flux_loops(rec%loops, pf, rec%rule)
            flux = pf
            if (size(extcur) > 0) flux = flux + matmul(rec%flux_m, extcur)
            call finalize_flux_signals(rec%loops, flux, rec%plasma%plasma%diamagnetic_flux)
        end if
        if (allocated(rec%segs)) then
            call rec%plasma%segrog(rec%segs, ps, rec%rule)
            seg = ps
            if (size(extcur) > 0) seg = seg + matmul(rec%seg_m, extcur)
            seg = seg * rec%segs%turn_scale
        end if
        if (allocated(rec%probes)) then
            call rec%plasma%bprobes(rec%probes, pp)
            probe = pp
            if (size(extcur) > 0) probe = probe + matmul(rec%probe_m, extcur)
            probe = probe * rec%probes%turn_scale
        end if
    end subroutine model_signals

    subroutine use_equilibrium(x, solved)
        !! Solve the equilibrium for x, unless the last solve was for the same
        !! equilibrium parameters, and load it into the plasma model.
        real(dp), intent(in) :: x(:)
        logical, intent(out), optional :: solved
        real(dp), allocatable :: values(:)
        type(vmec_boundary_t) :: vb
        logical :: ok

        values = eq_values(x)
        if (present(solved)) solved = .true.
        if (allocated(rec%eq_x)) then
            if (size(rec%eq_x) == size(values)) then
                if (all(rec%eq_x == values) .and. rec%eq%solved) return
            end if
        end if
        call rec%eq%solve(values, ok)
        if (.not. ok) then
            if (present(solved)) solved = .false.
            if (.not. present(solved)) error stop 'VMEC++ found no equilibrium'
            if (allocated(rec%eq_x)) deallocate(rec%eq_x)
            return
        end if
        call rec%eq%boundary(vb)
        call rec%plasma%enable_plasma_from_boundary(vb, int(rec%nphi, i32), int(rec%ntheta, i32))
        rec%eq_x = values
    end subroutine use_equilibrium

    real(dp) function signal_value(kind, index, flux, seg, probe)
        integer, intent(in) :: kind, index
        real(dp), intent(in) :: flux(:), seg(:), probe(:)
        select case (kind)
        case (KIND_FLUX)
            signal_value = flux(index)
        case (KIND_SEGROG)
            signal_value = seg(index)
        case default
            signal_value = probe(index)
        end select
    end function signal_value

    subroutine residual_at(x, r)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:)
        real(dp), allocatable :: flux(:), seg(:), probe(:)
        integer :: i, k
        logical :: ok

        call model_signals(x, flux, seg, probe, ok)
        rec%evaluations = rec%evaluations + 1
        if (.not. ok) then
            ! the line search backtracks from non-finite residuals
            r = ieee_value(1.0_dp, ieee_quiet_nan)
            write(*, '(A,I4,A)') 'evaluation ', rec%evaluations, ': no equilibrium, step rejected'
            return
        end if
        do i = 1, rec%n_meas
            r(i) = (signal_value(rec%meas_kind(i), rec%meas_index(i), flux, seg, probe) &
                - rec%meas_value(i)) / rec%meas_sigma(i)
        end do
        do i = 1, rec%n_fb
            r(rec%n_meas + i) = probe(rec%n_diag_probes + i) / rec%fb_sigma
        end do
        do i = 1, rec%n_prior
            k = rec%prior_param(i)
            r(rec%n_meas + rec%n_fb + i) = (x(k) - rec%prior(k)) / rec%prior_sigma(k)
        end do
        rec%history = [rec%history, sum(r**2)]
        write(*, '(A,I4,A,ES12.4,A,*(1X,A,"=",ES13.6))') 'evaluation ', rec%evaluations, &
            ': chi^2 = ', sum(r**2), ' ', (trim(rec%names(k)), x(k), k = 1, rec%n)
    end subroutine residual_at

    subroutine jacobian_at(x, jac)
        !! d r / d x (unscaled parameters), rows as in residual_at.
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: fr(:, :), sr(:, :), pr(:, :), fs(:, :), ss(:, :), prs(:, :)
        real(dp), allocatable :: cot(:, :), jeq(:, :), dphi(:), col(:), row(:)
        integer :: i, k, nb, ny, c, neq

        allocate(jac(n_rows(), rec%n))
        jac = 0.0_dp
        neq = size(rec%eq%names)
        if (neq > 0) then
            ! C = dS/dy at the current equilibrium (field then shape columns)
            call use_equilibrium(x)
            call rec%plasma%plasma_response(rec%loops, rec%segs, rec%probes, fr, sr, pr, &
                rec%rule, fs, ss, prs)
            nb = size(fr, 2)
            ny = nb + size(fs, 2)
            if (allocated(rec%loops)) then
                do c = 1, nb
                    call finalize_flux_signals(rec%loops, fr(:, c))
                end do
                do c = 1, size(fs, 2)
                    call finalize_flux_signals(rec%loops, fs(:, c))
                end do
            end if
            allocate(cot(rec%n_meas + rec%n_fb, ny))
            do i = 1, rec%n_meas
                select case (rec%meas_kind(i))
                case (KIND_FLUX)
                    row = [fr(rec%meas_index(i), :), fs(rec%meas_index(i), :)]
                case (KIND_SEGROG)
                    row = [sr(rec%meas_index(i), :), ss(rec%meas_index(i), :)] * &
                        rec%segs(rec%meas_index(i))%turn_scale
                case default
                    row = [pr(rec%meas_index(i), :), prs(rec%meas_index(i), :)] * &
                        rec%probes(rec%meas_index(i))%turn_scale
                end select
                cot(i, :) = row / rec%meas_sigma(i)
            end do
            do i = 1, rec%n_fb
                k = rec%n_diag_probes + i
                cot(rec%n_meas + i, :) = [pr(k, :), prs(k, :)] / rec%fb_sigma
            end do
            call rec%eq%jacobian(eq_values(x), cot, jeq)
            do k = 1, rec%n
                if (.not. rec%is_coil(k)) jac(:rec%n_meas + rec%n_fb, k) = jeq(:, rec%eq_index(k))
            end do
            ! idia = 1 loops: + phiedge * signgs
            do k = 1, rec%n
                if (rec%is_coil(k) .or. trim(upper(rec%names(k))) /= 'PHIEDGE') cycle
                if (.not. allocated(rec%loops)) cycle
                allocate(dphi(size(rec%loops)))
                dphi = 0.0_dp
                call finalize_flux_signals(rec%loops, dphi, &
                    sign(1.0_dp, rec%plasma%plasma%diamagnetic_flux * x(k)))
                do i = 1, rec%n_meas
                    if (rec%meas_kind(i) == KIND_FLUX) jac(i, k) = jac(i, k) + &
                        dphi(rec%meas_index(i)) / rec%meas_sigma(i)
                end do
                deallocate(dphi)
            end do
        end if
        ! coil currents: response matrices
        do k = 1, rec%n
            if (.not. rec%is_coil(k)) cycle
            if (allocated(rec%loops)) then
                col = rec%flux_m(:, rec%group(k))
                call finalize_flux_signals(rec%loops, col)
            end if
            do i = 1, rec%n_meas
                select case (rec%meas_kind(i))
                case (KIND_FLUX)
                    jac(i, k) = col(rec%meas_index(i))
                case (KIND_SEGROG)
                    jac(i, k) = rec%seg_m(rec%meas_index(i), rec%group(k)) * &
                        rec%segs(rec%meas_index(i))%turn_scale
                case default
                    jac(i, k) = rec%probe_m(rec%meas_index(i), rec%group(k)) * &
                        rec%probes(rec%meas_index(i))%turn_scale
                end select
                jac(i, k) = jac(i, k) / rec%meas_sigma(i)
            end do
            do i = 1, rec%n_fb
                jac(rec%n_meas + i, k) = rec%probe_m(rec%n_diag_probes + i, rec%group(k)) &
                    / rec%fb_sigma
            end do
        end do
        do i = 1, rec%n_prior
            k = rec%prior_param(i)
            jac(rec%n_meas + rec%n_fb + i, k) = 1.0_dp / rec%prior_sigma(k)
        end do
    end subroutine jacobian_at

    ! ---- fortopt callbacks on scaled parameters z = x / scale -------------------

    subroutine value_cb(z, residual, status)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: residual(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: x(:)

        x = parameters_of(z)
        if (allocated(cache_x)) then
            if (all(cache_x == x)) then
                residual = cache_r
                call status_set(status, FORTNUM_OK, '')
                return
            end if
        end if
        call residual_at(x, residual)
        cache_x = x
        cache_r = residual
        call status_set(status, FORTNUM_OK, '')
    end subroutine value_cb

    subroutine ensure_jacobian(z)
        real(dp), intent(in) :: z(:)
        real(dp), allocatable :: x(:), jx(:, :)
        integer :: k

        x = parameters_of(z)
        if (allocated(cache_jx)) then
            if (all(cache_jx == x)) return
        end if
        call jacobian_at(x, jx)
        do k = 1, rec%n
            jx(:, k) = jx(:, k) * rec%scale(k)
        end do
        cache_j = jx
        cache_jx = x
    end subroutine ensure_jacobian

    subroutine jvp_cb(z, direction, residual_dot, status)
        real(dp), intent(in) :: z(:), direction(:)
        real(dp), intent(out) :: residual_dot(:)
        type(fortnum_status_t), intent(out) :: status

        call ensure_jacobian(z)
        residual_dot = matmul(cache_j, direction)
        call status_set(status, FORTNUM_OK, '')
    end subroutine jvp_cb

    subroutine vjp_cb(z, residual_bar, gradient, status)
        real(dp), intent(in) :: z(:), residual_bar(:)
        real(dp), intent(out) :: gradient(:)
        type(fortnum_status_t), intent(out) :: status

        call ensure_jacobian(z)
        gradient = matmul(residual_bar, cache_j)
        call status_set(status, FORTNUM_OK, '')
    end subroutine vjp_cb

    subroutine fit(x, max_iterations, steps, message)
        !! Levenberg-Marquardt from x (updated in place).
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: max_iterations
        integer, intent(out) :: steps
        character(len=:), allocatable, intent(out) :: message
        type(least_squares_t) :: problem
        type(levenberg_marquardt_t) :: lm
        type(least_squares_options_t) :: options
        type(least_squares_result_t) :: result
        type(fortnum_status_t) :: status
        real(dp), allocatable :: z(:)

        call problem%initialize(rec%n, n_rows(), value_cb, jvp_cb, vjp_cb, status)
        options%max_iterations = max_iterations
        ! Stop at VMEC++'s noise floor: with FTOL ~ 1e-14 chi^2 is reproducible
        ! to ~1e-3, so smaller changes (and steps below 1e-4 of the parameter
        ! scales) carry no information.
        options%gradient_tolerance = 1.0e-10_dp
        options%step_tolerance = 1.0e-4_dp
        options%objective_tolerance = 1.0e-3_dp
        ! A rejected step raises the damping after one halving instead of halving
        ! a Gauss-Newton-like step many times: every trial is a VMEC++ solve.
        options%max_backtracking = 2
        options%max_damping_attempts = 6
        z = x / rec%scale
        call lm%minimize(problem, z, options, result, status)
        x = parameters_of(z)
        steps = result%accepted_steps
        message = trim(status%msg)
        if (status%code /= FORTNUM_OK) then
            if (steps > 0 .and. status%code == FORTNUM_CONVERGENCE_ERROR) then
                ! no decrease left at the noise floor of the equilibrium solves
                message = 'converged (no further decrease of chi^2)'
            else
                write(error_unit, '(A)') 'WARNING: Levenberg-Marquardt: '//message
            end if
        else if (result%state%converged) then
            message = 'converged'
        else
            message = 'iteration limit reached'
        end if
    end subroutine fit

    subroutine covariance(x, cov, jac, info)
        !! Linearized posterior covariance (J^T J)^-1 of the whitened problem.
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: cov(:, :), jac(:, :)
        integer, intent(out) :: info
        external :: dpotrf, dpotri
        integer :: i, j

        call jacobian_at(x, jac)
        cov = matmul(transpose(jac), jac)
        call dpotrf('U', rec%n, cov, rec%n, info)
        if (info /= 0) return
        call dpotri('U', rec%n, cov, rec%n, info)
        do j = 1, rec%n
            do i = j + 1, rec%n
                cov(i, j) = cov(j, i)
            end do
        end do
    end subroutine covariance

    subroutine synthesize(x_true, sigma_rel, sigma_flux, sigma_segrog, sigma_bprobe, seed, &
            values, sigmas, noise)
        !! Measurements of the model at x_true with Gaussian noise of standard
        !! deviation sigma_rel |S| + sigma_<kind> (deterministic seed); exact
        !! values (with those sigmas) if noise = .false.
        real(dp), intent(in) :: x_true(:), sigma_rel, sigma_flux, sigma_segrog, sigma_bprobe
        integer, intent(in) :: seed
        real(dp), allocatable, intent(out) :: values(:), sigmas(:)
        logical, intent(in), optional :: noise
        real(dp), allocatable :: flux(:), seg(:), probe(:)
        real(dp) :: exact, floor, u(2)
        integer, allocatable :: state(:)
        integer :: i, n

        call model_signals(x_true, flux, seg, probe)
        call random_seed(size=n)
        allocate(state(n))
        state = seed + 37 * [(i, i = 1, n)]
        call random_seed(put=state)
        allocate(values(rec%n_meas), sigmas(rec%n_meas))
        do i = 1, rec%n_meas
            exact = signal_value(rec%meas_kind(i), rec%meas_index(i), flux, seg, probe)
            select case (rec%meas_kind(i))
            case (KIND_FLUX)
                floor = sigma_flux
            case (KIND_SEGROG)
                floor = sigma_segrog
            case default
                floor = sigma_bprobe
            end select
            sigmas(i) = sigma_rel * abs(exact) + floor
            call random_number(u)
            u(1) = max(u(1), tiny(1.0_dp))
            values(i) = exact + sigmas(i) * sqrt(-2.0_dp * log(u(1))) * cos(2.0_dp * acos(-1.0_dp) * u(2))
            if (present(noise)) then
                if (.not. noise) values(i) = exact
            end if
        end do
    end subroutine synthesize

    pure function upper(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i
        out = text
        do i = 1, len(text)
            if (text(i:i) >= 'a' .and. text(i:i) <= 'z') out(i:i) = achar(iachar(text(i:i)) - 32)
        end do
    end function upper
end module tiago_reconstruction
