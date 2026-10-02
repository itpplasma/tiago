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
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
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
        ! Plasma responses (raw, before DIAGNO post-processing and turns) at the
        ! boundary geometry resp_geometry: the plasma part of the signals is
        ! exactly linear in the edge field y_B at fixed boundary shape, so while
        ! the boundary is unchanged S_plasma = R y_B replaces the sheet sums.
        real(dp), allocatable :: resp_geometry(:), resp_field(:)
        real(dp), allocatable :: r_flux(:, :), r_seg(:, :), r_probe(:, :)
        real(dp), allocatable :: s_flux(:, :), s_seg(:, :), s_probe(:, :)
    end type reconstruction_t

    type(reconstruction_t), public, target, save :: rec

    ! the last residual and Jacobian with their parameters
    real(dp), allocatable, save :: cache_rx(:), cache_r(:), cache_jx(:), cache_j(:, :)

    public :: n_rows, parameters_of, model_signals, residual_at, jacobian_at
    public :: fit, covariance, synthesize, signal_value

contains

    integer function n_rows()
        n_rows = rec%n_meas + rec%n_fb + rec%n_prior
    end function n_rows

    function parameters_of(z) result(x)
        !! Parameters from scaled parameters z = x / scale.
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
        real(dp), allocatable :: extcur(:), pf(:), ps(:), pp(:), yb(:)
        logical :: solved

        extcur = extcur_of(x)
        call use_equilibrium(x, solved)
        if (present(ok)) ok = solved
        if (.not. solved) then
            if (.not. present(ok)) error stop 'VMEC++ found no equilibrium'
            return
        end if
        allocate(flux(0), seg(0), probe(0))
        if (responses_current()) then
            yb = rec%eq%y(:size(rec%r_flux, 2))
            if (allocated(rec%loops)) pf = matmul(rec%r_flux, yb)
            if (allocated(rec%segs)) ps = matmul(rec%r_seg, yb)
            if (allocated(rec%probes)) pp = matmul(rec%r_probe, yb)
        else
            if (allocated(rec%loops)) call rec%plasma%flux_loops(rec%loops, pf, rec%rule)
            if (allocated(rec%segs)) call rec%plasma%segrog(rec%segs, ps, rec%rule)
            if (allocated(rec%probes)) call rec%plasma%bprobes(rec%probes, pp)
        end if
        ! (no matmul with a zero inner extent: gfortran's inlined matmul sizes
        ! the sum from it and overruns the result)
        if (allocated(rec%loops)) then
            flux = pf
            if (size(extcur) > 0) flux = flux + matmul(rec%flux_m, extcur)
            call finalize_flux_signals(rec%loops, flux, rec%plasma%plasma%diamagnetic_flux)
        end if
        if (allocated(rec%segs)) then
            seg = ps
            if (size(extcur) > 0) seg = seg + matmul(rec%seg_m, extcur)
            seg = seg * rec%segs%turn_scale
        end if
        if (allocated(rec%probes)) then
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
        !! Whitened residuals at x; a repeated x returns the last result without
        !! counting a new evaluation.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: r(:)
        real(dp), allocatable :: flux(:), seg(:), probe(:)
        integer :: i, k
        logical :: ok

        if (allocated(cache_rx)) then
            if (size(cache_rx) == size(x)) then
                if (all(cache_rx == x)) then
                    r = cache_r
                    return
                end if
            end if
        end if
        call model_signals(x, flux, seg, probe, ok)
        rec%evaluations = rec%evaluations + 1
        if (.not. ok) then
            ! the fit rejects a step to non-finite residuals
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
        cache_rx = x
        cache_r = r
        write(*, '(A,I4,A,ES12.4,A,*(1X,A,"=",ES13.6))') 'evaluation ', rec%evaluations, &
            ': chi^2 = ', sum(r**2), ' ', (trim(rec%names(k)), x(k), k = 1, rec%n)
    end subroutine residual_at

    logical function responses_current(with_shape)
        !! Field responses depend on geometry; shape responses also need the
        !! boundary field at which their geometry derivative was evaluated.
        logical, intent(in), optional :: with_shape
        integer :: nb

        responses_current = .false.
        if (.not. allocated(rec%resp_geometry)) return
        if (size(rec%resp_geometry) /= size(boundary_geometry())) return
        if (.not. all(rec%resp_geometry == boundary_geometry())) return
        if (present(with_shape)) then
            if (with_shape) then
                if (.not. allocated(rec%resp_field)) return
                nb = 2 * rec%eq%edge%mnmax_nyq
                if (size(rec%resp_field) /= nb) return
                if (.not. all(rec%resp_field == rec%eq%y(:nb))) return
            end if
        end if
        responses_current = .true.
    end function responses_current

    function boundary_geometry() result(g)
        !! Boundary rmnc, zmns of the current equilibrium (the edge part of y).
        real(dp), allocatable :: g(:)
        g = rec%eq%y(2 * rec%eq%edge%mnmax_nyq + 1:)
    end function boundary_geometry

    subroutine update_responses()
        !! Plasma field and shape responses at the current equilibrium's boundary.
        !! Without RBC/ZBS parameters the boundary is that of the VMEC input, so
        !! the shape responses do not enter the Jacobian and are left zero.
        integer :: k, ng
        logical :: shape

        shape = .false.
        do k = 1, rec%n
            shape = shape .or. index(upper(rec%names(k)), 'RBC(') == 1 .or. &
                index(upper(rec%names(k)), 'ZBS(') == 1
        end do
        if (responses_current(shape)) return
        if (shape) then
            call rec%plasma%plasma_response(rec%loops, rec%segs, rec%probes, rec%r_flux, &
                rec%r_seg, rec%r_probe, rec%rule, rec%s_flux, rec%s_seg, rec%s_probe)
        else
            call rec%plasma%plasma_response(rec%loops, rec%segs, rec%probes, rec%r_flux, &
                rec%r_seg, rec%r_probe, rec%rule)
            ng = rec%plasma%plasma%n_shape_columns()
            if (allocated(rec%s_flux)) deallocate(rec%s_flux)
            if (allocated(rec%s_seg)) deallocate(rec%s_seg)
            if (allocated(rec%s_probe)) deallocate(rec%s_probe)
            allocate(rec%s_flux(size(rec%r_flux, 1), ng), rec%s_seg(size(rec%r_seg, 1), ng), &
                rec%s_probe(size(rec%r_probe, 1), ng))
            rec%s_flux = 0.0_dp
            rec%s_seg = 0.0_dp
            rec%s_probe = 0.0_dp
        end if
        rec%resp_geometry = boundary_geometry()
        rec%resp_field = rec%eq%y(:2 * rec%eq%edge%mnmax_nyq)
    end subroutine update_responses

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
            call update_responses()
            fr = rec%r_flux
            sr = rec%r_seg
            pr = rec%r_probe
            fs = rec%s_flux
            ss = rec%s_seg
            prs = rec%s_probe
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

    subroutine jacobian_cached(x, jac)
        !! jacobian_at, reusing the last result at the same x (the fit ends at an
        !! accepted point whose Jacobian the covariance needs again).
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        if (allocated(cache_jx)) then
            if (size(cache_jx) == size(x)) then
                if (all(cache_jx == x)) then
                    jac = cache_j
                    return
                end if
            end if
        end if
        call jacobian_at(x, jac)
        cache_jx = x
        cache_j = jac
    end subroutine jacobian_cached

    subroutine fit(x, max_iterations, steps, message)
        !! Levenberg-Marquardt (Madsen, Nielsen and Tingleff 2004, with Marquardt
        !! scaling) from x, updated in place, on the scaled parameters x / scale.
        !! A trial point costs a VMEC++ solve, so the fit stops as soon as the
        !! Gauss-Newton model predicts a chi^2 decrease below the reproducibility
        !! of chi^2 (~1e-3: equilibria are solved to FTOL, hot-restarted), rather
        !! than probing the noise floor with further solves.
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: max_iterations
        integer, intent(out) :: steps
        character(len=:), allocatable, intent(out) :: message
        real(dp), parameter :: floor = 1.0e-3_dp
        real(dp), allocatable :: r(:), r_trial(:), jac(:, :), a(:, :), m(:, :), g(:), h(:), d(:)
        real(dp), allocatable :: x_trial(:)
        real(dp) :: chi2, chi2_trial, predicted, rho, mu, nu
        integer :: n, i, info
        external :: dposv

        n = size(x)
        allocate(r(n_rows()), r_trial(n_rows()))
        call residual_at(x, r)
        if (.not. all(ieee_is_finite(r))) then
            message = 'no equilibrium at the start values'
            steps = 0
            return
        end if
        chi2 = sum(r**2)
        call jacobian_cached(x, jac)
        jac = jac * spread(rec%scale, 1, size(jac, 1))
        a = matmul(transpose(jac), jac)
        g = matmul(transpose(jac), r)
        mu = 1.0e-3_dp
        nu = 2.0_dp
        steps = 0
        message = 'iteration limit reached'
        do while (steps < max_iterations)
            ! (A + mu diag(A)) h = -g
            d = [(max(a(i, i), tiny(1.0_dp)), i = 1, n)]
            m = a
            do i = 1, n
                m(i, i) = m(i, i) + mu * d(i)
            end do
            h = -g
            call dposv('U', n, 1, m, n, h, n, info)
            if (info /= 0) then
                mu = mu * nu
                nu = 2.0_dp * nu
                cycle
            end if
            ! chi2 decrease predicted by the linear model
            predicted = -2.0_dp * dot_product(h, g) - dot_product(h, matmul(a, h))
            if (predicted <= floor) then
                message = 'converged (predicted chi^2 decrease below the solver noise floor)'
                exit
            end if
            if (maxval(abs(h)) <= 1.0e-8_dp * (maxval(abs(x / rec%scale)) + 1.0e-8_dp)) then
                message = 'converged (step below 1e-8 of the parameters)'
                exit
            end if
            x_trial = x + h * rec%scale
            call residual_at(x_trial, r_trial)
            chi2_trial = sum(r_trial**2)
            rho = -1.0_dp
            if (ieee_is_finite(chi2_trial)) rho = (chi2 - chi2_trial) / predicted
            if (rho > 0.0_dp) then
                x = x_trial
                r = r_trial
                chi2 = chi2_trial
                steps = steps + 1
                call jacobian_cached(x, jac)
                jac = jac * spread(rec%scale, 1, size(jac, 1))
                a = matmul(transpose(jac), jac)
                g = matmul(transpose(jac), r)
                mu = mu * max(1.0_dp / 3.0_dp, 1.0_dp - (2.0_dp * rho - 1.0_dp)**3)
                nu = 2.0_dp
            else
                mu = mu * nu
                nu = 2.0_dp * nu
                if (mu > 1.0e12_dp) then
                    message = 'converged (no further decrease of chi^2)'
                    exit
                end if
            end if
        end do
    end subroutine fit

    subroutine covariance(x, cov, jac, info)
        !! Linearized posterior covariance (J^T J)^-1 of the whitened problem.
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: cov(:, :), jac(:, :)
        integer, intent(out) :: info
        external :: dpotrf, dpotri
        integer :: i, j

        call jacobian_cached(x, jac)
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
