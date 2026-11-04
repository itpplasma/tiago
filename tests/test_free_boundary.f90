program test_free_boundary
    !! Public CTH-like free boundary: complete signal derivatives and local UQ.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tiago_reconstruction, only: rec, jacobian_at, residual_at, invalidate_caches, &
        model_signals, signal_value, n_rows, fit, covariance
    use tiago_synthetic_case, only: setup_magnetics, place_probes, &
        synthesize_magnetics, synthesize_profiles
    use tiago_joint_model, only: joint_model_t
    use tiago_vacuum_forward, only: vacuum_solver_t
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    implicit none
    character(len=512) :: fixture, awk, work, mode
    character(len=16), parameter :: names(4) = [character(len=16) :: &
        'EXTCUR(1)', 'PHIEDGE', 'CURTOR', &
        'PRES_SCALE']
    type(joint_model_t) :: model
    type(vacuum_solver_t) :: coils
    type(laplace_t) :: lap
    type(rng_t) :: rng
    type(fortnum_status_t) :: st
    real(dp), allocatable :: x(:), truth(:), prior(:), sigma(:), theta(:), theta0(:)
    real(dp), allocatable :: r(:), rp(:), rm(:), jac(:, :), fd(:), cov(:, :)
    real(dp), allocatable :: flux(:), seg(:), probe(:), saved(:), y0(:)
    real(dp) :: h, err, initial, current_scale(2), log_extra
    integer :: i, k, npar, steps, info, unit, before
    logical :: ok, joint
    character(len=:), allocatable :: message

    call get_command_argument(1, fixture)
    call get_command_argument(2, awk)
    call get_command_argument(3, work)
    call get_command_argument(4, mode)
    joint = trim(mode) /= 'magnetic'
    npar = 4
    if (joint) npar = 3
    prior = [4700.0_dp, -0.035_dp, 43229.08092460368_dp, 432.29080924603676_dp]
    sigma = [100.0_dp, 0.001_dp, 2000.0_dp, 100.0_dp]
    call execute_command_line('mkdir -p "'//trim(work)//'" && awk '// &
        '-v R0=0.78 -v rs=0.5 -v nphi=3 -v nth=3 -v nprobe=12 '// &
        '-v dir="'//trim(work)//'" -f "'//trim(awk)//'"')
    call setup_magnetics(trim(fixture)//'/input.json', trim(work)//'/work', &
        trim(work), &
        names(:npar), prior(:npar), sigma(:npar), 1.0e-18_dp)
    if (.not. rec%eq%free_boundary) error stop 'fixture is not free-boundary'
    if (index(rec%eq%derivative_method, 'finite differences') == 0) &
        error stop 'missing FD provenance'
    rec%is_coil(1) = .true.
    rec%group(1) = 1
    call coils%init(trim(fixture)//'/coils.coils')
    rec%coils = coils
    call rec%coils%response(rec%loops, rec%segs, rec%probes, rec%flux_m, &
        rec%seg_m, rec%probe_m, rec%rule)
    ! TIAGO normalizes by the first filament current; restore published turns
    ! per ampere of VMEC EXTCUR (96 for the helical group, -36 for the vertical).
    current_scale = [96.0_dp, -36.0_dp]
    rec%flux_m = rec%flux_m*spread(current_scale, 1, size(rec%flux_m, 1))
    rec%seg_m = rec%seg_m*spread(current_scale, 1, size(rec%seg_m, 1))
    rec%probe_m = rec%probe_m*spread(current_scale, 1, size(rec%probe_m, 1))
    rec%extcur = [4700.0_dp, 1000.0_dp]
    rec%quiet = .true.
    call rng_seed(rng, 71_int64, st)
    if (.not. joint) then
        x = prior
        truth = prior + 0.4_dp*sigma
        call model_signals(truth, flux, seg, probe, ok)
        if (.not. ok) error stop 'planted equilibrium did not converge'
        do i = 1, rec%n_meas
            rec%meas_value(i) = signal_value(rec%meas_kind(i), rec%meas_index(i), &
                flux, seg, probe)
            rec%meas_sigma(i) = 0.003_dp*abs(rec%meas_value(i)) + 1.0e-6_dp
        end do
        allocate (r(n_rows()), rp(n_rows()), rm(n_rows()))
        call residual_at(x, r)
        call jacobian_at(x, jac)
        y0 = rec%eq%y
        before = rec%eq%solves
        call jacobian_at(x, cov)
        if (before /= rec%eq%solves) error stop 'FD state tangents were not cached'
        if (maxval(abs(y0 - rec%eq%y)) > 1.0e-12_dp) &
            error stop 'FD cache changed base state'
        call check_columns()
        call invalidate_caches()
        call rec%eq%vmec%set_input('hot_restart', 1.0_dp)
        call residual_at(x, r)
        initial = sum(r**2)
        call fit(x, 12, steps, message)
        call residual_at(x, r)
        if (sum(r**2) >= initial) error stop 'free-boundary fit did not improve'
        call covariance(x, cov, jac, info)
        if (info /= 0) error stop 'free-boundary covariance failed'
        if (any(abs(x - truth) > 2*sqrt([(cov(k, k), k=1, npar)]))) &
            error stop 'planted free-boundary parameters outside Laplace intervals'
        print '(A,ES12.4)', 'magnetic maximum parameter pull: ', &
            maxval(abs(x - truth)/sqrt([(cov(k, k), k=1, npar)]))
        open (newunit=unit, file=trim(work)//'/recovery.csv', status='replace')
        write (unit, '(A)') 'parameter,truth,fit,sd'
        do k = 1, npar
            write (unit, '(A,3(A,ES22.14))') trim(names(k)), ',', truth(k), ',', &
                x(k), ',', sqrt(cov(k, k))
        end do
        close (unit)
    else
        call model%prof%init([0.0_dp, 1.0_dp])
        associate (g => model%prof%greville)
            call model%prof%set_profile(P_NE, 1.0e19_dp*(1 - 0.4_dp*g), .true., &
                0.1_dp, 0.6_dp, 0.3_dp, 0.5_dp, 0.2_dp)
            call model%prof%set_profile(P_TE, 100.0_dp*(1 - 0.8_dp*g), .true., &
                0.1_dp, 0.6_dp, 0.3_dp, 0.5_dp, 0.2_dp)
            call model%prof%set_profile(P_TI, 80.0_dp*(1 - 0.8_dp*g), .true., &
                0.1_dp, 0.6_dp, 0.3_dp, 0.5_dp, 0.2_dp)
        end associate
        if (trim(mode) == 'joint') &
            call model%prof%set_nonthermal([0.0_dp, 1.0_dp], [2.0_dp, 8.0_dp])
        allocate (model%diags(0))
        call model%setup()
        theta = [rec%prior/rec%scale, [(0.0_dp, k=1, model%n_prof)]]
        allocate (r(model%n_residuals()))
        call model%residual(theta, r, ok)
        if (.not. ok) error stop 'joint kinetic equilibrium did not converge'
        call place_probes(model, [0.3_dp, 0.5_dp, 0.7_dp, 0.85_dp], &
            0.4_dp, 0.0_dp, 0.2_dp, 0.0_dp)
        call model%setup()
        deallocate (r)
        allocate (r(model%n_residuals()), rp(model%n_residuals()), &
            rm(model%n_residuals()))
        theta0 = theta
        truth = theta + [0.3_dp*sigma(:npar)/rec%scale, &
            [(0.2_dp*sin(real(k, dp)), k=1, model%n_prof)]]
        call model%residual(truth, r, ok)
        if (.not. ok) error stop 'planted joint equilibrium did not converge'
        call synthesize_magnetics(model, truth(:npar)*rec%scale, 0.003_dp, &
            1.0e-6_dp, 1.0e-7_dp, 1.0e-6_dp, rng, .false.)
        call synthesize_profiles(model, truth, 0.03_dp, [0.0_dp, 0.0_dp], rng, .false.)
        call model%residual(theta, r, ok)
        call model%jacobian(theta, jac)
        call check_columns()
        call rec%eq%vmec%set_input('hot_restart', 1.0_dp)
        theta = theta0
        call model%residual(theta, r, ok)
        initial = sum(r**2)
        call map_fit(model, theta, 12, steps, message)
        call model%residual(theta, r, ok)
        if (.not. ok .or. sum(r**2) >= initial) error stop 'joint fit did not improve'
        call laplace(model, theta, lap)
        if (lap%info /= 0 .or. .not. lap%proper) error stop 'joint Laplace failed'
        if (.not. ieee_is_finite(lap%log_evidence)) &
            error stop 'nonfinite free-boundary evidence'
        if (any(abs(theta - truth) > 2*sqrt([(lap%cov(k, k), k=1, model%n)]))) &
            error stop 'planted joint parameters outside Laplace intervals'
        print '(A,ES12.4)', 'joint maximum parameter pull: ', &
            maxval(abs(theta - truth)/sqrt([(lap%cov(k, k), k=1, model%n)]))
        print '(A,F14.5)', 'log evidence: ', lap%log_evidence
        open (newunit=unit, file=trim(work)//'/recovery.csv', status='replace')
        write (unit, '(A)') 'parameter,truth,fit,sd'
        do k = 1, model%n
            if (k <= npar) then
                write (unit, '(A,3(A,ES22.14))') trim(names(k)), ',', &
                    truth(k)*rec%scale(k), &
                    ',', theta(k)*rec%scale(k), ',', sqrt(lap%cov(k, k))*rec%scale(k)
            else
                write (unit, '(A,I0,3(A,ES22.14))') 'profile_latent_', k - npar, &
                    ',', truth(k), &
                    ',', theta(k), ',', sqrt(lap%cov(k, k))
            end if
        end do
        close (unit)
        if (trim(mode) == 'joint') then
            log_extra = lap%log_evidence
            call model%prof%set_nonthermal([0.0_dp, 1.0_dp], [2.0_dp, 8.0_dp], &
                enabled=.false.)
            theta = theta(:npar + model%prof%n_latent)
            call model%setup()
            if (size(theta) /= model%n) error stop 'thermal ablation latent layout'
            call map_fit(model, theta, 12, steps, message)
            call laplace(model, theta, lap)
            if (lap%info /= 0 .or. .not. lap%proper) &
                error stop 'thermal ablation Laplace failed'
            if (.not. ieee_is_finite(lap%log_evidence)) &
                error stop 'nonfinite ablation evidence'
            print '(A,F14.5)', 'same-data log evidence (extra minus thermal): ', &
                log_extra - lap%log_evidence
        end if
    end if
    print '(A)', 'free-boundary '//trim(mode)//': all checks passed'
contains
    subroutine check_columns()
        real(dp), allocatable :: cc(:, :), bc(:, :), uc(:, :), vc(:, :), zero(:, :)
        real(dp), allocatable :: jf(:, :), ff(:), fp(:), fm(:)
        integer :: col, nb, ng, j
        if (.not. joint) then
            call rec%eq%solve(x, ok)
            if (.not. ok) error stop 'spectral derivative base failed'
            ng = size(rec%eq%coef)
            nb = rec%eq%edge%mnmax_nyq*(rec%eq%ns - 1)
            allocate (cc(ng, 4), bc(nb, 4), uc(nb, 4), vc(nb, 4), zero(4, rec%eq%n_y()))
            cc = 0; bc = 0; uc = 0; vc = 0; zero = 0
            cc(:, 1) = [(sin(real(j, dp)), j=1, ng)]
            bc(:, 2) = [(cos(real(j, dp)), j=1, nb)]
            uc(:, 3) = bc(:, 2)
            vc(:, 4) = bc(:, 2)
            call rec%eq%jacobian(x, zero, jf, cc, bc, uc, vc)
        end if
        call rec%eq%vmec%set_input('hot_restart', 0.0_dp)
        do col = 1, size(jac, 2)
            if (joint) then
                saved = theta
                h = 2.0e-3_dp*max(abs(theta(col)), 0.1_dp)
                theta(col) = saved(col) + h
                call model%residual(theta, rp, ok)
                if (.not. ok) error stop 'positive joint FD failed'
                theta(col) = saved(col) - h
                call model%residual(theta, rm, ok)
                if (.not. ok) error stop 'negative joint FD failed'
                theta = saved
            else
                saved = x
                h = 2.0e-3_dp*abs(x(col))
                x(col) = saved(col) + h
                call residual_at(x, rp)
                fp = field_values()
                x(col) = saved(col) - h
                call residual_at(x, rm)
                fm = field_values()
                x = saved
                ff = (fp - fm)/(2*h)
                err = maxval(abs(jf(:, col) - ff))/max(maxval(abs(ff)), 1.0e-10_dp)
                print '(A,I3,A,ES12.4)', 'spectral column ', col, &
                    ' cold-resolve FD error: ', err
                if (err > 3.0e-3_dp) &
                    error stop 'free-boundary output spectrum derivative mismatch'
            end if
            fd = (rp - rm)/(2*h)
            if (.not. all(ieee_is_finite(fd))) error stop 'nonfinite FD reference'
            err = maxval(abs(jac(:, col) - fd))/max(maxval(abs(fd)), 1.0e-10_dp)
            print '(A,I3,A,ES12.4)', 'column ', col, &
                ' relative cold-resolve FD error: ', err
            if (err > 3.0e-3_dp) error stop 'free-boundary complete derivative mismatch'
        end do
    end subroutine check_columns

    function field_values() result(v)
        !! Output contractions for the real-space diagnostic derivative check.
        real(dp) :: v(4)
        real(dp), allocatable :: b(:)
        character(len=8), parameter :: fields(3) = ['bmnc    ', 'bsupumnc', 'bsupvmnc']
        integer :: j, nq, nb, ng
        ng = size(rec%eq%coef)
        nq = rec%eq%edge%mnmax_nyq
        nb = nq*(rec%eq%ns - 1)
        v(1) = dot_product(reshape(rec%eq%coef, [ng]), [(sin(real(j, dp)), j=1, ng)])
        do j = 1, 3
            b = rec%eq%vmec%wout(trim(fields(j)), nq*rec%eq%ns)
            v(j + 1) = dot_product(b(nq + 1:), [(cos(real(k, dp)), k=1, nb)])
        end do
    end function field_values
end program test_free_boundary
