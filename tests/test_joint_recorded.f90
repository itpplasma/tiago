program test_joint_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use joint_recorded_fixture, only: configure_joint, require
    use joint_signal_io, only: signal_roundtrip
    use joint_validation, only: add_noise, write_profiles, read_start, write_calibrations, &
                                calibration_start
    use tiago_joint_model, only: joint_model_t, diagnostic_slot_t
    use tiago_energy_consistency, only: energy_consistency_t
    use tiago_recorded_diagnostic, only: calibrated_diagnostic_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_point_probe, only: point_probe_t
    use tiago_soft_xray, only: soft_xray_t, pinhole_camera_t, SXR_BREMSSTRAHLUNG_MODEL
    use tiago_profile_tape, only: profile_tape_t
    use tiago_ece_recorded, only: recorded_ece_t
    use tiago_profile_diagnostic, only: plasma_sample_t, Q_TE
    use tiago_real_space_state, only: located_t
    use tiago_reconstruction, only: rec
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use fortnum_status, only: fortnum_status_t
    use fortnum_rng, only: rng_t, rng_seed
    implicit none
    type(joint_model_t)::model, control
    type(laplace_t)::post
    type(fortnum_status_t)::status
    real(dp), allocatable::truth(:), r(:), rp(:), rm(:), theta(:), direction(:), &
                            design(:, :), length(:), surrogate_r(:), exact_r(:), ablation(:), &
                            jdot(:), weight(:), bar(:)
    real(dp), allocatable :: candidate(:), scan(:), scores(:), starts(:, :)
    real(dp)::h, error, t0, t1, profile_error, best_chi, chi, truth_chi
    character(1024)::input, sensors, work, sources, start_path
    character(32)::rep_text
    character(:), allocatable::message
    integer::i, k, steps, unit, replica, chosen(1), j
    integer :: calibration_indices(4)
    integer(int64) :: first_tick, last_tick, tick_rate
    logical::ok
    call get_command_argument(1, input)
    call get_command_argument(2, sensors)
    call get_command_argument(3, work)
    call get_command_argument(4, sources)
    call configure_joint(model, trim(input), trim(sensors), trim(work), truth)
    call get_command_argument(5, rep_text)
    replica = 0
    if (len_trim(rep_text) > 0) read (rep_text, *) replica
    if (replica == -5) then
        call check_wrapped_ece(model, truth)
        stop
    end if
    if (replica == -4) then
        call check_grouped(model, truth)
        stop
    end if
    if (replica == -6) then
        call check_soft_xray(model, truth)
        stop
    end if
    if (replica == -2) then
        allocate (direction(model%n), jdot(model%n_residuals()))
        direction = [(sin(real(i, dp)), i=1, model%n)]
        direction(1) = .01_dp*direction(1)
        call model%jvp(truth, direction, jdot)
        call system_clock(first_tick, tick_rate)
        call cpu_time(t0)
        do i = 1, 3
            call model%jvp(truth, direction, jdot)
        end do
        call cpu_time(t1)
        call system_clock(last_tick)
        print '(A,3ES20.10)', 'joint JVP wall/CPU seconds per call and checksum=', &
            real(last_tick - first_tick, dp)/tick_rate/3, (t1 - t0)/3, sum(jdot)
        stop
    end if
    if (replica > 0) then
        call get_command_argument(6, start_path)
        allocate (theta(model%n))
        call read_start(trim(start_path), theta)
        call add_noise(model, replica)
        call map_fit(model, theta, 30, steps, message)
        call laplace(model, theta, post)
        call require(post%info == 0, 'replicated posterior covariance')
        call write_profiles(model, theta, post%cov, truth, trim(work)//'/profiles.csv', &
                            profile_error)
        call write_calibrations(model, theta, post%cov, truth, trim(work)//'/calibrations.csv')
        print *, 'replica', replica, 'steps', steps, 'profile error', profile_error, message
        call require(index(message, 'converged') == 1, 'replicated fit converges')
        stop
    end if
    allocate (r(model%n_residuals()), rp(model%n_residuals()), rm(model%n_residuals()))
    call model%residual(truth, r, ok)
    call require(ok, 'exact truth residual')
    truth_chi = sum(r**2)
    print '(A,ES14.6)', 'truth chi-square=', truth_chi
    theta = rec%meas_value
    call signal_roundtrip(trim(sources), 'recorded-mag', 'Wb', theta, rec%meas_sigma**2)
    rec%meas_value = theta
    do i = 1, 3
        select case (i)
        case (1)
            call signal_roundtrip(trim(sources), 'recorded-ts', 'V s', &
                                  model%diags(i)%data, [(model%diags(i)%covariance(k, k), &
                                                         k=1, size(model%diags(i)%data))])
        case (2)
            call signal_roundtrip(trim(sources), 'recorded-di', 'V', &
                                  model%diags(i)%data, [(model%diags(i)%covariance(k, k), &
                                                         k=1, size(model%diags(i)%data))])
        case (3)
            call signal_roundtrip(trim(sources), 'recorded-ece', 'V', &
                                  model%diags(i)%data, [(model%diags(i)%covariance(k, k), &
                                                         k=1, size(model%diags(i)%data))])
        end select
    end do
    call model%residual(truth, rp, ok)
    call require(maxval(abs(r - rp)) < 1.e-10_dp, 'source posterior parity')
    call require(maxval(abs(r(:model%n_mag))) < 1.e-6_dp, 'calibrated magnetic truth')
    allocate (direction(model%n))
    direction = [(sin(real(i, dp)), i=1, model%n)]
    direction(1) = .01_dp*direction(1)
    allocate (jdot(model%n_residuals()), weight(model%n_residuals()), bar(model%n))
    weight = [(cos(real(i, dp)), i=1, model%n_residuals())]
    call model%jvp(truth, direction, jdot)
    call model%vjp(truth, weight, bar)
    error = abs(dot_product(weight, jdot) - dot_product(direction, bar))
    print '(A,ES14.6)', 'joint adjoint absolute error=', error
    call require(error < 1.e-9_dp*max(1._dp, norm2(jdot)), 'joint adjoint identity')
    h = 1.e-3_dp
    call model%residual(truth + h*direction, rp, ok)
    call require(ok, 'positive joint perturbation')
    call model%residual(truth - h*direction, rm, ok)
    call require(ok, 'negative joint perturbation')
    error = norm2((rp - rm)/(2*h) - jdot)/norm2(jdot)
    print '(A,ES14.6)', 'full recorded joint derivative relative error=', error
    print *, 'mag errors', maxval(abs((rp(:model%n_mag) - rm(:model%n_mag))/(2*h) - &
                                      jdot(:model%n_mag)))
    k = model%n_mag + model%n_fb + model%n_rec_prior
    do i = 1, size(model%diags)
        print *, 'diag derivative', i, norm2((rp(k + 1:k + size(model%diags(i)%data)) - &
                                              rm(k + 1:k + size(model%diags(i)%data)))/(2*h) - &
                                             jdot(k + 1:k + size(model%diags(i)%data))), &
            norm2(jdot(k + 1:k + size(model%diags(i)%data)))
        k = k + size(model%diags(i)%data)
    end do
    call require(error < 5.e-4_dp, 'full recorded joint derivative')
    call check_options(model, truth)
    call check_energy(model, truth)
    if (replica == -3) stop
    theta = 0*truth; theta(1) = rec%prior(1)/rec%scale(1)
    call model%residual(theta, r, ok)
    print '(A,ES14.6)', 'initial chi-square=', sum(r**2)
    best_chi = huge(1._dp)
    allocate (candidate(model%n))
    allocate (scan(201), scores(201), starts(model%n, 201))
    do j = 1, size(scan)
        scan(j) = rec%prior(1)/rec%scale(1) + (real(j - 101, dp)/100)* &
                  rec%prior_sigma(1)/rec%scale(1)
        candidate = 0; candidate(1) = scan(j)
        call model%residual(candidate, r, ok)
        scores(j) = huge(1._dp)
        if (ok) then
            call calibration_start(model, candidate, r)
            scores(j) = sum(r**2)
        end if
        starts(:, j) = candidate
    end do
    do k = 1, 3
        chosen = minloc(scores)
        candidate = starts(:, chosen(1))
        call map_fit(model, candidate, 40, steps, message)
        call model%residual(candidate, r, ok)
        call require(ok, 'multistart final equilibrium')
        chi = sum(r**2)
        print *, 'start', k, 'accepted steps', steps, 'chi-square', chi, message
        if (chi < best_chi) then
            best_chi = chi; theta = candidate
        end if
        where (abs(scan - scan(chosen(1))) < .0005_dp) scores = huge(1._dp)
    end do
    call laplace(model, theta, post)
    call require(post%info == 0, 'joint covariance')
    call model%residual(theta, r, ok)
    print '(A,ES14.6)', 'fitted chi-square=', sum(r**2)
    call require(sum(r**2) < truth_chi + .01_dp, 'fit reaches the truth posterior density')
    print '(A,4F12.5)', 'injected calibration latents=', &
        truth(model%diags(1)%eta_offset + 1), truth(model%diags(2)%eta_offset + 6), &
        truth(model%diags(3)%eta_offset + 1), truth(model%mag_calibration_index + 3)
    print '(A,4F12.5)', 'recovered calibration latents=', &
        theta(model%diags(1)%eta_offset + 1), theta(model%diags(2)%eta_offset + 6), &
        theta(model%diags(3)%eta_offset + 1), theta(model%mag_calibration_index + 3)
    calibration_indices = [model%diags(1)%eta_offset + 1, model%diags(2)%eta_offset + 6, &
                           model%diags(3)%eta_offset + 1, model%mag_calibration_index + 3]
    do i = 1, size(calibration_indices)
        k = calibration_indices(i)
        call require(abs(theta(k) - truth(k)) <= 1.95996398454_dp*sqrt(post%cov(k, k)), &
                     'injected calibration lies in the posterior interval')
        call require(sqrt(post%cov(k, k)) < .9_dp, 'joint calibration contracts its prior')
    end do
    open (newunit=unit, file=trim(work)//'/fit.csv', status='replace')
    write (unit, '(A)') 'parameter,truth,map,sd'
    do i = 1, model%n
        write (unit, '(I0,3(",",ES22.14))') i, truth(i), theta(i), sqrt(post%cov(i, i))
    end do
    close (unit)
    call write_calibrations(model, theta, post%cov, truth, trim(work)//'/calibrations.csv')
    call write_profiles(model, theta, post%cov, truth, trim(work)//'/profiles-joint.csv', &
                        profile_error)
    print '(A,ES14.6)', 'joint relative profile RMSE=', profile_error
    allocate (design(model%n, 3), length(model%n))
    design(:, 1) = theta; design(:, 2) = theta + .02_dp*direction; design(:, 3) = theta - .02_dp*direction
    length = 1
    call model%train_surrogate(design, length, 1._dp, 1._dp, 1.e-6_dp, status)
    call require(status%code == 0, 'paired exact/surrogate training')
    allocate (exact_r(model%n_residuals()), surrogate_r(model%n_residuals()))
    call model%residual(design(:, 3), exact_r, ok)
    model%use_surrogate = .true.
    call model%residual(design(:, 3), surrogate_r, ok)
    call require(ok, 'surrogate evaluation')
    error = maxval(abs(exact_r - surrogate_r))
    print '(A,ES14.6)', 'surrogate training residual error=', error
    call require(error < 1.e-6_dp, 'autoregressive correction matches exact training value')
    call cpu_time(t0)
    do i = 1, 100
        call model%residual(theta, surrogate_r, ok)
    end do
    call cpu_time(t1)
    print '(A,ES14.6)', 'surrogate CPU seconds per call=', (t1 - t0)/100
    call model%jvp(theta, direction, jdot)
    call model%residual(theta + 1.e-5_dp*direction, rp, ok)
    call require(ok, 'positive surrogate derivative perturbation')
    call model%residual(theta - 1.e-5_dp*direction, rm, ok)
    call require(ok, 'negative surrogate derivative perturbation')
    error = norm2((rp - rm)/2.e-5_dp - jdot)/max(1._dp, norm2(jdot))
    print '(A,ES14.6)', 'surrogate derivative relative error=', error
    call require(error < 1.e-7_dp, 'surrogate finite differences')
    call model%residual(theta + .01_dp*direction, surrogate_r, ok)
    call require(ok, 'held-out surrogate query')
    model%use_surrogate = .false.
    call cpu_time(t0)
    call model%residual(theta + .01_dp*direction, exact_r, ok)
    call cpu_time(t1)
    call require(ok, 'held-out exact query')
    print '(A,ES14.6)', 'cold exact CPU seconds per call=', t1 - t0
    error = maxval(abs(exact_r - surrogate_r))
    print '(A,ES14.6)', 'held-out surrogate error in noise standard deviations=', error
    call require(error < .05_dp, 'held-out surrogate error budget')
    model%use_surrogate = .true.
    call model%residual(theta + direction, rp, ok)
    call require(.not. ok, 'surrogate refuses extrapolation')
    call model%update_covariances()
    model%use_surrogate = .true.
    call model%residual(theta, rp, ok)
    call require(.not. ok, 'covariance update invalidates surrogate')
    model%use_surrogate = .false.
    call model%residual(theta, exact_r, ok)
    call require(ok, 'exact final recheck')
    if (replica == -1) then
        control = model; control%eta_factor = 0; control%mag_calibration%factor = 0
        ablation = theta; ablation(model%n_rec + model%n_prof + 1:) = 0
        call map_fit(control, ablation, 40, steps, message)
        call laplace(control, ablation, post)
        call require(post%info == 0, 'fixed calibration covariance')
        call write_profiles(control, ablation, post%cov, truth, &
                            trim(work)//'/profiles-fixed-calibration.csv', profile_error)
        print '(A,ES14.6)', 'fixed calibration relative profile RMSE=', profile_error
        control = model
        control%mag_covariance = 1.e12_dp*control%mag_covariance
        do i = 2, 3
            control%diags(i)%covariance = 1.e12_dp*control%diags(i)%covariance
        end do
        call control%update_covariances()
        ablation = theta
        call map_fit(control, ablation, 40, steps, message)
        call laplace(control, ablation, post)
        call require(post%info == 0, 'uncoupled covariance')
        call write_profiles(control, ablation, post%cov, truth, &
                            trim(work)//'/profiles-thomson-only.csv', profile_error)
        print '(A,ES14.6)', 'Thomson-only relative profile RMSE=', profile_error
    end if
contains
    subroutine check_wrapped_ece(base, injected)
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: wrapped_model
        type(calibrated_diagnostic_t) :: wrapper, inner_wrapper
        type(calibration_prior_t) :: scalar_prior, grouped_prior
        class(profile_tape_t), allocatable :: tape
        type(recorded_ece_t) :: native
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        type(fortnum_status_t) :: code
        real(dp), allocatable :: z(:), eta(:), dx(:, :, :), direction(:)
        real(dp), allocatable :: r(:), rp(:), rm(:), tangent(:), bar(:), weights(:)
        real(dp), allocatable :: pdot(:), pbar(:), edot(:), ebar(:)
        real(dp) :: raw(4), predicted(4), expected(4), ydot(4), ybar(4), c(2)
        real(dp) :: h, fd_error, coarse_error, adjoint_error, tape_error
        integer :: mode, i, k, info, np
        logical :: ok
        call scalar_prior%init([1.05_dp, .02_dp], &
            reshape([.0004_dp, 0._dp, 0._dp, .0001_dp], [2, 2]), info)
        call require(info == 0, 'scalar ECE wrapper prior')
        call grouped_prior%init(log([1.02_dp, .98_dp]), &
            reshape([.0004_dp, .0002_dp, .0002_dp, .0005_dp], [2, 2]), info)
        call require(info == 0, 'grouped ECE wrapper prior')
        call base%predict_diagnostic(injected, 3, raw, sample, loc, eta, dx, native)
        call require(native%valid, 'bare ECE legacy tape parity')
        do mode = 1, 3
            wrapped_model = base
            if (mode == 1) then
                call wrapper%init(base%diags(3)%d, scalar_prior, info)
            else if (mode == 2) then
                call wrapper%init(base%diags(3)%d, grouped_prior, info, &
                    groups=[1, 1, 2, 2], log_gains=.true.)
            else
                call inner_wrapper%init(base%diags(3)%d, scalar_prior, info)
                call require(info == 0, 'inner ECE wrapper setup')
                call wrapper%init(inner_wrapper, grouped_prior, info, &
                    groups=[1, 1, 2, 2], log_gains=.true.)
            end if
            call require(info == 0, 'ECE wrapper setup')
            deallocate (wrapped_model%diags(3)%d, wrapped_model%diags(3)%eta_mean, &
                wrapped_model%diags(3)%eta_sigma)
            allocate (wrapped_model%diags(3)%d, source=wrapper)
            call wrapped_model%setup()
            wrapped_model%eta_factor(:base%n_eta, :base%n_eta) = base%eta_factor
            k = base%mag_calibration_index
            allocate (z(wrapped_model%n)); z = 0
            z(:k) = injected(:k)
            z(k + 1:k + 2) = [.3_dp, -.2_dp]
            if (mode == 3) z(k + 3:k + 4) = [.15_dp, -.1_dp]
            z(wrapped_model%mag_calibration_index + 1:) = injected(k + 1:)
            allocate (r(wrapped_model%n_residuals()))
            call wrapped_model%residual(z, r, ok)
            call require(ok, 'calibrated ECE joint preparation')
            call wrapped_model%predict_diagnostic(z, 3, predicted, sample, &
                loc, eta, dx, record=tape)
            call require(allocated(tape), 'wrapped ECE provider tape')
            call require(tape%valid, 'wrapped ECE valid tape')
            expected = raw
            if (mode /= 2) then
                expected = (1.05_dp + .02_dp*.3_dp)*raw + .02_dp - .01_dp*.2_dp
            end if
            if (mode > 1) then
                c = log([1.02_dp, .98_dp]) + [.02_dp*.3_dp, .01_dp*.3_dp - .02_dp*.2_dp]
                if (mode == 3) c = log([1.02_dp, .98_dp]) + &
                    [.02_dp*.15_dp, .01_dp*.15_dp - .02_dp*.1_dp]
                expected(:2) = exp(c(1))*expected(:2)
                expected(3:) = exp(c(2))*expected(3:)
            end if
            call require(maxval(abs(predicted - expected))/maxval(abs(expected)) < &
                1.e-12_dp, 'wrapped ECE independent affine/log gain oracle')
            np = wrapped_model%state%n_params() + wrapped_model%n_prof
            allocate (pdot(np), pbar(np), edot(size(eta)), ebar(size(eta)))
            pdot = [(sin(.2_dp*i), i=1, np)]
            edot = [(cos(.3_dp*i), i=1, size(eta))]
            ybar = [.3_dp, -.2_dp, .7_dp, -.1_dp]
            call tape%signal_jvp(pdot, edot, ydot, code)
            call require(code%code == 0, 'wrapped ECE tape tangent')
            call tape%signal_vjp(ybar, pbar, ebar, code)
            call require(code%code == 0, 'wrapped ECE tape adjoint')
            tape_error = abs(dot_product(ydot, ybar) - dot_product(pdot, pbar) - &
                dot_product(edot, ebar))/max(1._dp, norm2(ydot)*norm2(ybar))
            call require(tape_error < 1.e-10_dp, 'independent wrapped ECE tape adjoint')
            allocate (direction(wrapped_model%n), tangent(size(r)), bar(size(z)), &
                weights(size(r)), rp(size(r)), rm(size(r)))
            direction = [(sin(real(i, dp)), i=1, size(z))]
            direction(1) = .01_dp*direction(1)
            weights = [(cos(real(i, dp)), i=1, size(r))]
            call wrapped_model%jvp(z, direction, tangent)
            call wrapped_model%vjp(z, weights, bar)
            adjoint_error = abs(dot_product(tangent, weights) - &
                dot_product(direction, bar))/max(1._dp, norm2(tangent))
            h = 1.e-3_dp
            call wrapped_model%residual(z + h*direction, rp, ok)
            call require(ok, 'wrapped ECE positive perturbation')
            call wrapped_model%residual(z - h*direction, rm, ok)
            call require(ok, 'wrapped ECE negative perturbation')
            coarse_error = norm2((rp - rm)/(2*h) - tangent)/norm2(tangent)
            h = h/2
            call wrapped_model%residual(z + h*direction, rp, ok)
            call require(ok, 'wrapped ECE refined positive perturbation')
            call wrapped_model%residual(z - h*direction, rm, ok)
            call require(ok, 'wrapped ECE refined negative perturbation')
            fd_error = norm2((rp - rm)/(2*h) - tangent)/norm2(tangent)
            call require(fd_error < .35_dp*coarse_error + 1.e-7_dp, &
                'wrapped ECE finite-difference step convergence')
            call require(fd_error < 5.e-4_dp, &
                'wrapped ECE full moving-ray finite difference')
            call require(adjoint_error < 1.e-9_dp, 'wrapped ECE joint adjoint')
            print '(A,I0,4ES14.6)', 'wrapped ECE mode/coarse FD/fine FD/adjoints=', &
                mode, coarse_error, fd_error, adjoint_error, tape_error
            deallocate (z, r, rp, rm, tangent, bar, weights, direction, &
                pdot, pbar, edot, ebar)
        end do
        call check_wrapped_fidelity(base, injected)
    end subroutine

    subroutine check_wrapped_fidelity(base, injected)
        use tiago_multifidelity_profile, only: ar1_profile_t
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: bare, wrapped
        type(ar1_profile_t) :: prior
        type(calibration_prior_t) :: scale, gain
        type(calibrated_diagnostic_t) :: wrapper
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: loc(:)
        real(dp), allocatable :: z(:), wz(:), residual(:), eta(:), dx(:, :, :)
        real(dp) :: low(4), high(4), predicted(4), error, separation
        integer :: info, i, k, offset, n
        logical :: ok
        bare = base
        call scale%init([1.2_dp], reshape([.1_dp**2], [1, 1]), info)
        call require(info == 0, 'wrapped ECE AR1 scale prior')
        scale%record%record_id = 'synthetic-ECE-scale'
        call prior%init(bare%prof%greville, &
            1200*exp(-.6_dp*bare%prof%greville)/1.2_dp, scale, &
            log([100._dp, .12_dp, 50._dp, .1_dp]), [.2_dp, .2_dp, .2_dp, .2_dp], &
            'synthetic', 'constant and calibration tangent', info)
        call require(info == 0, 'wrapped ECE AR1 profile')
        call bare%prof%set_multifidelity_te(prior)
        call bare%setup()
        allocate (z(bare%n), residual(bare%n_residuals()))
        z = 0
        z(:bare%n_rec) = injected(:bare%n_rec)
        offset = bare%n_rec + bare%prof%te_ar1%offset
        n = bare%prof%te_ar1%n_params
        z(offset + 1:offset + n) = [(0.2_dp*sin(.4_dp*i), i=1, n)]
        call bare%residual(z, residual, ok)
        call require(ok, 'bare AR1 physical equilibrium')
        call bare%predict_diagnostic(z, 3, low, sample, loc, eta, dx)
        bare%diags(3)%low_fidelity_te = .false.
        call bare%predict_diagnostic(z, 3, high, sample, loc, eta, dx)
        separation = norm2(high - low)/norm2(low)
        call require(separation > 1.e-3_dp, 'AR1 low/high prediction control differs')
        bare%diags(3)%low_fidelity_te = .true.
        wrapped = bare
        call gain%init([1._dp], reshape([.02_dp**2], [1, 1]), info)
        call require(info == 0, 'AR1 identity calibration prior')
        call wrapper%init(bare%diags(3)%d, gain, info)
        call require(info == 0, 'AR1 ECE wrapper setup')
        deallocate (wrapped%diags(3)%d, wrapped%diags(3)%eta_mean, &
            wrapped%diags(3)%eta_sigma)
        allocate (wrapped%diags(3)%d, source=wrapper)
        wrapped%diags(3)%low_fidelity_te = .false.
        call wrapped%setup()
        wrapped%eta_factor(:bare%n_eta, :bare%n_eta) = bare%eta_factor
        allocate (wz(wrapped%n))
        k = bare%mag_calibration_index
        wz = 0
        wz(:k) = z(:k)
        wz(wrapped%mag_calibration_index + 1:) = z(k + 1:)
        deallocate (residual)
        allocate (residual(wrapped%n_residuals()))
        call wrapped%residual(wz, residual, ok)
        call require(ok, 'wrapped AR1 physical equilibrium')
        call wrapped%predict_diagnostic(wz, 3, predicted, sample, loc, eta, dx)
        error = norm2(predicted - low)/norm2(low)
        call require(error < 1.e-12_dp, 'active AR1 bare/wrapped ECE prediction parity')
        print '(A,2ES14.6)', 'wrapped AR1 parity error/low-high separation=', &
            error, separation
    end subroutine

    subroutine check_soft_xray(base, injected)
        !! Two pinhole cameras with one log gain each replace the ECE slot
        !! (emissivity-profile model, latent emissivity) and form a fourth slot
        !! (bremsstrahlung from the shared ne, Te). Oracles: the complete joint
        !! directional derivative against central differences of re-solved
        !! equilibria, the adjoint identity, and a local MAP fit from unit
        !! bremsstrahlung camera gains that recovers their ratio. Their common
        !! factor trades against the density and temperature scale
        !! (eps ~ ne^2 Te^(1/2)) and is set by the other diagnostics and priors.
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: fitted
        type(soft_xray_t) :: profile_sxr, brems_sxr
        type(pinhole_camera_t) :: cams(2)
        type(calibrated_diagnostic_t) :: wrapped(2)
        type(calibration_prior_t) :: prior
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: locations(:)
        type(diagnostic_slot_t), allocatable :: slots(:)
        real(dp), allocatable :: z(:), start(:), eta(:), dx(:, :, :), y(:), r(:), rp(:), rm(:)
        real(dp), allocatable :: tangent(:), weights(:), direction(:), bar(:), mean(:)
        real(dp) :: xa(3), xb(3), xu(3, 3), b, db(3), size0, h, fd_error, adjoint_error
        character(:), allocatable :: message
        integer :: code, i, j, k, old_eta, steps, first(2)
        logical :: ok

        fitted = base
        call fitted%state%eval([0._dp, 0._dp, 0._dp], xa, xu, b, db)
        call fitted%state%eval([1._dp, 0._dp, 0._dp], xb, xu, b, db)
        size0 = norm2(xb - xa)
        cams(1)%name = 'top'
        cams(1)%aperture = xa + [0._dp, 0._dp, 2*size0]
        allocate (cams(1)%detector(3, 5))
        do i = 1, 5
            cams(1)%detector(:, i) = cams(1)%aperture + &
                .05_dp*size0*[.15_dp*(i - 3), 0._dp, 1._dp]
        end do
        cams(1)%length = 4*size0
        cams(2)%name = 'outboard'
        cams(2)%aperture = xb + (xb - xa)
        allocate (cams(2)%detector(3, 4))
        do i = 1, 4
            cams(2)%detector(:, i) = cams(2)%aperture + .05_dp*size0* &
                ((xb - xa)/size0 + [0._dp, 0._dp, .2_dp*(i - 2.5_dp)])
        end do
        cams(2)%length = 3*size0
        cams(2)%filter%energy = [300._dp, 1500._dp]
        cams(2)%filter%transmission = [0._dp, .8_dp]
        call profile_sxr%set_cameras(cams, 40)
        call profile_sxr%set_emissivity_basis([0._dp, .5_dp, 1._dp])
        mean = 1.e4_dp*(1 - .9_dp*profile_sxr%emissivity%greville)
        call profile_sxr%set_emissivity_prior(mean, .5_dp, .5_dp, .3_dp, .6_dp, .2_dp)
        call brems_sxr%set_cameras(cams, 40)
        brems_sxr%model = SXR_BREMSSTRAHLUNG_MODEL
        call make_log_gains(prior)
        call wrapped(1)%init(profile_sxr, prior, code, groups=profile_sxr%camera, &
            log_gains=.true.)
        call require(code == 0, 'joint SXR profile wrapper')
        call wrapped(2)%init(brems_sxr, prior, code, groups=brems_sxr%camera, &
            log_gains=.true.)
        call require(code == 0, 'joint SXR bremsstrahlung wrapper')

        call move_alloc(fitted%diags, slots)
        allocate (fitted%diags(4))
        fitted%diags(1:2) = slots(1:2)
        do i = 3, 4
            allocate (fitted%diags(i)%d, source=wrapped(i - 2))
            j = wrapped(i - 2)%n_signals()
            allocate (fitted%diags(i)%data(j), fitted%diags(i)%covariance(j, j))
            fitted%diags(i)%data = 0
            fitted%diags(i)%covariance = 0
            do k = 1, j
                fitted%diags(i)%covariance(k, k) = 1
            end do
        end do
        old_eta = base%diags(3)%eta_offset - base%n_rec - base%n_prof
        call fitted%setup()
        fitted%eta_factor(:old_eta, :old_eta) = base%eta_factor(:old_eta, :old_eta)
        allocate (z(fitted%n)); z = 0
        k = fitted%diags(3)%eta_offset
        z(:k) = injected(:k)
        z(fitted%mag_calibration_index + 1:) = injected(base%mag_calibration_index + 1:)
        z(k + 1:k + profile_sxr%n_nuisance()) = [(.3_dp*sin(1.1_dp*i), &
            i=1, profile_sxr%n_nuisance())]
        first = [fitted%diags(3)%eta_offset + profile_sxr%n_nuisance(), &
            fitted%diags(4)%eta_offset]
        z(first(1) + 1:first(1) + 2) = [.5_dp, -.4_dp]
        z(first(2) + 1:first(2) + 2) = [.6_dp, -.7_dp]
        do i = 3, 4
            allocate (y(fitted%diags(i)%d%n_signals()))
            call fitted%predict_diagnostic(z, i, y, sample, locations, eta, dx)
            call require(all(y > 0), 'positive synthetic SXR signals')
            fitted%diags(i)%data = y
            do k = 1, size(y)
                fitted%diags(i)%covariance(k, k) = (.005_dp*y(k))**2
            end do
            deallocate (y)
        end do
        call fitted%update_covariances()

        allocate (r(fitted%n_residuals()), rp(fitted%n_residuals()), &
            rm(fitted%n_residuals()), tangent(fitted%n_residuals()), &
            weights(fitted%n_residuals()), direction(fitted%n), bar(fitted%n))
        direction = [(sin(real(i, dp)), i=1, fitted%n)]
        direction(1) = .01_dp*direction(1)
        weights = [(cos(real(i, dp)), i=1, size(weights))]
        call fitted%jvp(z, direction, tangent)
        call fitted%vjp(z, weights, bar)
        adjoint_error = abs(dot_product(weights, tangent) - &
            dot_product(direction, bar))/max(1._dp, norm2(tangent))
        h = 1.e-3_dp
        call fitted%residual(z + h*direction, rp, ok)
        call require(ok, 'joint SXR positive perturbation')
        call fitted%residual(z - h*direction, rm, ok)
        call require(ok, 'joint SXR negative perturbation')
        fd_error = norm2((rp - rm)/(2*h) - tangent)/norm2(tangent)
        print '(A,2ES14.6)', 'joint SXR FD/adjoint relative errors=', fd_error, adjoint_error
        call require(fd_error < 5.e-4_dp, 'joint SXR full directional oracle')
        call require(adjoint_error < 1.e-9_dp, 'joint SXR adjoint oracle')

        start = z
        start(first(2) + 1:first(2) + 2) = 0
        call map_fit(fitted, start, 20, steps, message)
        call fitted%residual(start, r, ok)
        call require(ok .and. index(message, 'converged') == 1, 'joint SXR MAP convergence')
        print '(A,2F12.5)', 'SXR bremsstrahlung camera log-gain latents injected=', &
            z(first(2) + 1:first(2) + 2)
        print '(A,2F12.5)', 'SXR bremsstrahlung camera log-gain latents recovered=', &
            start(first(2) + 1:first(2) + 2)
        call require(abs((start(first(2) + 1) - start(first(2) + 2)) - &
            (z(first(2) + 1) - z(first(2) + 2))) < .01_dp, 'SXR camera gain ratio recovered')
    end subroutine check_soft_xray

    subroutine make_log_gains(prior)
        !! Camera log gains log(1e5) and log(2e5), standard deviation 0.3.
        type(calibration_prior_t), intent(out) :: prior
        integer :: code
        call prior%init(log([1.e5_dp, 2.e5_dp]), reshape([.09_dp, 0._dp, 0._dp, .09_dp], &
            [2, 2]), code)
        call require(code == 0, 'SXR log-gain prior')
    end subroutine make_log_gains

    subroutine check_grouped(base, injected)
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: fitted, fixed
        type(calibrated_diagnostic_t) :: grouped
        type(calibration_prior_t) :: prior
        type(point_probe_t) :: probe
        type(rng_t) :: rng
        type(fortnum_status_t) :: rng_status
        type(plasma_sample_t) :: sample
        type(located_t), allocatable :: locations(:)
        real(dp), allocatable :: z(:), start(:), baseline(:), eta(:), dx(:, :, :)
        real(dp), allocatable :: prediction(:), r(:), rp(:), rm(:), tangent(:)
        real(dp), allocatable :: direction(:), weights(:), bar(:)
        real(dp) :: covariance(2, 2), factor(2, 2), h, fd_error, adjoint_error
        real(dp) :: fitted_chi, fixed_chi, xu(3, 3), b, db(3)
        real(dp) :: truth_gain(2), fitted_gain(2), fixed_gain(2)
        character(:), allocatable :: message
        integer :: code, i, k, old_eta, steps, gain_unit
        logical :: ok
        fitted = base
        factor = 0
        do i = 1, 2
            factor(i, i) = .15_dp
        end do
        factor(2, 1) = .1_dp
        covariance = matmul(factor, transpose(factor))
        call prior%init([0._dp, 0._dp], covariance, code)
        call require(code == 0, 'joint grouped prior')
        allocate (probe%x0(3, 4))
        probe%quantity = [Q_TE, Q_TE, Q_TE, Q_TE]
        probe%aligned = .true.
        do i = 1, 4
            call fitted%state%eval([.2_dp*i, 0._dp, 0._dp], probe%x0(:, i), xu, b, db)
        end do
        call grouped%init(probe, prior, code, groups=[1, 1, 2, 2], &
            log_gains=.true.)
        call require(code == 0, 'joint grouped wrapper')
        deallocate (fitted%diags(3)%d, fitted%diags(3)%eta_sigma, &
            fitted%diags(3)%eta_mean)
        allocate (fitted%diags(3)%d, source=grouped)
        old_eta = base%diags(3)%eta_offset - base%n_rec - base%n_prof
        call fitted%setup()
        fitted%eta_factor(:old_eta, :old_eta) = base%eta_factor(:old_eta, :old_eta)
        allocate (z(fitted%n)); z = 0
        k = fitted%diags(3)%eta_offset
        z(:k) = injected(:k)
        z(fitted%mag_calibration_index + 1:) = &
            injected(base%mag_calibration_index + 1:)
        z(k + 1:k + 2) = [.8_dp, -.6_dp]
        allocate (prediction(4))
        call fitted%predict_diagnostic(z, 3, prediction, sample, locations, eta, dx)
        fitted%diags(3)%data = prediction
        fitted%diags(3)%covariance = 0
        do i = 1, 4
            fitted%diags(3)%covariance(i, i) = (.005_dp*prediction(i))**2
        end do
        call fitted%update_covariances()
        allocate (r(fitted%n_residuals()), rp(fitted%n_residuals()), &
            rm(fitted%n_residuals()), tangent(fitted%n_residuals()), &
            weights(fitted%n_residuals()), direction(fitted%n), bar(fitted%n))
        direction = [(sin(real(i, dp)), i=1, fitted%n)]
        direction(1) = .01_dp*direction(1)
        weights = [(cos(real(i, dp)), i=1, size(weights))]
        call fitted%jvp(z, direction, tangent)
        call fitted%vjp(z, weights, bar)
        adjoint_error = abs(dot_product(weights, tangent) - &
            dot_product(direction, bar))/max(1._dp, norm2(tangent))
        h = 1.e-3_dp
        call fitted%residual(z + h*direction, rp, ok)
        call require(ok, 'joint grouped positive perturbation')
        call fitted%residual(z - h*direction, rm, ok)
        call require(ok, 'joint grouped negative perturbation')
        fd_error = norm2((rp - rm)/(2*h) - tangent)/norm2(tangent)
        call require(fd_error < 5.e-4_dp, 'joint grouped full directional oracle')
        call require(adjoint_error < 1.e-9_dp, 'joint grouped adjoint oracle')
        start = z
        call rng_seed(rng, 20261105_int64, rng_status)
        call fitted%draw_prior(rng, start, ok)
        call require(ok, 'grouped calibration included in joint prior sampling')
        start = z
        start(k + 1:k + 2) = 0
        fixed = fitted
        fixed%eta_factor(old_eta + 1:, :) = 0
        baseline = start
        call map_fit(fixed, baseline, 20, steps, message)
        call fixed%residual(baseline, r, ok)
        call require(ok, 'incorrect fixed calibration fit')
        fixed_chi = sum(r**2)
        call map_fit(fitted, start, 20, steps, message)
        call fitted%residual(start, r, ok)
        call require(ok, 'trusted-anchor grouped fit')
        fitted_chi = sum(r**2)
        call require(index(message, 'converged') == 1, 'grouped joint MAP convergence')
        call require(fitted_chi < fixed_chi, &
            'anchors prefer fitted correlated calibration')
        call require(maxval(abs(start(k + 1:k + 2) - z(k + 1:k + 2))) < .05_dp, &
            'anchored grouped calibration recovery in prior coordinates')
        print '(A,2ES14.6)', 'grouped joint FD/adjoint relative errors=', &
            fd_error, adjoint_error
        print '(A,2ES14.6)', 'grouped/fixed-calibration posterior chi-square=', &
            fitted_chi, fixed_chi
        print '(A,2F12.5)', 'grouped injected latents=', z(k + 1:k + 2)
        print '(A,2F12.5)', 'grouped recovered latents=', start(k + 1:k + 2)
        truth_gain = prior%physical(z(k + 1:k + 2))
        fitted_gain = prior%physical(start(k + 1:k + 2))
        fixed_gain = prior%physical(baseline(k + 1:k + 2))
        open (newunit=gain_unit, file=trim(work)//'/grouped-gains.csv', status='replace')
        write (gain_unit, '(A)') 'group,injected_log_gain,fitted_log_gain,fixed_log_gain'
        do i = 1, 2
            write (gain_unit, '(I0,3(",",ES25.17))') &
                i, truth_gain(i), fitted_gain(i), fixed_gain(i)
        end do
        close (gain_unit)
    end subroutine
    subroutine check_energy(base, injected)
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: energy_model
        type(energy_consistency_t) :: energy
        real(dp), allocatable :: eta(:), dx(:, :, :), point(:), tangent(:), product(:), high(:), low(:)
        integer :: old_eta_end, i, k
        logical :: valid
        real(dp) :: relative
        call base%unpack(injected, eta=eta, d=1)
        call base%diags(1)%d%sample_points(eta, energy%points, dx)
        allocate (energy%volume_weights(size(energy%points, 2)))
        energy%volume_weights = 1._dp/size(energy%points, 2)
        energy_model = base
        energy_model%diags = [base%diags, diagnostic_slot_t()]
        allocate (energy_model%diags(4)%d, source=energy)
        energy_model%diags(4)%data = [0._dp]
        allocate (energy_model%diags(4)%covariance(1, 1))
        energy_model%diags(4)%covariance = 1.e4_dp
        energy_model%diags(4)%eta_mean = [-2.e-5_dp, 1._dp, 1.42_dp, 0._dp]
        energy_model%diags(4)%eta_sigma = [1.e-6_dp, .01_dp, .01_dp, 1.e-6_dp]
        call energy_model%setup()
        energy_model%eta_factor(:base%n_eta, :base%n_eta) = base%eta_factor
        allocate (point(energy_model%n), tangent(energy_model%n), &
                  product(energy_model%n_residuals()), high(energy_model%n_residuals()), &
                  low(energy_model%n_residuals()))
        point = 0
        old_eta_end = base%n_rec + base%n_prof + base%n_eta
        point(:old_eta_end) = injected(:old_eta_end)
        k = energy_model%mag_calibration_index
        point(k + 1:k + 3) = injected(base%mag_calibration_index + 1:base%mag_calibration_index + 3)
        tangent = [(sin(real(i, dp)), i=1, energy_model%n)]
        tangent(1) = .01_dp*tangent(1)
        call energy_model%jvp(point, tangent, product)
        call energy_model%residual(point + 1.e-4_dp*tangent, high, valid)
        call require(valid, 'positive energy consistency perturbation')
        call energy_model%residual(point - 1.e-4_dp*tangent, low, valid)
        call require(valid, 'negative energy consistency perturbation')
        k = energy_model%n_mag + energy_model%n_fb + energy_model%n_rec_prior + &
            energy_model%diags(4)%row_offset + 1
        relative = abs((high(k) - low(k))/2.e-4_dp - product(k))/max(1._dp, abs(product(k)))
        print '(A,ES14.6)', 'joint energy consistency derivative relative error=', relative
        call require(relative < 5.e-4_dp, 'joint energy consistency derivative')
    end subroutine

    subroutine check_options(base, injected)
        type(joint_model_t), intent(in) :: base
        real(dp), intent(in) :: injected(:)
        type(joint_model_t) :: robust
        real(dp), allocatable :: point(:), tangent(:), product(:), high(:), low(:), &
                                 protected(:), discrepancy(:)
        integer :: n, i
        logical :: valid
        real(dp) :: relative
        robust = base
        n = size(robust%diags(1)%data)
        deallocate (robust%diags(1)%discrepancy_basis)
        allocate (robust%diags(1)%discrepancy_basis(n, 1), &
                  robust%diags(1)%protected_jacobian(n, 1))
        do i = 1, n
            robust%diags(1)%discrepancy_basis(i, 1) = &
                cos(real(i, dp))*robust%diags(1)%chol(i, i)
            robust%diags(1)%protected_jacobian(i, 1) = &
                sin(real(i, dp))*robust%diags(1)%chol(i, i)
        end do
        robust%diags(1)%student_nu = 4
        robust%diags(1)%eta_mean(1) = .03_dp
        robust%mag_student_nu = 5
        call robust%setup()
        ! Updating the noise metric must retain the discrepancy constraint.
        robust%diags(1)%covariance(1, 1) = 2*robust%diags(1)%covariance(1, 1)
        call robust%update_covariances()
        protected = robust%diags(1)%protected_jacobian(:, 1)
        discrepancy = robust%diags(1)%discrepancy_basis(:, 1)
        do i = 1, n
            protected(i) = protected(i)/robust%diags(1)%chol(i, i)
            discrepancy(i) = discrepancy(i)/robust%diags(1)%chol(i, i)
        end do
        call require(abs(dot_product(protected, discrepancy)) < 1.e-13_dp, &
                     'joint covariance update preserves protected discrepancy')
        point = [injected, .2_dp]
        tangent = [(sin(real(i, dp)), i=1, robust%n)]
        tangent(1) = .01_dp*tangent(1)
        point = point + .003_dp*tangent
        allocate (product(robust%n_residuals()), high(robust%n_residuals()), &
                  low(robust%n_residuals()))
        call robust%jvp(point, tangent, product)
        call robust%residual(point + 1.e-4_dp*tangent, high, valid)
        call require(valid, 'positive robust joint perturbation')
        call robust%residual(point - 1.e-4_dp*tangent, low, valid)
        call require(valid, 'negative robust joint perturbation')
        relative = norm2((high - low)/2.e-4_dp - product)/norm2(product)
        print '(A,ES14.6)', 'robust/discrepancy joint derivative relative error=', relative
        call require(relative < 5.e-4_dp, 'robust/discrepancy joint derivative')
    end subroutine
end program
