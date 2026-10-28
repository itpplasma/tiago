program test_joint_recorded
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use joint_recorded_fixture, only: configure_joint, require
    use joint_signal_io, only: signal_roundtrip
    use joint_validation, only: add_noise, write_profiles, read_start, write_calibrations, &
                                calibration_start
    use tiago_joint_model, only: joint_model_t, diagnostic_slot_t
    use tiago_energy_consistency, only: energy_consistency_t
    use tiago_reconstruction, only: rec
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    use fortnum_status, only: fortnum_status_t
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
