program test_calibration
    !! Log-scale calibration latents of the joint model on pinned LI383:
    !!  1. the Student-t prior row of the campaign latent: r^2 / 2 equals
    !!     (nu + 1) / 2 log(1 + z^2 / nu), the log_normalization change from
    !!     the Gaussian equals the quadrature of exp(-r^2 / 2), and the prior
    !!     row's derivative agrees with central differences also in the series
    !!     branch near z = 0;
    !!  2. calibration layout: Thomson with three lasers (campaign latent with
    !!     Student-t prior, two Helmert contrasts, sum of offsets zero), an
    !!     interferometer with one shared systematic, diamagnetic loops only;
    !!     Te signals are not scaled;
    !!  3. the joint Jacobian columns of the calibration latents, a profile
    !!     latent and CURTOR agree with central differences of re-solved
    !!     equilibria;
    !!  4. noise-free data at a truth with C = 1.8: the residual at the truth
    !!     vanishes in the data rows, while C = 1 leaves a residual.
    !! Usage: test_calibration <VMEC input> <sensors.awk> <workdir>
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use tiago_reconstruction, only: rec
    use tiago_joint_model, only: joint_model_t
    use tiago_synthetic_case, only: setup_magnetics, place_instruments, synthesize_magnetics, &
        synthesize_profiles
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    implicit none
    character(len=512) :: input, awk, work
    type(joint_model_t), target :: model
    type(rng_t) :: rng
    type(fortnum_status_t) :: st
    real(dp), allocatable :: theta0(:), theta_true(:), theta(:), r(:), rp(:), rm(:), jac(:, :), &
        colfd(:), b(:), c(:)
    real(dp) :: h, err, log_c, total, zq, nu, ndata
    integer :: i, k, col, fails, nrow
    integer, allocatable :: cols(:)
    logical :: ok

    fails = 0
    call get_command_argument(1, input)
    call get_command_argument(2, awk)
    call get_command_argument(3, work)

    ! --- 2. layout
    call execute_command_line('mkdir -p "'//trim(work)//'" && awk -v R0=1.42 -v rs=0.85 -v nphi=4 '// &
        '-v nth=4 -v nprobe=12 -v dir="'//trim(work)//'" -f "'//trim(awk)//'"')
    call setup_magnetics(trim(input), trim(work)//'/work', trim(work), &
        [character(len=64) :: 'CURTOR', 'PHIEDGE'], [-174250.0_dp, 0.514386_dp], &
        [2.0e4_dp, 0.02_dp], 1.0e-14_dp, 'diamagnetic')
    call expect(rec%n_meas == 4 .and. all(index(rec%meas_label, 'DIA_') == 1), &
        'magnetics reduced to the diamagnetic loops')
    call model%prof%init([0.0_dp, 0.2_dp, 0.45_dp, 0.7_dp, 1.0_dp])
    associate (g => model%prof%greville)
        call model%prof%set_profile(P_NE, 1.0e20_dp * (1.0_dp - 0.6_dp * g), .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_TE, 2.5e3_dp * (1.0_dp - g) + 50.0_dp, .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_TI, 2.0e3_dp * (1.0_dp - g) + 50.0_dp, .true., &
            0.3_dp, 0.4_dp, 0.2_dp, 0.6_dp, 0.2_dp)
        call model%prof%set_profile(P_ZEFF, 1.5_dp + 0.0_dp * g, .false.)
    end associate
    allocate(model%diags(0))
    call model%setup()
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n_prof)]]
    allocate(r(model%n_residuals()))
    call model%residual(theta0, r, ok)
    if (.not. ok) error stop 'no equilibrium at the prior mean'
    call place_instruments(model, [0.3_dp, 0.6_dp], 0.3_dp, 0.0_dp, 3, 0.1_dp, 0.5_dp, 4.0_dp, &
        0.15_dp, 2, 0.03_dp, .false., 0.0_dp)
    call model%setup()
    call expect(model%n_cal == 4 .and. model%diags(1)%n_laser == 3 .and. &
        model%diags(1)%laser_offset == model%diags(1)%cal_index .and. &
        model%diags(2)%cal_index == model%diags(1)%cal_index + 3 .and. model%n == &
        model%n_rec + model%n_prof + 4, 'calibration latent layout')
    theta0 = [rec%prior / rec%scale, [(0.0_dp, i = 1, model%n - model%n_rec)]]
    theta_true = theta0 + [0.03_dp, -0.02_dp, [(0.6_dp * sin(1.7_dp * i), i = 1, model%n_prof)], &
        (log(1.8_dp) - 0.1_dp) / 0.5_dp, 0.4_dp, -0.7_dp, 0.8_dp]
    call model%calibration(theta_true, 1, log_c, b)
    call expect(abs(log_c - log(1.8_dp)) < 1.0e-14_dp .and. abs(sum(b)) < 1.0e-14_dp, &
        'log C and zero-sum laser offsets')
    c = model%signal_scale(theta_true, 1)
    k = size(c) / 2
    call expect(all(c(k + 1:) == 1.0_dp) .and. abs(c(1) - exp(log_c + b(1))) < 1.0e-14_dp .and. &
        abs(c(k) - exp(log_c + b(3))) < 1.0e-14_dp, 'ne scaled per laser, Te not scaled')

    ! --- 4. noise-free synthetic data at the truth
    deallocate(r)
    allocate(r(model%n_residuals()), rp(model%n_residuals()), rm(model%n_residuals()))
    call rng_seed(rng, 1_int64, st)
    call model%residual(theta_true, r, ok)
    call synthesize_magnetics(model, theta_true(1:2) * rec%scale, 0.01_dp, 1.0e-4_dp, 1.0e-8_dp, &
        1.0e-7_dp, rng, .false.)
    call synthesize_profiles(model, theta_true, 0.05_dp, [0.0_dp, 0.0_dp], rng, .false., &
        [0.05_dp, 0.01_dp])
    call model%residual(theta_true, r, ok)
    ndata = model%n_mag + model%n_rec_prior + model%n_diag_rows
    nrow = model%n_mag + model%n_fb + model%n_rec_prior
    print '(A,ES10.2)', 'data residual at the truth: ', maxval(abs(r(nrow + 1:nrow + model%n_diag_rows)))
    call expect(maxval(abs(r(nrow + 1:nrow + model%n_diag_rows))) < 1.0e-9_dp, &
        'calibrated data reproduced at the truth')
    theta = theta_true
    theta(model%diags(1)%cal_index) = -0.1_dp / 0.5_dp
    call model%residual(theta, rp, ok)
    print '(A,ES10.2)', 'Thomson residual at C = 1: ', &
        maxval(abs(rp(nrow + 1:nrow + model%diags(1)%d%n_signals())))
    call expect(maxval(abs(rp(nrow + 1:nrow + model%diags(1)%d%n_signals()))) > 5.0_dp, &
        'C = 1 does not fit data made with C = 1.8')

    ! --- 3. Jacobian columns against central differences (cold solves)
    call rec%eq%vmec%set_input('hot_restart', 0.0_dp)
    theta = 0.5_dp * (theta0 + theta_true)
    theta(model%diags(1)%cal_index) = 1.3_dp   ! Student-t row away from the series
    call model%residual(theta, r, ok)
    call model%jacobian(theta, jac)
    cols = [1, model%n_rec + 2, model%n - 3, model%n - 2, model%n - 1, model%n]
    h = 1.0e-4_dp
    do k = 1, size(cols)
        col = cols(k)
        theta(col) = theta(col) + h
        call model%residual(theta, rp, ok)
        theta(col) = theta(col) - 2.0_dp * h
        call model%residual(theta, rm, ok)
        theta(col) = theta(col) + h
        colfd = (rp - rm) / (2.0_dp * h)
        err = maxval(abs(colfd - jac(:, col))) / max(maxval(abs(colfd)), 1.0e-12_dp)
        print '(A,I3,A,ES10.2,A,ES10.2)', 'column ', col, ': max |J - FD| / max |FD| = ', err, &
            '   max |FD| = ', maxval(abs(colfd))
        call expect(err < 2.0e-4_dp, 'Jacobian column vs central differences')
    end do

    ! --- 1. Student-t prior row through the model (cheap: same equilibrium)
    nu = model%diags(1)%cal_nu
    k = model%diags(1)%cal_index
    nrow = model%n_mag + model%n_fb + model%n_rec_prior + model%n_diag_rows + k - model%n_rec
    err = 0.0_dp
    do i = -6, 6
        theta(k) = 0.7_dp * i + 2.0e-3_dp
        call model%residual(theta, r, ok)
        err = max(err, abs(r(nrow)**2 / 2.0_dp - 0.5_dp * (nu + 1.0_dp) * &
            log(1.0_dp + theta(k)**2 / nu)) / max(1.0_dp, r(nrow)**2))
        call model%jacobian(theta, jac)
        h = 1.0e-6_dp
        theta(k) = theta(k) + h
        call model%residual(theta, rp, ok)
        theta(k) = theta(k) - 2.0_dp * h
        call model%residual(theta, rm, ok)
        theta(k) = theta(k) + h
        err = max(err, abs((rp(nrow) - rm(nrow)) / (2.0_dp * h) - jac(nrow, k)))
        if (sign(1.0_dp, r(nrow)) /= sign(1.0_dp, theta(k))) err = huge(1.0_dp)
    end do
    print '(A,ES10.2)', 'Student-t prior row: max error of r^2/2 and dr/dz ', err
    call expect(err < 1.0e-7_dp, 'Student-t prior row and its derivative')
    total = 0.0_dp
    do i = -400000, 400000
        zq = 5.0e-3_dp * i
        total = total + 5.0e-3_dp * (1.0_dp + zq * zq / nu)**(-0.5_dp * (nu + 1.0_dp))
    end do
    log_c = model%log_normalization()
    model%diags(1)%cal_nu = 0.0_dp
    log_c = log_c - model%log_normalization()
    model%diags(1)%cal_nu = nu
    print '(A,2F14.10)', 'Student-t normalization change, quadrature: ', log_c, &
        0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(total)
    call expect(abs(log_c - (0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(total))) < 1.0e-8_dp, &
        'Student-t normalization')
    if (fails > 0) error stop 1
    print '(A)', 'calibration: all checks passed'
contains
    subroutine expect(cond, what)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: what
        if (.not. cond) then
            fails = fails + 1
            write(error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine expect
end program test_calibration
