program test_profile_posterior
    !! Profile posterior in a fixed G-EQDSK equilibrium (analytic limited
    !! Cerfon-Freidberg case, tests/analytic): the exact Jacobian against
    !! central differences, and recovery of planted profile, alignment and
    !! calibration latents from noisy point measurements within four Laplace
    !! standard deviations.
    !! Usage: test_profile_posterior <geqdsk>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    use fortnum_status, only: fortnum_status_t
    use tiago_geqdsk_plasma, only: read_geqdsk_plasma
    use tiago_point_probe, only: point_probe_t
    use tiago_profile_diagnostic, only: Q_NE, Q_TE, Q_TI, plasma_sample_t
    use tiago_kinetic_profiles, only: P_NE, P_TE, P_TI, P_ZEFF
    use tiago_profile_posterior, only: profile_posterior_t
    use tiago_inference_engines, only: map_fit, laplace, laplace_t
    implicit none
    type(profile_posterior_t) :: model
    type(point_probe_t) :: thomson, cxrs
    type(rng_t) :: rng
    type(fortnum_status_t) :: status
    type(laplace_t) :: lap
    type(plasma_sample_t) :: sample
    character(len=1024) :: path
    character(len=:), allocatable :: message
    real(dp), allocatable :: truth(:), theta(:), r(:), rp(:), rm(:), jac(:, :), y(:), &
        eta(:), dx(:, :, :), z(:), g(:)
    real(dp) :: c, h, err, sd
    integer :: ierr, i, k, steps, n
    logical :: ok
    real(dp), parameter :: breaks(6) = [0.0_dp, 0.3_dp, 0.6_dp, 0.85_dp, 0.95_dp, 1.0_dp]

    call get_command_argument(1, path)
    call read_geqdsk_plasma(trim(path), model%eq, ierr, message)
    if (ierr /= 0) error stop message

    call model%prof%init(breaks, 6.0_dp)
    g = model%prof%greville
    call model%prof%set_profile(P_NE, 5.0e19_dp * (1 - 0.7_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call model%prof%set_profile(P_TE, 3.0e3_dp * (1 - 0.9_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call model%prof%set_profile(P_TI, 2.5e3_dp * (1 - 0.9_dp * g**2), .true., &
        sigma_f=0.5_dp, l_core=0.3_dp, l_edge=0.05_dp, s0=0.85_dp, sw=0.05_dp)
    call model%prof%set_profile(P_ZEFF, [(1.5_dp, i = 1, size(g))], .false.)

    ! Vertical Thomson chord at R = 1.70 m, n_e and T_e at 12 volumes.
    n = 12
    allocate(thomson%x0(3, 2 * n), thomson%quantity(2 * n))
    do i = 1, n
        thomson%x0(:, i) = [1.70_dp, 0.0_dp, 0.05_dp + 0.80_dp * (i - 1) / (n - 1)]
        thomson%x0(:, n + i) = thomson%x0(:, i)
    end do
    thomson%quantity = [(Q_NE, i = 1, n), (Q_TE, i = 1, n)]
    thomson%aligned = .false.
    ! Outboard-midplane T_i at 8 points with a calibration.
    allocate(cxrs%x0(3, 8), cxrs%quantity(8))
    do i = 1, 8
        cxrs%x0(:, i) = [1.80_dp + 0.38_dp * (i - 1) / 7, 0.0_dp, 0.02_dp]
    end do
    cxrs%quantity = Q_TI
    cxrs%aligned = .false.

    allocate(model%slots(2))
    allocate(model%slots(1)%d, source=thomson)
    allocate(model%slots(2)%d, source=cxrs)
    model%slots(1)%eta_sigma = [0.01_dp, 0.01_dp]
    model%slots(2)%eta_sigma = [0.01_dp, 0.01_dp]
    model%slots(2)%cal_sigma = 0.05_dp
    do k = 1, 2
        allocate(model%slots(k)%data(model%slots(k)%d%n_signals()), source=1.0_dp)
        allocate(model%slots(k)%sigma(model%slots(k)%d%n_signals()), source=1.0_dp)
    end do
    call model%setup()

    ! Planted latents and noisy data (3 % noise).
    call rng_seed(rng, 7_8, status)
    allocate(truth(model%n))
    do i = 1, model%n
        call rng_normal(rng, truth(i))
        truth(i) = 0.7_dp * truth(i)
    end do
    do k = 1, 2
        call model%predict(truth, k, y, sample, eta, dx, c)
        do i = 1, size(y)
            model%slots(k)%sigma(i) = 0.03_dp * abs(c * y(i))
            call rng_normal(rng, h)
            model%slots(k)%data(i) = c * y(i) + model%slots(k)%sigma(i) * h
        end do
    end do

    ! Exact Jacobian against central differences at a point between prior and truth.
    theta = 0.5_dp * truth
    call model%jacobian(theta, jac)
    allocate(r(model%n_residuals()), rp(model%n_residuals()), rm(model%n_residuals()))
    err = 0
    do i = 1, model%n
        h = 1.0e-6_dp
        z = theta
        z(i) = theta(i) + h
        call model%residual(z, rp, ok)
        z(i) = theta(i) - h
        call model%residual(z, rm, ok)
        err = max(err, maxval(abs((rp - rm) / (2 * h) - jac(:, i))) / max(maxval(abs(jac(:, i))), 1.0_dp))
    end do
    print '(A,ES10.3)', 'max relative Jacobian error ', err
    if (err > 1.0e-5_dp) error stop 'Jacobian differs from central differences'

    theta = 0
    call map_fit(model, theta, 100, steps, message)
    call laplace(model, theta, lap)
    print '(A,I0,A,ES10.3)', 'MAP steps ', steps, ', chi2 ', lap%chi2
    do i = 1, model%n
        sd = sqrt(lap%cov(i, i))
        if (abs(theta(i) - truth(i)) > 4 * sd) then
            print '(A,I0,3ES12.4)', 'latent, MAP, truth, sd: ', i, theta(i), truth(i), sd
            error stop 'planted latent not recovered'
        end if
    end do
    print '(A)', 'tiago_profile_posterior: OK'
end program test_profile_posterior
