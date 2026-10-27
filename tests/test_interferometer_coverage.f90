program test_interferometer_coverage
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tiago_dispersion_interferometer, only: dispersion_interferometer_t
    use tiago_profile_diagnostic, only: plasma_sample_t, N_QUANTITIES, Q_NE
    use tiago_sight_geometry, only: chord_t
    use tiago_calibration, only: calibration_prior_t
    implicit none
    integer, parameter :: nc = 48, np = 9, ns = 400
    real(dp), parameter :: sigma = 0.01_dp
    type(dispersion_interferometer_t) :: model
    type(calibration_prior_t) :: prior
    type(chord_t) :: chords(nc)
    type(plasma_sample_t) :: smp
    real(dp), allocatable :: x(:, :), dx(:, :, :), qdot(:, :)
    real(dp) :: mean(8), cov(8, 8), sd(8), neutral(nc), env(nc), base(nc)
    real(dp) :: truth(np), fit(np), data(2*nc), y(2*nc), j(2*nc, np)
    real(dp) :: hess(np, np), post(np, np), rhs(np, 1), delta(np), dev(np)
    real(dp) :: count1(np), count2(np), pulls(np), pull2(np), contraction(np)
    real(dp) :: etadot(8), density_dir(2*nc), error, value, trial(np), trial_y(2*nc)
    integer :: status, a, b, seedsize, k, s, it, info, ls
    integer, allocatable :: seed(:)
    call random_seed(size=seedsize)
    allocate (seed(seedsize))
    seed = [(731 + 37*k, k=1, seedsize)]
    call random_seed(put=seed)
    mean = [1.0e-5_dp, 1.1_dp, 0.9_dp, 0.17_dp, -0.23_dp, 0.3_dp, 1.0_dp, 0.2_dp]
    sd = [1.0e-8_dp, 0.03_dp, 0.03_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.05_dp, 0.1_dp]
    cov = 0
    do k = 1, 8
        cov(k, k) = sd(k)**2
    end do
    cov(4, 5) = 0.002_dp
    cov(5, 4) = cov(4, 5)
    call prior%init(mean, cov, status)
    if (status /= 0) error stop 'coverage prior'
    do k = 1, nc
        chords(k) = chord_t('synthetic', [3.2_dp, 0.0_dp, 0.0_dp], &
                            [0.0_dp, 0.0_dp, 1.0_dp], 1.0_dp, 1, 4)
        base(k) = 1.0e20_dp*(0.2_dp + 4.0_dp*k/nc)
        neutral(k) = 1.0e-7_dp*(1 + sin(0.71_dp*k))
        env(k) = cos(0.23_dp*k)
    end do
    call model%init(chords, prior, neutral, env, status)
    if (status /= 0) error stop 'coverage model'
    call model%sample_points(spread(0.0_dp, 1, 8), x, dx)
    smp%n = size(x, 2)
    allocate (smp%q(N_QUANTITIES, smp%n), qdot(N_QUANTITIES, smp%n))
    smp%q = 0
    count1 = 0
    count2 = 0
    pulls = 0
    pull2 = 0
    contraction = 0
    do s = 1, ns
        do a = 1, np
            truth(a) = normal()
        end do
        call predict(truth, y)
        do k = 1, 2*nc
            data(k) = y(k) + sigma*normal()
        end do
        fit = 0
        do it = 1, 30
            call predict(fit, y)
            qdot = 0
            do k = 1, smp%n
                b = mod(model%integral%owner(k) - 1, nc) + 1
                qdot(Q_NE, k) = 0.03_dp*base(b)
            end do
            etadot = 0
            call model%jvp(smp, fit(2:), qdot, etadot, density_dir)
            j(:, 1) = density_dir/sigma
            qdot = 0
            do a = 1, 8
                etadot = 0
                etadot(a) = 1
                call model%jvp(smp, fit(2:), qdot, etadot, j(:, a + 1))
            end do
            j(:, 2:) = j(:, 2:)/sigma
            hess = matmul(transpose(j), j)
            do a = 1, np
                hess(a, a) = hess(a, a) + 1
            end do
            post = hess
            rhs(:, 1) = matmul(transpose(j), (data - y)/sigma) - fit
            call dposv('L', np, 1, post, np, rhs, np, info)
            if (info /= 0) error stop 'coverage solve'
            delta = rhs(:, 1)
            if (maxval(abs(delta)) < 1.0e-8_dp) exit
            value = sum(((data - y)/sigma)**2) + sum(fit**2)
            do ls = 0, 15
                trial = fit + delta*0.5_dp**ls
                call predict(trial, trial_y)
                error = sum(((data - trial_y)/sigma)**2) + sum(trial**2)
                if (error <= value) exit
            end do
            if (ls > 15) then
                if (maxval(abs(delta)) < 1.0e-6_dp) exit
                error stop 'coverage line search'
            end if
            fit = trial
        end do
        if (it > 30) error stop 'coverage did not converge'
        post = hess
        call dpotrf('L', np, post, np, info)
        if (info /= 0) error stop 'coverage factor'
        call dpotri('L', np, post, np, info)
        if (info /= 0) error stop 'coverage inverse'
        do a = 1, np
            dev(a) = sqrt(post(a, a))
        end do
        where (abs(fit - truth) < dev) count1 = count1 + 1
        where (abs(fit - truth) < 2*dev) count2 = count2 + 1
        pulls = pulls + (fit - truth)/dev
        pull2 = pull2 + ((fit - truth)/dev)**2
        contraction = contraction + dev
    end do
    print '(A)', 'prior predictive Gaussian noise, 400 independent realizations'
    print '(A)', 'parameter,68pct,95pct,pull_mean,pull_rms,posterior_sd/prior_sd'
    do a = 1, np
        print '(I0,5(A,F9.5))', a, ',', count1(a)/ns, ',', count2(a)/ns, &
            ',', pulls(a)/ns, ',', sqrt(pull2(a)/ns), ',', contraction(a)/ns
        ! Four binomial standard deviations for each marginal; no pooled claim.
        if (abs(count1(a)/ns - 0.682689_dp) > 0.095_dp) error stop '68% coverage'
        if (abs(count2(a)/ns - 0.954500_dp) > 0.045_dp) error stop '95% coverage'
        if (abs(pulls(a)/ns) > 0.2_dp) error stop 'biased coverage pull'
    end do
    if (contraction(1)/ns > 0.1_dp) error stop 'density not recovered'
    if (maxval(contraction(5:6))/ns > 0.1_dp) error stop 'I/Q offsets not recovered'
contains
    subroutine predict(z, prediction)
        real(dp), intent(in) :: z(np)
        real(dp), intent(out) :: prediction(2*nc)
        integer :: p, c
        do p = 1, smp%n
            c = mod(model%integral%owner(p) - 1, nc) + 1
            smp%q(Q_NE, p) = (1 + 0.03_dp*z(1))*base(c)
        end do
        call model%predict(smp, z(2:), prediction)
    end subroutine

    real(dp) function normal()
        real(dp) :: u(2)
        call random_number(u)
        normal = sqrt(-2*log(max(u(1), tiny(1.0_dp))))*cos(2*acos(-1.0_dp)*u(2))
    end function
end program
