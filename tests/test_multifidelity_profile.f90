program test_multifidelity_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tiago_multifidelity_profile, only: ar1_profile_t
    use tiago_kinetic_profiles, only: kinetic_profiles_t, P_TE, N_PROFILES
    use tiago_calibration, only: calibration_prior_t
    use fortnum_status, only: fortnum_status_t
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    implicit none
    type(kinetic_profiles_t) :: profiles, low_profiles
    type(ar1_profile_t) :: prior
    type(calibration_prior_t) :: calibration
    type(fortnum_status_t) :: status
    type(rng_t) :: rng
    integer, parameter :: trials = 96, nobs = 13
    real(dp) :: s(nobs), data(2*nobs), truth_data(2*nobs), noise, logz, ctrue
    real(dp), allocatable :: truth(:), theta(:), covariance(:, :), gradient(:), hessian(:, :)
    real(dp), allocatable :: lo(:), hi(:), lj(:, :), hj(:, :), dt(:), dm(:), dc(:, :)
    real(dp) :: cwidth, dwidth, cpull(trials), dpull(trials), h, value, fp, fm
    integer :: i, j, code, n, o, nc, cindex, hindex, unit, covered(2)
    character(1024) :: output
    external :: dposv, dpotrf
    call profiles%init([0._dp, .4_dp, 1._dp])
    call calibration%init([1.15_dp], reshape([.2_dp**2], [1, 1]), code)
    call require(code == 0, 'calibration prior')
    calibration%record%record_id = 'synthetic-temperature-scale'
    call prior%init(profiles%greville, 1000 - 600*profiles%greville, calibration, &
        log([100._dp, .25_dp, 70._dp, .2_dp]), [.25_dp, .2_dp, .25_dp, .2_dp], &
        'synthetic lognormal hyperpriors', 'constant and calibration tangent', code)
    call require(code == 0, 'AR1 prior setup')
    call profiles%set_multifidelity_te(prior)
    low_profiles = profiles
    low_profiles%low_fidelity_te = .true.
    n = profiles%n_latent; nc = profiles%ncoef; o = profiles%te_ar1%offset
    cindex = o + prior%scale_index(); hindex = o + prior%hyper_index()
    allocate (truth(n), theta(n), covariance(n, n), gradient(n), hessian(n, n), &
        lo(nc), hi(nc), lj(nc, prior%n_params), hj(nc, prior%n_params), &
        dt(nc), dm(nc), dc(nc, n))
    s = [(real(i - 1, dp)/(nobs - 1), i=1, nobs)]
    truth = 0
    truth(:nc) = [(0.4_dp*sin(.7_dp*i), i=1, nc)]
    truth(nc + 1:nc + prior%n_shape) = &
        [(0.65_dp*cos(.7_dp*i), i=1, prior%n_shape)]
    truth(cindex) = .65_dp
    truth(hindex + 1:hindex + 4) = [.3_dp, -.3_dp, .2_dp, -.2_dp]
    call prior%coefficients(truth, lo, hi, status, lj, hj)
    call require(status%code == 0, 'physical coefficients')
    dt = hi - prior%scale(truth)*lo
    call require(abs(sum(dt)) < 1.e-10_dp, 'zero constant discrepancy')
    call require(abs(dot_product(prior%mean_low, dt)) < 1.e-7_dp, 'protected calibration tangent')
    h = 1.e-5_dp
    do j = 1, n
        theta = truth; theta(j) = theta(j) + h
        call prior%coefficients(theta, lo, dm, status)
        theta(j) = theta(j) - 2*h
        call prior%coefficients(theta, lo, hi, status)
        call require(maxval(abs((dm - hi)/(2*h) - hj(:, j))) < 2.e-6_dp, &
            'coefficient derivatives including calibration and hyperparameters')
    end do
    call profile_derivatives(profiles, truth)
    call profile_derivatives(low_profiles, truth)
    call single_fidelity()
    call forward(truth, truth_data)
    call rng_seed(rng, 20261030_int64, status)
    output = 'multifidelity-profile.csv'
    if (command_argument_count() > 0) call get_command_argument(1, output)
    open (newunit=unit, file=trim(output), status='replace')
    write (unit, '(A)') 'trial,C,C_sd,delta_eV,delta_sd_eV,log_evidence,C_true,delta_true_eV'
    covered = 0
    do i = 1, trials
        do j = 1, 2*nobs
            call rng_normal(rng, noise)
            data(j) = truth_data(j) + 10*noise
        end do
        call fit(data, theta, covariance, logz)
        ctrue = prior%scale(truth)
        cwidth = calibration%factor(1, 1)*sqrt(covariance(cindex, cindex))
        cpull(i) = (prior%scale(theta) - ctrue)/cwidth
        call prior%coefficients(theta, lo, hi, status, lj, hj)
        dm = hi - prior%scale(theta)*lo
        dc = hj - prior%scale(theta)*lj
        dc(:, cindex) = dc(:, cindex) - calibration%factor(1, 1)*lo
        j = maxloc(abs(dt), dim=1)
        dwidth = sqrt(dot_product(dc(j, :), matmul(covariance, dc(j, :))))
        dpull(i) = (dm(j) - dt(j))/dwidth
        if (abs(cpull(i)) < 1.95996398454_dp) covered(1) = covered(1) + 1
        if (abs(dpull(i)) < 1.95996398454_dp) covered(2) = covered(2) + 1
        if (i == 1) then
            call require(abs(cpull(i)) < 1.95996398454_dp, 'planted calibration in Laplace interval')
            call require(abs(dpull(i)) < 1.95996398454_dp, 'planted discrepancy in Laplace interval')
            call require(maxval(abs(theta(hindex + 1:hindex + 4))) > .01_dp, &
                'hyperparameters inferred with posterior')
        end if
        write (unit, '(I0,7(",",ES20.12))') &
            i, prior%scale(theta), cwidth, dm(j), dwidth, logz, ctrue, dt(j)
    end do
    close (unit)
    print '(A,2I5,A,I0)', '95 percent calibration/discrepancy coverage=', covered, ' / ', trials
    print '(A,2F12.6)', 'mean pulls=', sum(cpull)/trials, sum(dpull)/trials
    print '(A,2F12.6)', 'RMS pulls=', sqrt(sum(cpull**2)/trials), sqrt(sum(dpull**2)/trials)
    call require(all(covered >= 85), 'repeated-noise 95 percent coverage lower binomial bound')
    call require(sqrt(sum(cpull**2)/trials) > .5_dp, 'calibration interval scale')
    call require(sqrt(sum(dpull**2)/trials) > .5_dp, 'discrepancy interval scale')
    call require(abs(sum(cpull)/trials) < .4_dp, 'calibration bias')
    call require(abs(sum(dpull)/trials) < .4_dp, 'discrepancy bias')
contains
    subroutine require(ok, text)
        logical, intent(in) :: ok
        character(*), intent(in) :: text
        if (.not. ok) then
            print *, text
            error stop 1
        end if
    end subroutine

    subroutine forward(z, y, jac)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: y(:)
        real(dp), intent(out), optional :: jac(:, :)
        real(dp) :: f(N_PROFILES, nobs), df(N_PROFILES, nobs), fb(N_PROFILES, nobs)
        real(dp) :: zb(size(z)), sb(nobs)
        integer :: k
        call profiles%evaluate(z, s, f, df)
        y(:nobs) = f(P_TE, :)
        call low_profiles%evaluate(z, s, f, df)
        y(nobs + 1:) = f(P_TE, :)
        if (.not. present(jac)) return
        do k = 1, nobs
            fb = 0; fb(P_TE, k) = 1
            zb = 0
            call profiles%vjp(z, s, fb, zb, sb)
            jac(k, :) = zb
            zb = 0
            call low_profiles%vjp(z, s, fb, zb, sb)
            jac(nobs + k, :) = zb
        end do
    end subroutine

    subroutine objective(z, observations, f, g, fisher)
        real(dp), intent(in) :: z(:), observations(:)
        real(dp), intent(out) :: f, g(:), fisher(:, :)
        real(dp) :: y(2*nobs), jac(2*nobs, size(z)), r(2*nobs)
        integer :: k
        call forward(z, y, jac)
        r = (y - observations)/10
        jac = jac/10
        f = (sum(r**2) + sum(z**2))/2
        g = matmul(r, jac) + z
        fisher = matmul(transpose(jac), jac)
        do k = 1, size(z)
            fisher(k, k) = fisher(k, k) + 1
        end do
    end subroutine

    subroutine fit(observations, z, cov, log_evidence)
        real(dp), intent(in) :: observations(:)
        real(dp), intent(out) :: z(:), cov(:, :), log_evidence
        real(dp) :: g(size(z)), a(size(z), size(z)), f, trial_f, trial_g(size(z))
        real(dp) :: trial_a(size(z), size(z)), step(size(z)), m(size(z), size(z))
        real(dp) :: zp(size(z)), gm(size(z)), gp(size(z)), fraction, eps, ld
        integer :: k, it, ls, info
        z = 0
        do it = 1, 100
            call objective(z, observations, f, g, a)
            m = a; step = -g
            call dposv('U', size(z), 1, m, size(z), step, size(z), info)
            call require(info == 0, 'MAP step')
            if (maxval(abs(step)) < 2.e-7_dp) exit
            fraction = 1
            do ls = 1, 20
                call objective(z + fraction*step, observations, trial_f, trial_g, trial_a)
                if (trial_f < f) exit
                fraction = fraction/2
            end do
            call require(ls <= 20, 'MAP line search')
            z = z + fraction*step
        end do
        call require(it <= 100, 'MAP convergence')
        call objective(z, observations, f, g, a)
        call require(maxval(abs(g)) < 1.e-4_dp, 'MAP score')
        do k = 1, size(z)
            eps = 2.e-5_dp*max(1._dp, abs(z(k)))
            zp = z; zp(k) = zp(k) + eps
            call objective(zp, observations, trial_f, gp, trial_a)
            zp(k) = zp(k) - 2*eps
            call objective(zp, observations, trial_f, gm, trial_a)
            a(:, k) = (gp - gm)/(2*eps)
        end do
        a = (a + transpose(a))/2
        m = a
        call dpotrf('U', size(z), m, size(z), info)
        call require(info == 0, 'observed posterior Hessian positive')
        ld = 2*sum([(log(m(k, k)), k=1, size(z))])
        log_evidence = -f - 2*nobs*log(10._dp) - nobs*log(2*acos(-1._dp)) - ld/2
        cov = 0
        do k = 1, size(z)
            cov(k, k) = 1
        end do
        call dposv('U', size(z), size(z), a, size(z), cov, size(z), info)
        call require(info == 0, 'Laplace covariance')
    end subroutine

    subroutine profile_derivatives(p, z)
        type(kinetic_profiles_t), intent(in) :: p
        real(dp), intent(in) :: z(:)
        real(dp) :: zb(size(z)), sb(4), fb(N_PROFILES, 4), f(N_PROFILES, 4), df(N_PROFILES, 4)
        real(dp) :: fp(N_PROFILES, 4), fm(N_PROFILES, 4), points(4), d(size(z)), eps
        integer :: k
        points = [.13_dp, .4_dp, .7_dp, 1.1_dp]; eps = 1.e-6_dp
        d = [(cos(.3_dp*k), k=1, size(z))]
        fb = 0; fb(P_TE, :) = [.2_dp, -.3_dp, .4_dp, .1_dp]
        zb = 0
        call p%vjp(z, points, fb, zb, sb)
        call p%evaluate(z + eps*d, points, fp, df)
        call p%evaluate(z - eps*d, points, fm, df)
        call require(abs(sum(fb*(fp - fm))/(2*eps) - dot_product(zb, d)) < 1.e-5_dp, &
            'profile VJP including hyperparameters')
        call p%evaluate(z, points + eps, fp, df)
        call p%evaluate(z, points - eps, fm, df)
        call require(abs(sum(fb*(fp - fm))/(2*eps) - sum(sb)) < 1.e-5_dp, 'profile spatial VJP')
        call p%evaluate(z, points, f, df)
        call require(maxval(abs((fp(P_TE, :) - fm(P_TE, :))/(2*eps) - df(P_TE, :))) &
            < 1.e-5_dp, 'physical profile slope')
    end subroutine

    subroutine single_fidelity()
        use fortnum_ar1, only: rbf_covariance, covariance_factor
        type(ar1_profile_t) :: p
        real(dp), allocatable :: z(:), low(:), high(:), k(:, :), dk(:, :, :), l(:, :), dl(:, :, :)
        integer :: code
        call p%init(profiles%greville, prior%mean_low, calibration, prior%log_mean, prior%log_sd, &
            'synthetic', 'synthetic', code, discrepancy=.false., fixed_scale=.true.)
        call require(code == 0, 'single fidelity setup')
        allocate (z(p%n_params), low(p%n), high(p%n), k(p%n, p%n), dk(p%n, p%n, 2), &
            l(p%n, p%n), dl(p%n, p%n, 2))
        z = .3_dp
        z(p%hyper_index() + 1:p%hyper_index() + 2) = 0
        call p%coefficients(z, low, high, status)
        call require(status%code == 0, 'single fidelity coefficients')
        call rbf_covariance(p%nodes, p%nodes, exp(p%log_mean(1)), exp(p%log_mean(2)), k, dk, status)
        call covariance_factor(k, dk, l, dl, status)
        call require(maxval(abs(low - p%mean_low - matmul(l, z(:p%n)))) < 1.e-12_dp, &
            'single fidelity GP exact coefficients')
        call require(maxval(abs(high - calibration%mean(1)*low)) < 1.e-12_dp, &
            'fixed scale zero discrepancy reduction')
        call require(p%n_params == p%n + 2, 'no unused discrepancy or calibration parameters')
    end subroutine
end program
