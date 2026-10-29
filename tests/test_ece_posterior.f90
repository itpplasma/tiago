program test_ece_posterior
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use ece_recorded_fixture, only: configure
    use tiago_ece_state_medium, only: state_medium_t
    use tiago_ece_radiometer, only: ece_t, ece_channel_t
    use tiago_ece_recorded, only: recorded_ece_t, record_ece
    use tiago_ece_receiver, only: ece_response_channel_t
    use tiago_calibration, only: calibration_prior_t
    use tiago_kinetic_profiles, only: P_TE
    use fortnum_status, only: fortnum_status_t
    use fortnum_rng, only: rng_t, rng_seed, rng_normal
    implicit none
    integer, parameter :: nc = 6, nd = 3, trials = 128
    type(state_medium_t) :: base
    type(ece_t) :: diagnostic
    type(recorded_ece_t) :: truth
    type(fortnum_status_t) :: status
    type(rng_t) :: rng
    real(dp), allocatable :: p0(:), pd(:, :), ed(:, :), zd(:)
    real(dp) :: theta(nd), start(nd), true_theta(nd), noise, sd(nd), pull(nd), &
        mean_pull(nd), mean_square(nd), error(nd), estimate(nd), &
        reported_sd(nd), truth_signal(nc), truth_variance(nc), gain, te, &
        data(nc, trials), estimates(nd, trials), widths(nd, trials)
    integer :: i, j, count(nd), unit, ns, np
    character(1024) :: output
    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
    call configure(base, diagnostic, nc, 48, 1.0_dp, .true.)
    call setup_response()
    ns = base%state%n_params()
    np = base%n_params()
    allocate (p0(ns), pd(np, nd), ed(5*nc, nd), zd(base%profiles%n_latent))
    call base%state%get_params(p0)
    pd = 0
    ed = 0
    zd = 0
    call base%profiles%fit_latent(P_TE, base%profiles%greville, &
        1000*exp(-.5_dp*base%profiles%greville + &
        .1_dp), zd)
    zd = zd/norm2(zd)
    pd(ns + 1:, 1) = zd
    do i = 1, nc
        ed(3*i - 2, 2) = 1/sqrt(real(nc, dp))
    end do
    ! The final block contains contravariant toroidal field coefficients.
    pd(ns - 7:ns, 3) = .001_dp*p0(ns - 7:ns)
    true_theta = [.6_dp, 1.0_dp, .8_dp]
    call forward(true_theta, truth, status)
    if (status%code /= 0) print *, status%msg
    call require(status%code == 0, 'posterior truth')
    truth_signal = truth%voltage
    truth_variance = truth%variance
    call rng_seed(rng, 20261028_int64, status)
    call require(status%code == 0, 'posterior seed')
    output = 'ece-posterior.csv'
    if (command_argument_count() > 0) call get_command_argument(1, output)
    open (newunit=unit, file=trim(output), status='replace')
    write (unit, &
        '(A)') 'trial,Te_latent,gain_latent,B_shift_latent,sd_Te,sd_gain,sd_shift'
    count = 0
    mean_pull = 0
    mean_square = 0
    estimate = 0
    reported_sd = 0
    do i = 1, trials
        do j = 1, nc
            call rng_normal(rng, noise)
            data(j, i) = truth_signal(j) + sqrt(truth_variance(j) + &
                .002_dp**2)*noise
        end do
    end do
    ! Starting every trial at the maximum for the averaged data, itself found
    ! from the prior mean, saves the far-field Fisher steps of a cold start.
    call fit(sum(data, 2)/trials, [(0.0_dp, j=1, nd)], start, sd)
    !$omp parallel do default(shared) private(i, theta, sd, error, pull) &
    !$omp reduction(+:count, mean_pull, mean_square, estimate, reported_sd)
    do i = 1, trials
        call fit(data(:, i), start, theta, sd)
        estimates(:, i) = theta
        widths(:, i) = sd
        error = theta - true_theta
        pull = error/sd
        where (abs(pull) <= 1.959963984540054_dp) count = count + 1
        mean_pull = mean_pull + pull/trials
        mean_square = mean_square + pull**2/trials
        estimate = estimate + theta/trials
        reported_sd = reported_sd + sd/trials
    end do
    !$omp end parallel do
    do i = 1, trials
        write (unit, '(I0,6(",",ES18.10))') i, estimates(:, i), widths(:, i)
    end do
    close (unit)
    print '(A,3F10.5)', 'truth latents (Te, gain, B shift)=', true_theta
    print '(A,3F10.5)', 'mean recovered latents=', estimate
    print '(A,3F10.5)', 'mean posterior sd=', reported_sd
    print '(A,3F10.5)', 'mean pulls=', mean_pull
    print '(A,3F10.5)', 'rms pulls=', sqrt(mean_square)
    print '(A,3I5,A,I0)', '95% interval coverage counts=', count, ' / ', trials
    gain = 1 + .1_dp*true_theta(2)/sqrt(real(nc, dp))
    te = .001_dp*true_theta(3)
    print '(A,F10.6,A,F10.6)', 'true gain multiplier=', gain, &
        ' true relative field shift=', te
    call require(all(count >= 115), 'at least 90% repeated-noise posterior coverage')
    call require(all(abs(estimate - true_theta) < reported_sd), &
        'mean truth within posterior sd')
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine

    subroutine fit(observations, initial, theta, sd)
        real(dp), intent(in) :: observations(nc), initial(nd)
        real(dp), intent(out) :: theta(nd), sd(nd)
        real(dp) :: covariance(nd, nd), fisher(nd, nd), gradient(nd), &
            next_gradient(nd), next_fisher(nd, nd), step(nd, 1), &
            inverse(nd, nd), value, next_value, fraction
        type(recorded_ece_t) :: tape
        integer :: iteration, j, k, info
        theta = initial
        call objective(theta, observations, value, gradient, fisher)
        do iteration = 1, 24
            covariance = fisher
            step(:, 1) = -gradient
            call dposv('L', nd, 1, covariance, nd, step, nd, info)
            call require(info == 0, 'posterior positive Fisher matrix')
            if (maxval(abs(step)) < 2.0e-6_dp) exit
            ! Ray/coordinate roundoff can exceed the likelihood decrement
            ! near the root; require the score to decrease there instead.
            if (abs(dot_product(gradient, step(:, 1))) < 1.0e-8_dp) then
                call objective(theta + step(:, 1), observations, next_value, &
                    next_gradient, next_fisher)
                if (norm2(next_gradient) < norm2(gradient)) then
                    theta = theta + step(:, 1)
                    value = next_value
                    gradient = next_gradient
                    fisher = next_fisher
                    cycle
                end if
            end if
            fraction = 1
            do k = 1, 16
                call evaluate(theta + fraction*step(:, 1), observations, tape, &
                    next_value)
                if (next_value <= value) exit
                fraction = fraction/2
            end do
            call require(k <= 16, 'posterior line search')
            theta = theta + fraction*step(:, 1)
            value = next_value
            call derivatives(tape, theta, observations, gradient, fisher)
        end do
        call require(iteration <= 24, 'posterior convergence')
        inverse = 0
        do j = 1, nd
            inverse(j, j) = 1
        end do
        call dposv('L', nd, nd, fisher, nd, inverse, nd, info)
        call require(info == 0, 'posterior covariance')
        sd = [(sqrt(inverse(j, j)), j=1, nd)]
    end subroutine

    subroutine setup_response()
        type(ece_channel_t) :: channels(nc)
        type(ece_response_channel_t) :: response(nc)
        type(calibration_prior_t) :: prior, walls
        integer :: j, code
        channels = diagnostic%channels
        do j = 1, nc
            channels(j)%frequency = 2*1.7588200107721631e11_dp* &
                .3_dp*(3.50_dp + .028_dp*(j - &
                1))/(2*acos(-1.0_dp))
            response(j)%rays = reshape([j], [1, 1])
            response(j)%frequency = [channels(j)%frequency]
            response(j)%weights = [1.0e8_dp]
            response(j)%response = [1.0_dp]
            response(j)%antenna = [1.0_dp]
            response(j)%duration = .01_dp
        end do
        prior = diagnostic%calibration
        walls = diagnostic%wall_prior
        call diagnostic%init(channels, response, prior, walls, code)
        call require(code == 0, 'posterior response')
    end subroutine

    subroutine forward(t, tape, st)
        real(dp), intent(in) :: t(nd)
        type(recorded_ece_t), intent(out) :: tape
        type(fortnum_status_t), intent(out) :: st
        type(state_medium_t) :: medium
        medium = base
        call medium%state%set_params(p0 + matmul(pd(:ns, :), t))
        medium%z = base%z + matmul(pd(ns + 1:, :), t)
        call record_ece(diagnostic, medium, matmul(ed, t), tape, st)
    end subroutine

    subroutine objective(t, data, result, g, f)
        real(dp), intent(in) :: t(nd), data(nc)
        real(dp), intent(out) :: result
        real(dp), intent(out), optional :: g(nd), f(nd, nd)
        type(recorded_ece_t) :: tape
        call evaluate(t, data, tape, result)
        if (present(g)) call derivatives(tape, t, data, g, f)
    end subroutine

    subroutine evaluate(t, data, tape, result)
        real(dp), intent(in) :: t(nd), data(nc)
        type(recorded_ece_t), intent(out) :: tape
        real(dp), intent(out) :: result
        type(fortnum_status_t) :: st
        real(dp) :: v(nc), r(nc)
        call forward(t, tape, st)
        if (st%code /= 0) print *, st%msg
        call require(st%code == 0, 'posterior forward')
        v = tape%variance + .002_dp**2
        r = tape%voltage - data
        ! The normalization is essential for parameter-dependent noise.
        result = .5_dp*(sum(r**2/v + log(v)) + sum(t**2))
    end subroutine

    subroutine derivatives(tape, t, data, g, f)
        type(recorded_ece_t), intent(in) :: tape
        real(dp), intent(in) :: t(nd), data(nc)
        real(dp), intent(out) :: g(nd)
        real(dp), intent(out), optional :: f(nd, nd)
        type(fortnum_status_t) :: st
        real(dp) :: v(nc), r(nc), jv(nc, nd), jn(nc, nd)
        integer :: j, k
        v = tape%variance + .002_dp**2
        r = tape%voltage - data
        do j = 1, nd
            call tape%jvp(pd(:, j), ed(:, j), jv(:, j), jn(:, j), st)
            call require(st%code == 0, 'posterior sensitivity')
        end do
        g = matmul(r/v, jv) + .5_dp*matmul((1 - r**2/v)/v, jn) + t
        if (.not. present(f)) return
        do j = 1, nd
            do k = 1, nd
                f(j, k) = sum(jv(:, j)*jv(:, k)/v) + &
                    .5_dp*sum(jn(:, j)*jn(:, k)/v**2)
            end do
            f(j, j) = f(j, j) + 1
        end do
    end subroutine
end program
